"""Run steps 1--3 from source inputs, with stage checkpoints and paired outputs."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
from typing import TYPE_CHECKING, Callable, cast

from .docking_config import finite, integer, load_config, object_value
from .profiles import chain_ranges, peptide_profile, receptor_profile

SRC = Path(__file__).resolve().parents[1]
PEPTIDE_STEP = "0001-cyclic-peptide-preparation"
RECEPTOR_STEP = "0002-receptor-gprotein-preparation"
DOCKING_STEP = "0003-binding-site-selection"
Record = dict[str, object]
PREPARATION_MODULES: dict[str, set[str]] = {
    PEPTIDE_STEP: {"cleprep/states.py", "cleprep/conformers.py", "cleprep/common/peptides.py",
                   "cleprep/common/conformer_qc.py", "cleprep/common/geometry.py", "cleprep/common/io.py",
                   "cleprep/common/protein.py", "production/profiles.py", "production/cli.py"},
    RECEPTOR_STEP: {"cleprep/selection.py", "cleprep/repair.py", "cleprep/receptor_qc.py", "cleprep/tails.py",
                    "cleprep/chemistry.py", "cleprep/orientation.py", "cleprep/pdbqt.py", "cleprep/common/geometry.py",
                    "cleprep/common/io.py", "cleprep/common/protein.py", "production/receptor_models.py",
                    "production/profiles.py", "production/cli.py"},
}

if TYPE_CHECKING:
    from cleprep.workflow import WorkflowConfig


@dataclass(frozen=True)
class PeptideSettings:
    input_dir: Path
    profile: Path
    states: tuple[str, ...]
    attempts_per_state: int
    keep: int
    seed_batches: int
    seed: int
    threads: int


@dataclass(frozen=True)
class ReceptorSettings:
    input: Path
    profile: Path
    forcefield: Path
    transform: Path
    repair_seeds: tuple[int, ...]
    models: int
    loops: int
    keep: int
    ph: float
    seed: int
    cpu_threads: int


@dataclass(frozen=True)
class PipelineConfig:
    peptides: PeptideSettings
    receptor: ReceptorSettings
    docking_profile: Path


def load_pipeline_config(path: Path) -> PipelineConfig:
    path = path.expanduser().resolve()
    data = object_value(json.loads(path.read_text()), {"version", "peptides", "receptor", "docking_profile"},
                        {"version", "peptides", "receptor", "docking_profile"})
    if type(data["version"]) is not int or data["version"] != 1:
        raise ValueError("Pipeline configuration requires version=1")

    def resource(value: object, directory: bool = False) -> Path:
        if not isinstance(value, str) or not value:
            raise ValueError("Resource paths must be nonempty strings")
        result = Path(value).expanduser()
        result = (result if result.is_absolute() else path.parent / result).resolve()
        if not (result.is_dir() if directory else result.is_file()):
            raise FileNotFoundError(result)
        return result

    peptide = object_value(data["peptides"], set(PeptideSettings.__dataclass_fields__), {"input_dir", "profile"})
    preset = resource(peptide["profile"])
    known = peptide_profile(preset)["states"]
    states = peptide.get("states", list(known))
    if (not isinstance(states, list) or not states or any(not isinstance(state, str) or state not in known for state in states)
            or len(set(states)) != len(states)):
        raise ValueError("states must be unique state names from the peptide profile")
    peptides = PeptideSettings(resource(peptide["input_dir"], True), preset, tuple(states),
                               integer(peptide.get("attempts_per_state", 500), 1), integer(peptide.get("keep", 40), 1),
                               integer(peptide.get("seed_batches", 5), 1), integer(peptide.get("seed", 20260910), 0),
                               integer(peptide.get("threads", 2), 1))
    receptor = object_value(data["receptor"], set(ReceptorSettings.__dataclass_fields__),
                            {"input", "profile", "forcefield", "transform", "repair_seeds"})
    seeds = receptor["repair_seeds"]
    if not isinstance(seeds, list) or not seeds:
        raise ValueError("repair_seeds must be a nonempty list of unique negative integers")
    repair_seeds = tuple(integer(seed, -2147483647, -1) for seed in seeds)
    if len(set(repair_seeds)) != len(repair_seeds):
        raise ValueError("Duplicate receptor repair seeds")
    receptors = ReceptorSettings(resource(receptor["input"]), resource(receptor["profile"]),
                                 resource(receptor["forcefield"]), resource(receptor["transform"]), repair_seeds,
                                 integer(receptor.get("models", 1), 1), integer(receptor.get("loops", 4), 1),
                                 integer(receptor.get("keep", 3), 1), finite(receptor.get("ph", 7.4), 0, 14),
                                 integer(receptor.get("seed", 20260910), 0), integer(receptor.get("cpu_threads", 2), 1))
    return PipelineConfig(peptides, receptors, resource(data["docking_profile"]))


def check_inputs(config: PipelineConfig, limit_jobs: int = 0) -> Record:
    """Validate source chemistry, source identity, frame, force field and license."""
    from contextlib import redirect_stdout
    from io import StringIO
    import modeller
    from openmm.app import ForceField
    from cleprep.common.io import read_transform
    from cleprep.states import prepare_parent
    from cleprep.workflow import sha256

    if limit_jobs < 0:
        raise ValueError("limit_jobs must be >=0")
    peptides = peptide_profile(config.peptides.profile)
    states = [mol for name, spec in peptides["parents"].items()
              for mol in prepare_parent(config.peptides.input_dir / f"{name}.sdf",
                                        expected_residues=spec["residues"], expected_prolines=spec["prolines"],
                                        expected_sulfoxides=spec["sulfoxides"], expected_charge=spec["charge"])]
    if {mol.GetProp("state_id") for mol in states} != set(peptides["states"]):
        raise ValueError("Source peptide states differ from the peptide profile")
    for state in config.peptides.states:
        families = 2 ** peptides["parents"][peptides["states"][state]["parent"]]["prolines"]
        if config.peptides.keep % families or config.peptides.attempts_per_state < families * config.peptides.seed_batches:
            raise ValueError(f"{state}: keep must be divisible by {families}; attempts must cover every family and seed batch")
    receptor = receptor_profile(config.receptor.profile)
    if receptor.get("source_sha256") and sha256(config.receptor.input) != receptor["source_sha256"]:
        raise ValueError("Source receptor mmCIF differs from its profile")
    if (any(spec["name"] in {"CLR", "CHL", "CHL1"} for spec in receptor.get("ligands", {}).values())
            or receptor.get("chemistry", {}).get("expected_cholesterol_count", 0) != 0):
        raise ValueError("Full pipeline requires a cholesterol-free receptor profile")
    read_transform(config.receptor.transform)
    half_width = json.loads(config.receptor.transform.read_text()).get("pdbtm_membrane_half_width_angstrom")
    finite(half_width, 1e-10)
    ForceField(str(config.receptor.forcefield))
    with redirect_stdout(StringIO()):
        modeller.Environ(rand_seed=config.receptor.repair_seeds[0])
    docking = load_config(config.docking_profile, check_paths=False)
    if docking.existing_runs:
        raise ValueError("Full pipeline docking profile must not reuse existing_runs")
    if config.receptor.keep > len(docking.receptors):
        raise ValueError("Receptor keep exceeds the number of configured receptor IDs")
    if not set(config.peptides.states) <= {ligand.id for ligand in docking.ligands}:
        raise ValueError("Docking profile does not cover every requested peptide state")
    chains = tuple(receptor["chains"])
    if any(item.protein_chains != chains for item in docking.receptors):
        raise ValueError("Docking protein chains differ from the receptor construct")
    selected = {item.id for item in docking.receptors[:config.receptor.keep]}
    count = sum(len(selected.intersection(search.receptors)) for search in docking.searches)
    count *= len(config.peptides.states) * config.peptides.keep * len(docking.docking.seeds)
    if not count:
        raise ValueError("No search covers the requested receptor ensemble")
    return {"steps": [PEPTIDE_STEP, RECEPTOR_STEP, DOCKING_STEP], "states": list(config.peptides.states),
            "representatives_per_state": config.peptides.keep, "receptor_model_limit": config.receptor.keep,
            "docking_job_count_upper_bound": min(count, limit_jobs) if limit_jobs else count,
            "searches": [search.id for search in docking.searches], "cholesterol_retained": False,
            "versions": {**{name: version(name) for name in ("rdkit", "vina", "meeko", "openmm", "propka", "numpy", "scipy")},
                         "modeller": modeller.__version__}}


def run_script(step: str, name: str, arguments: list[str], working_directory: Path) -> None:
    command = [sys.executable, "-B", "-u", str(SRC / step / f"{name}.py"), *arguments]
    subprocess.run(command, cwd=working_directory, check=True)


def prepare_peptides(config: PeptideSettings, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=False)
    states, conformers = directory / "states", directory / "conformers"
    run_script(PEPTIDE_STEP, "prepare_states", ["--input-dir", str(config.input_dir), "--output-dir", str(states),
               "--profile", str(config.profile)], directory)
    run_script(PEPTIDE_STEP, "generate_conformers", ["--input-dir", str(states), "--output-dir", str(conformers),
               "--profile", str(config.profile), "--states", *config.states,
               "--attempts-per-state", str(config.attempts_per_state), "--keep", str(config.keep),
               "--seed-batches", str(config.seed_batches), "--seed", str(config.seed), "--threads", str(config.threads)], directory)


def prepare_receptors(config: ReceptorSettings, directory: Path) -> None:
    run_script(RECEPTOR_STEP, "prepare_models", ["--input", str(config.input), "--output-dir", str(directory),
               "--profile", str(config.profile), "--forcefield", str(config.forcefield), "--transform", str(config.transform),
               "--repair-seeds", *(str(seed) for seed in config.repair_seeds), "--models", str(config.models),
               "--loops", str(config.loops), "--keep", str(config.keep), "--ph", str(config.ph),
               "--seed", str(config.seed), "--cpu-threads", str(config.cpu_threads)], directory.parent)


def docking_config(config: PipelineConfig, peptides: Path, receptors: Path) -> WorkflowConfig:
    from cleprep.workflow import LigandInput, ReceptorInput, WorkflowConfig
    from .receptor_models import PreparedModel

    blueprint = load_config(config.docking_profile, check_paths=False)
    prepared = cast(list[PreparedModel], json.loads((receptors / "models.json").read_text()))
    ranges = chain_ranges(receptor_profile(config.receptor.profile))
    models = tuple(ReceptorInput(blueprint.receptors[index].id, Path(model["pdbqt"]), Path(model["pdb"]),
                   Path(model["chemistry_dir"]) / "atom_identity.json", Path(model["chemistry_dir"]) / "topology_bonds.json",
                   tuple(ranges), blueprint.receptors[index].contact_chains,
                   sum(last - first + 1 for first, last in ranges.values()))
                   for index, model in enumerate(prepared))
    available = {item.id for item in models}
    searches = tuple(replace(search, receptors=tuple(item for item in search.receptors if item in available))
                     for search in blueprint.searches if available.intersection(search.receptors))
    selected = {item.id: item for item in blueprint.ligands}
    ligands = tuple(LigandInput(state, peptides / "conformers" / state / "representatives.sdf", selected[state].records)
                    for state in config.peptides.states)
    return WorkflowConfig(models, ligands, searches, blueprint.docking, blueprint.geometry, blueprint.selection)


def run_pipeline(config: PipelineConfig, output: Path, *, resume: bool = False, limit_jobs: int = 0) -> Record:
    from cleprep.workflow import WorkflowConfig, atomic_json, portable_value, run_workflow, sha256

    checked = check_inputs(config, limit_jobs)
    resources = [config.peptides.profile, config.receptor.input, config.receptor.profile,
                 config.receptor.forcefield, config.receptor.transform, config.docking_profile]
    resources.extend(config.peptides.input_dir / f"{parent}.sdf" for parent in peptide_profile(config.peptides.profile)["parents"])
    methods = [file for folder in (SRC / "cleprep", SRC / "production", SRC / PEPTIDE_STEP, SRC / RECEPTOR_STEP, SRC / DOCKING_STEP)
               for file in sorted(folder.rglob("*.py"))]
    contract: Record = {"config": portable_value(asdict(config)), "limit_jobs": limit_jobs,
                        "input_sha256": {str(file): sha256(file) for file in resources},
                        "method_sha256": {str(file.relative_to(SRC)): sha256(file) for file in methods},
                        "versions": checked["versions"]}
    output = output.expanduser().resolve()
    record = output / "pipeline_manifest.json"
    method_digests = cast(dict[str, str], contract["method_sha256"])

    def preparation_methods(step: str, methods: dict[str, str]) -> dict[str, str]:
        return {name: digest for name, digest in methods.items()
                if name.startswith(step + "/") or name in PREPARATION_MODULES[step]}

    if resume:
        manifest = cast(Record, json.loads(record.read_text()))
        previous = cast(Record, manifest["contract"])
        if any(previous[key] != contract[key] for key in ("config", "limit_jobs", "input_sha256", "versions")):
            raise ValueError("Pipeline inputs, settings or software changed; use a new output directory")
        for step, state in cast(dict[str, Record], manifest["preparations"]).items():
            if state["status"] == "complete":
                original = state.get("method_sha256", preparation_methods(step, cast(dict[str, str], previous["method_sha256"])))
                if original != preparation_methods(step, method_digests):
                    raise ValueError(f"Completed preparation methods changed: {step}; use a new output directory")
        manifest["resume_method_sha256"] = method_digests
    else:
        output.mkdir(parents=True, exist_ok=False)
        manifest = {"version": 1, "contract": contract, "preparations": {}, "status": "running"}
    atomic_json(record, manifest)
    preparations = cast(dict[str, Record], manifest["preparations"])

    def preparation(step: str, action: Callable[[Path], None]) -> Path:
        state = preparations.setdefault(step, {"status": "pending", "attempts": []})
        if state["status"] == "complete":
            for filename, digest in cast(dict[str, str], state["output_sha256"]).items():
                if sha256(Path(filename)) != digest:
                    raise ValueError(f"Completed preparation output changed: {filename}")
            print(f"{step}: 已完成且校验通过。", flush=True)
            return Path(str(state["directory"]))
        attempts = cast(list[str], state["attempts"])
        directory = output / step / f"attempt-{len(attempts) + 1:03d}"
        directory.parent.mkdir(parents=True, exist_ok=True)
        attempts.append(str(directory))
        state.update(status="running", directory=str(directory))
        atomic_json(record, manifest)
        print(f"开始 {step}", flush=True)
        action(directory)
        files = sorted(file for file in directory.rglob("*") if file.is_file())
        state.update(status="complete", output_sha256={str(file): sha256(file) for file in files},
                     method_sha256=preparation_methods(step, method_digests))
        atomic_json(record, manifest)
        return directory

    try:
        peptides = preparation(PEPTIDE_STEP, lambda directory: prepare_peptides(config.peptides, directory))
        receptors = preparation(RECEPTOR_STEP, lambda directory: prepare_receptors(config.receptor, directory))
        workflow = docking_config(config, peptides, receptors)
        path = output / "docking_config.json"
        generated = cast(Record, portable_value(asdict(workflow)))
        generated["version"] = 1
        for row in cast(list[Record], generated["ligands"]):
            if row["records"] is None:
                row["records"] = "all"
        atomic_json(path, generated)
        docking_output = output / DOCKING_STEP
        print(f"开始 {DOCKING_STEP}", flush=True)
        summary = run_workflow(load_config(path), docking_output,
                               resume=resume and docking_output.exists(), limit_jobs=limit_jobs)
        if contract["input_sha256"] != {str(file): sha256(file) for file in resources}:
            raise ValueError("Source inputs changed during the pipeline")
        summary.update(peptide_dir=str(peptides), receptor_dir=str(receptors), docking_dir=str(docking_output))
        manifest["status"] = ("partial" if summary["failed_job_count"] else
                              "complete" if summary["exported_model_count"] else "no_candidates")
        manifest["summary"] = summary
        atomic_json(output / "summary.json", summary)
        atomic_json(record, manifest)
        return summary
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError, KeyboardInterrupt) as error:
        manifest.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed", error=str(error))
        atomic_json(record, manifest)
        raise
