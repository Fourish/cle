"""Serial docking batches, resumable records and geometry-filtered complexes."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from importlib.metadata import version
import hashlib
import json
from pathlib import Path
from typing import Callable, cast

import numpy as np
from rdkit import Chem

from .common.io import write_json
from .common.protein import require_no_cholesterol
from .complexes import (GeometryPolicy, HeavyReceptor, PoseAudit, Record, audit_pose,
                        check_receptor_frame, export_complex, heavy_graph,
                        read_receptor, same_cluster)
from .docking import dock_ligand, parse_index_map, read_ligand, restore_heavy_pose


@dataclass(frozen=True)
class ReceptorInput:
    id: str
    pdbqt: Path
    pdb: Path
    atom_identity: Path
    topology_bonds: Path
    protein_chains: tuple[str, ...]
    contact_chains: tuple[str, ...]
    expected_residues: int | None = None


@dataclass(frozen=True)
class LigandInput:
    id: str
    sdf: Path
    records: tuple[int, ...] | None = None


@dataclass(frozen=True)
class SearchInput:
    id: str
    receptors: tuple[str, ...]
    center: tuple[float, float, float]
    size: tuple[float, float, float]


@dataclass(frozen=True)
class DockingSettings:
    seeds: tuple[int, ...] = (20260910,)
    cpu: int = 1
    exhaustiveness: int = 8
    max_evals: int = 100000
    n_poses: int = 10


@dataclass(frozen=True)
class SelectionSettings:
    backbone_rmsd: float = 3.0
    contact_jaccard: float = 0.5
    max_representatives_per_group: int = 3
    ligand_chain: str = "P"


@dataclass(frozen=True)
class WorkflowConfig:
    receptors: tuple[ReceptorInput, ...]
    ligands: tuple[LigandInput, ...]
    searches: tuple[SearchInput, ...]
    docking: DockingSettings = DockingSettings()
    geometry: GeometryPolicy = GeometryPolicy()
    selection: SelectionSettings = SelectionSettings()
    existing_runs: tuple[Path, ...] = ()


@dataclass(frozen=True)
class DockJob:
    id: str
    receptor: ReceptorInput
    ligand: LigandInput
    search: SearchInput
    record: int
    seed: int

    def parameters(self, settings: DockingSettings) -> Record:
        return {"center_A": list(self.search.center), "size_A": list(self.search.size),
                "seed": self.seed, "cpu": settings.cpu, "exhaustiveness": settings.exhaustiveness,
                "max_evals": settings.max_evals, "n_poses": settings.n_poses,
                "spacing_A": 0.375, "min_rmsd_A": 1.0, "energy_range_kcal_mol": 100.0,
                "rigid_macrocycles": True, "flexible_amides": False, "scoring_function": "vina"}


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def portable_value(value: object) -> object:
    if isinstance(value, Path):
        return str(value.resolve())
    if isinstance(value, dict):
        return {str(key): portable_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [portable_value(item) for item in value]
    return value


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".writing")
    write_json(temporary, payload)
    temporary.replace(path)


def prepare_jobs(config: WorkflowConfig, limit_jobs: int = 0) -> tuple[list[DockJob], dict[str, HeavyReceptor]]:
    """Check paired receptors and stable source atom order before any search."""
    receptors: dict[str, HeavyReceptor] = {}
    for item in config.receptors:
        require_no_cholesterol(item.pdb.read_text().splitlines(), str(item.pdb))
        require_no_cholesterol(item.pdbqt.read_text().splitlines(), str(item.pdbqt))
        receptor = read_receptor(item.pdb, item.atom_identity, item.topology_bonds,
                                  protein_chains=item.protein_chains, contact_chains=item.contact_chains,
                                  expected_residues=item.expected_residues)
        check_receptor_frame(receptor, item.pdbqt)
        if config.selection.ligand_chain in {line[21] for line in receptor.lines}:
            raise ValueError("Configured ligand chain is already present in a receptor")
        receptors[item.id] = receptor
    records: dict[str, tuple[int, ...]] = {}
    for ligand in config.ligands:
        supplier = Chem.SDMolSupplier(str(ligand.sdf), removeHs=False)
        chosen = ligand.records if ligand.records is not None else tuple(range(len(supplier)))
        if not chosen:
            raise ValueError(f"Ligand has no selected records: {ligand.id}")
        signatures: set[tuple[object, ...]] = set()
        for index in chosen:
            source = read_ligand(ligand.sdf, index)
            signatures.add(heavy_graph(Chem.RemoveAllHs(source)))
        if len(signatures) != 1:
            raise ValueError(f"Ligand records change heavy chemistry, stereochemistry or source atom order: {ligand.id}")
        records[ligand.id] = chosen
    lookup = {item.id: item for item in config.receptors}
    jobs: list[DockJob] = []
    for search in config.searches:
        for receptor_id in search.receptors:
            for ligand in config.ligands:
                for record in records[ligand.id]:
                    for seed in config.docking.seeds:
                        jobs.append(DockJob(f"J{len(jobs) + 1:06d}", lookup[receptor_id], ligand, search, record, seed))
    if limit_jobs < 0 or not jobs:
        raise ValueError("Require nonempty jobs and limit_jobs >=0")
    return jobs[:limit_jobs] if limit_jobs else jobs, receptors


def input_contract(config: WorkflowConfig, jobs: list[DockJob], limit_jobs: int) -> Record:
    paths = {item.sdf for item in config.ligands}
    paths.update(path for item in config.receptors for path in (item.pdbqt, item.pdb, item.atom_identity, item.topology_bonds))
    paths.update(run / name for run in config.existing_runs for name in ("ligand.pdbqt", "poses.pdbqt", "poses_heavy.sdf", "summary.json"))
    files = {str(path.resolve()): sha256(path) for path in sorted(paths)}
    src = Path(__file__).resolve().parent
    methods = {name: sha256(src / name) for name in ("workflow.py", "complexes.py", "docking.py",
               "common/conformer_qc.py", "common/geometry.py", "common/peptides.py", "common/protein.py", "common/topology.py")}
    return {"config": portable_value(asdict(config)), "input_sha256": files, "method_sha256": methods,
            "versions": {name: version(name) for name in ("rdkit", "vina", "meeko", "numpy", "scipy")},
            "limit_jobs": limit_jobs, "job_count": len(jobs)}


def run_digests(directory: Path) -> dict[str, str]:
    return {name: sha256(directory / name) for name in ("ligand.pdbqt", "poses.pdbqt", "poses_heavy.sdf", "summary.json")}


def matching_run(job: DockJob, directory: Path, settings: DockingSettings) -> bool:
    summary = cast(Record, json.loads((directory / "summary.json").read_text()))
    return (Path(str(summary["receptor"])).resolve() == job.receptor.pdbqt.resolve()
            and Path(str(summary["ligand"])).resolve() == job.ligand.sdf.resolve()
            and summary["source_record_index0"] == job.record and summary["parameters"] == job.parameters(settings))


def check_source_conformer(job: DockJob, directory: Path) -> None:
    """Verify the cached ligand preparation against the currently supplied SDF."""
    source = read_ligand(job.ligand.sdf, job.record)
    text = (directory / "ligand.pdbqt").read_text()
    restored = restore_heavy_pose(source, text, parse_index_map(text))
    reference = Chem.RemoveAllHs(source)
    if heavy_graph(restored) != heavy_graph(reference) or np.max(np.abs(
            restored.GetConformer().GetPositions() - reference.GetConformer().GetPositions())) > 0.00051:
        raise ValueError("Cached ligand preparation differs from the selected source conformer")


@dataclass
class Candidate:
    job: DockJob
    pose_index1: int
    score: float
    audit: PoseAudit


def screen_completed(config: WorkflowConfig, jobs: list[DockJob], states: list[Record],
                     receptors: dict[str, HeavyReceptor], output: Path) -> Record:
    """Audit all poses; cluster within one chemical/receptor/search condition."""
    output.mkdir(parents=True, exist_ok=False)
    reports: list[Record] = []
    groups: dict[tuple[str, str, str], list[Candidate]] = defaultdict(list)
    sources: dict[tuple[str, int], Chem.Mol] = {}
    for job, state in zip(jobs, states, strict=True):
        if state["status"] != "complete":
            continue
        directory = Path(str(state["run_dir"]))
        summary = cast(Record, json.loads((directory / "summary.json").read_text()))
        poses = Chem.SDMolSupplier(str(directory / "poses_heavy.sdf"), removeHs=False)
        if len(poses) != summary["pose_count"]:
            raise ValueError("Docking pose count differs from its recorded summary")
        key = job.ligand.id, job.record
        if key not in sources:
            sources[key] = read_ligand(job.ligand.sdf, job.record)
        source = sources[key]
        for index, pose in enumerate(poses, 1):
            if pose is None:
                audit = PoseAudit(None, [], [], {"passed": False, "issues": ["unreadable_SDF_pose"]})
                score = None
            else:
                score = float(pose.GetProp("vina_score_kcal_mol"))
                row = cast(list[Record], summary["poses"])[index - 1]
                if not np.isfinite(score) or score != float(row["score_kcal_mol"]):
                    raise ValueError("Pose score is nonfinite or differs from the docking summary")
                audit = audit_pose(pose, source, receptors[job.receptor.id], config.geometry)
            reports.append({"job_id": job.id, "ligand_id": job.ligand.id, "receptor_id": job.receptor.id,
                            "search_id": job.search.id, "source_record_index0": job.record, "seed": job.seed,
                            "pose_index1": index, "score_kcal_mol": score, "run_dir": str(directory), **audit.report})
            if audit.report["passed"] and score is not None:
                groups[(job.ligand.id, job.receptor.id, job.search.id)].append(Candidate(job, index, score, audit))
    write_json(output / "poses_qc.json", reports)
    cluster_rows: list[Record] = []
    models: list[Record] = []
    for key, candidates in sorted(groups.items()):
        clusters: list[list[Candidate]] = []
        for candidate in sorted(candidates, key=lambda c: (c.score, c.job.id, c.pose_index1)):
            for cluster in clusters:
                if same_cluster(candidate.audit, cluster[0].audit, rmsd_cutoff=config.selection.backbone_rmsd,
                                contact_jaccard=config.selection.contact_jaccard):
                    cluster.append(candidate)
                    break
            else:
                clusters.append([candidate])
        clusters.sort(key=lambda c: (-len({p.job.seed for p in c}), -len({p.job.record for p in c}), c[0].score))
        for rank, cluster in enumerate(clusters, 1):
            representative = cluster[0]
            cluster_id = "__".join((*key, f"C{rank:03d}"))
            row: Record = {"cluster_id": cluster_id, "ligand_id": key[0], "receptor_id": key[1], "search_id": key[2],
                           "members": [{"job_id": c.job.id, "pose_index1": c.pose_index1} for c in cluster],
                           "independent_seeds": sorted({c.job.seed for c in cluster}),
                           "source_records": sorted({c.job.record for c in cluster}), "score_kcal_mol": representative.score,
                           "selected": rank <= config.selection.max_representatives_per_group}
            cluster_rows.append(row)
            if row["selected"]:
                job = representative.job
                directory = output / "models" / cluster_id
                manifest = export_complex(directory, representative.audit, sources[(job.ligand.id, job.record)],
                                          receptors[job.receptor.id], {**row, "job_id": job.id,
                                          "source_ligand": str(job.ligand.sdf), "source_record_index0": job.record,
                                          "source_receptor": str(job.receptor.pdb), "source_run_dir": str(next(
                                              state["run_dir"] for state in states if state["job_id"] == job.id)),
                                          "pose_index1": representative.pose_index1, "seed": job.seed},
                                          ligand_chain=config.selection.ligand_chain, policy=config.geometry)
                models.append({"cluster_id": cluster_id, "directory": str(directory),
                               "score_kcal_mol": representative.score, "geometry_passed": True,
                               "complex_pdb": str(directory / "complex_heavy.pdb"),
                               "peptide_sdf": str(directory / "peptide_heavy.sdf")})
    write_json(output / "clusters.json", cluster_rows)
    summary: Record = {"pose_count": len(reports), "geometry_passed_count": sum(bool(row["passed"]) for row in reports),
                       "exported_model_count": len(models), "models": models,
                       "rejection_counts": dict(Counter(str(issue) for row in reports for issue in cast(list[str], row["issues"]))),
                       "geometry_policy": asdict(config.geometry), "selection_policy": asdict(config.selection),
                       "scope": "Geometry-filtered heavy candidates; search IDs are provenance, no binding-site assignment or cross-model validation"}
    write_json(output / "summary.json", summary)
    return summary


def run_workflow(config: WorkflowConfig, output: Path, *, resume: bool = False,
                 screen_only: bool = False, limit_jobs: int = 0,
                 progress: Callable[[str], None] = print) -> Record:
    output = output.resolve()
    if output.exists() and not resume:
        raise FileExistsError(f"Output exists; use a new directory or --resume: {output}")
    if resume and not (output / "run_manifest.json").is_file():
        raise ValueError("Resume requires this workflow's run_manifest.json")
    jobs, receptors = prepare_jobs(config, limit_jobs)
    contract = input_contract(config, jobs, limit_jobs)
    if resume:
        manifest = cast(Record, json.loads((output / "run_manifest.json").read_text()))
        if manifest["contract"] != contract:
            raise ValueError("Inputs, settings, method or software changed; resume requires an identical contract")
        states = cast(list[Record], manifest["jobs"])
        for state in states:
            if state["status"] == "complete" and state["output_sha256"] != run_digests(Path(str(state["run_dir"]))):
                raise ValueError(f"Completed docking outputs changed: {state['job_id']}")
        if manifest.get("status") == "complete":
            for name, digest in cast(dict[str, str], manifest["delivery_sha256"]).items():
                if sha256(output / name) != digest:
                    raise ValueError(f"Completed delivery changed: {name}")
            progress("已完成且校验通过，直接复用交付模型。")
            return cast(Record, manifest["summary"])
    else:
        output.mkdir(parents=True, exist_ok=False)
        states = [{"job_id": job.id, "ligand_id": job.ligand.id, "receptor_id": job.receptor.id,
                   "search_id": job.search.id, "record": job.record, "seed": job.seed,
                   "status": "pending", "attempts": []} for job in jobs]
        manifest = {"version": 1, "contract": contract, "jobs": states, "status": "running"}
    atomic_json(output / "run_manifest.json", manifest)
    for index, (job, state) in enumerate(zip(jobs, states, strict=True), 1):
        if state["status"] == "complete":
            continue
        progress(f"[{index}/{len(jobs)}] {job.ligand.id} record={job.record} {job.receptor.id}/{job.search.id} seed={job.seed}")
        attempts = cast(list[str], state["attempts"])
        directory = output / "runs" / job.id / f"attempt-{len(attempts) + 1:03d}"
        state["status"] = "running"
        atomic_json(output / "run_manifest.json", manifest)
        try:
            matches = [run for run in config.existing_runs if matching_run(job, run, config.docking)]
            if len(matches) > 1:
                raise ValueError("Multiple cached runs match one job")
            if matches:
                directory = matches[0]
                check_source_conformer(job, directory)
                state["reused_existing_run"] = True
            elif screen_only:
                raise ValueError("No matching existing run; --screen-only does not perform docking")
            else:
                attempts.append(str(directory))
                atomic_json(output / "run_manifest.json", manifest)
                dock_ligand(job.receptor.pdbqt, job.ligand.sdf, directory, list(job.search.center), list(job.search.size),
                            record=job.record, seed=job.seed, cpu=config.docking.cpu,
                            exhaustiveness=config.docking.exhaustiveness, max_evals=config.docking.max_evals,
                            n_poses=config.docking.n_poses)
                state["reused_existing_run"] = False
            state.update({"status": "complete", "run_dir": str(directory), "output_sha256": run_digests(directory)})
            state.pop("error", None)
        except (ValueError, OSError, RuntimeError) as error:
            state.update({"status": "failed", "error": f"{type(error).__name__}: {error}"})
            progress(f"任务失败，已记录：{error}")
        atomic_json(output / "run_manifest.json", manifest)
    # New attempts preserve interrupted/failed screening files rather than overwrite.
    attempts_dir = output / "screening"
    number = 1
    while (attempts_dir / f"attempt-{number:03d}").exists():
        number += 1
    delivery = attempts_dir / f"attempt-{number:03d}"
    progress("检查对接姿态、完整复合物几何并导出代表模型。")
    summary = screen_completed(config, jobs, states, receptors, delivery)
    if contract != input_contract(config, jobs, limit_jobs):
        raise ValueError("Frozen input changed during this workflow")
    failed = sum(state["status"] != "complete" for state in states)
    summary.update({"job_count": len(jobs), "failed_job_count": failed, "delivery_dir": str(delivery)})
    atomic_json(output / "summary.json", summary)
    manifest.update({"status": "partial" if failed else "complete", "summary": summary,
                     "delivery_sha256": {str(path.relative_to(output)): sha256(path) for path in sorted(delivery.rglob("*")) if path.is_file()}})
    cast(dict[str, str], manifest["delivery_sha256"])["summary.json"] = sha256(output / "summary.json")
    atomic_json(output / "run_manifest.json", manifest)
    return summary
