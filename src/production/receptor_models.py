"""Prepare and screen receptor ensembles from an experimental mmCIF."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Callable, TypedDict, cast

from cleprep.chemistry import prepare
from cleprep.common.io import read_fasta, read_transform, write_json
from cleprep.orientation import orient
from cleprep.pdbqt import export_receptor
from cleprep.receptor_qc import QCSummary, inspect_structure, read_structure
from cleprep.repair import repair
from cleprep.selection import select
from cleprep.tails import TailRegion, refine_tails
from .profiles import ReceptorProfile, chain_ranges, chain_starts


class PreparedModel(TypedDict):
    id: str
    pdb: str
    pdbqt: str
    chemistry_dir: str
    source: str
    dope_score: float | None


@dataclass(frozen=True)
class Candidate:
    source: Path
    refined: Path
    dope_score: float | None

    @property
    def score(self) -> float:
        return self.dope_score if self.dope_score is not None else float("inf")


def prepare_models(source: Path, output: Path, forcefield: Path, transform: Path,
                   profile: ReceptorProfile, *, repair_seeds: tuple[int, ...],
                   models: int = 1, loops: int = 4, keep: int = 3,
                   ph: float = 7.4, seed: int = 20260910, cpu_threads: int = 2,
                   progress: Callable[[str], None] = print) -> list[PreparedModel]:
    """Repair candidates, reject failed geometry, then parameterize and orient."""
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    selection = output / "selection"
    targets = {chain: (spec["label"], spec["first"], spec["last"])
               for chain, spec in profile["chains"].items()}
    ligands = {label: (spec["name"], spec["heavy_atoms"])
               for label, spec in profile.get("ligands", {}).items()}
    select(source, selection, targets, ligands, profile.get("source_sha256"))
    template = selection / "complex_template.pdb"
    fasta = selection / "target_sequences.fasta"
    sequences = read_fasta(fasta)
    ranges = chain_ranges(profile)
    rotation, translation = read_transform(transform)
    half_width = float(json.loads(transform.read_text())["pdbtm_membrane_half_width_angstrom"])
    regions = [TailRegion(spec["chain"], spec["first"], spec["last"], spec["critical_start"],
                          spec.get("entry_phi_count", 1),
                          tuple(spec.get("second_hinges", [(0, "psi")])),
                          tuple(spec.get("third_hinges", [])))
               for spec in profile.get("tails", [])]
    candidates: list[Candidate] = []
    screening: list[dict[str, object]] = []
    qc = output / "qc"
    qc.mkdir()
    for batch, repair_seed in enumerate(repair_seeds, 1):
        progress(f"受体补建 [{batch}/{len(repair_seeds)}] seed={repair_seed}")
        directory = output / "repair" / f"batch-{batch:03d}"
        results = repair(template, fasta, directory, chain_starts(profile), repair_seed, models, loops)
        for result in results:
            row: dict[str, object] = {"batch": batch, "seed": repair_seed, **result, "accepted": False}
            screening.append(row)
            if result["failed"] is not None:
                continue
            repaired = directory / (Path(result["filename"]).stem + "_auth.pdb")
            refined = qc / f"batch-{batch:03d}_{repaired.stem}.pdb"
            try:
                if regions:
                    report = refine_tails(repaired, refined, template, sequences, ranges,
                                          regions, rotation, translation, half_width)
                    summary = cast(QCSummary, report["summary"])
                else:
                    refined.write_bytes(repaired.read_bytes())
                    qc_report = inspect_structure(read_structure(refined), read_structure(template), sequences, ranges)
                    write_json(refined.with_suffix(".json"), qc_report)
                    summary = qc_report["summary"]
                if summary["has_blocking_geometry_issues"]:
                    row["rejection"] = "blocking_geometry_issues"
                    continue
            except ValueError as error:
                row["rejection"] = str(error)
                continue
            row["accepted"] = True
            row["refined"] = str(refined)
            candidates.append(Candidate(repaired, refined, result["dope_score"]))
    write_json(output / "candidate_screening.json", screening)
    if not candidates:
        raise RuntimeError("No repaired receptor passed geometry screening; see candidate_screening.json and qc/")
    chosen = sorted(candidates, key=lambda item: (item.score, str(item.source)))[:keep]
    progress(f"受体几何通过 {len(candidates)} 个候选，准备 {len(chosen)} 个模型。")
    chemistry = profile.get("chemistry", {})
    prepared: list[PreparedModel] = []
    for number, candidate in enumerate(chosen, 1):
        model_id = f"M{number:03d}"
        directory = output / "models" / model_id
        chemistry_dir = directory / "chemistry"
        progress(f"受体化学处理与定向 [{number}/{len(chosen)}] {model_id}")
        prepare(candidate.refined, chemistry_dir, forcefield, ranges,
                context_residues=frozenset(chemistry.get("context_residues", [])),
                expected_context_atoms=chemistry.get("expected_context_atoms"),
                expected_cholesterol_count=0, ph=ph, seed=seed, cpu_threads=cpu_threads,
                progress=progress)
        pdb = directory / "receptor.pdb"
        pdbqt = directory / "receptor.pdbqt"
        orient(chemistry_dir / "protein_cholesterol.pdb", pdb, rotation, translation)
        export_receptor(pdb, chemistry_dir, pdbqt, keep_cholesterol=False,
                        chain_starts=chain_starts(profile),
                        expected_heavy_atoms=profile.get("expected_protein_heavy_atoms"))
        prepared.append({"id": model_id, "pdb": str(pdb), "pdbqt": str(pdbqt),
                         "chemistry_dir": str(chemistry_dir), "source": str(candidate.source),
                         "dope_score": candidate.dope_score})
    write_json(output / "models.json", prepared)
    return prepared
