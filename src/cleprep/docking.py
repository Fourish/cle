"""Single-conformer Vina docking with a fixed macrocycle and mapped SDF export."""

from __future__ import annotations

from importlib.metadata import version
import json
from pathlib import Path
import time

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from vina import Vina

from .common.conformer_qc import audit_conformer
from .common.io import Json, write_json
from .common.peptides import peptide_backbone, stereo_labels, write_sdf
from .common.protein import require_no_cholesterol


def read_ligand(path: Path, record: int = 0) -> Chem.Mol:
    """Read one zero-based SDF record; require explicit H and 3D coordinates."""
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, strictParsing=True)
    if not 0 <= record < len(supplier):
        raise ValueError(f"SDF record {record} is outside 0..{len(supplier) - 1}")
    mol: Chem.Mol | None = supplier[record]
    if mol is None or mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        raise ValueError("Ligand must be a readable 3D SDF record")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Ligand must be one connected molecule")
    if any(atom.GetNumImplicitHs() or atom.GetNumExplicitHs() for atom in mol.GetAtoms()):
        raise ValueError("Ligand requires hydrogen atoms with explicit coordinates")
    return mol


def parse_index_map(pdbqt: str) -> dict[int, int]:
    """Map PDBQT serial1 to original SDF index0 using Meeko INDEX MAP pairs."""
    mapping: dict[int, int] = {}
    for line in pdbqt.splitlines():
        if line.startswith("REMARK INDEX MAP "):
            tokens = [int(token) for token in line.split()[3:]]
            if len(tokens) % 2:
                raise ValueError("Incomplete Meeko INDEX MAP pair")
            for source1, serial1 in zip(tokens[::2], tokens[1::2], strict=True):
                if source1 < 1 or serial1 < 1 or serial1 in mapping:
                    raise ValueError("Invalid or duplicate Meeko INDEX MAP entry")
                mapping[serial1] = source1 - 1
    if not mapping or len(set(mapping.values())) != len(mapping):
        raise ValueError("Missing or nonunique Meeko INDEX MAP")
    return mapping


def restore_heavy_pose(source: Chem.Mol, block: str, mapping: dict[int, int]) -> Chem.Mol:
    """Copy the source graph, remove all H and replace mapped heavy coordinates."""
    coordinates: dict[int, list[float]] = {}
    for line in block.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            source_index = mapping[int(line[6:11])]
            if source_index in coordinates or not 0 <= source_index < source.GetNumAtoms():
                raise ValueError("Invalid or duplicate output atom")
            coordinates[source_index] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    mol = Chem.Mol(source)
    for atom in mol.GetAtoms():
        atom.SetIntProp("source_atom_index0", atom.GetIdx())
    mol = Chem.RemoveAllHs(mol)
    indices = [atom.GetIntProp("source_atom_index0") for atom in mol.GetAtoms()]
    if not set(indices).issubset(coordinates):
        raise ValueError("Docking output is missing source heavy atoms")
    for index, source_index in enumerate(indices):
        mol.GetConformer().SetAtomPosition(index, coordinates[source_index])
    # Source metadata may refer to removed H or contain obsolete conformer energies.
    for prop in list(mol.GetPropNames(includePrivate=True, includeComputed=False)):
        mol.ClearProp(prop)
    mol.SetProp("pose_atom_to_source_index0", json.dumps(indices))
    return mol


def ligand_geometry_issues(mol: Chem.Mol, reference: Chem.Mol) -> list[str]:
    """Reuse peptide QC and preserve every input amide's cis/trans assignment."""
    backbone, amides, _ = peptide_backbone(reference)
    issues = audit_conformer(mol, 0, stereo_labels(Chem.Mol(reference)), backbone, amides)
    for atoms in amides:
        cis = abs(rdMolTransforms.GetDihedralDeg(mol.GetConformer(), *atoms)) < 90
        reference_cis = abs(rdMolTransforms.GetDihedralDeg(reference.GetConformer(), *atoms)) < 90
        if cis != reference_cis:
            issues.append(f"amide_cis_trans_changed:{atoms}")
    return issues


def prepare_ligand(source: Chem.Mol) -> str:
    """Keep ring bonds and amides fixed, retain normal sidechain torsions."""
    heavy = Chem.RemoveAllHs(source)
    issues = ligand_geometry_issues(heavy, heavy)
    if issues:
        raise ValueError(f"Input ligand geometry failed: {issues}")
    preparation = MoleculePreparation(rigid_macrocycles=True, flexible_amides=False, add_index_map=True)
    setups = preparation.prepare(source)
    if len(setups) != 1:
        raise ValueError("Expected one fixed-macrocycle Meeko preparation")
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
    if not ok:
        raise ValueError(str(error))
    mapping = parse_index_map(pdbqt)
    for line in pdbqt.splitlines():
        if line.startswith("BRANCH "):
            a, b = [mapping[int(token)] for token in line.split()[1:3]]
            bond = source.GetBondBetweenAtoms(a, b)
            if bond is None or bond.IsInRing():
                raise ValueError("Meeko left a ring bond flexible")
    restored = restore_heavy_pose(source, pdbqt, mapping)
    error_A = np.max(np.abs(restored.GetConformer().GetPositions() - heavy.GetConformer().GetPositions()))
    if error_A > 0.00051 or ligand_geometry_issues(restored, heavy):
        raise ValueError("Ligand PDBQT roundtrip changed coordinates or chemistry")
    return pdbqt


def pose_blocks(pdbqt: str) -> list[str]:
    """Split Vina MODEL records, also accepting a single unwrapped pose."""
    if not any(line.startswith("MODEL ") for line in pdbqt.splitlines()):
        return [pdbqt] if any(line.startswith(("ATOM  ", "HETATM")) for line in pdbqt.splitlines()) else []
    blocks: list[str] = []
    current: list[str] = []
    for line in pdbqt.splitlines():
        if line.startswith("MODEL "):
            current = [line]
        elif current:
            current.append(line)
            if line.startswith("ENDMDL"):
                blocks.append("\n".join(current) + "\n")
                current = []
    return blocks


def dock_ligand(receptor: Path, ligand: Path, output_dir: Path,
                center: list[float], size: list[float], *, record: int = 0,
                seed: int = 20260910, cpu: int = 1, exhaustiveness: int = 8,
                max_evals: int = 100000, n_poses: int = 10) -> dict[str, Json]:
    """Prepare, dock and export one conformer. All coordinates/box lengths are A."""
    if output_dir.exists():
        raise FileExistsError(f"Output directory exists: {output_dir}")
    if not receptor.is_file():
        raise FileNotFoundError(receptor)
    require_no_cholesterol(receptor.read_text().splitlines(), str(receptor))
    if (len(center) != 3 or len(size) != 3 or not np.isfinite(center + size).all()
            or min(size) <= 0):
        raise ValueError("Box requires three finite center values and three positive full lengths")
    if not 1 <= seed <= 2147483647 or min(cpu, exhaustiveness, n_poses) < 1 or max_evals < 0:
        raise ValueError("Require seed in 1..2147483647, positive counts and max_evals >= 0")
    started = time.perf_counter()
    source = read_ligand(ligand, record)
    reference = Chem.RemoveAllHs(source)
    pdbqt = prepare_ligand(source)
    mapping = parse_index_map(pdbqt)
    engine = Vina(sf_name="vina", cpu=cpu, seed=seed, verbosity=1)
    engine.set_receptor(str(receptor))
    engine.set_ligand_from_string(pdbqt)
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "ligand.pdbqt").write_text(pdbqt, encoding="ascii")
    engine.compute_vina_maps(center=center, box_size=size, spacing=0.375)
    engine.dock(exhaustiveness=exhaustiveness, n_poses=n_poses, min_rmsd=1.0, max_evals=max_evals)
    raw = engine.poses(n_poses=n_poses, energy_range=100.0)
    (output_dir / "poses.pdbqt").write_text(raw, encoding="ascii")
    blocks = pose_blocks(raw)
    energies = np.asarray(engine.energies(n_poses=n_poses, energy_range=100.0), dtype=np.float64)
    if not blocks or energies.shape != (len(blocks), 5) or not np.isfinite(energies).all():
        raise ValueError("Empty poses, nonfinite scores or pose/score count mismatch")
    poses: list[Chem.Mol] = []
    rows: list[Json] = []
    source_name = source.GetProp("_Name") if source.HasProp("_Name") else ligand.stem
    for index, block in enumerate(blocks, start=1):
        mol = restore_heavy_pose(source, block, mapping)
        issues = ligand_geometry_issues(mol, reference)
        score = float(energies[index - 1, 0])
        mol.SetProp("_Name", f"{source_name}_pose_{index:02d}")
        mol.SetProp("source_conformer", source_name)
        if source.HasProp("state_id"):
            mol.SetProp("state_id", source.GetProp("state_id"))
        mol.SetIntProp("pose_index1", index)
        mol.SetDoubleProp("vina_score_kcal_mol", score)
        mol.SetBoolProp("ligand_geometry_passed", not issues)
        mol.SetProp("ligand_geometry_issues", json.dumps(issues))
        poses.append(mol)
        rows.append({"pose_index1": index, "score_kcal_mol": score,
                     "ligand_geometry_passed": not issues, "issues": list(issues)})
    write_sdf(output_dir / "poses_heavy.sdf", poses)
    summary: dict[str, Json] = {
        "receptor": str(receptor.resolve()), "ligand": str(ligand.resolve()),
        "source_record_index0": record, "source_conformer": source_name,
        "cholesterol_retained": False,
        "parameters": {"center_A": list(center), "size_A": list(size), "seed": seed, "cpu": cpu,
                       "exhaustiveness": exhaustiveness, "max_evals": max_evals, "n_poses": n_poses,
                       "spacing_A": 0.375, "min_rmsd_A": 1.0, "energy_range_kcal_mol": 100.0,
                       "rigid_macrocycles": True, "flexible_amides": False, "scoring_function": "vina"},
        "versions": {name: version(name) for name in ("vina", "meeko", "rdkit", "numpy")},
        "pose_count": len(poses), "ligand_geometry_passed_count": sum(pose.GetBoolProp("ligand_geometry_passed") for pose in poses),
        "elapsed_seconds": time.perf_counter() - started, "poses": rows,
    }
    write_json(output_dir / "summary.json", summary)
    return summary
