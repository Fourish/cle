"""Generate, optimize and filter 3D cyclic-peptide conformers from explicit-H states."""

from __future__ import annotations

from collections import Counter
from itertools import product
from typing import Callable

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem, DistanceGeometry
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolTransforms

from .common.conformer_qc import audit_conformer
from .common.geometry import rigid_fit
from .common.peptides import Omega, peptide_backbone, stereo_labels, sulfoxide_atoms


def make_bounds(mol: Chem.Mol, amides: list[Omega], pro_amides: list[Omega],
                family: str) -> NDArray[np.float64]:
    """Guide cis/trans sampling with CA-CA bounds; keep non-Pro amides trans."""
    bounds = rdDistGeom.GetMoleculeBoundsMatrix(
        mol, set15bounds=False, doTriangleSmoothing=False,
        useMacrocycle14config=False, forceTransAmides=False,
    )
    cis_amides = {atoms for atoms, tag in zip(pro_amides, family) if tag == "C"}
    for atoms in amides:
        low, high = sorted((atoms[0], atoms[3]))
        bounds[high, low], bounds[low, high] = (
            (2.65, 3.15) if atoms in cis_amides else (3.55, 4.10)
        )
    if not DistanceGeometry.DoTriangleSmoothing(bounds):
        raise ValueError(f"Inconsistent distance bounds for {family}")
    return bounds


def aligned_rmsds(first: Chem.Mol, second: Chem.Mol, backbone: list[int]) -> tuple[float, float]:
    """Fit fixed atom pairs by translation and proper rotation, never reflection."""
    a = np.asarray(first.GetConformer().GetPositions(), dtype=np.float64)
    b = np.asarray(second.GetConformer().GetPositions(), dtype=np.float64)
    rotation, translation = rigid_fit(a[backbone], b[backbone])
    delta = a @ rotation + translation - b
    heavy = [atom.GetIdx() for atom in first.GetAtoms() if atom.GetAtomicNum() > 1]
    return (float(np.sqrt(np.mean(np.sum(delta[backbone] ** 2, axis=1)))),
            float(np.sqrt(np.mean(np.sum(delta[heavy] ** 2, axis=1)))))


def select_representatives(candidates: list[Chem.Mol], backbone: list[int],
                           families: list[str], keep: int) -> list[Chem.Mol]:
    """Select energy-ordered, geometrically distinct conformers within each family."""
    selected: list[Chem.Mol] = []
    quota = keep // len(families)
    ordered = sorted(candidates, key=lambda mol: (mol.GetDoubleProp("MMFF94s_energy_kcal_mol"),
                                                  mol.GetProp("_Name")))
    for family in families:
        group: list[Chem.Mol] = []
        for mol in ordered:
            if mol.GetProp("observed_x_pro_family") != family:
                continue
            if all(bb >= 1.0 or heavy >= 1.75
                   for bb, heavy in (aligned_rmsds(mol, other, backbone) for other in group)):
                group.append(mol)
                if len(group) == quota:
                    break
        selected.extend(group)
    return sorted(selected, key=lambda mol: (mol.GetDoubleProp("MMFF94s_energy_kcal_mol"),
                                             mol.GetProp("_Name")))


def generate(molecule: Chem.Mol, state_id: str, attempts: int, keep: int,
             seed: int, seed_batches: int, threads: int,
             progress: Callable[[str], None] | None = None) -> tuple[list[Chem.Mol], list[Chem.Mol], dict[str, object]]:
    validate_state(molecule)
    backbone, amides, pro_amides = peptide_backbone(molecule)
    family_count = 2 ** len(pro_amides)
    if (min(attempts, keep, seed_batches, threads) < 1 or seed < 0
            or keep % family_count or attempts < family_count * seed_batches
            or seed + (family_count - 1) * 1000 + seed_batches - 1 > 2147483647):
        raise ValueError("Invalid counts, family quota or RDKit random seed")
    expected_cip = stereo_labels(molecule)
    families = ["".join(tags) or "noPro" for tags in product("TC", repeat=len(pro_amides))]
    base, remainder = divmod(attempts, len(families) * seed_batches)
    properties = rdForceFieldHelpers.MMFFGetMoleculeProperties(molecule, mmffVariant="MMFF94s")
    accepted: list[Chem.Mol] = []
    rejected: Counter[str] = Counter()
    for family_index, family in enumerate(families):
        for batch_index in range(seed_batches):
            ordinal = family_index * seed_batches + batch_index
            requested = base + int(ordinal < remainder)
            batch_seed = seed + family_index * 1000 + batch_index
            mol = Chem.Mol(molecule)
            mol.RemoveAllConformers()
            params = rdDistGeom.ETKDGv3()
            params.randomSeed = batch_seed
            params.numThreads = threads
            params.maxIterations = 200
            params.timeout = 0
            params.pruneRmsThresh = -1.0
            params.enforceChirality = True
            params.forceTransAmides = False
            params.useMacrocycleTorsions = True
            params.useMacrocycle14config = False
            params.SetBoundsMat(make_bounds(mol, amides, pro_amides, family))
            returned_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=requested, params=params)
            actual_ids = {conf.GetId() for conf in mol.GetConformers()}
            conf_ids = [i for i in returned_ids if i >= 0 and i in actual_ids]
            batch_accepted = 0
            for conf_id in conf_ids:
                constrained = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, properties, confId=conf_id)
                cis_amides = {atoms for atoms, tag in zip(pro_amides, family) if tag == "C"}
                for atoms in amides:
                    low, high = (-15.0, 15.0) if atoms in cis_amides else (165.0, 195.0)
                    constrained.MMFFAddTorsionConstraint(*atoms, False, low, high, 10.0)
                constrained.Minimize(maxIts=1000)
                free = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, properties, confId=conf_id)
                status = free.Minimize(maxIts=2000)
                if status != 0:
                    status = free.Minimize(maxIts=2000)
                energy = float(free.CalcEnergy())
                issues = audit_conformer(mol, conf_id, expected_cip, backbone, amides)
                if status != 0:
                    issues.append("unconstrained_MMFF_not_converged")
                if not np.isfinite(energy):
                    issues.append("nonfinite_energy")
                if any(abs(rdMolTransforms.GetDihedralDeg(mol.GetConformer(conf_id), *a)) < 90
                       for a in amides if a not in pro_amides):
                    issues.append("nonpro_cis_outside_sampling_model")
                if issues:
                    rejected.update(issue.split(":", 1)[0] for issue in issues)
                    continue
                observed = "".join("C" if abs(rdMolTransforms.GetDihedralDeg(
                    mol.GetConformer(conf_id), *a)) < 90 else "T" for a in pro_amides) or "noPro"
                single = Chem.Mol(mol)
                single.RemoveAllConformers()
                single.AddConformer(Chem.Conformer(mol.GetConformer(conf_id)), assignId=True)
                single.GetConformer().Set3D(True)
                single.SetProp("_Name", f"{state_id}_{family}_{batch_seed}_{conf_id:03d}")
                single.SetProp("coordinate_status", "3D; unconstrained MMFF94s optimized; geometry checked")
                single.SetProp("target_x_pro_family", family)
                single.SetProp("observed_x_pro_family", observed)
                single.SetIntProp("random_seed", batch_seed)
                single.SetDoubleProp("MMFF94s_energy_kcal_mol", energy)
                single.SetProp("geometry_qc", "PASS")
                accepted.append(single)
                batch_accepted += 1
            if progress is not None:
                progress(f"{state_id} {family} seed={batch_seed}: "
                         f"embedded={len(conf_ids)}/{requested}, accepted={batch_accepted}")
    representatives = select_representatives(accepted, backbone, families, keep)
    covered = Counter(mol.GetProp("observed_x_pro_family") for mol in representatives)
    summary: dict[str, object] = {
        "state_id": state_id, "requested": attempts, "accepted": len(accepted),
        "representatives_requested": keep, "representatives": len(representatives),
        "accepted_by_family": dict(Counter(m.GetProp("observed_x_pro_family") for m in accepted)),
        "representatives_by_family": dict(covered),
        "missing_families": [family for family in families if not covered[family]],
        "rejection_reasons": dict(rejected),
    }
    return accepted, representatives, summary


def validate_state(mol: Chem.Mol, *, expected_residues: int | None = None,
                   expected_prolines: int | None = None, expected_charge: int | None = None,
                   expected_sulfur_cip: str | None = None) -> None:
    if len(Chem.GetMolFrags(mol)) != 1 or (expected_charge is not None and Chem.GetFormalCharge(mol) != expected_charge):
        raise ValueError("State connectivity or formal charge mismatch")
    if any(a.GetNumImplicitHs() or a.GetNumExplicitHs() for a in mol.GetAtoms() if a.GetAtomicNum() != 1):
        raise ValueError("All hydrogens must be explicit atoms")
    cip = stereo_labels(mol)
    if any(info.specified == Chem.StereoSpecified.Unspecified for info in Chem.FindPotentialStereo(mol)):
        raise ValueError("Unresolved stereochemistry")
    sulfurs = sulfoxide_atoms(mol)
    actual = "".join(cip.get(i, "?") for i in sulfurs) or "none"
    if expected_sulfur_cip is not None and actual != expected_sulfur_cip:
        raise ValueError("Sulfur configuration does not match the selected state")
    _, amides, pro_amides = peptide_backbone(mol)
    if ((expected_residues is not None and len(amides) != expected_residues)
            or (expected_prolines is not None and len(pro_amides) != expected_prolines)):
        raise ValueError("Peptide-unit or Pro count mismatch")
    if not rdForceFieldHelpers.MMFFHasAllMoleculeParams(mol):
        raise ValueError("Incomplete MMFF94s parameter coverage")
