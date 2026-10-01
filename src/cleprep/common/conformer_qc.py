"""Screen cyclic-peptide coordinates, 3D chirality and geometry in Angstrom."""

from __future__ import annotations

from itertools import combinations
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from .geometry import segment_distance
from .peptides import Omega, stereo_labels


def expected_bond_length(bond: Chem.Bond) -> float:
    table = Chem.GetPeriodicTable()
    radii = (table.GetRcovalent(bond.GetBeginAtom().GetAtomicNum())
             + table.GetRcovalent(bond.GetEndAtom().GetAtomicNum()))
    order = bond.GetBondTypeAsDouble()
    correction = 0.92 if bond.GetIsAromatic() else 0.80 if order >= 3 else 0.88 if order >= 2 else 1.0
    return float(radii * correction)


def audit_conformer(mol: Chem.Mol, conf_id: int, expected_cip: dict[int, str],
                    backbone_indices: list[int], peptide_omegas: list[Omega]) -> list[str]:
    """Return rejection reasons; thresholds are the original geometric screens."""
    indices = set(expected_cip) | set(backbone_indices) | {i for omega in peptide_omegas for i in omega}
    if any(i < 0 or i >= mol.GetNumAtoms() for i in indices):
        raise ValueError("QC atom indices must be within the molecule")
    conformer = mol.GetConformer(conf_id)
    xyz = np.asarray(conformer.GetPositions(), dtype=np.float64)
    if not np.isfinite(xyz).all():
        return ["nonfinite_coordinates"]
    issues: list[str] = []
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if len(heavy) >= 4:
        singular = np.linalg.svd(xyz[heavy] - xyz[heavy].mean(axis=0), compute_uv=False)
        if singular[-1] / math.sqrt(len(heavy)) <= 0.03:
            issues.append("nearly_planar_coordinates")

    measured = Chem.Mol(mol)
    Chem.RemoveStereochemistry(measured)
    measured.GetConformer(conf_id).Set3D(True)
    Chem.AssignStereochemistryFrom3D(measured, confId=conf_id, replaceExistingTags=True)
    actual_cip = stereo_labels(measured)
    for i, expected in expected_cip.items():
        if actual_cip.get(i) != expected:
            issues.append(f"3D_CIP_mismatch:{i} expected={expected} actual={actual_cip.get(i)}")

    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ratio = float(np.linalg.norm(xyz[a] - xyz[b])) / expected_bond_length(bond)
        if not 0.75 <= ratio <= 1.25:
            issues.append(f"abnormal_bond:{a}-{b} ratio={ratio:.3f}")

    topological = Chem.GetDistanceMatrix(mol)
    table = Chem.GetPeriodicTable()
    for a, b in combinations(heavy, 2):
        if topological[a, b] <= 2:
            continue
        radii = table.GetRvdw(mol.GetAtomWithIdx(a).GetAtomicNum()) + table.GetRvdw(mol.GetAtomWithIdx(b).GetAtomicNum())
        if np.linalg.norm(xyz[a] - xyz[b]) / radii < 0.55:
            issues.append(f"nonbonded_overlap:{a}-{b}")

    for omega in peptide_omegas:
        degrees = float(rdMolTransforms.GetDihedralDeg(conformer, *omega))
        if not math.isfinite(degrees) or min(abs(degrees), abs(180 - abs(degrees))) > 35.0:
            issues.append(f"twisted_amide:{omega}")

    backbone = set(backbone_indices)
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
             if b.GetBeginAtomIdx() in backbone and b.GetEndAtomIdx() in backbone]
    for (a, b), (c, d) in combinations(bonds, 2):
        if min(topological[i, j] for i in (a, b) for j in (c, d)) <= 2:
            continue
        separation, s, t = segment_distance(xyz[a], xyz[b], xyz[c], xyz[d])
        if separation < 0.55 and 0.05 < s < 0.95 and 0.05 < t < 0.95:
            issues.append(f"backbone_near_intersection:{a}-{b}/{c}-{d}")
    return issues
