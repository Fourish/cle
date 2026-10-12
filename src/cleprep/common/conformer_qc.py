"""Screen cyclic-peptide coordinates, 3D chirality and geometry in Angstrom."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from .geometry import segment_distance
from .peptides import Omega, stereo_labels


@dataclass(frozen=True)
class QCThresholds:
    """Geometry tolerances; chirality and severe clashes remain mandatory."""

    bond_ratio_min: float = 0.75
    bond_ratio_max: float = 1.25
    amide_deviation_degrees: float = 35.0
    peptide_cn_min_angstrom: float = 1.20
    peptide_cn_max_angstrom: float = 1.60


STRICT_QC = QCThresholds()
RELAXED_QC = QCThresholds(
    bond_ratio_min=0.70, bond_ratio_max=1.30, amide_deviation_degrees=45.0,
    peptide_cn_min_angstrom=1.15, peptide_cn_max_angstrom=1.75,
)
QC_PROFILES = {"strict": STRICT_QC, "relaxed": RELAXED_QC}


def expected_bond_length(bond: Chem.Bond) -> float:
    table = Chem.GetPeriodicTable()
    radii = (table.GetRcovalent(bond.GetBeginAtom().GetAtomicNum())
             + table.GetRcovalent(bond.GetEndAtom().GetAtomicNum()))
    order = bond.GetBondTypeAsDouble()
    correction = 0.92 if bond.GetIsAromatic() else 0.80 if order >= 3 else 0.88 if order >= 2 else 1.0
    return float(radii * correction)


def audit_conformer(mol: Chem.Mol, conf_id: int, expected_cip: dict[int, str],
                    backbone_indices: list[int], peptide_omegas: list[Omega],
                    *, thresholds: QCThresholds = STRICT_QC) -> list[str]:
    """Return rejection reasons, preserving original strict defaults."""
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
        if not thresholds.bond_ratio_min <= ratio <= thresholds.bond_ratio_max:
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
        if (not math.isfinite(degrees)
                or min(abs(degrees), abs(180 - abs(degrees))) > thresholds.amide_deviation_degrees):
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
