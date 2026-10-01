"""Screen standard-protein repair candidates; report only geometry failures.

Distances are in Angstrom. Nonprotein bonds are inferred from the fixed template.
Hydrogens are excluded. These geometric screens do not establish MD convergence.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Mapping, TypeAlias, TypedDict

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .common.protein import AA_ONE
from .common.geometry import Array, dihedral, segment_distance
from .common.io import Json

Vector: TypeAlias = Array
ResidueKey: TypeAlias = tuple[str, int, str]
AtomKey: TypeAlias = tuple[str, int, str, str, str]
ChiralCenter: TypeAlias = tuple[str, tuple[str, str, str]]
Ring: TypeAlias = tuple[str, ...]


class QCSummary(TypedDict):
    has_blocking_geometry_issues: bool
    protein_residues: int
    new_residues: int
    heavy_atoms: int


class QCReport(TypedDict):
    summary: QCSummary
    issues: dict[str, list[str] | list[dict[str, Json]]]

SIDE_BONDS: dict[str, str] = {
    "ALA": "CA-CB", "ARG": "CA-CB CB-CG CG-CD CD-NE NE-CZ CZ-NH1 CZ-NH2",
    "ASN": "CA-CB CB-CG CG-OD1 CG-ND2", "ASP": "CA-CB CB-CG CG-OD1 CG-OD2",
    "CYS": "CA-CB CB-SG", "GLN": "CA-CB CB-CG CG-CD CD-OE1 CD-NE2",
    "GLU": "CA-CB CB-CG CG-CD CD-OE1 CD-OE2", "GLY": "",
    "HIS": "CA-CB CB-CG CG-ND1 ND1-CE1 CE1-NE2 NE2-CD2 CD2-CG",
    "ILE": "CA-CB CB-CG1 CB-CG2 CG1-CD1", "LEU": "CA-CB CB-CG CG-CD1 CG-CD2",
    "LYS": "CA-CB CB-CG CG-CD CD-CE CE-NZ", "MET": "CA-CB CB-CG CG-SD SD-CE",
    "PHE": "CA-CB CB-CG CG-CD1 CG-CD2 CD1-CE1 CD2-CE2 CE1-CZ CE2-CZ",
    "PRO": "CA-CB CB-CG CG-CD CD-N", "SER": "CA-CB CB-OG",
    "THR": "CA-CB CB-OG1 CB-CG2",
    "TRP": "CA-CB CB-CG CG-CD1 CG-CD2 CD1-NE1 NE1-CE2 CE2-CD2 CD2-CE3 CE3-CZ3 CZ3-CH2 CH2-CZ2 CZ2-CE2",
    "TYR": "CA-CB CB-CG CG-CD1 CG-CD2 CD1-CE1 CD2-CE2 CE1-CZ CE2-CZ CZ-OH",
    "VAL": "CA-CB CB-CG1 CB-CG2",
}
STANDARD_BONDS: dict[str, tuple[tuple[str, str], ...]] = {
    name: tuple(tuple(item.split("-")) for item in ("N-CA CA-C C-O " + side).split())
    for name, side in SIDE_BONDS.items()
}
for _histidine in ("HSD", "HSE", "HSP"):
    STANDARD_BONDS[_histidine] = STANDARD_BONDS["HIS"]
EXPECTED_ATOMS: dict[str, frozenset[str]] = {
    name: frozenset(atom for pair in bonds for atom in pair)
    for name, bonds in STANDARD_BONDS.items()
}
COVALENT_RADII: dict[str, float] = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "F": 0.57, "CL": 1.02}
AROMATIC_RINGS: dict[str, tuple[Ring, ...]] = {
    "PHE": (("CG", "CD1", "CE1", "CZ", "CE2", "CD2"),),
    "TYR": (("CG", "CD1", "CE1", "CZ", "CE2", "CD2"),),
    "HIS": (("CG", "ND1", "CE1", "NE2", "CD2"),),
    "TRP": (("CG", "CD1", "NE1", "CE2", "CD2"),
            ("CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3")),
}
PLANAR_SIDECHAIN_GROUPS: dict[str, tuple[str, ...]] = {
    "ASN": ("CG", "CB", "OD1", "ND2"), "GLN": ("CD", "CG", "OE1", "NE2"),
    "ASP": ("CG", "CB", "OD1", "OD2"), "GLU": ("CD", "CG", "OE1", "OE2"),
    "ARG": ("CZ", "NE", "NH1", "NH2"),
}


@dataclass(frozen=True)
class AtomRecord:
    index: int
    chain: str
    number: int
    insertion: str
    residue: str
    name: str
    element: str
    xyz: Vector

    @property
    def residue_key(self) -> ResidueKey:
        return self.chain, self.number, self.insertion

    @property
    def key(self) -> AtomKey:
        return self.chain, self.number, self.insertion, self.residue, self.name

    @property
    def protein(self) -> bool:
        return self.residue in AA_ONE

    def label(self) -> str:
        return f"{self.chain}:{self.number}{self.insertion}:{self.residue}:{self.name}"


@dataclass
class StructureData:
    atoms: list[AtomRecord]
    residues: dict[ResidueKey, dict[str, AtomRecord]]
    read_issues: list[dict[str, Json]]


def residue_label(key: ResidueKey) -> str:
    return f"{key[0]}:{key[1]}{key[2]}"


def read_structure(path: Path) -> StructureData:
    """Read one model and choose highest occupancy for an alternate atom name."""
    structure: gemmi.Structure = gemmi.read_structure(str(path))
    issues: list[dict[str, Json]] = []
    if len(structure) != 1:
        issues.append({"kind": "model_count", "count": len(structure), "used": 1})
    if not structure:
        raise ValueError(f"No coordinate model in {path}")
    atoms: list[AtomRecord] = []
    residues: dict[ResidueKey, dict[str, AtomRecord]] = {}
    for chain in structure[0]:
        for residue in chain:
            insertion: str = residue.seqid.icode.strip()
            key: ResidueKey = (chain.name, residue.seqid.num, insertion)
            if key in residues:
                issues.append({"kind": "duplicate_residue_identity", "residue": residue_label(key)})
                continue
            selected: dict[str, gemmi.Atom] = {}
            counts: dict[str, int] = {}
            for atom in residue:
                if atom.element.is_hydrogen:
                    continue
                name: str = atom.name.strip()
                counts[name] = counts.get(name, 0) + 1
                previous: gemmi.Atom | None = selected.get(name)
                if previous is None or atom.occ > previous.occ:
                    selected[name] = atom
            atom_map: dict[str, AtomRecord] = {}
            for name, atom in selected.items():
                record = AtomRecord(len(atoms), chain.name, residue.seqid.num, insertion,
                                    residue.name, name, atom.element.name.upper(),
                                    np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=np.float64))
                if not np.isfinite(record.xyz).all():
                    raise ValueError(f"Nonfinite coordinates: {record.label()}")
                atoms.append(record)
                atom_map[name] = record
                if counts[name] > 1 or atom.altloc not in ("\x00", " "):
                    issues.append({"kind": "alternate_or_duplicate_atom", "atom": record.label(),
                                   "number_of_records": counts[name], "selected_altloc": atom.altloc,
                                   "selection": "highest occupancy, first record breaks ties"})
            if atom_map:
                residues[key] = atom_map
    return StructureData(atoms, residues, issues)


def distance(a: AtomRecord, b: AtomRecord) -> float:
    return float(np.linalg.norm(a.xyz - b.xyz))


def graph_neighbors_within_three(graph: list[set[int]], origin: int) -> set[int]:
    visited: set[int] = {origin}
    frontier: set[int] = {origin}
    for _ in range(3):
        frontier = set().union(*(graph[node] for node in frontier)) - visited if frontier else set()
        visited.update(frontier)
    return visited


def chiral_centers(residue_name: str) -> tuple[ChiralCenter, ...]:
    """Named atom ordering defines geometry without CIP R/S conventions."""
    if residue_name not in AA_ONE or residue_name == "GLY":
        return ()
    alpha: ChiralCenter = ("CA", ("N", "C", "CB"))
    if residue_name == "THR":
        return alpha, ("CB", ("CA", "OG1", "CG2"))
    if residue_name == "ILE":
        return alpha, ("CB", ("CA", "CG1", "CG2"))
    return (alpha,)


def signed_triple_product(atoms: Mapping[str, AtomRecord], center: str,
                          neighbors: tuple[str, str, str]) -> float | None:
    if center not in atoms or any(name not in atoms for name in neighbors):
        return None
    origin: Vector = atoms[center].xyz
    first: Vector = atoms[neighbors[0]].xyz - origin
    second: Vector = atoms[neighbors[1]].xyz - origin
    third: Vector = atoms[neighbors[2]].xyz - origin
    return float(np.dot(np.cross(first, second), third))


def canonical_residue_name(name: str) -> str:
    return "HIS" if name in {"HSD", "HSE", "HSP"} else name


def protein_bond_range(residue_name: str, name1: str, name2: str,
                       element1: str, element2: str) -> tuple[str, float, float]:
    """Broad operational limits by standard protein bond class."""
    pair: frozenset[str] = frozenset((name1, name2))
    residue: str = canonical_residue_name(residue_name)
    if pair == {"C", "O"}:
        return "backbone_carbonyl", 1.10, 1.35
    if pair == {"C", "OXT"}:
        return "terminal_carboxyl", 1.10, 1.45
    if pair == {"N", "CA"}:
        return "backbone_n_ca", 1.25, 1.65
    if pair == {"CA", "C"}:
        return "backbone_ca_c", 1.35, 1.70
    if (residue == "ASN" and pair == {"CG", "OD1"}) or (residue == "GLN" and pair == {"CD", "OE1"}):
        return "sidechain_carbonyl", 1.10, 1.35
    if (residue == "ASP" and name1 == "CG" and name2 in {"OD1", "OD2"}) or (
            residue == "GLU" and name1 == "CD" and name2 in {"OE1", "OE2"}):
        return "sidechain_carboxyl", 1.10, 1.40
    if residue in {"PHE", "TYR", "HIS", "TRP"} and not pair & {"CA", "CB", "OH", "N", "C", "O", "OXT"}:
        return "aromatic_ring", 1.20, 1.55
    if residue == "ARG" and "CZ" in pair:
        return "guanidinium_cn", 1.20, 1.50
    elements: frozenset[str] = frozenset((element1, element2))
    if elements == {"C", "S"}:
        return "carbon_sulfur", 1.65, 1.98
    if elements == {"C", "O"}:
        return "carbon_oxygen_single", 1.20, 1.60
    if elements == {"C", "N"}:
        return "carbon_nitrogen", 1.20, 1.65
    if elements == {"C"}:
        return "carbon_carbon_single", 1.35, 1.70
    return "other_standard_heavy_bond", 1.10, 2.00



def check_sequence(data: StructureData, targets: Mapping[str, str],
                   chain_ranges: Mapping[str, tuple[int, int]], allowed_ligands: set[str]) -> list[str]:
    expected = {(chain, chain_ranges[chain][0] + i, ""): aa
                for chain, sequence in targets.items() for i, aa in enumerate(sequence)}
    observed = {key: atoms for key, atoms in data.residues.items()
                if next(iter(atoms.values())).protein}
    issues = [f"missing_residue:{residue_label(k)}" for k in sorted(set(expected) - set(observed))]
    issues += [f"unexpected_residue:{residue_label(k)}" for k in sorted(set(observed) - set(expected))]
    for chain, sequence in targets.items():
        start, end = chain_ranges[chain]
        if len(sequence) != end - start + 1:
            issues.append(f"target_length:{chain} {len(sequence)}")
    for key, atoms in data.residues.items():
        first = next(iter(atoms.values()))
        label = residue_label(key)
        if not first.protein:
            if first.residue not in allowed_ligands:
                issues.append(f"unsupported_residue:{label} {first.residue}")
            continue
        if key in expected and AA_ONE[first.residue] != expected[key]:
            issues.append(f"sequence_mismatch:{label} {first.residue}/{expected[key]}")
        absent = EXPECTED_ATOMS[first.residue] - set(atoms)
        extra = set(atoms) - EXPECTED_ATOMS[first.residue] - {"OXT"}
        if absent or extra:
            issues.append(f"atom_coverage:{label} missing={sorted(absent)} extra={sorted(extra)}")
    return issues


def check_backbone(data: StructureData) -> list[str]:
    issues: list[str] = []
    for chain in sorted({a.chain for a in data.atoms if a.protein}):
        keys = sorted(k for k, v in data.residues.items() if k[0] == chain and next(iter(v.values())).protein)
        for left_key, right_key in zip(keys, keys[1:]):
            label = f"{residue_label(left_key)}/{residue_label(right_key)}"
            if right_key[1] != left_key[1] + 1 or left_key[2] or right_key[2]:
                issues.append(f"backbone_gap:{label}")
                continue
            left, right = data.residues[left_key], data.residues[right_key]
            if not {"CA", "C"} <= left.keys() or not {"N", "CA"} <= right.keys():
                issues.append(f"missing_backbone_atoms:{label}")
                continue
            if not 2.7 <= distance(left["CA"], right["CA"]) <= 4.1:
                issues.append(f"CA_distance:{label}")
            omega = dihedral(left["CA"].xyz, left["C"].xyz, right["N"].xyz, right["CA"].xyz)
            if omega is None:
                issues.append(f"undefined_omega:{label}")
            else:
                deviation = 180 - abs(omega)
                if right["CA"].residue == "PRO":
                    deviation = min(abs(omega), deviation)
                if deviation > 30.0:
                    issues.append(f"peptide_omega:{label} {omega:.2f} deg")
    return issues


def build_bond_graph(data: StructureData, template: StructureData) -> list[set[int]]:
    graph: list[set[int]] = [set() for _ in data.atoms]

    def connect(a: AtomRecord, b: AtomRecord) -> None:
        graph[a.index].add(b.index)
        graph[b.index].add(a.index)

    for key, atoms in data.residues.items():
        first = next(iter(atoms.values()))
        if first.protein:
            for a, b in STANDARD_BONDS[first.residue] + (("C", "OXT"),):
                if a in atoms and b in atoms:
                    connect(atoms[a], atoms[b])
            following = data.residues.get((key[0], key[1] + 1, ""), {})
            if not key[2] and "C" in atoms and "N" in following and distance(atoms["C"], following["N"]) <= 1.8:
                connect(atoms["C"], following["N"])
        elif not first.protein:
            reference = template.residues.get(key, {})
            for a, b in combinations(sorted(set(atoms) & set(reference)), 2):
                limit = COVALENT_RADII.get(reference[a].element, 0.7) + COVALENT_RADII.get(reference[b].element, 0.7) + 0.40
                if 0.7 <= distance(reference[a], reference[b]) <= limit:
                    connect(atoms[a], atoms[b])
    return graph


def check_clashes(data: StructureData, graph: list[set[int]]) -> list[str]:
    xyz = np.array([a.xyz for a in data.atoms], dtype=np.float64)
    pairs = sorted(cKDTree(xyz).query_pairs(2.0)) if data.atoms else []
    excluded: dict[int, set[int]] = {}
    issues: list[str] = []
    for i, j in pairs:
        if i not in excluded:
            excluded[i] = graph_neighbors_within_three(graph, i)
        if j not in excluded[i] and distance(data.atoms[i], data.atoms[j]) < 2.0:
            issues.append(f"clash:{data.atoms[i].label()}/{data.atoms[j].label()}")
    return issues


def check_chirality(data: StructureData, template: StructureData) -> list[str]:
    reference_values: dict[tuple[str, str], list[float]] = {}
    corresponding: dict[tuple[ResidueKey, str], float] = {}
    for key, atoms in template.residues.items():
        name = canonical_residue_name(next(iter(atoms.values())).residue)
        for center, neighbors in chiral_centers(name):
            value = signed_triple_product(atoms, center, neighbors)
            if value is not None:
                reference_values.setdefault((name, center), []).append(value)
                corresponding[key, center] = value
    signs: dict[tuple[str, str], int] = {}
    for identity, values in reference_values.items():
        nonplanar = [v for v in values if abs(v) >= 0.50]
        if nonplanar:
            signs[identity] = 1 if sum(v > 0 for v in nonplanar) > len(nonplanar) / 2 else -1
    issues: list[str] = []
    for label, structure in (("template", template), ("model", data)):
        for key, atoms in structure.residues.items():
            name = canonical_residue_name(next(iter(atoms.values())).residue)
            for center, neighbors in chiral_centers(name):
                value = signed_triple_product(atoms, center, neighbors)
                if value is None:
                    if label == "model":
                        issues.append(f"missing_chiral_center:{residue_label(key)}:{center}")
                    continue
                reference = corresponding.get((key, center))
                sign = signs.get((name, center))
                if sign is None and reference is not None and abs(reference) >= 0.50:
                    sign = 1 if reference > 0 else -1
                if (abs(value) < 0.50 or (sign is not None and value * sign < 0)
                        or (label == "model" and (sign is None or (reference is not None and value * reference < 0)))):
                    issues.append(f"chirality:{label}:{residue_label(key)}:{center} volume={value:.3f}")
    return issues


def check_bond_lengths(data: StructureData) -> list[str]:
    issues: list[str] = []

    def check(a: AtomRecord, b: AtomRecord, lower: float, upper: float) -> None:
        value = distance(a, b)
        if not lower <= value <= upper:
            issues.append(f"bond_length:{a.label()}/{b.label()} {value:.3f} A")

    for key, atoms in data.residues.items():
        name = next(iter(atoms.values())).residue
        if name not in STANDARD_BONDS:
            continue
        for left, right in STANDARD_BONDS[name] + (("C", "OXT"),):
            if left in atoms and right in atoms:
                _, lower, upper = protein_bond_range(name, left, right, atoms[left].element, atoms[right].element)
                check(atoms[left], atoms[right], lower, upper)
        following = data.residues.get((key[0], key[1] + 1, ""), {})
        if not key[2] and "C" in atoms and "N" in following and following["N"].protein:
            check(atoms["C"], following["N"], 1.15, 1.55)
    return issues


def check_local_geometry(data: StructureData, graph: list[set[int]]) -> list[str]:
    issues: list[str] = []
    for atom in data.atoms:
        visited, frontier = {atom.index}, {atom.index}
        for steps in range(1, 4):
            frontier = set().union(*(graph[i] for i in frontier)) - visited if frontier else set()
            visited.update(frontier)
            if steps >= 2:
                for i in sorted(frontier):
                    if i > atom.index and distance(atom, data.atoms[i]) < 1.80:
                        issues.append(f"local_graph_collapse:{atom.label()}/{data.atoms[i].label()}")
    for key, atoms in data.residues.items():
        name = canonical_residue_name(next(iter(atoms.values())).residue)
        rings = AROMATIC_RINGS.get(name, ())
        groups = [tuple(dict.fromkeys(a for ring in rings for a in ring))] if rings else []
        if name in PLANAR_SIDECHAIN_GROUPS:
            groups.append(PLANAR_SIDECHAIN_GROUPS[name])
        for group in groups:
            if any(a not in atoms for a in group):
                issues.append(f"missing_rigid_group:{residue_label(key)}")
                continue
            xyz = np.array([atoms[a].xyz for a in group], dtype=np.float64)
            centered = xyz - xyz.mean(axis=0)
            _, _, vt = np.linalg.svd(centered, full_matrices=False)
            if np.max(np.abs(centered @ vt[-1])) > 0.15:
                issues.append(f"nonplanar_rigid_group:{residue_label(key)}")
        for ring in rings:
            if any(a not in atoms for a in ring):
                continue
            for i, center in enumerate(ring):
                left = atoms[ring[(i - 1) % len(ring)]].xyz - atoms[center].xyz
                right = atoms[ring[(i + 1) % len(ring)]].xyz - atoms[center].xyz
                angle = float(np.degrees(np.arctan2(np.linalg.norm(np.cross(left, right)), np.dot(left, right))))
                if not 95 <= angle <= 145:
                    issues.append(f"aromatic_angle:{residue_label(key)}:{center} {angle:.2f} deg")
    return issues


def screen_segments(data: StructureData, new_residues: set[ResidueKey]) -> list[str]:
    segments: list[tuple[AtomRecord, AtomRecord]] = []
    neighborhood = {(c, n + delta, ins) for c, n, ins in new_residues for delta in (-1, 0, 1)}
    for key, atoms in sorted(data.residues.items()):
        following = data.residues.get((key[0], key[1] + 1, ""), {})
        if (not key[2] and "CA" in atoms and "CA" in following and atoms["CA"].protein
                and following["CA"].protein and 2.7 <= distance(atoms["CA"], following["CA"]) <= 4.1):
            segments.append((atoms["CA"], following["CA"]))
    issues: list[str] = []
    for (a, b), (c, d) in combinations(segments, 2):
        if not any(atom.residue_key in neighborhood for atom in (a, b, c, d)):
            continue
        if a.chain == c.chain and min(abs(x.number - y.number) for x in (a, b) for y in (c, d)) <= 1:
            continue
        if np.linalg.norm((a.xyz + b.xyz - c.xyz - d.xyz) / 2) > 5.1:
            continue
        separation, _, _ = segment_distance(a.xyz, b.xyz, c.xyz, d.xyz)
        if separation < 1.0:
            issues.append(f"CA_segments_near:{a.label()}-{b.label()}/{c.label()}-{d.label()}")
    return issues


def inspect_structure(model: StructureData, template: StructureData,
                      target_sequences: Mapping[str, str],
                      chain_ranges: Mapping[str, tuple[int, int]]) -> QCReport:
    if not chain_ranges or set(target_sequences) != set(chain_ranges):
        raise ValueError("FASTA chains must match the supplied chain ranges")
    model_protein = {a.residue_key for a in model.atoms if a.protein}
    reference_protein = {a.residue_key for a in template.atoms if a.protein}
    graph = build_bond_graph(model, template)
    allowed_ligands = {a.residue for a in template.atoms if not a.protein}
    issues: dict[str, list[str] | list[dict[str, Json]]] = {
        "read": model.read_issues + template.read_issues,
        "sequence": check_sequence(model, target_sequences, chain_ranges, allowed_ligands),
        "backbone": check_backbone(model),
        "bond_lengths": check_bond_lengths(model),
        "clashes": check_clashes(model, graph),
        "chirality": check_chirality(model, template),
        "local_geometry": check_local_geometry(model, graph),
        "segments": screen_segments(model, model_protein - reference_protein),
    }
    return {"summary": {"has_blocking_geometry_issues": any(bool(v) for v in issues.values()),
                        "protein_residues": len(model_protein), "new_residues": len(model_protein - reference_protein),
                        "heavy_atoms": len(model.atoms)}, "issues": issues}
