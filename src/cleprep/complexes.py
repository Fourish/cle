"""Mapped heavy-complex geometry screens and chemistry-preserving PDB export."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
from typing import TypeAlias, cast

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

from .common.geometry import Array, segment_triangle_intersections
from .common.io import write_json
from .common.peptides import peptide_backbone, stereo_labels, write_sdf
from .common.protein import require_no_cholesterol
from .common.topology import read_topology
from .docking import ligand_geometry_issues

Record: TypeAlias = dict[str, object]
AtomKey: TypeAlias = tuple[str, str, str]


@dataclass(frozen=True)
class GeometryPolicy:
    external_vdw_ratio: float = 0.75
    internal_vdw_ratio: float = 0.75
    minimum_external_distance: float = 1.8
    contact_distance: float = 4.0
    minimum_contact_atoms: int = 8
    minimum_contact_residues: int = 3


@dataclass
class HeavyReceptor:
    path: Path
    lines: list[str]
    keys: list[str]
    coordinates: Array
    radii: Array
    tree: cKDTree
    heavy_bonds: list[tuple[int, int]]
    backbone_bonds: list[tuple[int, int]]
    protein_chains: tuple[str, ...]
    contact_chains: tuple[str, ...]


def read_receptor(pdb: Path, atom_identity: Path, topology_bonds: Path, *,
                  protein_chains: tuple[str, ...], contact_chains: tuple[str, ...],
                  expected_residues: int | None = None) -> HeavyReceptor:
    """Read all retained heavy atoms and declared bonds, without guessing bonds."""
    lines: list[str] = []
    keys: list[str] = []
    identities: dict[AtomKey, int] = {}
    coordinates: list[list[float]] = []
    radii: list[float] = []
    residues: dict[tuple[str, str], set[str]] = {}
    table = Chem.GetPeriodicTable()
    models = 0
    for line in pdb.read_text().splitlines():
        if line.startswith("MODEL "):
            models += 1
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        symbol = line[76:78].strip()
        if symbol in {"H", "D"}:
            continue
        if not symbol or line[16:17].strip() or line[26:27].strip():
            raise ValueError("Receptor needs explicit elements, unique atoms and no altloc/insertion codes")
        number = table.GetAtomicNumber(symbol)
        if number <= 1:
            raise ValueError(f"Invalid receptor heavy element: {symbol}")
        key = (line[21], str(int(line[22:26])), line[12:16].strip())
        if key in identities:
            raise ValueError(f"Duplicate receptor atom: {key}")
        identities[key] = len(lines)
        lines.append(line)
        keys.append(f"{key[0]}:{key[1]}:{line[17:20].strip()}:{key[2]}")
        coordinates.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        radii.append(table.GetRvdw(number))
        if key[0] in protein_chains:
            residues.setdefault(key[:2], set()).add(key[2])
    xyz = np.asarray(coordinates, dtype=np.float64)
    if models > 1 or not lines or not np.isfinite(xyz).all():
        raise ValueError("Receptor must be a single nonempty model with finite coordinates")
    if set(protein_chains) != {key[0] for key in residues} or not set(contact_chains) <= set(protein_chains):
        raise ValueError("Declared protein/contact chains are missing or inconsistent")
    if any(not {"N", "CA", "C"} <= names for names in residues.values()):
        raise ValueError("Protein backbone is incomplete")
    if expected_residues is not None and len(residues) != expected_residues:
        raise ValueError(f"Expected {expected_residues} protein residues, observed {len(residues)}")
    rows, stored = read_topology(atom_identity, topology_bonds)
    source_to_heavy: dict[int, int] = {}
    matched: set[int] = set()
    for row in rows:
        key = (row["chain"], row["residue_id"], row["atom_name"])
        if key in identities:
            index = identities[key]
            if index in matched:
                raise ValueError(f"Duplicate topology identity: {key}")
            # CHARMM names such as CHL1 use their three-character PDB alias.
            if row["element"] != lines[index][76:78].strip() or row["residue_name"][:3] != lines[index][17:20].strip():
                raise ValueError(f"Coordinate/topology identity mismatch: {key}")
            source_to_heavy[int(row["index"])] = index
            matched.add(index)
    if set(source_to_heavy.values()) != set(range(len(lines))):
        raise ValueError("Topology identities do not cover every retained receptor heavy atom")
    bonds = sorted({tuple(sorted((source_to_heavy[a], source_to_heavy[b])))
                    for a, b in stored if a in source_to_heavy and b in source_to_heavy})
    if {i for pair in bonds for i in pair} != set(range(len(lines))):
        raise ValueError("Retained receptor heavy atoms have missing topology bonds")
    backbone = [(a, b) for a, b in bonds if lines[a][21] in protein_chains and lines[b][21] in protein_chains
                and lines[a][12:16].strip() in {"N", "CA", "C"} and lines[b][12:16].strip() in {"N", "CA", "C"}]
    for chain in protein_chains:
        atoms = {i for i, line in enumerate(lines) if line[21] == chain and line[12:16].strip() in {"N", "CA", "C"}}
        neighbors: dict[int, set[int]] = {i: set() for i in atoms}
        for a, b in backbone:
            if a in atoms and b in atoms:
                neighbors[a].add(b)
                neighbors[b].add(a)
        reached: set[int] = set()
        pending = [min(atoms)]
        while pending:
            index = pending.pop()
            if index not in reached:
                reached.add(index)
                pending.extend(neighbors[index] - reached)
        if reached != atoms:
            raise ValueError(f"Declared backbone topology is disconnected in protein chain {chain}")
    return HeavyReceptor(pdb, lines, keys, xyz, np.asarray(radii, dtype=np.float64),
                         cKDTree(xyz), bonds, backbone, protein_chains, contact_chains)


def check_receptor_frame(receptor: HeavyReceptor, pdbqt: Path) -> None:
    """Require identical heavy components and rounded coordinates in PDB/PDBQT."""
    elements = {"A": "C", "NA": "N", "OA": "O", "SA": "S"}
    actual: Counter[tuple[str, float, float, float]] = Counter()
    for line in pdbqt.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        ad_type = line.split()[-1]
        if ad_type in {"H", "HD", "HS"}:
            continue
        symbol = elements.get(ad_type, ad_type)
        actual[(symbol, *(round(float(line[a:b]), 3) for a, b in ((30, 38), (38, 46), (46, 54))))] += 1
    expected = Counter((line[76:78].strip(), *(round(float(v), 3) for v in xyz))
                       for line, xyz in zip(receptor.lines, receptor.coordinates, strict=True))
    if actual != expected:
        raise ValueError("Receptor PDB and PDBQT have different heavy components or coordinate frames")


def heavy_graph(mol: Chem.Mol) -> tuple[object, ...]:
    """Index-sensitive chemistry signature for a mapped heavy molecule."""
    atoms = tuple((a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(), a.GetIsAromatic()) for a in mol.GetAtoms())
    bonds = tuple(sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                          b.GetBondTypeAsDouble(), b.GetIsAromatic()) for b in mol.GetBonds()))
    return atoms, bonds, tuple(sorted(stereo_labels(Chem.Mol(mol)).items()))


def mapped_pose(pose: Chem.Mol, source: Chem.Mol) -> tuple[Chem.Mol, list[int]]:
    """Validate complete mapping and restore original heavy-atom order."""
    indices = [a.GetIdx() for a in source.GetAtoms() if a.GetAtomicNum() > 1]
    if not pose.HasProp("pose_atom_to_source_index0"):
        raise ValueError("Missing pose-to-source atom mapping")
    mapping = json.loads(pose.GetProp("pose_atom_to_source_index0"))
    if (not isinstance(mapping, list) or any(type(i) is not int for i in mapping)
            or len(mapping) != pose.GetNumAtoms() or len(set(mapping)) != len(mapping) or set(mapping) != set(indices)):
        raise ValueError("Incomplete or duplicated pose-to-source atom mapping")
    reverse = {index: local for local, index in enumerate(mapping)}
    ordered = Chem.RenumberAtoms(pose, [reverse[i] for i in indices])
    reference = Chem.RemoveAllHs(source)
    if heavy_graph(ordered) != heavy_graph(reference):
        raise ValueError("Mapped atom elements, charges, bond orders or declared stereochemistry changed")
    ordered.SetProp("pose_atom_to_source_index0", json.dumps(indices))
    return ordered, indices


def interface_metrics(xyz: Array, atomic_numbers: list[int], receptor: HeavyReceptor,
                      policy: GeometryPolicy) -> Record:
    """Measure contacts and clashes against every retained complex component."""
    table = Chem.GetPeriodicTable()
    radii = np.asarray([table.GetRvdw(number) for number in atomic_numbers], dtype=np.float64)
    distances, nearest = receptor.tree.query(xyz, k=1)
    # Include every pair that could violate a cutoff, even with larger atoms.
    cutoff = max(5.0, policy.contact_distance, policy.minimum_external_distance,
                 max(0.8, policy.external_vdw_ratio) * float(radii.max() + receptor.radii.max()))
    contacts: set[str] = set()
    contacting: set[int] = set()
    severe: list[Record] = []
    moderate: list[Record] = []
    failures: list[Record] = []
    sensitivity_count = 0
    minimum_ratio = float("inf")
    for i, targets in enumerate(receptor.tree.query_ball_point(xyz, cutoff)):
        if not targets:
            continue
        d = np.linalg.norm(receptor.coordinates[targets] - xyz[i], axis=1)
        ratios = d / (radii[i] + receptor.radii[targets])
        minimum_ratio = min(minimum_ratio, float(ratios.min()))
        for j, distance, ratio in zip(targets, d, ratios, strict=True):
            residue = ":".join(receptor.keys[j].split(":")[:3])
            if distance <= policy.contact_distance:
                contacts.add(residue)
                if receptor.lines[j][21] in receptor.contact_chains:
                    contacting.add(i)
            pair: Record = {"pose_heavy_index0": i, "receptor_atom": receptor.keys[j],
                            "distance_A": float(distance), "vdw_ratio": float(ratio)}
            if ratio < 0.65 or distance < 1.8:
                severe.append(pair)
            if ratio < 0.75:
                moderate.append(pair)
            if ratio < policy.external_vdw_ratio or distance < policy.minimum_external_distance:
                failures.append(pair)
            sensitivity_count += int(ratio < 0.8)
    if not np.isfinite(minimum_ratio):
        minimum_ratio = float(np.min(distances / (radii + receptor.radii[nearest])))
    contact_residues = sorted(key for key in contacts if key.split(":")[0] in receptor.contact_chains)
    return {"minimum_complex_distance_A": float(distances.min()), "minimum_vdw_ratio": minimum_ratio,
            "severe_overlap_count": len(severe), "overlap_count_vdw_ratio_lt0p75": len(moderate),
            "overlap_count_vdw_ratio_lt0p8": sensitivity_count, "external_failure_pair_count": len(failures),
            "external_failure_pairs_first20": failures[:20], "contact_residues": sorted(contacts),
            "target_contact_residues": contact_residues, "target_contact_atom_count": len(contacting),
            "external_geometry_passed": not failures,
            "meaningful_target_contacts": len(contact_residues) >= policy.minimum_contact_residues
                and len(contacting) >= policy.minimum_contact_atoms}


def internal_steric_metrics(molecule: Chem.Mol, policy: GeometryPolicy) -> Record:
    """Exclude graph distances <=3, including bonded and 1-4 pairs."""
    xyz = np.asarray(molecule.GetConformer().GetPositions(), dtype=np.float64)
    first, second = np.where(np.triu(Chem.GetDistanceMatrix(molecule) > 3.0, 1))
    radii = np.asarray([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in molecule.GetAtoms()])
    distances = np.linalg.norm(xyz[first] - xyz[second], axis=1)
    ratios = distances / (radii[first] + radii[second])
    return {"internal_nonbonded_minimum_vdw_ratio": float(ratios.min()) if len(ratios) else None,
            "internal_nonbonded_pairs_below_0p65": int(np.sum(ratios < 0.65)),
            "internal_nonbonded_pairs_below_0p75": int(np.sum(ratios < 0.75)),
            "internal_nonbonded_pairs_below_0p8": int(np.sum(ratios < 0.8)),
            "internal_geometry_passed": not bool(np.any(ratios < policy.internal_vdw_ratio)),
            "internal_pair_exclusion": "Shortest topological distance <=3"}


def threading_diagnostic(xyz: Array, ordered_ring: list[int], receptor: HeavyReceptor) -> Record:
    """Conservative backbone/ring centroid-fan penetration diagnostic."""
    ring = xyz[ordered_ring]
    center = ring.mean(axis=0)
    radius = float(np.max(np.linalg.norm(ring - center, axis=1))) + 2.0
    bonds = receptor.backbone_bonds
    starts = receptor.coordinates[[a for a, _ in bonds]]
    ends = receptor.coordinates[[b for _, b in bonds]]
    selected = np.flatnonzero(np.sum(((starts + ends) * 0.5 - center) ** 2, axis=1) <= radius * radius)
    hits: list[Record] = []
    seen: list[Array] = []
    for i in range(len(ring)):
        for j, point in segment_triangle_intersections(starts[selected], ends[selected], center, ring[i], ring[(i + 1) % len(ring)]):
            if any(float(np.linalg.norm(point - previous)) < 0.01 for previous in seen):
                continue
            seen.append(point)
            a, b = bonds[int(selected[j])]
            hits.append({"receptor_bond": [receptor.keys[a], receptor.keys[b]],
                         "intersection_xyz_A": point.tolist(), "fan_triangle_index0": i})
    return {"potential_backbone_threading": bool(hits), "ring_fan_intersections": hits,
            "scope": "Backbone segments intersecting a nonplanar ring centroid fan; conservative flag, not a topological proof"}


@dataclass
class PoseAudit:
    molecule: Chem.Mol | None
    source_indices: list[int]
    backbone: list[int]
    report: Record


def audit_pose(pose: Chem.Mol, source: Chem.Mol, receptor: HeavyReceptor,
               policy: GeometryPolicy) -> PoseAudit:
    issues: list[str] = []
    try:
        ordered, mapping = mapped_pose(pose, source)
        if ordered.GetNumConformers() != 1 or not ordered.GetConformer().Is3D():
            raise ValueError("Pose requires one 3D conformer")
        xyz = np.asarray(ordered.GetConformer().GetPositions(), dtype=np.float64)
        if not np.isfinite(xyz).all():
            raise ValueError("Nonfinite pose coordinates")
        reference = Chem.RemoveAllHs(source)
        backbone, _, _ = peptide_backbone(reference)
        issues.extend(ligand_geometry_issues(ordered, reference))
    except ValueError as error:
        return PoseAudit(None, [], [], {"passed": False, "issues": [str(error)]})
    external = interface_metrics(xyz, [a.GetAtomicNum() for a in ordered.GetAtoms()], receptor, policy)
    internal = internal_steric_metrics(ordered, policy)
    threading = threading_diagnostic(xyz, backbone, receptor)
    if not external["external_geometry_passed"]:
        issues.append("external_complex_overlap")
    if not external["meaningful_target_contacts"]:
        issues.append("insufficient_target_contacts")
    if not internal["internal_geometry_passed"]:
        issues.append("internal_nonbonded_overlap")
    if threading["potential_backbone_threading"]:
        issues.append("potential_backbone_threading")
    return PoseAudit(ordered, mapping, backbone, {"passed": not issues, "issues": issues,
                     "source_heavy_graph_matches": True, "interface": external,
                     "internal": internal, "threading": threading})


def same_cluster(first: PoseAudit, second: PoseAudit, *, rmsd_cutoff: float,
                 contact_jaccard: float) -> bool:
    if first.molecule is None or second.molecule is None or first.backbone != second.backbone:
        return False
    xyz1 = first.molecule.GetConformer().GetPositions()[first.backbone]
    xyz2 = second.molecule.GetConformer().GetPositions()[second.backbone]
    rmsd = float(np.sqrt(np.mean(np.sum((xyz1 - xyz2) ** 2, axis=1))))
    contacts1 = set(cast(list[str], cast(Record, first.report["interface"])["target_contact_residues"]))
    contacts2 = set(cast(list[str], cast(Record, second.report["interface"])["target_contact_residues"]))
    union = contacts1 | contacts2
    similarity = len(contacts1 & contacts2) / len(union) if union else 0.0
    return rmsd <= rmsd_cutoff and similarity >= contact_jaccard


def export_complex(output: Path, audit: PoseAudit, source: Chem.Mol,
                   receptor: HeavyReceptor, provenance: Record, *, ligand_chain: str,
                   policy: GeometryPolicy = GeometryPolicy()) -> Record:
    """Write a heavy PDB, SDF and complete source mapping for an accepted pose."""
    if not audit.report["passed"] or audit.molecule is None:
        raise ValueError("Only passing poses can be exported")
    require_no_cholesterol(receptor.lines, str(receptor.path))
    if len(ligand_chain) != 1 or ligand_chain in {line[21] for line in receptor.lines}:
        raise ValueError("Ligand chain must be one character unused by the receptor")
    molecule = audit.molecule
    receptor_chains = list(dict.fromkeys(line[21] for line in receptor.lines))
    if len(receptor.lines) + molecule.GetNumAtoms() + len(receptor_chains) + 1 > 99999:
        raise ValueError("Complex exceeds PDB atom serial capacity")
    xyz = molecule.GetConformer().GetPositions()
    if np.max(np.abs(xyz)) >= 1000:
        raise ValueError("Ligand coordinates exceed supported PDB field width")
    rounded = Chem.Mol(molecule)
    for i, position in enumerate(np.round(xyz, 3)):
        rounded.GetConformer().SetAtomPosition(i, position)
    pdb_audit = audit_pose(rounded, source, receptor, policy)
    if not pdb_audit.report["passed"]:
        raise ValueError("PDB coordinate rounding changes geometry acceptance")
    output.mkdir(parents=True, exist_ok=False)
    lines: list[str] = []
    serial = 1
    receptor_serials = [0] * len(receptor.lines)
    chain_order = {chain: i for i, chain in enumerate(receptor_chains)}
    receptor_order = sorted(range(len(receptor.lines)),
                            key=lambda i: (chain_order[receptor.lines[i][21]], int(receptor.lines[i][22:26])))
    previous: str | None = None
    for index in receptor_order:
        line = receptor.lines[index]
        if previous is not None and line[21] != previous[21]:
            lines.append(f"TER   {serial:5d}      {previous[17:27]}")
            serial += 1
        receptor_serials[index] = serial
        lines.append(line[:6] + f"{serial:5d}" + line[11:])
        serial += 1
        previous = line
    if previous is not None:
        lines.append(f"TER   {serial:5d}      {previous[17:27]}")
        serial += 1
    mapping: list[Record] = []
    named = all(all(source.GetAtomWithIdx(index).HasProp(prop) for prop in ("residue_id", "residue_name", "atom_name"))
                for index in audit.source_indices)
    names_seen: set[tuple[int, str]] = set()
    labels: list[tuple[int, str, str, str]] = []
    for i, source_index in enumerate(audit.source_indices):
        atom = source.GetAtomWithIdx(source_index)
        symbol = molecule.GetAtomWithIdx(i).GetSymbol()
        resid = int(atom.GetProp("residue_id")) if named else 1
        resname = atom.GetProp("residue_name") if named else "LIG"
        name = atom.GetProp("atom_name") if named else f"{symbol}{i + 1}"
        if not (1 <= resid <= 9999 and 1 <= len(resname) <= 3 and 1 <= len(name) <= 4) or (resid, name) in names_seen:
            raise ValueError("Source atom/residue names are invalid for unique PDB export")
        names_seen.add((resid, name))
        labels.append((resid, resname, name, symbol))
    ligand_serials = [0] * molecule.GetNumAtoms()
    # Oxidation atoms can be appended after other residues in the source SDF.
    # Regroup PDB records without changing the SDF order or any coordinates.
    for i in sorted(range(molecule.GetNumAtoms()), key=lambda i: labels[i][0]):
        resid, resname, name, symbol = labels[i]
        source_index = audit.source_indices[i]
        ligand_serials[i] = serial
        field = f" {name:<3}" if len(name) < 4 else name
        x, y, z = xyz[i]
        lines.append(f"HETATM{serial:5d} {field:4s} {resname:3s} {ligand_chain}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {symbol:>2s}  ")
        mapping.append({"pose_heavy_index0": i, "source_atom_index0": source_index, "complex_serial1": serial,
                        "chain": ligand_chain, "resid": resid, "resname": resname, "name": name, "element": symbol})
        serial += 1
    lines.append(f"TER   {serial:5d}      {lines[-1][17:27]}")
    mapping.sort(key=lambda row: cast(int, row["pose_heavy_index0"]))
    complex_bonds = [(receptor_serials[a], receptor_serials[b]) for a, b in receptor.heavy_bonds]
    complex_bonds += [(ligand_serials[b.GetBeginAtomIdx()], ligand_serials[b.GetEndAtomIdx()])
                     for b in molecule.GetBonds()]
    lines += [f"CONECT{a:5d}{b:5d}" for a, b in complex_bonds]
    pdb = output / "complex_heavy.pdb"
    pdb.write_text("REMARK Geometry-filtered heavy complex; ligand chemistry is defined by companion SDF.\n"
                   + "\n".join(lines) + "\nEND\n", encoding="ascii")
    write_sdf(output / "peptide_heavy.sdf", [molecule])
    manifest = {**provenance, "geometry": audit.report, "pdb_geometry": pdb_audit.report,
                "protein_chains": list(receptor.protein_chains), "cholesterol_retained": False,
                "ligand_chain": ligand_chain, "ligand_naming": "source_residues" if named else "single_LIG_residue",
                "pdb_atom_order": "grouped_by_chain_and_residue",
                "receptor_heavy_atom_count": len(receptor.lines), "ligand_heavy_atom_count": molecule.GetNumAtoms(),
                "source_mapping": mapping, "receptor_source_mapping": [
                    {"complex_serial1": receptor_serials[i], "source_serial1": int(line[6:11]), "atom_key": key}
                    for i, (line, key) in enumerate(zip(receptor.lines, receptor.keys, strict=True))],
                "heavy_bonds_serial1": [list(pair) for pair in complex_bonds],
                "interpretation": "Static geometry-filtered candidate; hydrogens, complete parameters and membrane preparation required before MD"}
    write_json(output / "manifest.json", manifest)
    return manifest
