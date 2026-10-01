"""Select explicit author/label chain ranges and ligands from an mmCIF."""

from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
from pathlib import Path

import gemmi

from typing import Mapping

from .common.io import write_fasta, write_json

Row = dict[str, str]


def rows(block: gemmi.cif.Block, category: str) -> list[Row]:
    table = block.find_mmcif_category(category)
    names = [tag[len(category):] for tag in table.tags]
    return [dict(zip(names, row, strict=True)) for row in table]


def pdb_atom_line(atom: Row, serial: int) -> str:
    name: str = atom["auth_atom_id"]
    element: str = atom["type_symbol"]
    atom_field: str = f" {name:<3}" if len(name) < 4 and len(element) == 1 else f"{name:<4}"
    altloc: str = " " if atom["label_alt_id"] in {".", "?"} else atom["label_alt_id"]
    insertion: str = " " if atom["pdbx_PDB_ins_code"] in {".", "?"} else atom["pdbx_PDB_ins_code"]
    return (
        f'{atom["group_PDB"]:<6}{serial:5d} {atom_field}{altloc}{atom["auth_comp_id"]:>3} '
        f'{atom["auth_asym_id"]}{int(atom["auth_seq_id"]):4d}{insertion}   '
        f'{float(atom["Cartn_x"]):8.3f}{float(atom["Cartn_y"]):8.3f}{float(atom["Cartn_z"]):8.3f}'
        f'{float(atom["occupancy"]):6.2f}{float(atom["B_iso_or_equiv"]):6.2f}          {element:>2}  '
    )


def write_pdb(path: Path, groups: list[tuple[bool, list[Row]]], bonds: list[Row]) -> None:
    lines: list[str] = ["REMARK Observed template; original author numbering; gaps unrepaired"]
    ligand_keys: dict[tuple[str, str], int] = {}
    serial: int = 1
    for is_protein, atoms in groups:
        if not atoms:
            continue
        for atom in atoms:
            lines.append(pdb_atom_line(atom, serial))
            if not is_protein:
                ligand_keys[(atom["label_asym_id"], atom["label_atom_id"])] = serial
            serial += 1
        if is_protein:
            last: Row = atoms[-1]
            lines.append(f'TER   {serial:5d}      {last["auth_comp_id"]:>3} {last["auth_asym_id"]}{int(last["auth_seq_id"]):4d}')
            serial += 1
    adjacency: dict[int, set[int]] = defaultdict(set)
    for bond in bonds:
        first: int | None = ligand_keys.get((bond["label_asym_id"], bond["atom_id_1"]))
        second: int | None = ligand_keys.get((bond["label_asym_id"], bond["atom_id_2"]))
        if first is not None and second is not None:
            adjacency[first].add(second)
            adjacency[second].add(first)
    for first, neighbors in sorted(adjacency.items()):
        ordered: list[int] = sorted(neighbors)
        for offset in range(0, len(ordered), 4):
            lines.append(f"CONECT{first:5d}" + "".join(f"{second:5d}" for second in ordered[offset:offset + 4]))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def select(source: Path, output: Path, targets: Mapping[str, tuple[str, int, int]],
           ligands_to_keep: Mapping[str, tuple[str, int | None]],
           expected_sha256: str | None = None) -> dict[str, object]:
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError("Source mmCIF checksum does not match the requested profile")
    if not targets or any(len(chain) != 1 or first > last for chain, (_, first, last) in targets.items()):
        raise ValueError("PDB output needs one-character author chains and nonempty residue ranges")
    block = gemmi.cif.read_file(str(source)).sole_block()
    atoms = rows(block, "_atom_site.")
    schemes = rows(block, "_pdbx_poly_seq_scheme.")
    if {a["pdbx_PDB_model_num"] for a in atoms} != {"1"} or any(a["label_alt_id"] != "." for a in atoms):
        raise ValueError("Expected one model without alternate locations")
    grouped: dict[tuple[str, str], list[Row]] = defaultdict(list)
    for atom in atoms:
        grouped[atom["label_asym_id"], atom["label_seq_id"]].append(atom)
    sequences: dict[str, str] = {}
    descriptions: dict[str, str] = {}
    target_residues = 0
    missing_residues = 0
    groups: list[tuple[bool, list[Row]]] = []
    for chain, (label, start, end) in targets.items():
        selected = sorted((r for r in schemes if r["asym_id"] == label
                           and start <= int(r["pdb_seq_num"]) <= end),
                          key=lambda r: int(r["pdb_seq_num"]))
        if [int(r["pdb_seq_num"]) for r in selected] != list(range(start, end + 1)):
            raise ValueError(f"Incomplete deposited sequence for {chain}")
        sequence = "".join(gemmi.find_tabulated_residue(r["mon_id"]).one_letter_code for r in selected)
        sequences[chain] = sequence
        descriptions[chain] = f"auth={start}-{end} label={label}"
        target_residues += len(selected)
        for residue in selected:
            observed = grouped.get((label, residue["seq_id"]), [])
            if any(a["auth_asym_id"] != chain or a["auth_seq_id"] != residue["pdb_seq_num"]
                   or a["label_comp_id"] != residue["mon_id"] for a in observed):
                raise ValueError("Atom/sequence identity mismatch")
            missing_residues += int(not observed)
        groups.append((True, [a for a in atoms if a["label_asym_id"] == label
                              and start <= int(a["auth_seq_id"]) <= end]))
    ligands = [a for a in atoms if a["label_asym_id"] in ligands_to_keep]
    observed_ligands = Counter((a["label_asym_id"], a["auth_comp_id"]) for a in ligands)
    for label, (name, expected_atoms) in ligands_to_keep.items():
        count = observed_ligands[label, name]
        if not count or (expected_atoms is not None and count != expected_atoms):
            raise ValueError(f"Ligand identity/count mismatch: {label}:{name}")
    if any(name != ligands_to_keep[label][0] for label, name in observed_ligands):
        raise ValueError("Unexpected ligand identity")
    groups.append((False, ligands))
    bonds: list[Row] = []
    for label, (name, _) in ligands_to_keep.items():
        bonds.extend({**r, "label_asym_id": label} for r in rows(block, "_chem_comp_bond.") if r["comp_id"] == name)
    output.mkdir(parents=True, exist_ok=False)
    write_pdb(output / "complex_template.pdb", groups, bonds)
    write_fasta(output / "target_sequences.fasta", sequences, descriptions)
    summary: dict[str, object] = {
        "target_residues": target_residues,
        "missing_residues": missing_residues,
        "protein_atoms": sum(len(group) for protein, group in groups if protein),
        "ligand_heavy_atoms": len(ligands),
    }
    write_json(output / "selection.json", summary)
    return summary
