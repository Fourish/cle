"""Shared cyclic-peptide chemistry helpers (zero-based atom indices)."""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

Omega = tuple[int, int, int, int]


def read_molecule(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, strictParsing=True)
    if len(supplier) != 1 or supplier[0] is None:
        raise ValueError(f"Expected exactly one readable molecule: {path}")
    return supplier[0]


def stereo_labels(mol: Chem.Mol) -> dict[int, str]:
    """Use modern CIP; legacy labels can be reversed for these sulfoxides."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            atom.ClearProp("_CIPCode")
    rdCIPLabeler.AssignCIPLabels(mol)
    return {atom.GetIdx(): atom.GetProp("_CIPCode")
            for atom in mol.GetAtoms() if atom.HasProp("_CIPCode")}


def sulfoxide_atoms(mol: Chem.Mol) -> list[int]:
    return [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16
            and any(bond.GetOtherAtom(atom).GetAtomicNum() == 8
                    and bond.GetBondType() == Chem.BondType.DOUBLE
                    for bond in atom.GetBonds())]


def peptide_backbone(mol: Chem.Mol) -> tuple[list[int], list[Omega], list[Omega]]:
    """Follow the single C(=O)-N cycle and return zero-based backbone indices."""
    query = Chem.MolFromSmarts("[NX3]-[CX4]-[CX3](=[OX1])")
    cores = list(mol.GetSubstructMatches(query))
    if not cores or len({core[1] for core in cores}) != len(cores):
        raise ValueError("Missing or ambiguous peptide backbone")
    nitrogen = {core[0]: i for i, core in enumerate(cores)}
    successor: dict[int, int] = {}
    for i, (_, _, carbon, _) in enumerate(cores):
        following = [nitrogen[a.GetIdx()] for a in mol.GetAtomWithIdx(carbon).GetNeighbors()
                     if a.GetIdx() in nitrogen]
        if len(following) != 1:
            raise ValueError("Peptide backbone is open or branched")
        successor[i] = following[0]
    if len(set(successor.values())) != len(cores):
        raise ValueError("Peptide backbone has duplicate connections")
    start = min(range(len(cores)), key=lambda i: cores[i][1])
    order = [start]
    while successor[order[-1]] not in order:
        order.append(successor[order[-1]])
    if successor[order[-1]] != start or len(order) != len(cores):
        raise ValueError("Backbone is not one cycle containing all peptide units")
    backbone = [atom for i in order for atom in cores[i][:3]]
    amides: list[Omega] = []
    pro_amides: list[Omega] = []
    for i in order:
        n, ca, _, _ = cores[successor[i]]
        omega = (cores[i][1], cores[i][2], n, ca)
        amides.append(omega)
        if mol.GetRingInfo().AreAtomsInSameRingOfSize(n, ca, 5):
            pro_amides.append(omega)
    return backbone, amides, pro_amides


def write_sdf(path: Path, molecules: list[Chem.Mol]) -> None:
    with Chem.SDWriter(str(path)) as writer:
        writer.SetForceV3000(True)
        for mol in molecules:
            writer.write(mol)


