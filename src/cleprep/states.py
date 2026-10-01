"""Prepare cyclic-peptide states while preserving specified stereochemistry."""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors

from itertools import product

from .common.peptides import peptide_backbone, read_molecule, stereo_labels, sulfoxide_atoms


def prepare_parent(path: Path, *, expected_residues: int | None = None,
                   expected_prolines: int | None = None, expected_sulfoxides: int | None = None,
                   expected_charge: int | None = None) -> list[Chem.Mol]:
    parent = path.stem
    source = read_molecule(path)
    charge = Chem.GetFormalCharge(source)
    if len(Chem.GetMolFrags(source)) != 1 or (expected_charge is not None and charge != expected_charge):
        raise ValueError(f"{parent}: connectivity or formal charge mismatch")
    source_cip = stereo_labels(source)
    sulfurs = sulfoxide_atoms(source)
    if expected_sulfoxides is not None and len(sulfurs) != expected_sulfoxides:
        raise ValueError(f"{parent}: unexpected number of sulfoxides")
    unspecified = [i for i in sulfurs if i not in source_cip]
    states: list[Chem.Mol] = []
    for labels in product("RS", repeat=len(unspecified)):
        mol = Chem.Mol(source)
        for index, label in zip(unspecified, labels, strict=True):
            sulfur = mol.GetAtomWithIdx(index)
            sulfur.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            if stereo_labels(mol).get(index) != label:
                sulfur.InvertChirality()
        sulfur_cip = "".join(stereo_labels(mol)[i] for i in sulfurs)
        state_id = f"{parent}_S-{sulfur_cip}" if unspecified else parent
        mol = Chem.AddHs(mol, addCoords=False)
        cip = stereo_labels(mol)
        if any(cip.get(i) != label for i, label in source_cip.items()):
            raise ValueError(f"{state_id}: source stereochemistry changed")
        if any(cip.get(i) != label for i, label in zip(unspecified, labels, strict=True)):
            raise ValueError(f"{state_id}: sulfur CIP assignment failed")
        if any(info.specified == Chem.StereoSpecified.Unspecified
               for info in Chem.FindPotentialStereo(mol)):
            raise ValueError(f"{state_id}: unresolved stereochemistry")
        if (rdMolDescriptors.CalcMolFormula(mol) != rdMolDescriptors.CalcMolFormula(source)
                or Chem.GetFormalCharge(mol) != charge):
            raise ValueError(f"{state_id}: formula or charge changed after adding hydrogens")
        _, amides, pro_amides = peptide_backbone(mol)
        if ((expected_residues is not None and len(amides) != expected_residues)
                or (expected_prolines is not None and len(pro_amides) != expected_prolines)):
            raise ValueError(f"{state_id}: unexpected peptide-unit or Pro count")
        mol.SetProp("_Name", state_id)
        mol.SetProp("state_id", state_id)
        mol.SetProp("sulfur_cip", sulfur_cip or "none")
        mol.SetProp("coordinate_status", "2D chemical state; not a 3D conformer")
        rdDepictor.Compute2DCoords(mol)
        states.append(mol)
    return states
