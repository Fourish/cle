"""Build five HighFold3 monomer inputs from this repository's peptide SDFs.

Run with the cle environment. The CCD JSON is exported from the deployed AF3
dictionary so that standard atom names and stereochemistry are checked against
the actual inference installation, not inferred from a sequence label.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cleprep.common.peptides import peptide_backbone, stereo_labels, sulfoxide_atoms, write_sdf
from cleprep.states import prepare_parent


STANDARD = {"MET": "M", "LEU": "L", "ILE": "I", "PRO": "P",
            "PHE": "F", "VAL": "V", "TRP": "W"}
MODIFIED = {"SME": "R", "MHO": "S"}


def implicit_hydrogens(mol):
    for atom in mol.GetAtoms():
        if not (atom.GetIsAromatic() and atom.GetSymbol() == "N"):
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(False)
    Chem.SanitizeMol(mol)
    return mol


def read_component(data):
    mol = Chem.RWMol()
    names = data["_chem_comp_atom.atom_id"]
    for name, element, charge in zip(names, data["_chem_comp_atom.type_symbol"],
                                    data["_chem_comp_atom.charge"], strict=True):
        atom = Chem.Atom(element.title())
        atom.SetFormalCharge(int(charge))
        atom.SetProp("atom_name", name)
        mol.AddAtom(atom)
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE}
    for a, b, order in zip(data["_chem_comp_bond.atom_id_1"],
                           data["_chem_comp_bond.atom_id_2"],
                           data["_chem_comp_bond.value_order"], strict=True):
        mol.AddBond(names.index(a), names.index(b), orders[order])
    conformer = Chem.Conformer(len(names))
    for i, xyz in enumerate(zip(*(data[f"_chem_comp_atom.pdbx_model_Cartn_{a}_ideal"]
                                  for a in "xyz"), strict=True)):
        conformer.SetAtomPosition(i, [float(v) for v in xyz])
    mol.AddConformer(conformer)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    return implicit_hydrogens(Chem.RemoveHs(mol))


def without_oxt(mol):
    edit = Chem.RWMol(mol)
    edit.RemoveAtom(next(a.GetIdx() for a in edit.GetAtoms()
                         if a.GetProp("atom_name") == "OXT"))
    return implicit_hydrogens(edit.GetMol())


def build(ccd_path, output, root, center_modification=False, seeds=(1,)):
    if not seeds or len(set(seeds)) != len(seeds) or any(seed < 0 for seed in seeds):
        raise ValueError("Seeds must be unique nonnegative integers")
    output.mkdir(parents=True, exist_ok=True)
    (output / "inputs").mkdir(exist_ok=True)
    (output / "states").mkdir(exist_ok=True)
    data = json.loads(ccd_path.read_text())
    templates = {k: read_component(data[k]) for k in (*STANDARD, *MODIFIED)}
    for code, sulfur_cip in MODIFIED.items():
        mol = templates[code]
        sulfur = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "S")
        if stereo_labels(mol).get(sulfur) != sulfur_cip:
            raise ValueError(f"{code}: CCD sulfur stereochemistry differs from expected {sulfur_cip}")
    fragments = {k: without_oxt(v) for k, v in templates.items()}
    profile = json.loads((root / "src/profiles/peptides.json").read_text())
    manifest = {"ccd_sha256": hashlib.sha256(ccd_path.read_bytes()).hexdigest(), "states": []}
    for parent, settings in profile["parents"].items():
        source = root / "data/inputs/peptides" / (parent + ".sdf")
        states = prepare_parent(source, **{"expected_" + k: v for k, v in settings.items()})
        for state in states:
            name = state.GetProp("_Name")
            mol = implicit_hydrogens(Chem.RemoveHs(state))
            bb, omegas, pro_omegas = peptide_backbone(mol)
            edit = Chem.RWMol(mol)
            for i in range(0, len(bb), 3):
                edit.RemoveBond(bb[i + 2], bb[(i + 3) % len(bb)])
            maps = []
            parts = Chem.GetMolFrags(edit.GetMol(), asMols=True, sanitizeFrags=False,
                                    fragsMolAtomMapping=maps)
            parts = sorted(zip(parts, maps), key=lambda p: bb.index(next(a for a in p[1] if a in bb)) // 3)
            rotation = 0
            if center_modification:
                modified = [i for i, (part, _) in enumerate(parts) if sulfoxide_atoms(part)]
                if modified:
                    if len(modified) != 1:
                        raise ValueError("Only a single sulfoxide residue is supported")
                    rotation = (modified[0] - len(parts) // 2) % len(parts)
                    parts = parts[rotation:] + parts[:rotation]
            codes, atom_mapping, assembly = [], [], None
            links, offset = [], 0
            for rid, (part, indices) in enumerate(parts, 1):
                part = implicit_hydrogens(part)
                matches = [(code, part.GetSubstructMatch(ref, useChirality=True))
                           for code, ref in fragments.items() if ref.GetNumAtoms() == part.GetNumAtoms()]
                matches = [(code, match) for code, match in matches if match]
                if len(matches) != 1:
                    raise ValueError(f"{name}:{rid}: ambiguous or unsupported residue {matches}")
                code, match = matches[0]
                ref = fragments[code]
                codes.append(code)
                by_name = {a.GetProp("atom_name"): a.GetIdx() for a in ref.GetAtoms()}
                links.append((offset + by_name["N"], offset + by_name["C"]))
                for a, idx in zip(ref.GetAtoms(), match, strict=True):
                    atom_mapping.append({"source_heavy_atom": indices[idx], "residue_id": rid,
                                         "atom_name": a.GetProp("atom_name"), "component": code})
                assembly = Chem.Mol(ref) if assembly is None else Chem.CombineMols(assembly, ref)
                offset += ref.GetNumAtoms()
            assembly = Chem.RWMol(assembly)
            for i, (_, carbon) in enumerate(links):
                assembly.AddBond(carbon, links[(i + 1) % len(links)][0], Chem.BondType.SINGLE)
            assembly = implicit_hydrogens(assembly.GetMol())
            if Chem.MolToSmiles(assembly) != Chem.MolToSmiles(mol):
                raise ValueError(f"{name}: full cyclic graph/stereochemistry round trip failed")
            sequence = "".join(STANDARD.get(code, "M") for code in codes)
            modifications = [{"ptmPosition": i, "ptmType": code} for i, code in enumerate(codes, 1) if code in MODIFIED]
            job = {"name": name, "dialect": "alphafold3", "version": 2, "modelSeeds": list(seeds),
                   "sequences": [{"protein": {"id": "A", "sequence": sequence,
                       "modifications": modifications, "unpairedMsa": "", "pairedMsa": "", "templates": []}}]}
            (output / "inputs" / f"{name}.json").write_text(json.dumps(job, indent=2) + "\n")
            write_sdf(output / "states" / f"{name}.sdf", [mol])
            manifest["states"].append({"state": name, "source": str(source.relative_to(root)),
                "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "sequence": sequence, "components": codes, "modifications": modifications,
                "first_residue_in_source_order": rotation + 1,
                "model_seeds": list(seeds),
                "sulfur_cip": state.GetProp("sulfur_cip"),
                "formula": rdMolDescriptors.CalcMolFormula(mol), "heavy_atoms": mol.GetNumAtoms(),
                "smiles": Chem.MolToSmiles(mol), "cip": stereo_labels(mol),
                "atom_mapping": atom_mapping, "backbone": bb,
                "amide_omegas": omegas, "proline_omegas": pro_omegas,
                "cyclic_graph_and_stereochemistry_roundtrip": True})
            print(name, sequence, codes, rdMolDescriptors.CalcMolFormula(mol), "graph/stereo PASS")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ccd-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--center-modification", action="store_true",
                        help="Rotate ring numbering to place the modified residue in the middle")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1])
    args = parser.parse_args()
    build(args.ccd_json, args.output, Path(__file__).resolve().parents[2], args.center_modification, args.seeds)
