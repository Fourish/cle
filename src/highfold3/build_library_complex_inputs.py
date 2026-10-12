"""Extend an existing four-protein HighFold3 baseline with audited cyclic SDFs."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from build_inputs import STANDARD, MODIFIED, implicit_hydrogens, read_component, without_oxt
from cleprep.common.peptides import peptide_backbone, read_molecule, stereo_labels, write_sdf


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(row, source, fragments, seeds):
    """Match actual residue graphs, then verify the reassembled cyclic molecule."""
    name = row["state"]
    mol = implicit_hydrogens(Chem.RemoveHs(read_molecule(source)))
    assert Chem.MolToSmiles(mol) == row["canonical_isomeric_smiles"], name
    assert rdMolDescriptors.CalcMolFormula(mol) == row["neutral_formula_from_graph"], name
    bb, omegas, pro_omegas = peptide_backbone(mol)
    edit = Chem.RWMol(mol)
    for i in range(0, len(bb), 3):
        edit.RemoveBond(bb[i + 2], bb[(i + 3) % len(bb)])
    maps = []
    parts = Chem.GetMolFrags(edit.GetMol(), asMols=True, sanitizeFrags=False,
                            fragsMolAtomMapping=maps)
    parts = sorted(zip(parts, maps),
                   key=lambda p: bb.index(next(a for a in p[1] if a in bb)) // 3)
    matches = []
    for part, indices in parts:
        part = implicit_hydrogens(part)
        found = [(code, part.GetSubstructMatch(ref, useChirality=True))
                 for code, ref in fragments.items() if ref.GetNumAtoms() == part.GetNumAtoms()]
        found = [(code, match) for code, match in found if match]
        if len(found) != 1:
            raise ValueError(f"{name}: unsupported or ambiguous residue {found}")
        code, match = found[0]
        matches.append((code, match, indices))
    modified = [i for i, (code, _, _) in enumerate(matches) if code not in STANDARD]
    n = len(matches)
    rotation = 0
    if modified:
        # The single-PTM case exactly retains the preceding baseline's rule.
        # For multiple PTMs, maximize distance from numbering ends, then center
        # the group on floor(n/2). No prediction outcome enters this choice.
        def score(r):
            positions = [(i - r) % n for i in modified]
            return (min(min(p, n - 1 - p) for p in positions),
                    -sum((p - n // 2) ** 2 for p in positions), -r)
        rotation = max(range(n), key=score)
        matches = matches[rotation:] + matches[:rotation]
    codes, mapping, links, assembly, offset = [], [], [], None, 0
    for rid, (code, match, indices) in enumerate(matches, 1):
        ref = fragments[code]
        codes.append(code)
        by_name = {a.GetProp("atom_name"): a.GetIdx() for a in ref.GetAtoms()}
        links.append((offset + by_name["N"], offset + by_name["C"]))
        for atom, idx in zip(ref.GetAtoms(), match, strict=True):
            mapping.append({"source_heavy_atom": indices[idx], "residue_id": rid,
                            "source_residue_id": (rotation + rid - 1) % n + 1,
                            "atom_name": atom.GetProp("atom_name"), "component": code})
        assembly = Chem.Mol(ref) if assembly is None else Chem.CombineMols(assembly, ref)
        offset += ref.GetNumAtoms()
    assembly = Chem.RWMol(assembly)
    for i, (_, carbon) in enumerate(links):
        assembly.AddBond(carbon, links[(i + 1) % n][0], Chem.BondType.SINGLE)
    assembly = implicit_hydrogens(assembly.GetMol())
    assert Chem.MolToSmiles(assembly) == Chem.MolToSmiles(mol), name
    assert sorted(a["source_heavy_atom"] for a in mapping) == list(range(mol.GetNumAtoms()))
    source_codes = [codes[(i - rotation) % n] for i in range(n)]
    source_sequence = [STANDARD.get(c, "MetO2" if c == "OMT" else "MetO") for c in source_codes]
    assert source_sequence == row["graph_residue_sequence_in_package_numbering"], (name, source_sequence)
    source_sulfur = {str(i + 1): MODIFIED[c] for i, c in enumerate(source_codes) if c in MODIFIED}
    assert source_sulfur == row["sulfoxide_cip_by_package_residue"], name
    state = {
        "state": name, "parent": row["package_parent_name"],
        "source": str(source), "source_sha256": sha256(source),
        "sequence": "".join(STANDARD.get(c, "M") for c in codes),
        "components": codes,
        "modifications": [{"ptmPosition": i, "ptmType": c}
                          for i, c in enumerate(codes, 1) if c not in STANDARD],
        "first_residue_in_source_order": rotation + 1, "model_seeds": seeds,
        "sulfur_cip": source_sulfur,
        "formula": rdMolDescriptors.CalcMolFormula(mol), "heavy_atoms": mol.GetNumAtoms(),
        "smiles": Chem.MolToSmiles(mol), "cip": stereo_labels(mol),
        "atom_mapping": mapping, "backbone": bb, "amide_omegas": omegas,
        "proline_omegas": pro_omegas, "cyclic_graph_and_stereochemistry_roundtrip": True,
        "package_graph_and_stereo_match": True, "peptide_chain": "A",
        "source_coordinate_type": row["coordinate_type"], "stereo_status": row["stereo_status"],
    }
    return mol, state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--ccd-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    args = parser.parse_args()
    assert args.seeds and len(set(args.seeds)) == len(args.seeds) and min(args.seeds) >= 0
    out = args.output
    assert not (out / "manifest.json").exists(), "Refusing to replace an existing prediction baseline"
    for folder in ("source", "states", "inputs", "feature_only_inputs"):
        (out / folder).mkdir(parents=True, exist_ok=True)
    prior = json.loads((args.baseline / "manifest.json").read_text())
    catalog_path = args.package / "structure_check_PMC10416992/audit.json"
    catalog = json.loads(catalog_path.read_text())
    excluded = {s["state"] for s in prior["states"]}
    rows = [s for s in catalog["states"] if s["state"] not in excluded]
    assert len(rows) == 25 and len({s["package_parent_name"] for s in rows}) == 14
    ccd = json.loads(args.ccd_json.read_text())
    templates = {code: read_component(ccd[code]) for code in (*STANDARD, *MODIFIED, "OMT")}
    for code, label in MODIFIED.items():
        sulfur = next(a.GetIdx() for a in templates[code].GetAtoms() if a.GetSymbol() == "S")
        assert stereo_labels(templates[code])[sulfur] == label
    fragments = {code: without_oxt(mol) for code, mol in templates.items()}
    prior_input = args.baseline / "inputs" / (prior["states"][0]["state"] + ".json")
    proteins = [s for s in json.loads(prior_input.read_text())["sequences"] if s["protein"]["id"] != "A"]
    expected = {p["id"]: p["sequence"] for p in prior["protein_chains"]}
    assert {p["protein"]["id"]: p["protein"]["sequence"] for p in proteins} == expected
    protein_payload_hash = hashlib.sha256(json.dumps(proteins, sort_keys=True).encode()).hexdigest()
    manifest = {k: copy.deepcopy(prior[k]) for k in (
        "peptide_chain", "receptor_chain", "chain_order", "source_to_model_chain",
        "fasta_source", "fasta_sha256", "protein_chains", "templates_policy", "protein_data_ready")}
    assert sha256(Path(manifest["fasta_source"])) == manifest["fasta_sha256"]
    manifest.update(
        ccd_sha256=sha256(args.ccd_json), prior_baseline=str(args.baseline.resolve()),
        prior_manifest_sha256=sha256(args.baseline / "manifest.json"),
        source_catalog_sha256=sha256(catalog_path), excluded_completed_states=sorted(excluded),
        protein_payload_sha256=protein_payload_hash,
        numbering_policy="Maximize modified-residue distance from numbering ends, then center at floor(n/2); keep source order if unmodified",
        states=[])
    shutil.copyfile(catalog_path, out / "source/chemical_catalog_audit.json")
    shutil.copyfile(args.package / "paper_2024_other_ligands/compounds.json", out / "source/compounds.json")
    for filename in ("TAS2R14_Gprotein.fasta", "prepared_protein_nochol.cif"):
        shutil.copyfile(args.baseline / "source" / filename, out / "source" / filename)
    shutil.copyfile(args.baseline / "deployment.json", out / "deployment.json")
    for row in rows:
        source = args.package / row["sdf"]
        assert sha256(source) == catalog["input_sha256"][row["sdf"]], source
        mol, state = encode(row, source, fragments, args.seeds)
        name = state["state"]
        write_sdf(out / "states" / (name + ".sdf"), [mol])
        peptide = {"protein": {"id": "A", "sequence": state["sequence"],
                   "modifications": state["modifications"], "unpairedMsa": "", "pairedMsa": "", "templates": []}}
        job = {"name": name, "dialect": "alphafold3", "version": 2, "modelSeeds": args.seeds,
               "sequences": [peptide, *proteins]}
        assert [s["protein"]["id"] for s in job["sequences"]] == manifest["chain_order"]
        path = out / "inputs" / (name + ".json")
        path.write_text(json.dumps(job, indent=2) + "\n")
        state["input_sha256"] = sha256(path)
        small = copy.deepcopy(job)
        for chain in small["sequences"]:
            chain["protein"].update(unpairedMsa="", pairedMsa="", templates=[])
        (out / "feature_only_inputs" / path.name).write_text(json.dumps(small, indent=2) + "\n")
        manifest["states"].append(state)
        print(name, state["sequence"], state["components"], state["formula"], "graph/stereo PASS", flush=True)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
