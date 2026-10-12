"""Reuse validated cyclic peptide chemistry with the current four protein chains."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil

from Bio import SeqIO
import gemmi
from rdkit import Chem

from build_inputs import implicit_hydrogens
from cleprep.common.peptides import read_molecule


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def msa_rows(text: str) -> int:
    return text.count("\n>") + int(text.startswith(">"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--monomers", type=Path, required=True)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protein-data", type=Path,
                        help="Downloaded per-chain AF3 data pipeline output directory")
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    for folder in ("source", "states", "protein_queries", "inputs"):
        (output / folder).mkdir(exist_ok=True)
    fasta = args.package / "receptor/TAS2R14_Gprotein.fasta"
    records = list(SeqIO.parse(fasta, "fasta"))
    assert [(r.id, len(r)) for r in records] == [("R", 316), ("A", 249), ("B", 340), ("C", 71)]
    shutil.copyfile(fasta, output / "source/TAS2R14_Gprotein.fasta")
    shutil.copyfile(args.package / "manifest.json", output / "source/package_manifest.json")
    manifest = copy.deepcopy(json.loads((args.monomers / "manifest.json").read_text()))
    manifest.update(
        monomer_baseline=str(args.monomers.resolve()),
        monomer_manifest_sha256=sha256(args.monomers / "manifest.json"),
        peptide_chain="A", receptor_chain="R", chain_order=["A", "B", "C", "G", "R"],
        source_to_model_chain={"R": "R", "A": "G", "B": "B", "C": "C"},
        fasta_source=str(fasta), fasta_sha256=sha256(fasta),
        protein_chains=[],
        templates_policy="Automatic AF3 search; no docked pose or custom structure template",
    )
    proteins = []
    for record in records:
        name = "tas2r14_" + record.id.lower()
        protein = {"id": record.id, "sequence": str(record.seq)}
        query = {"name": name, "dialect": "alphafold3", "version": 2,
                 "modelSeeds": [1], "sequences": [{"protein": protein}]}
        (output / "protein_queries" / (name + ".json")).write_text(json.dumps(query, indent=2) + "\n")
        model_id = manifest["source_to_model_chain"][record.id]
        info = {"id": model_id, "source_id": record.id, "length": len(record), "sequence": str(record.seq),
                "description": record.description,
                "author_residue_offset": 1 if record.id == "R" else 0}
        if args.protein_data:
            path = args.protein_data / name / (name + "_data.json")
            cached = json.loads(path.read_text())
            protein = cached["sequences"][0]["protein"]
            assert protein["id"] == record.id and protein["sequence"] == str(record.seq)
            for key in ("unpairedMsa", "pairedMsa", "templates"):
                assert protein[key] is not None
            info.update(data_file=str(path), data_sha256=sha256(path),
                        unpaired_msa_rows=msa_rows(protein["unpairedMsa"]),
                        paired_msa_rows=msa_rows(protein["pairedMsa"]),
                        templates=len(protein["templates"]),
                        template_alignments=[{
                            "pdb": gemmi.cif.read_string(t["mmcif"]).sole_block().name,
                            "aligned_query_residues": len(set(t["queryIndices"])),
                            "query_coverage": len(set(t["queryIndices"])) / len(record),
                        } for t in protein["templates"]])
        protein = copy.deepcopy(protein)
        protein["id"] = model_id
        manifest["protein_chains"].append(info)
        proteins.append({"protein": protein})
    for state in manifest["states"]:
        name = state["state"]
        package_sdf = args.package / "ligands" / (name + ".sdf")
        reference = implicit_hydrogens(Chem.RemoveHs(read_molecule(package_sdf)))
        assert Chem.MolToSmiles(reference) == state["smiles"], name
        state.update(package_sdf=str(package_sdf), package_sdf_sha256=sha256(package_sdf),
                     package_graph_and_stereo_match=True, peptide_chain="A")
        shutil.copyfile(args.monomers / "states" / (name + ".sdf"), output / "states" / (name + ".sdf"))
        if args.protein_data:
            job = json.loads((args.monomers / "inputs" / (name + ".json")).read_text())
            job["sequences"][0]["protein"]["id"] = "A"
            job["sequences"].extend(copy.deepcopy(proteins))
            job["sequences"].sort(key=lambda item: item["protein"]["id"])
            assert [s["protein"]["id"] for s in job["sequences"]] == manifest["chain_order"]
            (output / "inputs" / (name + ".json")).write_text(json.dumps(job, indent=2) + "\n")
            print(name, "seeds", job["modelSeeds"], "chains", manifest["chain_order"])
    manifest["protein_data_ready"] = args.protein_data is not None
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Protein queries prepared; all five chemical graphs match the referenced package.")


if __name__ == "__main__":
    main()
