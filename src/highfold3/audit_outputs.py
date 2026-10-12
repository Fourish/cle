"""Audit every raw HighFold3 sample against the source cyclic peptide graph."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import sys

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cleprep.common.conformer_qc import QC_PROFILES, audit_conformer
from cleprep.common.peptides import read_molecule, stereo_labels, write_sdf


def audit(sample, state, source, chain_id="A", *, qc_profile="strict"):
    thresholds = QC_PROFILES[qc_profile]
    block = gemmi.cif.read(str(sample)).sole_block()
    columns = {k: list(block.find_values("_atom_site." + k)) for k in (
        "label_asym_id", "label_seq_id", "label_atom_id", "label_comp_id",
        "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv")}
    atoms = {}
    confidences = []
    for i, chain in enumerate(columns["label_asym_id"]):
        if chain != chain_id:
            continue
        key = (int(columns["label_seq_id"][i]), columns["label_atom_id"][i])
        if key in atoms:
            raise ValueError(f"Duplicate atom in {sample}: {key}")
        atoms[key] = (columns["label_comp_id"][i], np.array(
            [float(columns["Cartn_" + a][i]) for a in "xyz"]))
        confidences.append(float(columns["B_iso_or_equiv"][i]))
    mapping = state["atom_mapping"]
    expected = {(a["residue_id"], a["atom_name"]): a for a in mapping}
    missing = sorted(set(expected) - set(atoms))
    extra = sorted(set(atoms) - set(expected))
    wrong_components = [key for key in set(expected) & set(atoms)
                        if expected[key]["component"] != atoms[key][0]]
    report = {"state": state["state"], "sample": sample.parent.name, "cif": str(sample),
              "expected_heavy_atoms": state["heavy_atoms"], "observed_heavy_atoms": len(atoms),
              "missing_atoms": missing, "extra_atoms": extra, "wrong_components": wrong_components,
              "mean_atom_plddt": float(np.mean(confidences)),
              "qc_profile": qc_profile, "qc_thresholds": asdict(thresholds)}
    if missing:
        report.update(accepted=False, issues=["missing_atoms"])
        return report, None
    mol = Chem.Mol(source)
    mol.RemoveAllConformers()
    conformer = Chem.Conformer(mol.GetNumAtoms())
    conformer.Set3D(True)
    for key, a in expected.items():
        conformer.SetAtomPosition(a["source_heavy_atom"], atoms[key][1])
    mol.AddConformer(conformer)
    expected_cip = {int(k): v for k, v in state["cip"].items()}
    issues = audit_conformer(mol, 0, expected_cip, state["backbone"], state["amide_omegas"],
                             thresholds=thresholds)
    if extra:
        issues.append("extra_atoms")
    if wrong_components:
        issues.append("wrong_components")
    measured = Chem.Mol(mol)
    Chem.RemoveStereochemistry(measured)
    Chem.AssignStereochemistryFrom3D(measured, replaceExistingTags=True)
    actual = stereo_labels(measured)
    report["chirality_mismatches"] = {str(i): {"expected": v, "actual": actual.get(i)}
                                      for i, v in expected_cip.items() if actual.get(i) != v}
    nres = len(state["components"])
    connections = set()
    connection_fields = [list(block.find_values("_struct_conn." + field)) for field in (
        "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_label_atom_id",
        "ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_label_atom_id")]
    for c1, r1, a1, c2, r2, a2 in zip(*connection_fields, strict=True):
        if c1 == c2 == chain_id:
            connections.add(frozenset(((int(r1), a1), (int(r2), a2))))
    peptide_bonds = {frozenset(((r, "C"), (r % nres + 1, "N")))
                     for r in range(1, nres + 1)}
    report["all_peptide_bonds_declared_in_cif"] = peptide_bonds <= connections
    if not report["all_peptide_bonds_declared_in_cif"]:
        issues.append("missing_peptide_bond_metadata")
    distances = [float(np.linalg.norm(atoms[(r, "C")][1] - atoms[(r % nres + 1, "N")][1]))
                 for r in range(1, nres + 1)]
    report["peptide_cn_angstrom"] = distances
    report["closure_cn_angstrom"] = distances[-1]
    report["all_peptide_cn_pass"] = all(
        thresholds.peptide_cn_min_angstrom <= d <= thresholds.peptide_cn_max_angstrom
        for d in distances)
    report["amide_omega_degrees"] = [float(rdMolTransforms.GetDihedralDeg(conformer, *q))
                                     for q in state["amide_omegas"]]
    report["amide_plane_deviations_degrees"] = [min(abs(d), abs(180 - abs(d)))
                                                for d in report["amide_omega_degrees"]]
    report["proline_omega_degrees"] = [float(rdMolTransforms.GetDihedralDeg(conformer, *q))
                                         for q in state["proline_omegas"]]
    summary_path = sample.parent / "summary_confidences.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        report["ranking_score"] = summary.get("ranking_score")
        report["ptm"] = summary.get("ptm")
    if not report["all_peptide_cn_pass"]:
        issues.append("peptide_cn_out_of_range")
    report.update(accepted=not issues, issues=issues)
    mol.SetProp("_Name", state["state"] + "__" + sample.parent.name)
    mol.SetProp("highfold3_raw_geometry_pass", str(report["accepted"]))
    mol.SetProp("highfold3_qc_profile", qc_profile)
    mol.SetProp("source_cif", str(sample))
    return report, mol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples-per-seed", type=int, default=5)
    parser.add_argument("--qc-profile", choices=QC_PROFILES, default="strict",
                        help="Relaxed broadens three geometry tolerances; keeps chirality and clash checks.")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.baseline / "manifest.json").read_text())
    reports, summaries = [], []
    for state in manifest["states"]:
        name = state["state"]
        job = json.loads((args.baseline / "inputs" / (name + ".json")).read_text())
        expected_samples = len(job["modelSeeds"]) * args.samples_per_seed
        source = read_molecule(args.baseline / "states" / (name + ".sdf"))
        directory = args.predictions / name.lower()
        samples = sorted(directory.glob("seed-*_sample-*/model.cif"))
        accepted, rejected, state_reports = [], [], []
        for sample in samples:
            report, mol = audit(sample, state, source, qc_profile=args.qc_profile)
            reports.append(report)
            state_reports.append(report)
            if mol is not None:
                (accepted if report["accepted"] else rejected).append(mol)
        if accepted:
            write_sdf(args.output / (name + "__accepted.sdf"), accepted)
            passed_reports = [r for r in state_reports if r["accepted"]]
            best_idx = max(range(len(passed_reports)), key=lambda i: (
                passed_reports[i].get("ranking_score") or 0.0,
                passed_reports[i]["mean_atom_plddt"]))
            best = passed_reports[best_idx]
            selected_dir = args.output / "selected"
            selected_dir.mkdir(exist_ok=True)
            shutil.copyfile(best["cif"], selected_dir / (name + ".cif"))
            write_sdf(selected_dir / (name + ".sdf"), [accepted[best_idx]])
            (selected_dir / (name + ".json")).write_text(json.dumps(best, indent=2) + "\n")
        else:
            (args.output / (name + "__accepted.sdf")).unlink(missing_ok=True)
            for suffix in (".cif", ".sdf", ".json"):
                (args.output / "selected" / (name + suffix)).unlink(missing_ok=True)
        if rejected:
            write_sdf(args.output / (name + "__rejected.sdf"), rejected)
        else:
            (args.output / (name + "__rejected.sdf")).unlink(missing_ok=True)
        summary = {"state": name, "samples": len(samples), "expected_samples": expected_samples,
                   "accepted": len(accepted),
                   "selected_sample": best["sample"] if accepted else None,
                   "closure_cn_angstrom": [r.get("closure_cn_angstrom") for r in state_reports],
                   "samples_with_chirality_mismatch": sum(bool(r.get("chirality_mismatches")) for r in state_reports)}
        summaries.append(summary)
        print(json.dumps(summary))
    (args.output / "audit.json").write_text(json.dumps({"qc_profile": args.qc_profile,
        "qc_thresholds": asdict(QC_PROFILES[args.qc_profile]),
        "summary": summaries, "samples": reports}, indent=2) + "\n")
    if not all(s["samples"] == s["expected_samples"] and s["accepted"] >= 1 for s in summaries):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
