"""Screen raw five-chain predictions and report the peptide–TAS2R14 interface."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path
import shutil

import gemmi
import numpy as np
from rdkit import Chem

from audit_outputs import audit
from cleprep.common.conformer_qc import QC_PROFILES
from cleprep.common.geometry import segment_triangle_intersections
from cleprep.common.peptides import read_molecule, write_sdf


def atom_table(path: Path) -> dict[str, np.ndarray]:
    block = gemmi.cif.read(str(path)).sole_block()
    names = ("label_asym_id", "label_seq_id", "label_comp_id", "label_atom_id",
             "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv")
    result = {key: np.asarray(list(block.find_values("_atom_site." + key))) for key in names}
    result["resid"] = result["label_seq_id"].astype(int)
    result["xyz"] = np.column_stack([result["Cartn_" + a].astype(float) for a in "xyz"])
    result["plddt"] = result["B_iso_or_equiv"].astype(float)
    return result


def fit_rmsd(mobile: np.ndarray, reference: np.ndarray) -> float:
    x = mobile - mobile.mean(axis=0)
    y = reference - reference.mean(axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    aligned = x @ u @ correction @ vt
    return float(np.sqrt(np.mean(np.sum((aligned - y) ** 2, axis=1))))


def backbone_threading(atoms: dict, manifest: dict) -> dict:
    """Reuse the repository's conservative ring-fan intersection diagnostic."""
    indices = {(str(c), int(r), str(a)): i for i, (c, r, a) in enumerate(zip(
        atoms["label_asym_id"], atoms["resid"], atoms["label_atom_id"], strict=True))}
    peptide = manifest["peptide_chain"]
    resids = sorted({r for c, r, a in indices if c == peptide})
    ring = atoms["xyz"][[indices[(peptide, r, a)] for r in resids for a in ("N", "CA", "C")]]
    center = ring.mean(axis=0)
    bonds = []
    for protein in manifest["protein_chains"]:
        chain = protein["id"]
        keys = [(chain, r, a) for r in range(1, protein["length"] + 1) for a in ("N", "CA", "C")]
        bonds.extend(zip(keys[:-1], keys[1:]))
    starts = atoms["xyz"][[indices[a] for a, _ in bonds]]
    ends = atoms["xyz"][[indices[b] for _, b in bonds]]
    radius = float(np.max(np.linalg.norm(ring - center, axis=1))) + 2.0
    selected = np.flatnonzero(np.sum(((starts + ends) * 0.5 - center) ** 2, axis=1) <= radius ** 2)
    hits, seen = [], []
    for i in range(len(ring)):
        for j, point in segment_triangle_intersections(
                starts[selected], ends[selected], center, ring[i], ring[(i + 1) % len(ring)]):
            if any(np.linalg.norm(point - previous) < 0.01 for previous in seen):
                continue
            seen.append(point)
            a, b = bonds[int(selected[j])]
            hits.append({"protein_bond": [list(a), list(b)], "intersection_xyz_A": point.tolist(),
                         "fan_triangle_index0": i})
    return {"potential_backbone_threading": bool(hits), "ring_fan_intersections": hits,
            "scope": "Protein backbone segments intersecting a nonplanar cyclic N-CA-C centroid fan; conservative flag, not a topological proof"}


def interface_metrics(sample: Path, manifest: dict, reference: dict | None) -> dict:
    atoms = atom_table(sample)
    conf = json.loads((sample.parent / "confidences.json").read_text())
    summary = json.loads((sample.parent / "summary_confidences.json").read_text())
    token_chain = np.asarray(conf["token_chain_ids"])
    token_resid = np.asarray(conf["token_res_ids"])
    # Deployed AF3 assigns asym IDs in chain appearance order; its pair matrix
    # uses sorted numeric asym IDs, not lexicographically sorted chain names.
    chain_order = list(dict.fromkeys(conf["token_chain_ids"]))
    if chain_order != manifest["chain_order"]:
        raise ValueError(f"Unexpected chain order: {chain_order}")
    peptide, receptor = manifest["peptide_chain"], manifest["receptor_chain"]
    pi, ri = chain_order.index(peptide), chain_order.index(receptor)
    pt, rt = np.flatnonzero(token_chain == peptide), np.flatnonzero(token_chain == receptor)
    pae = np.asarray(conf["pae"], dtype=float)
    pair_pae = np.concatenate((pae[np.ix_(pt, rt)].ravel(), pae[np.ix_(rt, pt)].ravel()))
    result = {
        "confidence_chain_order": chain_order,
        "receptor_peptide_iptm": summary["chain_pair_iptm"][pi][ri],
        "receptor_peptide_mean_pae": float(pair_pae.mean()),
        "receptor_peptide_min_pae": float(pair_pae.min()),
        "complex_iptm": summary["iptm"], "complex_ranking_score": summary["ranking_score"],
        "af3_has_clash": bool(summary["has_clash"]),
        "peptide_other_chain_iptm": {c: summary["chain_pair_iptm"][pi][i]
                                     for i, c in enumerate(chain_order) if c != peptide},
    }
    pm = (atoms["label_asym_id"] == peptide) & (atoms["type_symbol"] != "H")
    peptide_xyz = atoms["xyz"][pm]
    vdw = Chem.GetPeriodicTable()
    radii = np.array([vdw.GetRvdw(str(e).title()) for e in atoms["type_symbol"]])
    contacts, overlaps, chain_checks = {}, [], {}
    mobile_ca, reference_ca = [], []
    receptor_contact_resids = []
    for protein in manifest["protein_chains"]:
        chain = protein["id"]
        ca = (atoms["label_asym_id"] == chain) & (atoms["label_atom_id"] == "CA")
        sequence = "".join(gemmi.find_tabulated_residue(str(n)).one_letter_code
                           for n in atoms["label_comp_id"][ca])
        chain_checks[chain] = {
            "sequence_matches": sequence == protein["sequence"],
            "residue_numbers_match": np.array_equal(atoms["resid"][ca], np.arange(1, protein["length"] + 1)),
            "ca_count": int(ca.sum()), "mean_atom_plddt": float(atoms["plddt"][atoms["label_asym_id"] == chain].mean()),
        }
        cm = (atoms["label_asym_id"] == chain) & (atoms["type_symbol"] != "H")
        xyz = atoms["xyz"][cm]
        distances = np.linalg.norm(peptide_xyz[:, None, :] - xyz[None, :, :], axis=-1)
        close = distances <= 4.5
        resids = sorted(set(int(r) for r in atoms["resid"][cm][close.any(axis=0)]))
        labels = []
        for rid in resids:
            mask = atoms["resid"][cm] == rid
            labels.append({"query_resid": rid, "author_resid": rid + protein["author_residue_offset"],
                           "resname": str(atoms["label_comp_id"][cm][mask][0]),
                           "min_distance": float(distances[:, mask].min())})
        contacts[chain] = {"min_distance": float(distances.min()),
                           "atom_pairs_within_4_5": int(close.sum()),
                           "peptide_atoms_within_4_5": int(close.any(axis=1).sum()),
                           "residues": labels}
        if chain == receptor:
            receptor_contact_resids = resids
        ratios = distances / (radii[pm][:, None] + radii[cm][None, :])
        for i, j in zip(*np.where(ratios < 0.55), strict=True):
            overlaps.append({"chain": chain, "query_resid": int(atoms["resid"][cm][j]),
                             "atom": str(atoms["label_atom_id"][cm][j]),
                             "peptide_resid": int(atoms["resid"][pm][i]),
                             "peptide_atom": str(atoms["label_atom_id"][pm][i]),
                             "distance": float(distances[i, j]), "vdw_ratio": float(ratios[i, j])})
        if reference is not None:
            rm = (reference["label_asym_id"] == protein["source_id"]) & (reference["label_atom_id"] == "CA")
            if int(rm.sum()) == protein["length"] and int(ca.sum()) == protein["length"]:
                assert np.array_equal(reference["resid"][rm], np.arange(1, protein["length"] + 1) + protein["author_residue_offset"])
                chain_checks[chain]["ca_rmsd_to_prepared_reference"] = fit_rmsd(atoms["xyz"][ca], reference["xyz"][rm])
                mobile_ca.append(atoms["xyz"][ca])
                reference_ca.append(reference["xyz"][rm])
    if len(mobile_ca) == len(manifest["protein_chains"]):
        result["protein_assembly_ca_rmsd_to_prepared_reference"] = fit_rmsd(
            np.concatenate(mobile_ca), np.concatenate(reference_ca))
    contact_tokens = np.flatnonzero((token_chain == receptor) & np.isin(token_resid, receptor_contact_resids))
    result["receptor_contact_mean_pae"] = (float(np.mean(np.concatenate((
        pae[np.ix_(pt, contact_tokens)].ravel(), pae[np.ix_(contact_tokens, pt)].ravel()))))
        if len(contact_tokens) else None)
    result.update(chain_checks=chain_checks, peptide_protein_contacts=contacts,
                  peptide_protein_overlaps=overlaps,
                  threading=backbone_threading(atoms, manifest),
                  receptor_contact_residue_count=len(contacts[receptor]["residues"]),
                  nearest_protein_chain=min(contacts, key=lambda c: contacts[c]["min_distance"]),
                  protein_sequences_pass=all(x["sequence_matches"] and x["residue_numbers_match"] for x in chain_checks.values()))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--samples-per-seed", type=int, default=5)
    parser.add_argument("--qc-profile", choices=QC_PROFILES, default="strict",
                        help="Relaxed broadens three geometry tolerances; keeps chirality and clash checks.")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    selected = args.output / "selected"
    selected.mkdir(exist_ok=True)
    manifest = json.loads((args.baseline / "manifest.json").read_text())
    reference = atom_table(args.reference) if args.reference else None
    summaries, reports = [], []
    for state in manifest["states"]:
        name = state["state"]
        source = read_molecule(args.baseline / "states" / (name + ".sdf"))
        job = json.loads((args.baseline / "inputs" / (name + ".json")).read_text())
        samples = sorted((args.predictions / name.lower()).glob("seed-*_sample-*/model.cif"))
        local_reports, accepted_mols, rejected_mols, candidates = [], [], [], []
        for path in samples:
            report, mol = audit(path, state, source, chain_id=manifest["peptide_chain"],
                                qc_profile=args.qc_profile)
            report["peptide_geometry_accepted"] = report["accepted"]
            report.update(interface_metrics(path, manifest, reference))
            if not report["protein_sequences_pass"]:
                report["issues"].append("protein_sequence_or_residue_mismatch")
            if report["peptide_protein_overlaps"]:
                report["issues"].append("peptide_protein_overlap")
            if report["af3_has_clash"]:
                report["issues"].append("af3_clash_flag")
            if report["threading"]["potential_backbone_threading"]:
                report["issues"].append("potential_backbone_threading")
            if not report["peptide_protein_contacts"]["R"]["residues"]:
                report["issues"].append("no_receptor_contact_within_4.5A")
            report["accepted"] = not report["issues"]
            if mol is not None:
                mol.SetProp("highfold3_complex_geometry_pass", str(report["accepted"]))
                (accepted_mols if report["accepted"] else rejected_mols).append(mol)
                if report["accepted"]:
                    candidates.append((report, mol))
            local_reports.append(report)
        for suffix, molecules in (("accepted", accepted_mols), ("rejected", rejected_mols)):
            path = args.output / f"{name}__{suffix}.sdf"
            if molecules:
                write_sdf(path, molecules)
            else:
                path.unlink(missing_ok=True)
        best = None
        if candidates:
            best, mol = max(candidates, key=lambda pair: (
                pair[0]["receptor_peptide_iptm"],
                -pair[0]["receptor_contact_mean_pae"], pair[0]["mean_atom_plddt"]))
            shutil.copyfile(best["cif"], selected / (name + ".cif"))
            write_sdf(selected / (name + ".sdf"), [mol])
            (selected / (name + ".json")).write_text(json.dumps(best, indent=2) + "\n")
        else:
            for suffix in ("cif", "sdf", "json"):
                (selected / (name + "." + suffix)).unlink(missing_ok=True)
        item = {"state": name, "samples": len(samples),
                "expected_samples": len(job["modelSeeds"]) * args.samples_per_seed,
                "peptide_geometry_accepted": sum(r["peptide_geometry_accepted"] for r in local_reports),
                "complex_geometry_accepted": len(candidates),
                "selected_sample": best["sample"] if best else None,
                "selected_receptor_peptide_iptm": best["receptor_peptide_iptm"] if best else None,
                "selected_receptor_contact_mean_pae": best["receptor_contact_mean_pae"] if best else None}
        summaries.append(item)
        reports.extend(local_reports)
        print(json.dumps(item), flush=True)
    result = {"selection": "Pass peptide chemistry, protein sequence, no severe peptide/protein overlap, AF3 clash flag or potential backbone threading, and receptor contact <=4.5 A; rank by peptide/receptor pair ipTM, then contact PAE, then peptide pLDDT.",
              "reference": str(args.reference) if args.reference else None,
              "qc_profile": args.qc_profile, "qc_thresholds": asdict(QC_PROFILES[args.qc_profile]),
              "summary": summaries, "samples": reports}
    (args.output / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    fields = ["state", "sample", "accepted", "peptide_geometry_accepted", "mean_atom_plddt",
              "closure_cn_angstrom", "receptor_peptide_iptm", "receptor_contact_mean_pae",
              "receptor_peptide_mean_pae", "receptor_contact_residue_count", "nearest_protein_chain",
              "complex_iptm", "complex_ranking_score",
              "protein_assembly_ca_rmsd_to_prepared_reference", "issues"]
    with (args.output / "metrics.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(reports)
    if not all(s["samples"] == s["expected_samples"] and s["complex_geometry_accepted"] for s in summaries):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
