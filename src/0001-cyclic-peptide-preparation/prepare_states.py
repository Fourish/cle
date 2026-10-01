"""Prepare cyclic-peptide chemical states using an explicit production preset."""

import argparse
import sys
from pathlib import Path

# Allow direct execution from any working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


from production.cli import add_profile, require_new_output
from production.profiles import peptide_profile


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    add_profile(parser, "peptides.json")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from cleprep.common.io import write_json
    from cleprep.common.peptides import stereo_labels, write_sdf
    from cleprep.states import prepare_parent

    profile = peptide_profile(args.profile)
    require_new_output(parser, args.output_dir, "--output-dir")
    states = [mol for name, spec in profile["parents"].items()
              for mol in prepare_parent(args.input_dir / f"{name}.sdf", expected_residues=spec["residues"],
                                        expected_prolines=spec["prolines"], expected_sulfoxides=spec["sulfoxides"],
                                        expected_charge=spec["charge"])]
    if {mol.GetProp("state_id") for mol in states} != set(profile["states"]):
        parser.error("Generated state names differ from the production preset")
    args.output_dir.mkdir(parents=True)
    summary: list[dict[str, object]] = []
    for mol in states:
        state = mol.GetProp("state_id")
        write_sdf(args.output_dir / f"{state}.sdf", [mol])
        summary.append({"state_id": state, "formula": rdMolDescriptors.CalcMolFormula(mol),
                        "atoms": mol.GetNumAtoms(), "heavy_atoms": mol.GetNumHeavyAtoms(),
                        "formal_charge": Chem.GetFormalCharge(mol), "cip": stereo_labels(mol)})
        print(f"{state}: saved 2D state with {mol.GetNumAtoms()} atoms including H", flush=True)
    write_json(args.output_dir / "summary.json", summary)


if __name__ == "__main__":
    main()
