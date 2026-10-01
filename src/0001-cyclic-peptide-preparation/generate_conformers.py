"""Generate cyclic-peptide conformers with the selected production preset."""

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
    parser.add_argument("--states", nargs="+", help="State names from the preset; default: all")
    parser.add_argument("--attempts-per-state", type=int, default=500)
    parser.add_argument("--keep", type=int, default=40)
    parser.add_argument("--seed-batches", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--threads", type=int, default=2)
    add_profile(parser, "peptides.json")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from rdkit import Chem
    from cleprep.conformers import generate, validate_state
    from cleprep.common.io import write_json
    from cleprep.common.peptides import read_molecule, write_sdf

    profile = peptide_profile(args.profile)
    states = list(dict.fromkeys(args.states or profile["states"]))
    if any(state not in profile["states"] for state in states):
        parser.error("--states must come from the selected preset")
    if min(args.attempts_per_state, args.keep, args.seed_batches, args.threads) < 1 or args.seed < 0:
        parser.error("Counts must be positive and the seed must be nonnegative")
    require_new_output(parser, args.output_dir, "--output-dir")
    inputs: dict[str, Chem.Mol] = {}
    for state in states:
        spec = profile["states"][state]
        parent = profile["parents"][spec["parent"]]
        mol = read_molecule(args.input_dir / f"{state}.sdf")
        if not mol.HasProp("state_id") or mol.GetProp("state_id") != state:
            parser.error(f"{state}: SDF state_id does not match the filename")
        validate_state(mol, expected_residues=parent["residues"], expected_prolines=parent["prolines"],
                       expected_charge=parent["charge"], expected_sulfur_cip=spec["sulfur_cip"])
        families = 2 ** parent["prolines"]
        if args.keep % families or args.attempts_per_state < families * args.seed_batches:
            parser.error(f"{state}: --keep must be a multiple of {families}; attempts must cover every family and batch")
        inputs[state] = mol
    args.output_dir.mkdir(parents=True)
    incomplete = False
    for state, mol in inputs.items():
        state_seed = args.seed + profile["states"][state]["seed_offset"]
        valid, representatives, summary = generate(mol, state, args.attempts_per_state, args.keep,
                                                  state_seed, args.seed_batches, args.threads,
                                                  progress=lambda message: print(message, flush=True))
        folder = args.output_dir / state
        folder.mkdir()
        write_sdf(folder / "representatives.sdf", representatives)
        summary.update(seed=args.seed, state_seed=state_seed, seed_batches=args.seed_batches, threads=args.threads)
        write_json(folder / "summary.json", summary)
        incomplete |= bool(summary["missing_families"])
        print(f"{state}: accepted={len(valid)}, representatives={len(representatives)}/{args.keep}, "
              f"missing_families={summary['missing_families']}", flush=True)
    if incomplete:
        raise SystemExit("Some cis/trans families have no valid representative; see summary.json")


if __name__ == "__main__":
    main()
