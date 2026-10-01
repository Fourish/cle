"""Prepare receptor microstates and explicit H using a production preset."""

import argparse
import sys
import json
from pathlib import Path

# Allow direct execution from any working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from production.cli import add_profile
from production.profiles import chain_ranges, receptor_profile


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--source-cif", type=Path, help="Required when cholesterol is present")
    parser.add_argument("--forcefield", required=True, type=Path)
    parser.add_argument("--toppar", type=Path, help="CHARMM cholesterol stream archive; required for cholesterol")
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--cpu-threads", type=int, default=2)
    add_profile(parser, "8xql.json")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from cleprep.chemistry import prepare

    profile = receptor_profile(args.profile)
    chemistry = profile.get("chemistry", {})
    summary = prepare(args.input, args.output_dir, args.forcefield, chain_ranges(profile),
                      source_cif=args.source_cif, toppar=args.toppar,
                      context_residues=frozenset(chemistry.get("context_residues", [])),
                      expected_context_atoms=chemistry.get("expected_context_atoms"),
                      expected_cholesterol_count=chemistry.get("expected_cholesterol_count"),
                      cholesterol_chain=chemistry.get("cholesterol_chain", "L"),
                      cholesterol_residue_id=chemistry.get("cholesterol_residue_id", "1"),
                      ph=args.ph, seed=args.seed, cpu_threads=args.cpu_threads,
                      progress=lambda message: print(message, flush=True))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
