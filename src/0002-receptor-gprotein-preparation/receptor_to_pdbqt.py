"""Export a cholesterol-free rigid Vina receptor from a production preset."""

import argparse
import sys
import json
from pathlib import Path

# Allow direct execution from any working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from production.cli import add_profile
from production.profiles import chain_starts, receptor_profile


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--chemistry-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    add_profile(parser, "8xql.json")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from cleprep.pdbqt import export_receptor

    profile = receptor_profile(args.profile)
    expected = profile.get("expected_protein_heavy_atoms")
    summary = export_receptor(args.input, args.chemistry_dir, args.output,
                              keep_cholesterol=False,
                              chain_starts=chain_starts(profile), expected_heavy_atoms=expected)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
