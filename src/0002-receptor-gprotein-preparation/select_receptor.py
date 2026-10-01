"""Select a receptor construct from mmCIF using a production preset."""

import argparse
import sys
import json
from pathlib import Path

# Allow direct execution from any working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from production.cli import add_profile
from production.profiles import receptor_profile


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    add_profile(parser, "8xql.json")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from cleprep.selection import select

    profile = receptor_profile(args.profile)
    targets = {chain: (spec["label"], spec["first"], spec["last"]) for chain, spec in profile["chains"].items()}
    ligands = {label: (spec["name"], spec["heavy_atoms"]) for label, spec in profile.get("ligands", {}).items()}
    print(json.dumps(select(args.input, args.output_dir, targets, ligands, profile.get("source_sha256"))))


if __name__ == "__main__":
    main()
