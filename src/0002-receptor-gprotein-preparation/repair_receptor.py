"""Run restricted MODELLER repair with chain numbering from a production preset."""

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
    parser.add_argument("--seed", type=int, default=-20260910)
    parser.add_argument("--models", type=int, default=1)
    parser.add_argument("--loops", type=int, default=4)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    add_profile(parser, "8xql.json")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from cleprep.repair import repair

    profile = receptor_profile(args.profile)
    results = repair(args.template, args.fasta, args.output_dir, chain_starts(profile), args.seed, args.models, args.loops)
    print(json.dumps({"output_dir": str(args.output_dir), "successful": sum(r["failed"] is None for r in results)}))


if __name__ == "__main__":
    main()
