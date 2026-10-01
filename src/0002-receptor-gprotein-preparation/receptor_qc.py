"""Screen a repaired receptor against its template, sequence and production preset."""

import argparse
import sys
import json
from pathlib import Path

# Allow direct execution from any working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from production.cli import add_profile, require_new_output
from production.profiles import chain_ranges, receptor_profile


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    add_profile(parser, "8xql.json")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from cleprep.common.io import read_fasta, write_json
    from cleprep.receptor_qc import inspect_structure, read_structure

    profile = receptor_profile(args.profile)
    require_new_output(parser, args.output, "--output")
    report = inspect_structure(read_structure(args.model), read_structure(args.template),
                               read_fasta(args.fasta), chain_ranges(profile))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, report)
    print(json.dumps(report["summary"]))
    if report["summary"]["has_blocking_geometry_issues"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
