"""Refine unresolved tails using hinge regions from a production preset."""

import argparse
import sys
import json
from pathlib import Path
from typing import cast

# Allow direct execution from any working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from production.cli import add_profile
from production.profiles import chain_ranges, receptor_profile


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--transform", required=True, type=Path)
    add_profile(parser, "8xql.json")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from cleprep.common.io import read_fasta, read_transform
    from cleprep.receptor_qc import QCSummary
    from cleprep.tails import TailRegion, refine_tails

    profile = receptor_profile(args.profile)
    regions = [TailRegion(spec["chain"], spec["first"], spec["last"], spec["critical_start"],
                          spec.get("entry_phi_count", 1),
                          tuple((offset, kind) for offset, kind in spec.get("second_hinges", [(0, "psi")])),
                          tuple((offset, kind) for offset, kind in spec.get("third_hinges", [])))
               for spec in profile.get("tails", [])]
    rotation, translation = read_transform(args.transform)
    half_width = float(json.loads(args.transform.read_text())["pdbtm_membrane_half_width_angstrom"])
    report = refine_tails(args.input, args.output, args.template, read_fasta(args.fasta),
                          chain_ranges(profile), regions, rotation, translation, half_width)
    summary = cast(QCSummary, report["summary"])
    print(json.dumps(summary))
    if summary["has_blocking_geometry_issues"]:
        raise SystemExit("Tail model failed geometry screening; see the JSON report")


if __name__ == "__main__":
    main()
