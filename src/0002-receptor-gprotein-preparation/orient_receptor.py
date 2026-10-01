"""Apply a proper rigid transform to all atoms of a PDB file."""

import argparse
import sys
from pathlib import Path

# Allow direct execution from any working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--transform", required=True, type=Path)
    parser.add_argument("--inverse", action="store_true")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from cleprep.common.io import read_transform
    from cleprep.orientation import orient

    rotation, translation = read_transform(args.transform)
    count = orient(args.input, args.output, rotation, translation, args.inverse)
    print(f"Transformed {count} atoms: {args.output}")


if __name__ == "__main__":
    main()
