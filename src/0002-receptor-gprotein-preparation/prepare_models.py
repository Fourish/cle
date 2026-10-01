"""Prepare a cholesterol-free receptor ensemble from mmCIF through PDBQT."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from production.cli import add_profile, require_new_output
from production.profiles import receptor_profile


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--forcefield", required=True, type=Path)
    parser.add_argument("--transform", required=True, type=Path)
    parser.add_argument("--repair-seeds", nargs="+", type=int, default=[-20260910, -20260911, -20260912])
    parser.add_argument("--models", type=int, default=1)
    parser.add_argument("--loops", type=int, default=4)
    parser.add_argument("--keep", type=int, default=3)
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--cpu-threads", type=int, default=2)
    add_profile(parser, "8xql.json")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if (any(seed >= 0 for seed in args.repair_seeds) or len(set(args.repair_seeds)) != len(args.repair_seeds)
            or min(args.models, args.loops, args.keep, args.cpu_threads) < 1 or args.seed < 0
            or not 0 <= args.ph <= 14):
        parser.error("Require unique negative repair seeds, positive counts, seed >=0 and pH in [0,14]")
    require_new_output(parser, args.output_dir, "--output-dir")
    from production.receptor_models import prepare_models

    try:
        models = prepare_models(args.input.resolve(), args.output_dir, args.forcefield.resolve(),
                                args.transform.resolve(), receptor_profile(args.profile),
                                repair_seeds=tuple(args.repair_seeds), models=args.models,
                                loops=args.loops, keep=args.keep, ph=args.ph, seed=args.seed,
                                cpu_threads=args.cpu_threads)
    except (ValueError, OSError, RuntimeError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(json.dumps({"prepared_model_count": len(models), "output_dir": str(args.output_dir.resolve())}))


if __name__ == "__main__":
    main()
