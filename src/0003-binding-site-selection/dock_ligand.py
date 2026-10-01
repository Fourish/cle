"""Dock one cyclic-peptide SDF conformer into a prepared rigid Vina receptor."""

import argparse
import sys
from pathlib import Path

# Allow direct execution from any working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receptor", required=True, type=Path, help="Prepared cholesterol-free receptor PDBQT")
    parser.add_argument("--ligand", required=True, type=Path, help="Explicit-H 3D cyclic-peptide SDF")
    parser.add_argument("--output-dir", required=True, type=Path, help="New output directory")
    parser.add_argument("--center", required=True, type=float, nargs=3, metavar=("X", "Y", "Z"), help="Box center in Angstrom")
    parser.add_argument("--size", required=True, type=float, nargs=3, metavar=("LX", "LY", "LZ"), help="Full box side lengths in Angstrom")
    parser.add_argument("--record", type=int, default=0, help="Zero-based SDF record (default: 0)")
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--cpu", type=int, default=1)
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--max-evals", type=int, default=100000, help="Per-trajectory search budget; 0 uses Vina's automatic budget")
    parser.add_argument("--n-poses", type=int, default=10)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from cleprep.docking import dock_ligand

    try:
        summary = dock_ligand(args.receptor, args.ligand, args.output_dir, args.center, args.size,
                              record=args.record, seed=args.seed, cpu=args.cpu,
                              exhaustiveness=args.exhaustiveness, max_evals=args.max_evals, n_poses=args.n_poses)
    except (ValueError, OSError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(f"poses={summary['pose_count']}, ligand_geometry_passed={summary['ligand_geometry_passed_count']}, "
          f"output={args.output_dir}", flush=True)
    if summary["ligand_geometry_passed_count"] == 0:
        raise SystemExit("No poses passed ligand geometry checks; see summary.json")


if __name__ == "__main__":
    main()
