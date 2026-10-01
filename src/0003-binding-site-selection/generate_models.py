"""Run complete steps 1--3 from 2D peptides and receptor mmCIF to heavy complexes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path, help="Full pipeline profile; paths relative to this JSON")
    parser.add_argument("--output-dir", required=True, type=Path, help="Output directory for steps 1--3")
    parser.add_argument("--resume", action="store_true", help="Verify and reuse completed preparation and docking stages")
    parser.add_argument("--plan-only", action="store_true", help="Validate source inputs and dependencies without computation")
    parser.add_argument("--limit-jobs", type=int, default=0, help="Limit docking tasks only; 0 runs the whole plan")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    from production.pipeline import check_inputs, load_pipeline_config, run_pipeline

    try:
        config = load_pipeline_config(args.config)
        if args.plan_only:
            print(json.dumps(check_inputs(config, args.limit_jobs), ensure_ascii=False, indent=2))
            return
        summary = run_pipeline(config, args.output_dir, resume=args.resume, limit_jobs=args.limit_jobs)
    except KeyboardInterrupt:
        parser.exit(130, "Interrupted; use --resume with the same configuration and output directory.\n")
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if summary["failed_job_count"]:
        parser.exit(1, "Some docking tasks failed; see run_manifest.json and use --resume.\n")
    if not summary["exported_model_count"]:
        parser.exit(2, "No poses passed geometry screening; all rejection reports were saved.\n")


if __name__ == "__main__":
    main()
