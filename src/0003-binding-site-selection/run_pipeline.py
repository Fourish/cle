"""Batch prepared conformers through docking to geometry-filtered heavy complexes."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path, help="Workflow JSON; data paths are relative to this file")
    parser.add_argument("--output-dir", required=True, type=Path, help="New output directory; --resume reuses a matching run")
    parser.add_argument("--resume", action="store_true", help="Verify and resume this workflow's completed/incomplete tasks")
    parser.add_argument("--screen-only", action="store_true", help="Use matching existing_runs; perform no new docking")
    parser.add_argument("--limit-jobs", type=int, default=0, help="Run the first N planned tasks; 0 runs the whole plan")
    parser.add_argument("--plan-only", action="store_true", help="Check inputs and report task count without writing or searching")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    from production.docking_config import load_config
    from cleprep.workflow import prepare_jobs, run_workflow

    try:
        config = load_config(args.config)
        if args.plan_only:
            jobs, _ = prepare_jobs(config, args.limit_jobs)
            print(json.dumps({"job_count": len(jobs), "first_job": jobs[0].id,
                              "searches": sorted({job.search.id for job in jobs})}))
            return
        summary = run_workflow(config, args.output_dir, resume=args.resume,
                               screen_only=args.screen_only, limit_jobs=args.limit_jobs)
    except (ValueError, OSError, RuntimeError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if summary["failed_job_count"]:
        parser.exit(1, "Some docking tasks failed; see run_manifest.json and use --resume.\n")
    if summary["exported_model_count"] == 0:
        parser.exit(2, "No poses passed complex geometry screening; all QC results were saved.\n")


if __name__ == "__main__":
    main()
