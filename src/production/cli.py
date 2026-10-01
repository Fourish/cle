"""Small command-line helpers shared by the production steps."""

from __future__ import annotations

import argparse
from pathlib import Path


def add_profile(parser: argparse.ArgumentParser, filename: str) -> None:
    """Resolve bundled presets relative to src, independent of the launch directory."""
    parser.add_argument(
        "--profile", type=Path,
        default=Path(__file__).resolve().parents[1] / "profiles" / filename,
        help=f"Production preset JSON (default: profiles/{filename})",
    )


def require_new_output(parser: argparse.ArgumentParser, path: Path, option: str) -> None:
    """Reject an existing output before the entrypoint starts processing inputs."""
    if path.exists():
        label = "Output directory" if option == "--output-dir" else "Output"
        parser.error(f"{label} exists; select a new {option}")
