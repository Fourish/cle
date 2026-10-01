"""Apply an externally supplied rigid transform to every atom in a PDB file."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from .common.geometry import Array, checked_transform


def orient(source: Path, output: Path, rotation: Array, translation: Array, inverse: bool = False) -> int:
    rotation, translation = checked_transform(rotation, translation)
    if inverse:
        translation = -rotation.T @ translation
        rotation = rotation.T
    lines: list[str] = []
    count = 0
    for line in source.read_text().splitlines(keepends=True):
        if line.startswith("ANISOU"):
            raise ValueError("Anisotropic displacement tensors are not supported")
        if line.startswith(("ATOM  ", "HETATM")):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=np.float64)
            xyz = rotation @ xyz + translation
            if not np.isfinite(xyz).all():
                raise ValueError("Nonfinite atom coordinates")
            coordinates = "".join(f"{value:8.3f}" for value in xyz)
            if len(coordinates) != 24:
                raise ValueError("Coordinates exceed PDB field width")
            line = line[:30] + coordinates + line[54:]
            count += 1
        lines.append(line)
    if not count:
        raise ValueError("No PDB atoms found")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        handle.write("".join(lines))
    return count
