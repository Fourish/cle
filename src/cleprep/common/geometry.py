"""Rigid transforms and geometric measurements using float64 coordinates."""

from __future__ import annotations

import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer
from numpy.typing import ArrayLike, NDArray
from typing import TypeAlias

Array: TypeAlias = NDArray[np.float64]


def checked_transform(rotation: ArrayLike, translation: ArrayLike) -> tuple[Array, Array]:
    rotation = np.asarray(rotation, dtype=np.float64)
    translation = np.asarray(translation, dtype=np.float64)
    if (rotation.shape != (3, 3) or translation.shape != (3,)
            or not np.isfinite(rotation).all() or not np.isfinite(translation).all()
            or not np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-7)
            or not np.isclose(np.linalg.det(rotation), 1.0)):
        raise ValueError("Expected a finite proper rigid transform in Angstrom")
    return rotation, translation


def rigid_fit(mobile: Array, target: Array) -> tuple[Array, Array]:
    """Return proper row-vector rotation and translation; never reflect."""
    fitter = SVDSuperimposer()
    fitter.set(np.asarray(target, dtype=np.float64), np.asarray(mobile, dtype=np.float64))
    fitter.run()
    return fitter.get_rotran()


def dihedral(a: Array, b: Array, c: Array, d: Array) -> float | None:
    axis = c - b
    if np.linalg.norm(axis) < 1e-10:
        return None
    axis /= np.linalg.norm(axis)
    v = a - b - np.dot(a - b, axis) * axis
    w = d - c - np.dot(d - c, axis) * axis
    if min(np.linalg.norm(v), np.linalg.norm(w)) < 1e-10:
        return None
    return float(np.degrees(np.arctan2(np.dot(np.cross(axis, v), w), np.dot(v, w))))


def segment_distance(a: Array, b: Array, c: Array, d: Array) -> tuple[float, float, float]:
    """Return minimum finite-segment distance and closest-point fractions."""
    u, v, w = b - a, d - c, a - c
    uu, uv, vv, uw, vw = (float(np.dot(x, y)) for x, y in
                          ((u, u), (u, v), (v, v), (u, w), (v, w)))
    if uu <= 1e-12 or vv <= 1e-12:
        return float("inf"), 0.0, 0.0
    candidates = [(s, float(np.clip((vw + s * uv) / vv, 0.0, 1.0))) for s in (0.0, 1.0)]
    candidates += [(float(np.clip((t * uv - uw) / uu, 0.0, 1.0)), t) for t in (0.0, 1.0)]
    determinant = uu * vv - uv * uv
    if determinant > 1e-12:
        s, t = (uv * vw - vv * uw) / determinant, (uu * vw - uv * uw) / determinant
        if 0.0 <= s <= 1.0 and 0.0 <= t <= 1.0:
            candidates.append((s, t))
    return min((float(np.linalg.norm(w + s * u - t * v)), s, t) for s, t in candidates)


def segment_triangle_intersections(starts: Array, ends: Array, first: Array,
                                   second: Array, third: Array) -> list[tuple[int, Array]]:
    """Intersect finite segments with a triangle, excluding segment endpoints.

    These are the original Möller–Trumbore tolerances used by the ring-fan
    diagnostic. Coplanar segments have no unique intersection and are ignored.
    """
    directions = ends - starts
    edge1, edge2 = second - first, third - first
    h = np.cross(directions, edge2)
    determinant = h @ edge1
    valid = np.abs(determinant) > 1e-10
    inverse = np.zeros(len(starts), dtype=np.float64)
    inverse[valid] = 1.0 / determinant[valid]
    s = starts - first
    u = inverse * np.einsum("ij,ij->i", s, h)
    q = np.cross(s, edge1)
    v = inverse * np.einsum("ij,ij->i", directions, q)
    t = inverse * (q @ edge2)
    hit = valid & (u >= -1e-7) & (v >= -1e-7) & (u + v <= 1 + 1e-7) & (t > 1e-6) & (t < 1 - 1e-6)
    return [(int(index), starts[index] + t[index] * directions[index]) for index in np.flatnonzero(hit)]
