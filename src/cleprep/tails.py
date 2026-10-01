"""Redirect specified unresolved tails outside the membrane using phi/psi hinges.

Use the original two-hinge 10-degree grid and three-hinge 20-degree fallback.
This geometric construction leaves experimental atoms and peptide omega fixed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from .common.geometry import Array, checked_transform, dihedral
from .common.io import Json, write_json
from .receptor_qc import StructureData, build_bond_graph, graph_neighbors_within_three, inspect_structure, read_structure

Indices: TypeAlias = NDArray[np.int64]
BACKBONE_NAMES: set[str] = {"N", "CA", "C", "O", "OXT"}
CHARGED_NAMES: set[str] = {"ARG", "LYS", "ASP", "GLU", "HIS", "HSD", "HSE", "HSP"}


@dataclass(frozen=True)
class TailRegion:
    chain: str
    first: int
    last: int
    critical_start: int
    entry_phi_count: int = 1
    second_hinges: tuple[tuple[int, str], ...] = ((0, "psi"),)
    third_hinges: tuple[tuple[int, str], ...] = ()

    def __post_init__(self) -> None:
        if not self.first <= self.critical_start <= self.last or not 1 <= self.entry_phi_count <= self.last - self.first + 1:
            raise ValueError("Invalid tail boundaries or entry phi count")
        if not self.second_hinges or any(kind not in {"phi", "psi"} or not 0 <= offset <= self.last - self.first
                                        for offset, kind in self.second_hinges + self.third_hinges):
            raise ValueError("Tail hinges must be phi/psi at offsets within the tail")


def rotate_points(coordinates: Array, indices: Indices, start: Array, end: Array, degrees: float) -> Array:
    axis = end - start
    axis /= np.linalg.norm(axis)
    rotation = Rotation.from_rotvec(axis * degrees, degrees=True)
    updated = coordinates.copy()
    updated[indices] = rotation.apply(coordinates[indices] - start) + start
    return updated


def refine_tail(data: StructureData, template: StructureData, coordinates: Array,
                region: TailRegion, membrane_rotation: Array,
                membrane_translation: Array, intracellular_bound: float) -> tuple[Array, dict[str, Json]]:
    chain, first, critical_start = region.chain, region.first, region.critical_start
    tail = np.array([a.index for a in data.atoms if a.chain == chain and first <= a.number <= region.last], dtype=np.int64)
    if not len(tail):
        raise ValueError(f"Tail not present: {chain}:{first}-{region.last}")
    tail_set = set(int(i) for i in tail)
    critical = np.array([a.index for a in data.atoms if a.index in tail_set and
                        ((a.number >= critical_start and a.name in BACKBONE_NAMES) or a.residue in CHARGED_NAMES)], dtype=np.int64)
    graph = build_bond_graph(data, template)
    excluded = {int(i): graph_neighbors_within_three(graph, int(i)) for i in tail}
    static = [a.index for a in data.atoms if a.index not in tail_set]
    tree = cKDTree(coordinates[static])
    internal_pairs = np.array([(int(i), int(j)) for k, i in enumerate(tail) for j in tail[k + 1:]
                               if int(j) not in excluded[int(i)]], dtype=np.int64).reshape(-1, 2)
    short_pairs = np.array(sorted({tuple(sorted((int(i), j))) for i in tail for j in excluded[int(i)]
                                  if j != int(i) and j not in graph[int(i)]}), dtype=np.int64).reshape(-1, 2)

    def hinge(number: int, kind: str) -> tuple[Indices, int, int]:
        atoms = data.residues[(chain, number, "")]
        names = ("N", "CA") if kind == "phi" else ("CA", "C")
        moving = np.array([a.index for a in data.atoms if a.index in tail_set and
                           (a.number > number or (a.number == number and
                            (a.name not in {"N", "CA"} if kind == "phi" else a.name in {"O", "OXT"})))], dtype=np.int64)
        return moving, atoms[names[0]].index, atoms[names[1]].index

    def rotate(xyz: Array, axis: tuple[Indices, int, int], angle: int) -> Array:
        moving, a, b = axis
        return rotate_points(xyz, moving, xyz[a], xyz[b], float(angle))

    def phi(xyz: Array, number: int) -> float | None:
        previous = data.residues[(chain, number - 1, "")]
        current = data.residues[(chain, number, "")]
        return dihedral(xyz[previous["C"].index], xyz[current["N"].index],
                        xyz[current["CA"].index], xyz[current["C"].index])

    def entry_allowed(xyz: Array) -> bool:
        for number in range(first, first + region.entry_phi_count):
            value = phi(xyz, number)
            if value is None or not -170 <= value <= -30:
                return False
        return True

    def sterics(xyz: Array) -> tuple[int, float]:
        lengths = np.linalg.norm(xyz[internal_pairs[:, 0]] - xyz[internal_pairs[:, 1]], axis=1)
        clashes = int(np.sum(lengths < 2.0))
        penalty = float(np.sum(np.maximum(0.0, 2.8 - lengths) ** 2))
        for i, position in zip(tail, xyz[tail], strict=True):
            for local in tree.query_ball_point(position, 3.0):
                j = static[local]
                if j in excluded[int(i)]:
                    continue
                length = float(np.linalg.norm(position - xyz[j]))
                clashes += int(length < 2.0)
                penalty += max(0.0, 2.8 - length) ** 2
        short_lengths = np.linalg.norm(xyz[short_pairs[:, 0]] - xyz[short_pairs[:, 1]], axis=1)
        return clashes + int(np.sum(short_lengths < 1.8)), penalty

    sampled = 0
    minimum_excess = float("inf")
    best_score: tuple[int, float, float] | None = None
    best_coordinates: Array | None = None
    best_record: dict[str, Json] = {}

    def consider(xyz: Array, hinges: list[str], angles: list[int]) -> None:
        nonlocal sampled, minimum_excess, best_score, best_coordinates, best_record
        sampled += 1
        z_max = float(np.max((xyz[critical] @ membrane_rotation.T + membrane_translation)[:, 2]))
        excess = max(0.0, z_max - intracellular_bound)
        minimum_excess = min(minimum_excess, excess)
        if excess > 0 or not entry_allowed(xyz):
            return
        clashes, penalty = sterics(xyz)
        score = (clashes, penalty, float(sum(min(a, 360 - a) ** 2 for a in angles)))
        if best_score is None or score < best_score:
            best_score, best_coordinates = score, xyz
            best_record = {"hinges": hinges, "rotation_degrees": angles, "critical_z_max_A": z_max}

    first_phi = hinge(first, "phi")
    options = [(first + offset, kind) for offset, kind in region.second_hinges]
    for number, kind in options:
        second = hinge(number, kind)
        for angle1 in range(0, 360, 10):
            xyz1 = rotate(coordinates, first_phi, angle1)
            for angle2 in range(0, 360, 10):
                consider(rotate(xyz1, second, angle2), [f"{first}:phi", f"{number}:{kind}"], [angle1, angle2])
    if best_score is None and region.third_hinges:
        first_psi = hinge(first, "psi")
        for offset, kind in region.third_hinges:
            number = first + offset
            third = hinge(number, kind)
            for angle1 in range(0, 360, 20):
                xyz1 = rotate(coordinates, first_phi, angle1)
                value = phi(xyz1, first)
                if value is None or not -170 <= value <= -30:
                    continue
                for angle2 in range(0, 360, 20):
                    xyz2 = rotate(xyz1, first_psi, angle2)
                    for angle3 in range(0, 360, 20):
                        consider(rotate(xyz2, third, angle3), [f"{first}:phi", f"{first}:psi", f"{number}:{kind}"],
                                 [angle1, angle2, angle3])
    if best_score is None or best_coordinates is None:
        raise ValueError(f"No membrane-excluded orientation for {chain}{first}; residual excess {minimum_excess:.3f} A")
    if best_score[0] != 0:
        raise ValueError(f"All membrane-excluded orientations for {chain}{first} have collisions: {best_score[0]}")
    return best_coordinates, {"chain": chain, "sampled": sampled, "boundary_A": intracellular_bound, **best_record}


def refine_tails(source: Path, output: Path, template_path: Path,
                 target_sequences: Mapping[str, str], chain_ranges: Mapping[str, tuple[int, int]],
                 regions: list[TailRegion], rotation: Array, translation: Array,
                 membrane_half_width: float) -> dict[str, object]:
    report_path = output.with_suffix(".json")
    if any(path.exists() for path in (output, report_path)):
        raise FileExistsError("Output or report exists")
    rotation, translation = checked_transform(rotation, translation)
    if not np.isfinite(membrane_half_width) or membrane_half_width <= 0:
        raise ValueError("Membrane half width must be finite and positive")
    boundary = -membrane_half_width - 0.10
    model, template = read_structure(source), read_structure(template_path)
    if model.read_issues or template.read_issues:
        raise ValueError("Input or template has ambiguous atom identities")
    if not regions or any(region.chain not in chain_ranges or not
                          chain_ranges[region.chain][0] <= region.first <= region.last <= chain_ranges[region.chain][1]
                          for region in regions):
        raise ValueError("Tail regions must lie within the supplied protein chains")
    if any(a.chain == r.chain and r.first <= a.number <= r.last for a in template.atoms for r in regions):
        raise ValueError("Tail rotation would move experimental template atoms")
    coordinates = np.array([a.xyz for a in model.atoms], dtype=np.float64)
    tail_reports: list[dict[str, Json]] = []
    for region in regions:
        coordinates, result = refine_tail(model, template, coordinates, region, rotation, translation, boundary)
        tail_reports.append(result)
    moved = {a.index for a in model.atoms if any(a.chain == r.chain and r.first <= a.number <= r.last for r in regions)}
    by_identity = {(a.chain, a.number, a.insertion, a.name): a for a in model.atoms}
    lines: list[str] = []
    for line in source.read_text().splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            key = line[21], int(line[22:26]), line[26].strip(), line[12:16].strip()
            atom = by_identity[key]
            if atom.index in moved:
                xyz = coordinates[atom.index]
                line = line[:30] + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}" + line[54:]
        lines.append(line)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        handle.write("".join(lines))
    reread = read_structure(output)
    by_key = {a.key: a for a in reread.atoms}
    moved_count = 0
    for atom in model.atoms:
        drift = float(np.linalg.norm(by_key[atom.key].xyz - atom.xyz))
        if atom.index not in moved and drift != 0:
            raise ValueError("An atom outside the specified tails moved")
        moved_count += int(atom.index in moved and drift > 0)
    report = inspect_structure(reread, template, dict(target_sequences), dict(chain_ranges))
    result = {"moved_atoms": moved_count, "tails": tail_reports, **report}
    write_json(report_path, result)
    return result
