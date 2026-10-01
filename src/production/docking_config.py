"""Read explicit docking workflow inputs; paths are relative to the JSON file."""

from __future__ import annotations

from dataclasses import fields
import json
import math
from pathlib import Path
import re
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from cleprep.workflow import WorkflowConfig

Record = dict[str, object]


def object_value(value: object, allowed: set[str], required: set[str] | frozenset[str] = frozenset()) -> Record:
    if not isinstance(value, dict) or not required <= value.keys() or value.keys() - allowed:
        raise ValueError(f"JSON object requires {sorted(required)} and allows only {sorted(allowed)}")
    return cast(Record, value)


def rows(value: object) -> list[Record]:
    if not isinstance(value, list) or not value or any(not isinstance(row, dict) for row in value):
        raise ValueError("Require a nonempty list of JSON objects")
    return cast(list[Record], value)


def identifier(value: object) -> str:
    if not isinstance(value, str) or len(value) > 64 or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value) or "__" in value:
        raise ValueError("IDs require letters/numbers/._-, without '__'")
    return value


def chain_list(value: object) -> tuple[str, ...]:
    if (not isinstance(value, list) or not value or any(not isinstance(item, str) or not re.fullmatch(r"[A-Za-z0-9]", item) for item in value)
            or len(set(value)) != len(value)):
        raise ValueError("Chains require a nonempty list of unique one-character letters/numbers")
    return tuple(value)


def integer(value: object, minimum: int, maximum: int | None = None) -> int:
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f"Integer must be in {minimum}..{maximum or 'unbounded'}")
    return value


def finite(value: object, minimum: float, maximum: float | None = None) -> float:
    if type(value) not in {int, float} or not math.isfinite(float(value)) or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f"Finite number must be in {minimum}..{maximum or 'unbounded'}")
    return float(value)


def vector(value: object, positive: bool = False) -> tuple[float, float, float]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError("Box vector requires three numbers")
    result = tuple(finite(item, -float("inf")) for item in value)
    if positive and min(result) <= 0:
        raise ValueError("Box side lengths must be positive")
    return cast(tuple[float, float, float], result)


def load_config(path: Path, *, check_paths: bool = True) -> WorkflowConfig:
    # Keep TUI/parser loading free of scientific imports.
    from cleprep.complexes import GeometryPolicy
    from cleprep.workflow import (DockingSettings, LigandInput, ReceptorInput,
                                 SearchInput, SelectionSettings, WorkflowConfig)

    path = path.expanduser().resolve()
    data = object_value(json.loads(path.read_text()),
                        {"version", "receptors", "ligands", "searches", "docking", "geometry", "selection", "existing_runs"},
                        {"version", "receptors", "ligands", "searches"})
    if type(data["version"]) is not int or data["version"] != 1:
        raise ValueError("Workflow configuration requires version=1")

    def file_path(value: object, directory: bool = False) -> Path:
        if not isinstance(value, str) or not value:
            raise ValueError("Paths must be nonempty strings")
        resolved = Path(value).expanduser()
        resolved = (resolved if resolved.is_absolute() else path.parent / resolved).resolve()
        if check_paths and not (resolved.is_dir() if directory else resolved.is_file()):
            raise FileNotFoundError(resolved)
        return resolved

    receptors: list[ReceptorInput] = []
    for row in rows(data["receptors"]):
        allowed = {field.name for field in fields(ReceptorInput)}
        object_value(row, allowed, allowed - {"expected_residues"})
        count = integer(row["expected_residues"], 1) if "expected_residues" in row else None
        receptors.append(ReceptorInput(identifier(row["id"]), file_path(row["pdbqt"]), file_path(row["pdb"]),
                         file_path(row["atom_identity"]), file_path(row["topology_bonds"]),
                         chain_list(row["protein_chains"]), chain_list(row["contact_chains"]), count))
    ligands: list[LigandInput] = []
    for row in rows(data["ligands"]):
        object_value(row, {"id", "sdf", "records"}, {"id", "sdf"})
        chosen = row.get("records", "all")
        if chosen == "all":
            indices = None
        elif isinstance(chosen, list) and chosen:
            indices = tuple(integer(item, 0) for item in chosen)
            if len(set(indices)) != len(indices):
                raise ValueError("Duplicate SDF records would inflate sampling support")
        else:
            raise ValueError("records requires 'all' or nonempty zero-based indices")
        ligands.append(LigandInput(identifier(row["id"]), file_path(row["sdf"]), indices))
    searches: list[SearchInput] = []
    receptor_ids = {item.id for item in receptors}
    for row in rows(data["searches"]):
        object_value(row, {"id", "receptors", "center", "size"}, {"id", "receptors", "center", "size"})
        selected = row["receptors"]
        if (not isinstance(selected, list) or not selected or any(not isinstance(item, str) for item in selected)
                or len(set(selected)) != len(selected) or not set(selected) <= receptor_ids):
            raise ValueError("Search references unknown/duplicate receptor IDs")
        searches.append(SearchInput(identifier(row["id"]), tuple(selected), vector(row["center"]), vector(row["size"], True)))
    for items in (receptors, ligands, searches):
        if len({item.id for item in items}) != len(items):
            raise ValueError("Duplicate input/search IDs")
    docking_row = object_value(data.get("docking", {}), {field.name for field in fields(DockingSettings)})
    seeds = docking_row.get("seeds", [20260910])
    if not isinstance(seeds, list) or not seeds:
        raise ValueError("seeds requires a nonempty integer list")
    seed_values = tuple(integer(item, 1, 2147483647) for item in seeds)
    if len(set(seed_values)) != len(seed_values):
        raise ValueError("Duplicate seeds would inflate independent support")
    docking = DockingSettings(seed_values, integer(docking_row.get("cpu", 1), 1),
                              integer(docking_row.get("exhaustiveness", 8), 1), integer(docking_row.get("max_evals", 100000), 0),
                              integer(docking_row.get("n_poses", 10), 1))
    geometry_row = object_value(data.get("geometry", {}), {field.name for field in fields(GeometryPolicy)})
    geometry = GeometryPolicy(**{key: integer(value, 0) if key.startswith("minimum_contact_") else
        finite(value, 1e-10, 1.0 if key.endswith("vdw_ratio") else None) for key, value in geometry_row.items()})
    selection_row = object_value(data.get("selection", {}), {field.name for field in fields(SelectionSettings)})
    chain = selection_row.get("ligand_chain", "P")
    chain_list([chain])
    selection = SelectionSettings(finite(selection_row.get("backbone_rmsd", 3.0), 1e-10),
                  finite(selection_row.get("contact_jaccard", 0.5), 0.0, 1.0),
                  integer(selection_row.get("max_representatives_per_group", 3), 1), cast(str, chain))
    existing = data.get("existing_runs", [])
    if not isinstance(existing, list):
        raise ValueError("existing_runs requires a list of run directories")
    directories = tuple(file_path(item, True) for item in existing)
    if len(set(directories)) != len(directories):
        raise ValueError("Duplicate existing run directory")
    return WorkflowConfig(tuple(receptors), tuple(ligands), tuple(searches), docking, geometry, selection, directories)
