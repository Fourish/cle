"""Normalize prepared atom identities and explicit bonds without inferring chemistry."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TypedDict, cast


class AtomIdentity(TypedDict):
    index: int
    chain: str
    residue_id: str
    residue_name: str
    atom_name: str
    element: str


def read_topology(identity_path: Path, bonds_path: Path) -> tuple[list[AtomIdentity], list[tuple[int, int]]]:
    """Accept current preparation fields and existing short-name topology records."""
    raw = json.loads(identity_path.read_text())
    if not isinstance(raw, list) or not raw:
        raise ValueError("Atom identities require a nonempty list")
    identities: list[AtomIdentity] = []
    aliases = {"residue_id": "resid", "residue_name": "resname", "atom_name": "name"}
    for item in raw:
        if not isinstance(item, dict) or type(item.get("index")) is not int or item["index"] < 0:
            raise ValueError("Atom identities require nonnegative integer indices")
        fields: dict[str, str] = {}
        for name in ("chain", "residue_id", "residue_name", "atom_name", "element"):
            alias = aliases.get(name, name)
            value = item.get(name, item.get(alias))
            if not isinstance(value, (str, int)) or isinstance(value, bool) or not str(value):
                raise ValueError(f"Missing or invalid atom identity field: {name}")
            if name != alias and name in item and alias in item and str(item[name]) != str(item[alias]):
                raise ValueError(f"Conflicting atom identity fields: {name}/{alias}")
            fields[name] = str(value)
        identities.append({"index": item["index"], "chain": fields["chain"], "residue_id": fields["residue_id"],
                           "residue_name": fields["residue_name"], "atom_name": fields["atom_name"],
                           "element": fields["element"]})
    indices = {row["index"] for row in identities}
    if len(indices) != len(identities):
        raise ValueError("Duplicate topology atom indices")
    raw_bonds = json.loads(bonds_path.read_text())
    if not isinstance(raw_bonds, list):
        raise ValueError("Topology bonds require a list")
    bonds: list[tuple[int, int]] = []
    for item in raw_bonds:
        if isinstance(item, dict):
            pair = [item.get("atom_index_0"), item.get("atom_index_1")]
        elif isinstance(item, list) and len(item) == 2:
            pair = item
        else:
            raise ValueError("Invalid declared topology bond")
        if any(type(index) is not int or index not in indices for index in pair) or pair[0] == pair[1]:
            raise ValueError("Invalid declared topology bond")
        bonds.append((cast(int, pair[0]), cast(int, pair[1])))
    return identities, bonds
