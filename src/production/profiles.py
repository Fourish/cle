"""Small production presets; this module is not part of the reusable library."""

from __future__ import annotations

import json
from pathlib import Path
from typing import NotRequired, TypedDict, cast


class ParentSpec(TypedDict):
    residues: int
    prolines: int
    sulfoxides: int
    charge: int


class StateSpec(TypedDict):
    parent: str
    sulfur_cip: str
    seed_offset: int


class PeptideProfile(TypedDict):
    parents: dict[str, ParentSpec]
    states: dict[str, StateSpec]


class ChainSpec(TypedDict):
    label: str
    first: int
    last: int


class LigandSpec(TypedDict):
    name: str
    heavy_atoms: int


class ChemistrySpec(TypedDict, total=False):
    context_residues: list[str]
    expected_context_atoms: int
    expected_cholesterol_count: int
    cholesterol_chain: str
    cholesterol_residue_id: str


class TailSpec(TypedDict):
    chain: str
    first: int
    last: int
    critical_start: int
    entry_phi_count: NotRequired[int]
    second_hinges: NotRequired[list[tuple[int, str]]]
    third_hinges: NotRequired[list[tuple[int, str]]]


class ReceptorProfile(TypedDict):
    chains: dict[str, ChainSpec]
    ligands: NotRequired[dict[str, LigandSpec]]
    source_sha256: NotRequired[str]
    expected_protein_heavy_atoms: NotRequired[int]
    chemistry: NotRequired[ChemistrySpec]
    tails: NotRequired[list[TailSpec]]


def peptide_profile(path: Path) -> PeptideProfile:
    profile = cast(PeptideProfile, json.loads(path.read_text()))
    if not profile["parents"] or not profile["states"]:
        raise ValueError("Peptide preset requires parents and states")
    if any(state["parent"] not in profile["parents"] for state in profile["states"].values()):
        raise ValueError("Every state must reference a parent in the preset")
    return profile


def receptor_profile(path: Path) -> ReceptorProfile:
    profile = cast(ReceptorProfile, json.loads(path.read_text()))
    if not profile["chains"] or any(len(chain) != 1 or spec["first"] > spec["last"]
                                    for chain, spec in profile["chains"].items()):
        raise ValueError("Receptor preset requires one-character chains and ascending ranges")
    return profile


def chain_ranges(profile: ReceptorProfile) -> dict[str, tuple[int, int]]:
    return {chain: (spec["first"], spec["last"]) for chain, spec in profile["chains"].items()}


def chain_starts(profile: ReceptorProfile) -> dict[str, int]:
    return {chain: spec["first"] for chain, spec in profile["chains"].items()}
