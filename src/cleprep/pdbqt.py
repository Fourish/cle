"""Export a prepared protein/cholesterol complex as a rigid Vina receptor."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.resources import files
import json
from pathlib import Path
from typing import Mapping, TypeAlias, TypedDict

import numpy as np
import openmm
from openmm import XmlSerializer, unit
from meeko.atomtyper import AtomTyper
from meeko.molsetup import RDKitMoleculeSetup
from rdkit import Chem

from .common.geometry import Array
from .common.io import Json, write_json
from .common.protein import AA_ONE

ARRAY: TypeAlias = Array
Key: TypeAlias = tuple[str, str, str]


class ResidueTemplate(TypedDict):
    smiles: str
    atom_name: list[str]


class TemplateLibrary(TypedDict):
    residue_templates: dict[str, ResidueTemplate]


class SelectedSidechain(TypedDict):
    chain: str
    residue_id: int
    selected_variant: str


class StateSelection(TypedDict):
    sidechains: list[SelectedSidechain]


class BaseIdentity(TypedDict):
    index: int
    chain: str
    residue_id: str
    residue_name: str
    atom_name: str
    element: str


@dataclass(frozen=True)
class Atom:
    index: int
    chain: str
    resid: str
    resname: str
    name: str
    element: str
    charge: float

    @property
    def key(self) -> Key:
        return self.chain, self.resid, self.name


def read_pdb(path: Path) -> dict[Key, tuple[str, str, ARRAY]]:
    records: dict[Key, tuple[str, str, ARRAY]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = line[21], line[22:26].strip(), line[12:16].strip()
        if key in records:
            raise ValueError(f"Duplicate atom {key} in {path}")
        records[key] = (line[17:20].strip(), line[76:78].strip(), np.array([
            float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float))
    return records


def load_base(pdb: Path, chemistry: Path) -> tuple[list[Atom], list[tuple[int, int]], ARRAY]:
    identities: list[BaseIdentity] = json.loads((chemistry / "atom_identity.json").read_text())
    system = XmlSerializer.deserialize((chemistry / "system.xml").read_text())
    nonbonded = [force for force in system.getForces() if isinstance(force, openmm.NonbondedForce)]
    if len(nonbonded) != 1 or system.getNumParticles() != len(identities):
        raise ValueError("System/identity count mismatch or missing per-atom charges")
    charges = [float(nonbonded[0].getParticleParameters(i)[0].value_in_unit(unit.elementary_charge))
               for i in range(len(identities))]
    atoms = [Atom(r["index"], r["chain"], r["residue_id"], r["residue_name"],
                  r["atom_name"], r["element"], charges[i]) for i, r in enumerate(identities)]
    records = read_pdb(pdb)
    if ([a.index for a in atoms] != list(range(len(atoms)))
            or len(records) != len(atoms) or len({a.key for a in atoms}) != len(atoms)):
        raise ValueError("Atom identities are duplicated, incomplete or not contiguous")
    for atom in atoms:
        name, symbol, _ = records[atom.key]
        if name != atom.resname[:3] or symbol != atom.element:
            raise ValueError(f"Coordinate/chemistry identity mismatch: {atom.key}")
    stored: list[dict[str, int | None]] = json.loads((chemistry / "topology_bonds.json").read_text())
    bonds = [(int(r["atom_index_0"]), int(r["atom_index_1"])) for r in stored]
    coords = np.array([records[a.key][2] for a in atoms], dtype=np.float64)
    if not np.isfinite(coords).all() or not np.isfinite(charges).all():
        raise ValueError("Nonfinite coordinates or charges")
    return atoms, bonds, coords


def chemical_graph(atoms: list[Atom], bonds: list[tuple[int, int]], coordinates: ARRAY,
                   chemistry: Path, chain_starts: Mapping[str, int] | None = None) -> Chem.Mol:
    """Enrich validated connectivity with explicit, named chemical templates."""
    template_path = files("meeko").joinpath("data/residue_chem_templates.json")
    template_library: TemplateLibrary = json.loads(template_path.read_text())
    template_data = template_library["residue_templates"]
    states_data: StateSelection = json.loads((chemistry / "selected_states.json").read_text())
    variants: dict[tuple[str, str], str] = {
        (r["chain"], str(r["residue_id"])): r["selected_variant"] for r in states_data["sidechains"]}
    grouped: dict[tuple[str, str], list[Atom]] = defaultdict(list)
    adjacent: dict[int, list[int]] = defaultdict(list)
    for atom in atoms:
        grouped[(atom.chain, atom.resid)].append(atom)
    for left, right in bonds:
        adjacent[left].append(right)
        adjacent[right].append(left)
    mol = Chem.RWMol()
    for atom in atoms:
        rd = Chem.Atom(atom.element)
        rd.SetNoImplicit(True)
        rd.SetProp("source_key", ":".join(atom.key))
        mol.AddAtom(rd)
    orders: dict[frozenset[int], Chem.BondType] = {}
    heavy_expected: set[frozenset[int]] = set()
    chl_text = (chemistry / "chl1_locked_rtf_excerpt.txt").read_text() if any(a.resname == "CHL1" for a in atoms) else ""
    starts = dict(chain_starts) if chain_starts is not None else {
        chain: min(int(a.resid) for a in atoms if a.chain == chain and a.resname in AA_ONE)
        for chain in dict.fromkeys(a.chain for a in atoms if a.resname in AA_ONE)}
    for (chain, resid), members in grouped.items():
        names = {a.name: a for a in members}
        heavy = {a.name: a for a in members if a.element != "H"}
        if members[0].resname in AA_ONE:
            variant = variants.get((chain, resid), members[0].resname)
            if variant in {"HIS", "CYS", "LYS", "GLU", "ASP"}:
                variant = members[0].resname
            if resid == str(starts[chain]):
                variant = "N" + variant
            elif "OXT" in names:
                variant = "C" + variant
            row = template_data[variant]
            options = Chem.SmilesParserParams()
            options.removeHs = False
            template = Chem.MolFromSmiles(row["smiles"], options)
            if template is None:
                raise ValueError(f"Invalid Meeko template {variant}")
            labels: list[str] = row["atom_name"]
            template_atoms = {labels[a.GetIdx()]: a for a in template.GetAtoms() if a.GetAtomicNum() != 1}
            if set(template_atoms) != set(heavy):
                raise ValueError(f"Template heavy coverage mismatch {chain}:{resid}:{variant}")
            for name, templ_atom in template_atoms.items():
                rd = mol.GetAtomWithIdx(heavy[name].index)
                rd.SetFormalCharge(templ_atom.GetFormalCharge())
                rd.SetIsAromatic(templ_atom.GetIsAromatic())
                template_h = sum(a.GetAtomicNum() == 1 for a in templ_atom.GetNeighbors())
                actual_h = sum(atoms[i].element == "H" for i in adjacent[heavy[name].index])
                if actual_h != template_h:
                    raise ValueError(f"Explicit microstate H mismatch {chain}:{resid}:{name}: {actual_h} != {template_h}")
            for bond in template.GetBonds():
                left, right = labels[bond.GetBeginAtomIdx()], labels[bond.GetEndAtomIdx()]
                if left in heavy and right in heavy:
                    key = frozenset((heavy[left].index, heavy[right].index))
                    orders[key] = bond.GetBondType()
                    heavy_expected.add(key)
        elif members[0].resname == "CHL1":
            for line in chl_text.splitlines():
                words = line.split("!", 1)[0].split()
                if not words or words[0] not in {"BOND", "DOUBLE"}:
                    continue
                order = Chem.BondType.DOUBLE if words[0] == "DOUBLE" else Chem.BondType.SINGLE
                for index in range(1, len(words), 2):
                    left, right = words[index:index + 2]
                    if left in heavy and right in heavy:
                        key = frozenset((heavy[left].index, heavy[right].index))
                        orders[key] = order
                        heavy_expected.add(key)
        else:
            raise ValueError(f"Unsupported component {chain}:{resid}")
    for left, right in bonds:
        key = frozenset((left, right))
        if atoms[left].element != "H" and atoms[right].element != "H":
            if atoms[left].chain == atoms[right].chain and atoms[left].resid != atoms[right].resid:
                if {atoms[left].name, atoms[right].name} != {"C", "N"}:
                    raise ValueError("Unexpected interresidue bond")
                heavy_expected.add(key)
            if key not in heavy_expected:
                raise ValueError(f"Heavy bond absent in chemical templates: {left}-{right}")
        mol.AddBond(left, right, orders.get(key, Chem.BondType.SINGLE))
    graph_edges = {frozenset((a, b)) for a, b in bonds}
    if not heavy_expected <= graph_edges:
        raise ValueError("Chemical template contains bonds absent from audited topology")
    result = mol.GetMol()
    Chem.SanitizeMol(result)
    if any(a.GetNumRadicalElectrons() for a in result.GetAtoms()):
        raise ValueError("Reconstructed chemistry unexpectedly contains radicals")
    if abs(Chem.GetFormalCharge(result) - sum(a.charge for a in atoms)) > 1e-5:
        raise ValueError("Complete complex charge mismatch")
    conformer = Chem.Conformer(len(atoms))
    conformer.Set3D(True)
    for atom in atoms:
        conformer.SetAtomPosition(atom.index, coordinates[atom.index])
    result.AddConformer(conformer)
    return result


def meeko_types(molecule: Chem.Mol) -> list[str]:
    """Call Meeko's atom typer on each actual connected component."""
    params: dict[str, list[dict[str, Json]]] = json.loads((files("meeko").joinpath("data/params/ad4_types.json")).read_text())
    fragments: list[tuple[int, ...]] = []
    pieces = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True, fragsMolAtomMapping=fragments)
    types = [""] * molecule.GetNumAtoms()
    for piece, indices in zip(pieces, fragments, strict=True):
        setup = RDKitMoleculeSetup.from_mol(piece, charge_model="zero", compute_charges=True)
        AtomTyper.type_everything(setup, params, "zero")
        for local, global_index in enumerate(indices):
            atype = setup.get_atom_type(local)
            if not isinstance(atype, str) or not atype:
                raise ValueError(f"Meeko did not type atom {global_index}")
            types[global_index] = atype
    if not all(types):
        raise ValueError("Incomplete AutoDock atom type coverage")
    return types


def pdb_line(atom: Atom, serial: int, xyz: ARRAY, charge: float | None = None, atype: str | None = None) -> str:
    record = "ATOM  " if atom.resname in AA_ONE else "HETATM"
    name = f" {atom.name:<3}" if len(atom.name) < 4 else atom.name
    left = f"{record}{serial:5d} {name:4s} {atom.resname[:3]:3s} {atom.chain:1s}{int(atom.resid):4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
    if charge is not None and atype is not None:
        return f"{left}  1.00  0.00    {charge:6.3f} {atype:<2s}\n"
    return f"{left}  1.00  0.00          {atom.element:>2s}  \n"


def write_receptor(atoms: list[Atom], bonds: list[tuple[int, int]], types: list[str],
                   coords: ARRAY, output: Path, keep_cholesterol: bool,
                   expected_heavy_atoms: int | None = None) -> dict[str, object]:
    selected = {a.index for a in atoms if a.resname in AA_ONE
                or (a.resname == "CHL1" and keep_cholesterol)}
    parent: dict[int, int] = {}
    for left, right in bonds:
        if atoms[left].element == "H":
            parent[left] = right
        elif atoms[right].element == "H":
            parent[right] = left
    visible = [i for i in sorted(selected) if types[i] != "H"]
    serial = {index: number + 1 for number, index in enumerate(visible)}
    charges = {i: atoms[i].charge for i in visible}
    for index in sorted(selected):
        if types[index] == "H":
            if atoms[index].element != "H" or atoms[parent[index]].element != "C" or parent[index] not in serial:
                raise ValueError("Only carbon-bound nonpolar hydrogens can be merged")
            charges[parent[index]] += atoms[index].charge
    heavy = sum(atoms[i].element != "H" for i in selected)
    if expected_heavy_atoms is not None and heavy != expected_heavy_atoms:
        raise ValueError(f"Receptor heavy-atom coverage {heavy} != {expected_heavy_atoms}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="ascii") as handle:
        handle.write("REMARK VINA ONLY. CHARMM charges; nonpolar H merged into bonded C.\n")
        handle.write("REMARK Rigid standard-protein receptor; explicit microstates retained.\n")
        handle.writelines(pdb_line(atoms[i], serial[i], coords[i], charges[i], types[i]) for i in visible)
    return {"all_atoms": len(selected), "heavy_atoms": heavy, "pdbqt_atoms": len(visible),
            "nonpolar_hydrogens_merged": len(selected) - len(visible),
            "charge_e_before_rounding": sum(charges.values()), "cholesterol_retained": keep_cholesterol,
            "scoring_function": "vina", "autodock_types": dict(Counter(types[i] for i in visible))}


def export_receptor(source: Path, chemistry: Path, output: Path, *, keep_cholesterol: bool = False,
                    chain_starts: Mapping[str, int] | None = None,
                    expected_heavy_atoms: int | None = None) -> dict[str, object]:
    """Export one rigid receptor using matching identities, bonds and CHARMM charges."""
    if output.suffix.lower() != ".pdbqt":
        raise ValueError("Output must end in .pdbqt")
    if any(path.exists() for path in (output, output.with_suffix(".json"))):
        raise FileExistsError("Output or report exists")
    atoms, bonds, coords = load_base(source, chemistry)
    molecule = chemical_graph(atoms, bonds, coords, chemistry, chain_starts)
    types = meeko_types(molecule)
    summary = write_receptor(atoms, bonds, types, coords, output, keep_cholesterol, expected_heavy_atoms)
    write_json(output.with_suffix(".json"), summary)
    return summary
