"""Prepare standard proteins and optional CHARMM cholesterol with fixed heavy atoms.

Context ligands remain separate from the energy model. No disulfides are inferred.
"""

from __future__ import annotations

import math
import os
import random
import re
import subprocess
import sys
import tarfile
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Mapping, TypeAlias, cast

import gemmi
import numpy as np
import openmm
from openmm import Platform, Vec3, XmlSerializer, unit
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile, Topology, element
from openmm.app.topology import Atom
from rdkit import Chem

from .common.io import write_json
from .common.protein import PROTEIN_NAMES
AtomKey: TypeAlias = tuple[str, str, str]
HydrogenVariant: TypeAlias = str | list[tuple[str, str]] | None


@dataclass(frozen=True)
class PkaSite:
    group: str
    residue_id: int
    chain: str
    pka: float


@dataclass
class StateChoice:
    chain: str
    residue_id: int
    residue_name: str
    pka: float | None
    selected_variant: str


def require_element(parent: ET.Element, query: str) -> ET.Element:
    found = parent.find(query)
    if found is None:
        raise ValueError(f"Missing force-field XML element: {query}")
    return found


def charmm_cholesterol_graph(xml_root: ET.Element, toppar: Path) -> tuple[dict[str, str], list[tuple[str, str, int]], list[tuple[str, str]], str]:
    """Use XML atom types and the locked RTF's explicit bond orders."""
    residue = require_element(xml_root, "Residues/Residue[@name='CHL1']")
    types = {node.attrib["name"]: node.attrib["element"] for node in xml_root.findall("AtomTypes/Type") if "element" in node.attrib}
    symbols = {node.attrib["name"]: types[node.attrib["type"]] for node in residue.findall("Atom")}
    with tarfile.open(toppar) as archive:
        stream = archive.extractfile("toppar/stream/lipid/toppar_all36_lipid_cholesterol.str")
        if stream is None:
            raise ValueError("CHARMM cholesterol stream is absent")
        text = stream.read().decode()
    block = re.split(r"(?m)^RESI\s+CHL1\s", text)[1]
    block = re.split(r"(?m)^(?:RESI|PRES)\s", block)[0]
    all_bonds: list[tuple[str, str, int]] = []
    for line in block.splitlines():
        words = line.split("!", 1)[0].split()
        if words and words[0] in {"BOND", "DOUBLE"}:
            order = 2 if words[0] == "DOUBLE" else 1
            for index in range(1, len(words), 2):
                all_bonds.append((words[index], words[index + 1], order))
    xml_edges = {frozenset((b.attrib["atomName1"], b.attrib["atomName2"])) for b in residue.findall("Bond")}
    if xml_edges != {frozenset((left, right)) for left, right, _ in all_bonds}:
        raise ValueError("CHL1 XML and locked RTF bond graphs disagree")
    heavy = {name: symbol for name, symbol in symbols.items() if symbol != "H"}
    heavy_bonds = [(left, right, order) for left, right, order in all_bonds if left in heavy and right in heavy]
    hydrogen_pairs: list[tuple[str, str]] = []
    for left, right, _ in all_bonds:
        if symbols[left] == "H":
            hydrogen_pairs.append((left, right))
        elif symbols[right] == "H":
            hydrogen_pairs.append((right, left))
    return heavy, heavy_bonds, hydrogen_pairs, block


def graph_molecule(symbols: dict[str, str], bonds: list[tuple[str, str, int]]) -> tuple[Chem.Mol, list[str]]:
    names = sorted(symbols)
    indexes = {name: index for index, name in enumerate(names)}
    mol = Chem.RWMol()
    for name in names:
        mol.AddAtom(Chem.Atom(symbols[name]))
    for left, right, order in bonds:
        mol.AddBond(indexes[left], indexes[right], Chem.BondType.DOUBLE if order == 2 else Chem.BondType.SINGLE)
    result = mol.GetMol()
    Chem.SanitizeMol(result)
    return result, names


def cholesterol_mapping(xml_root: ET.Element, out: Path, source_cif: Path, toppar: Path) -> tuple[dict[str, str], list[tuple[str, str, int]], list[tuple[str, str]]]:
    source = gemmi.cif.read_file(str(source_cif)).sole_block()
    symbols: dict[str, str] = {}
    for row in source.find(["_chem_comp_atom.comp_id", "_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol"]):
        if row[0] == "CLR" and row[2] != "H":
            symbols[row[1]] = row[2]
    source_bonds: list[tuple[str, str, int]] = []
    for row in source.find(["_chem_comp_bond.comp_id", "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2", "_chem_comp_bond.value_order"]):
        if row[0] == "CLR" and row[1] in symbols and row[2] in symbols:
            if row[3].lower() not in {"sing", "doub"}:
                raise ValueError(f"Unsupported CLR bond order {row[3]}")
            source_bonds.append((row[1], row[2], 2 if row[3].lower() == "doub" else 1))
    target_symbols, target_bonds, hydrogen_pairs, rtf_block = charmm_cholesterol_graph(xml_root, toppar)
    source_mol, source_names = graph_molecule(symbols, source_bonds)
    target_mol, target_names = graph_molecule(target_symbols, target_bonds)
    matches = source_mol.GetSubstructMatches(target_mol, uniquify=False, useChirality=False, maxMatches=1000)
    if not matches or len(source_names) != len(target_names):
        raise ValueError("CLR CCD and CHL1 are not heavy-atom graph isomorphic")
    mappings = [{source_names[source_index]: target_names[target_index] for target_index, source_index in enumerate(match)} for match in matches]
    mapping = min(mappings, key=lambda m: (-sum(left == right for left, right in m.items()), sorted(m.items())))
    mapped_edges = {(frozenset((mapping[left], mapping[right])), order) for left, right, order in source_bonds}
    if mapped_edges != {(frozenset((left, right)), order) for left, right, order in target_bonds}:
        raise ValueError("Mapped cholesterol bond orders disagree")
    (out / "chl1_locked_rtf_excerpt.txt").write_text("RESI CHL1 " + rtf_block)
    return mapping, target_bonds, hydrogen_pairs


def vec3(position: gemmi.Position) -> Vec3:
    return Vec3(position.x * 0.1, position.y * 0.1, position.z * 0.1)


def build_topology(input_path: Path, out: Path, mapping: dict[str, str],
                   chl_bonds: list[tuple[str, str, int]], chain_ranges: Mapping[str, tuple[int, int]], *,
                   context_residues: frozenset[str] = frozenset(), expected_context_atoms: int | None = None,
                   expected_cholesterol_count: int | None = None, cholesterol_chain: str = "L",
                   cholesterol_residue_id: str = "1") -> tuple[Topology, list[Vec3], dict[AtomKey, Vec3]]:
    structure = gemmi.read_structure(str(input_path))
    if len(structure) != 1 or not chain_ranges or any(start > end for start, end in chain_ranges.values()):
        raise ValueError("One coordinate model and nonempty, ascending chain ranges are required")
    topology = Topology()
    positions: list[Vec3] = []
    original: dict[AtomKey, Vec3] = {}
    residue_lookup: dict[tuple[str, int], gemmi.Residue] = {}
    cholesterol: list[gemmi.Residue] = []
    context_atoms = 0
    for chain in structure[0]:
        for residue in chain:
            if residue.name in PROTEIN_NAMES:
                key = chain.name, residue.seqid.num
                if key in residue_lookup or residue.seqid.icode.strip():
                    raise ValueError(f"Duplicate residue or insertion code: {key}")
                residue_lookup[key] = residue
            elif residue.name == "CLR":
                cholesterol.append(residue)
            elif residue.name in context_residues:
                context_atoms += sum(atom.element.name != "H" for atom in residue)
            else:
                raise ValueError(f"Unexpected residue {chain.name}:{residue.seqid}:{residue.name}")
    expected = {(chain, number) for chain, (start, end) in chain_ranges.items() for number in range(start, end + 1)}
    if set(residue_lookup) != expected:
        raise ValueError(f"Protein sequence boundaries differ: missing={sorted(expected-set(residue_lookup))}, extra={sorted(set(residue_lookup)-expected)}")
    if len(cholesterol) > 1 or (expected_cholesterol_count is not None and len(cholesterol) != expected_cholesterol_count):
        raise ValueError("This preparation supports zero or one cholesterol; count differs from expectation")
    if expected_context_atoms is not None and context_atoms != expected_context_atoms:
        raise ValueError(f"Context heavy-atom count {context_atoms} != {expected_context_atoms}")
    if cholesterol and cholesterol_chain in chain_ranges:
        raise ValueError("Cholesterol needs a chain distinct from protein chains")
    for chain_id, (start, end) in chain_ranges.items():
        chain = topology.addChain(chain_id)
        for number in range(start, end + 1):
            source_res = residue_lookup[(chain_id, number)]
            residue = topology.addResidue(source_res.name, chain, str(number))
            names: dict[str, Vec3] = {}
            for source_atom in source_res:
                if source_atom.element.name == "H":
                    continue
                if source_atom.name in names:
                    raise ValueError(f"Duplicate atom/alternate conformation {chain_id}:{number}:{source_atom.name}")
                name = source_atom.name
                atom_position = vec3(source_atom.pos)
                topology.addAtom(name, element.Element.getBySymbol(source_atom.element.name), residue)
                positions.append(atom_position)
                names[name] = atom_position
                original[(chain_id, str(number), name)] = atom_position
            if number == end and "OXT" not in names:
                c = np.asarray(names["C"], dtype=np.float64)
                vca = np.asarray(names["CA"], dtype=np.float64) - c
                vo = np.asarray(names["O"], dtype=np.float64) - c
                direction = -(vca / np.linalg.norm(vca) + vo / np.linalg.norm(vo))
                oxt = c + 0.125 * direction / np.linalg.norm(direction)
                topology.addAtom("OXT", element.oxygen, residue)
                positions.append(Vec3(*oxt))
    topology.createStandardBonds()
    # Disulfide construction requires separate, explicit chemistry support.
    with (out / "propka_input.pdb").open("w") as handle:
        PDBFile.writeFile(topology, positions * unit.nanometer, handle, keepIds=True)
    if cholesterol:
        chain = topology.addChain(cholesterol_chain)
        clr_residue = topology.addResidue("CHL1", chain, cholesterol_residue_id)
        chl_atoms: dict[str, Atom] = {}
        seen: set[str] = set()
        for source_atom in cholesterol[0]:
            if source_atom.element.name == "H":
                continue
            if source_atom.name in seen or source_atom.name not in mapping:
                raise ValueError(f"Duplicate or unmapped cholesterol atom: {source_atom.name}")
            target_name = mapping[source_atom.name]
            seen.add(source_atom.name)
            chl_atoms[target_name] = topology.addAtom(target_name, element.Element.getBySymbol(source_atom.element.name), clr_residue)
            position = vec3(source_atom.pos)
            positions.append(position)
            original[(cholesterol_chain, cholesterol_residue_id, target_name)] = position
        if seen != set(mapping):
            raise ValueError("CLR heavy atom set is incomplete")
        for left, right, order in chl_bonds:
            topology.addBond(chl_atoms[left], chl_atoms[right], order=order)
    context_lines = [line for line in input_path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))
                     and line[17:20].strip() in context_residues]
    if context_lines:
        (out / "context_reference_unparameterized.pdb").write_text(
            "REMARK Context ligands excluded from the energy model.\n" + "\n".join(context_lines) + "\nEND\n")
    return topology, positions, original


def run_propka(out: Path, ph: float) -> list[PkaSite]:
    command = [sys.executable, "-m", "propka", "--pH", str(ph), "propka_input.pdb"]
    with (out / "propka.log").open("w") as handle:
        subprocess.run(command, cwd=out, stdout=handle, stderr=subprocess.STDOUT, check=True)
    text = (out / "propka_input.pka").read_text()
    summary = text.split("SUMMARY OF THIS PREDICTION", 1)[1]
    sites: list[PkaSite] = []
    for line in summary.splitlines():
        match = re.match(r"\s*(ASP|GLU|CYS|LYS|HIS|TYR|ARG|N\+|C-)\s+(-?\d+)\s+(\S)\s+(-?[\d.]+)\s+(-?[\d.]+)", line)
        if match:
            sites.append(PkaSite(match[1], int(match[2]), match[3], float(match[4])))
    if not sites:
        raise ValueError("No PROPKA summary entries parsed")
    return sites


def choose_states(topology: Topology, sites: list[PkaSite], ph: float,
                  hydrogen_pairs: list[tuple[str, str]]) -> tuple[list[HydrogenVariant], list[StateChoice]]:
    lookup = {(s.chain, s.residue_id, s.group): s.pka for s in sites}
    variants: list[HydrogenVariant] = []
    choices: list[StateChoice] = []
    for residue in topology.residues():
        if residue.name == "CHL1":
            variants.append(hydrogen_pairs)
            continue
        pka = lookup.get((residue.chain.id, int(residue.id), residue.name))
        variant: str | None = None
        if residue.name in {"ASP", "GLU"}:
            variant = ("ASH" if residue.name == "ASP" else "GLH") if pka is not None and pka > ph else residue.name
        elif residue.name == "LYS":
            variant = "LYN" if pka is not None and pka < ph else "LYS"
        elif residue.name == "CYS":
            variant = "CYX" if pka is not None and pka < ph else "CYS"
        elif residue.name == "HIS":
            variant = "HIP" if pka is not None and pka >= ph else None
        if residue.name in {"ASP", "GLU", "CYS", "LYS", "HIS", "ARG", "TYR"}:
            selected = variant or ("auto_neutral_histidine" if residue.name == "HIS" else residue.name)
            choices.append(StateChoice(residue.chain.id, int(residue.id), residue.name, pka, selected))
        variants.append(variant)
    return variants, choices


def validate_topology(topology: Topology, positions: unit.Quantity, original: dict[AtomKey, Vec3],
                      system: openmm.System) -> dict[str, int | float]:
    xyz = np.asarray(positions.value_in_unit(unit.nanometer), dtype=np.float64)
    if not np.isfinite(xyz).all() or system.getNumParticles() != topology.getNumAtoms():
        raise ValueError("Nonfinite coordinates or particle count mismatch")
    current = {(a.residue.chain.id, a.residue.id, a.name): xyz[a.index] for a in topology.atoms()}
    if not original.keys() <= current.keys():
        raise ValueError("Original heavy atoms missing")
    drift = max(float(np.linalg.norm(current[key] - position) * 10) for key, position in original.items())
    if drift > 1e-6:
        raise ValueError(f"Original heavy atoms moved: {drift} A")
    peptide_bonds = 0
    for a, b in topology.bonds():
        if a.name == b.name == "SG":
            raise ValueError("Unexpected disulfide")
        if a.residue != b.residue:
            if {a.name, b.name} != {"C", "N"}:
                raise ValueError("Unexpected interresidue bond")
            peptide_bonds += 1
    expected_bonds = sum(max(0, sum(residue.name in PROTEIN_NAMES for residue in chain.residues()) - 1)
                         for chain in topology.chains())
    if peptide_bonds != expected_bonds:
        raise ValueError(f"Protein connectivity mismatch: {peptide_bonds} peptide bonds")
    forces = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
    if len(forces) != 1:
        raise ValueError("Per-atom force-field charges unavailable")
    charges = [float(forces[0].getParticleParameters(i)[0].value_in_unit(unit.elementary_charge))
               for i in range(system.getNumParticles())]
    if not np.isfinite(charges).all():
        raise ValueError("Nonfinite charges")
    return {"atoms": topology.getNumAtoms(), "charge_e": sum(charges), "heavy_max_displacement_A": drift}


def relax_hydrogens(modeller: Modeller, system: openmm.System) -> float:
    """Relax only added H atoms; retain physical masses in the exported System."""
    fixed_system = cast(openmm.System, XmlSerializer.deserialize(XmlSerializer.serialize(system)))
    heavy: list[int] = []
    for atom in modeller.topology.atoms():
        if atom.element != element.hydrogen:
            heavy.append(atom.index)
            fixed_system.setParticleMass(atom.index, 0.0)
    before = np.asarray(modeller.positions.value_in_unit(unit.nanometer), dtype=float)
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    context = openmm.Context(fixed_system, integrator, Platform.getPlatformByName("CPU"))
    context.setPositions(modeller.positions)
    openmm.LocalEnergyMinimizer.minimize(context, 10.0, 1000)
    state = context.getState(getPositions=True, getEnergy=True)
    after = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float)
    drift = float(np.max(np.linalg.norm(after[heavy] - before[heavy], axis=1)) * 10.0)
    if drift > 1e-6 or not np.isfinite(after).all():
        raise ValueError(f"Hydrogen-only optimization moved heavy atoms: {drift} A")
    modeller.positions = state.getPositions()
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
    del context, integrator
    if not math.isfinite(energy):
        raise ValueError("Hydrogen relaxation energy is not finite")
    return energy


def prepare(source: Path, output: Path, forcefield_path: Path,
            chain_ranges: Mapping[str, tuple[int, int]], *, source_cif: Path | None = None,
            toppar: Path | None = None, context_residues: frozenset[str] = frozenset(),
            expected_context_atoms: int | None = None, expected_cholesterol_count: int | None = None,
            cholesterol_chain: str = "L", cholesterol_residue_id: str = "1", ph: float = 7.4,
            seed: int = 20260910, cpu_threads: int = 2,
            progress: Callable[[str], None] | None = None) -> dict[str, int | float]:
    """Run PROPKA, add explicit H, and relax H only; write a matching chemistry bundle."""
    if not math.isfinite(ph) or not 0 <= ph <= 14 or cpu_threads < 1:
        raise ValueError("pH must be in [0, 14] and cpu_threads must be positive")
    output.mkdir(parents=True, exist_ok=False)
    previous_random = random.getstate()
    previous_threads = os.environ.get("OPENMM_CPU_THREADS")
    try:
        random.seed(seed)
        os.environ["OPENMM_CPU_THREADS"] = str(cpu_threads)
        mapping: dict[str, str] = {}
        chl_bonds: list[tuple[str, str, int]] = []
        hydrogen_pairs: list[tuple[str, str]] = []
        structure = gemmi.read_structure(str(source))
        if any(residue.name == "CLR" for model in structure for chain in model for residue in chain):
            if source_cif is None or toppar is None:
                raise ValueError("Cholesterol requires its source CCD mmCIF and the CHARMM toppar archive")
            mapping, chl_bonds, hydrogen_pairs = cholesterol_mapping(ET.parse(forcefield_path).getroot(), output, source_cif, toppar)
        topology, positions, original = build_topology(
            source, output, mapping, chl_bonds, chain_ranges, context_residues=context_residues,
            expected_context_atoms=expected_context_atoms, expected_cholesterol_count=expected_cholesterol_count,
            cholesterol_chain=cholesterol_chain, cholesterol_residue_id=cholesterol_residue_id)
        if progress is not None:
            progress("Prepared heavy atoms and optional cholesterol mapping")
        sites = run_propka(output, ph)
        variants, choices = choose_states(topology, sites, ph, hydrogen_pairs)
        forcefield = ForceField(str(forcefield_path))
        modeller = Modeller(topology, positions * unit.nanometer)
        if progress is not None:
            progress("Adding and relaxing hydrogens with heavy atoms fixed")
        actual = cast(list[HydrogenVariant], modeller.addHydrogens(
            forcefield, pH=ph, variants=variants, platform=Platform.getPlatformByName("CPU")))
        choice_lookup = {(c.chain, c.residue_id): c for c in choices}
        for residue in modeller.topology.residues():
            if residue.name == "HIS":
                variant = actual[residue.index]
                if not isinstance(variant, str):
                    raise ValueError("Histidine variant was not determined")
                choice_lookup[residue.chain.id, int(residue.id)].selected_variant = variant
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)
        energy = relax_hydrogens(modeller, system)
        summary = validate_topology(modeller.topology, modeller.positions, original, system)
        (output / "system.xml").write_text(XmlSerializer.serialize(system))
        with (output / "protein_cholesterol.pdb").open("w") as handle:
            PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)
        write_json(output / "atom_identity.json", [{"index": a.index, "chain": a.residue.chain.id,
                   "residue_id": a.residue.id, "residue_name": a.residue.name,
                   "atom_name": a.name, "element": a.element.symbol} for a in modeller.topology.atoms()])
        write_json(output / "topology_bonds.json", [{"atom_index_0": a.index, "atom_index_1": b.index}
                                                 for a, b in modeller.topology.bonds()])
        write_json(output / "selected_states.json", {"pH": ph, "sidechains": [asdict(c) for c in choices]})
        summary.update(pH=ph, seed=seed, cpu_threads=cpu_threads, hydrogen_relaxation_energy_kJ_mol=energy)
        write_json(output / "summary.json", summary)
        return summary
    finally:
        random.setstate(previous_random)
        if previous_threads is None:
            os.environ.pop("OPENMM_CPU_THREADS", None)
        else:
            os.environ["OPENMM_CPU_THREADS"] = previous_threads
