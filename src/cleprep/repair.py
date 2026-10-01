"""Repair missing standard-protein regions with restricted MODELLER sampling."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
from typing import Mapping, TypedDict

import gemmi
import modeller
from modeller import Alignment, Environ, Selection
from modeller.automodel import LoopModel, assess, autosched, refine

from .common.io import read_fasta, write_json

AtomKey = tuple[str, int, str]
ResidueKey = tuple[str, int]


@dataclass
class RepairInputs:
    existing: set[AtomKey]
    movable: set[ResidueKey]
    chain_starts: dict[str, int]
    ligand_chain: str | None
    ligand_start: int
    ligand_identities: dict[int, tuple[str, int, str]]


class ModelOutput(TypedDict):
    filename: str
    kind: str
    dope_score: float | None
    molpdf: float | None
    failed: str | None


def prepare_inputs(source: Path, fasta: Path, directory: Path,
                   chain_starts: Mapping[str, int]) -> RepairInputs:
    sequences = read_fasta(fasta)
    if not chain_starts or set(sequences) != set(chain_starts):
        raise ValueError("FASTA chains must match the supplied chain starts")
    structure = gemmi.read_structure(str(source))
    if len(structure) != 1:
        raise ValueError("Exactly one template model is required")
    existing: set[AtomKey] = set()
    observed: set[ResidueKey] = set()
    for chain in structure[0]:
        if chain.name not in chain_starts:
            continue
        for residue in chain:
            if residue.het_flag != 'A':
                continue
            observed.add((chain.name, residue.seqid.num))
            existing.update((chain.name, residue.seqid.num, atom.name) for atom in residue)
    missing = {(chain, start + i) for chain, start in chain_starts.items()
               for i in range(len(sequences[chain])) if (chain, start + i) not in observed}
    movable = missing | {(chain, number + delta) for chain, number in missing for delta in (-1, 1)
                          if (chain, number + delta) in observed}
    source_lines = source.read_text().splitlines()
    output_lines: list[str] = []
    for chain in chain_starts:
        output_lines.extend(line for line in source_lines if line.startswith("ATOM  ") and line[21] == chain)
        output_lines.append("TER")
    ligand_lines = [line for line in source_lines if line.startswith("HETATM")]
    ligand_keys = list(dict.fromkeys((line[21], int(line[22:26]), line[26]) for line in ligand_lines))
    ligand_chain = next((c for c in "LABCDEFGHIJKMNOPQRSTUVWXYZ" if c not in chain_starts), None) if ligand_keys else None
    if ligand_keys and ligand_chain is None:
        raise ValueError("No spare one-character chain for MODELLER BLK ligands")
    ligand_start = min(key[1] for key in ligand_keys) if ligand_keys else 1
    ligand_identities = {ligand_start + i: key for i, key in enumerate(ligand_keys)}
    virtual_ids = {key: number for number, key in ligand_identities.items()}
    for line in ligand_lines:
        key = line[21], int(line[22:26]), line[26]
        output_lines.append(line[:21] + f"{ligand_chain}{virtual_ids[key]:4d} " + line[27:])
    if ligand_lines:
        output_lines.append("TER")
    output_lines.append("END")
    (directory / "template.pdb").write_text("\n".join(output_lines) + "\n")
    template_segments = ["".join(letter if (chain, start + i) in observed else "-"
                                for i, letter in enumerate(sequences[chain])) for chain, start in chain_starts.items()]
    target_segments = [sequences[chain] for chain in chain_starts]
    first_chain = next(iter(chain_starts))
    first_residue = min(n for c, n in observed if c == first_chain)
    last_chain = list(chain_starts)[-1]
    last_residue = max(n for c, n in observed if c == last_chain)
    if ligand_lines:
        template_segments.append("." * len(ligand_keys))
        target_segments.append("." * len(ligand_keys))
        last_chain, last_residue = ligand_chain, max(ligand_identities)
    alignment = (f">P1;template\nstructureX:template:{first_residue}:{first_chain}:{last_residue}:{last_chain}::::\n"
                 + "/".join(template_segments) + "*\n>P1;repaired\nsequence:repaired::::::::\n"
                 + "/".join(target_segments) + "*\n")
    (directory / "alignment.ali").write_text(alignment)
    write_json(directory / "selection.json", {"missing_residues": sorted(missing), "movable_residues": sorted(movable)})
    return RepairInputs(existing, movable, dict(chain_starts), ligand_chain, ligand_start, ligand_identities)


class RepairModel(LoopModel):
    inputs: RepairInputs

    def default_patches(self, aln: Alignment) -> None:
        return None

    def special_patches(self, aln: Alignment) -> None:
        ids = list(self.inputs.chain_starts)
        starts = list(self.inputs.chain_starts.values())
        if self.inputs.ligand_chain is not None:
            ids.append(self.inputs.ligand_chain)
            starts.append(self.inputs.ligand_start)
        self.rename_segments(segment_ids=ids, renumber_residues=starts)

    def select_atoms(self) -> Selection:
        selected = Selection()
        for atom in self.atoms:
            chain, number = atom.residue.chain.name, int(atom.residue.num)
            if chain not in self.inputs.chain_starts:
                continue
            name = {"OT1": "O", "OT2": "OXT"}.get(atom.name, atom.name)
            if atom.residue.name == "ILE" and name == "CD":
                name = "CD1"
            if (chain, number) in self.inputs.movable or (chain, number, name) not in self.inputs.existing:
                selected.add(atom)
        return selected

    def select_loop_atoms(self) -> Selection:
        return self.select_atoms()


def normalized_pdb(source: Path, target: Path, inputs: RepairInputs) -> None:
    lines: list[str] = []
    for line in source.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            name = line[12:16].strip()
            replacement = {"OT1": "O", "OT2": "OXT"}.get(name, name)
            if line[17:20].strip() == "ILE" and name == "CD":
                replacement = "CD1"
            line = line[:12] + f" {replacement:<3}" + line[16:]
            if line[21] == inputs.ligand_chain:
                chain, number, insertion = inputs.ligand_identities[int(line[22:26])]
                line = line[:21] + f"{chain}{number:4d}{insertion}" + line[27:]
            lines.append(line)
        elif line.startswith("TER"):
            lines.append("TER")
    target.write_text("\n".join(lines + ["END"]) + "\n")


def repair(template: Path, fasta: Path, output: Path, chain_starts: Mapping[str, int],
           seed: int = -20260910, models: int = 1, loops: int = 4) -> list[ModelOutput]:
    if seed >= 0 or min(models, loops) < 1:
        raise ValueError("MODELLER needs a negative seed and positive model/loop counts")
    output.mkdir(parents=True, exist_ok=False)
    inputs = prepare_inputs(template, fasta, output, chain_starts)
    previous = Path.cwd()
    try:
        os.chdir(output)
        modeller.log.minimal()
        env = Environ(rand_seed=seed)
        env.io.hetatm = True
        env.io.atom_files_directory = ['.']
        model = RepairModel(env, alnfile="alignment.ali", knowns="template", sequence="repaired",
                            assess_methods=(assess.DOPE,), loop_assess_methods=(assess.DOPE,))
        model.inputs = inputs
        model.starting_model, model.ending_model = 1, models
        model.library_schedule, model.max_var_iterations, model.md_level = autosched.slow, 300, refine.fast
        model.loop.starting_model, model.loop.ending_model = 1, loops
        model.loop.md_level, model.loop.max_var_iterations = refine.fast, 300
        model.final_malign3d = False
        model.make()
        outputs: list[ModelOutput] = []
        for kind, records in (("comparative", model.outputs), ("loop", model.loop.outputs)):
            for record in records:
                filename = str(record["name"])
                failure = str(record["failure"]) if record["failure"] else None
                outputs.append({"filename": filename, "kind": kind, "failed": failure,
                                "dope_score": float(record["DOPE score"]) if "DOPE score" in record else None,
                                "molpdf": float(record["molpdf"]) if "molpdf" in record else None})
                if failure is None:
                    normalized_pdb(Path(filename), Path(filename).with_name(Path(filename).stem + "_auth.pdb"), inputs)
        write_json(Path("results.json"), {"seed": seed, "models": models, "loops_per_model": loops, "outputs": outputs})
    finally:
        os.chdir(previous)
    if not any(r['failed'] is None for r in outputs):
        raise RuntimeError("No model was generated successfully")
    return outputs
