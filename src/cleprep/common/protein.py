"""Protein residue names and coordinate component checks."""

from collections.abc import Iterable

from Bio.Data.PDBData import protein_letters_3to1

AA_ONE: dict[str, str] = dict(protein_letters_3to1)
PROTEIN_NAMES: frozenset[str] = frozenset(AA_ONE)
AA_ONE.update({"HSD": "H", "HSE": "H", "HSP": "H"})

CHOLESTEROL_NAMES: frozenset[str] = frozenset({"CLR", "CHL", "CHL1"})


def require_no_cholesterol(records: Iterable[str], source: str) -> None:
    """Reject pocket cholesterol in production PDB/PDBQT coordinate records."""
    residues = {(line[21:22], line[22:26].strip(), line[17:21].strip())
                for line in records if line.startswith(("ATOM  ", "HETATM"))
                and line[17:21].strip() in CHOLESTEROL_NAMES}
    if residues:
        raise ValueError(f"Pocket cholesterol remains in {source}: {sorted(residues)}. "
                         "Use a cholesterol-free receptor PDB/PDBQT pair before docking or export.")
