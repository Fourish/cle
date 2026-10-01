"""JSON, FASTA and transform-file I/O without project assumptions."""

from __future__ import annotations

import json
from io import StringIO
from pathlib import Path
from typing import Mapping, TypeAlias, TypedDict, cast

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .geometry import Array, checked_transform

Json: TypeAlias = None | bool | int | float | str | list["Json"] | dict[str, "Json"]


class TransformData(TypedDict):
    rotation: list[list[float]]
    translation: list[float]


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=True, indent=2, allow_nan=False) + "\n", encoding="ascii")


def read_fasta(path: Path) -> dict[str, str]:
    """Read plain FASTA IDs and uppercase sequences; reject empty/duplicate records."""
    with StringIO(path.read_text().lstrip()) as handle:
        records: dict[str, SeqRecord] = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    if not records or any(not name or not record.seq for name, record in records.items()):
        raise ValueError("FASTA requires nonempty IDs and sequences")
    return {name: str(record.seq).upper() for name, record in records.items()}


def write_fasta(path: Path, sequences: Mapping[str, str],
                descriptions: Mapping[str, str] | None = None) -> None:
    """Write one sequence line per record, preserving IDs, order and descriptions."""
    if not sequences or any(not name or any(c.isspace() for c in name) or not sequence
                            or any(c.isspace() for c in sequence) for name, sequence in sequences.items()):
        raise ValueError("FASTA requires nonempty IDs and sequences without whitespace")
    records = [SeqRecord(Seq(sequence), id=name, description=(descriptions or {}).get(name, ""))
               for name, sequence in sequences.items()]
    with path.open("w", encoding="ascii") as handle:
        SeqIO.write(records, handle, "fasta-2line")


def read_transform(path: Path) -> tuple[Array, Array]:
    data = cast(TransformData, json.loads(path.read_text()))
    return checked_transform(data["rotation"], data["translation"])
