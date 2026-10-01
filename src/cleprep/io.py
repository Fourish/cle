"""Compatibility imports; maintain implementations in cleprep.common.io."""

from .common.io import Json, TransformData, read_fasta, read_transform, write_fasta, write_json

__all__ = ["Json", "TransformData", "read_fasta", "read_transform", "write_fasta", "write_json"]
