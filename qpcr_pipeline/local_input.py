"""Load local sequence files into normalized records."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Mapping

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass(frozen=True, slots=True)
class LocalSequenceRecord:
    sequence_id: str
    sequence: str
    metadata: Mapping[str, object] = field(default_factory=dict)


def _normalize(record: SeqRecord) -> LocalSequenceRecord:
    return LocalSequenceRecord(
        sequence_id=record.id,
        sequence=str(record.seq).upper(),
        metadata={
            "name": record.name,
            "description": record.description,
            "annotations": dict(record.annotations),
            "dbxrefs": tuple(record.dbxrefs),
            "features": tuple(record.features),
        },
    )


def load_fasta(path: str | Path) -> tuple[LocalSequenceRecord, ...]:
    fasta_path = Path(path)
    with fasta_path.open(encoding="utf-8") as handle:
        return tuple(_normalize(record) for record in SeqIO.parse(handle, "fasta"))


def load_genbank(path: str | Path) -> tuple[LocalSequenceRecord, ...]:
    genbank_path = Path(path)
    with genbank_path.open(encoding="utf-8") as handle:
        return tuple(_normalize(record) for record in SeqIO.parse(handle, "genbank"))


def load_local_sequences(
    path: str | Path,
    file_format: Literal["fasta", "genbank"],
) -> tuple[LocalSequenceRecord, ...]:
    if file_format == "fasta":
        return load_fasta(path)
    if file_format == "genbank":
        return load_genbank(path)
    raise ValueError(f"Unsupported local sequence format: {file_format}")
