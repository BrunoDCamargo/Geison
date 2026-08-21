"""Load local sequence files into normalized records."""

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass(frozen=True, slots=True)
class LocalSequenceRecord:
    sequence_id: str
    sequence: str


def load_fasta(path: str | Path) -> tuple[LocalSequenceRecord, ...]:
    fasta_path = Path(path)
    return tuple(
        LocalSequenceRecord(sequence_id=record.id, sequence=str(record.seq).upper())
        for record in SeqIO.parse(fasta_path, "fasta")
    )
