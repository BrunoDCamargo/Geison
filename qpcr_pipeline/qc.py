"""Quality control for local target sequences."""

from dataclasses import dataclass
from enum import Enum

from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import EvaluationSet, TargetSequenceSet


class QCStatus(str, Enum):
    ACCEPTED = "ACCEPTED"
    REJECTED = "REJECTED"


@dataclass(frozen=True, slots=True)
class QCRecord:
    sequence_id: str
    status: QCStatus
    reason_codes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class QCResult:
    records: tuple[QCRecord, ...]
    target_sequence_set: TargetSequenceSet
    evaluation_set: EvaluationSet


_VALID_NUCLEOTIDES = frozenset("ACGTRYSWKMBDHVN")


def evaluate_sequences(
    records: tuple[LocalSequenceRecord, ...],
    *,
    min_length: int | None = None,
) -> QCResult:
    qc_records: list[QCRecord] = []
    accepted_ids: list[str] = []

    for record in records:
        invalid_nucleotide = any(base not in _VALID_NUCLEOTIDES for base in record.sequence)

        if invalid_nucleotide:
            qc_records.append(
                QCRecord(
                    sequence_id=record.sequence_id,
                    status=QCStatus.REJECTED,
                    reason_codes=("INVALID_NUCLEOTIDE",),
                )
            )
            continue

        if min_length is not None and len(record.sequence) < min_length:
            qc_records.append(
                QCRecord(
                    sequence_id=record.sequence_id,
                    status=QCStatus.REJECTED,
                    reason_codes=("TOO_SHORT",),
                )
            )
            continue

        qc_records.append(
            QCRecord(
                sequence_id=record.sequence_id,
                status=QCStatus.ACCEPTED,
                reason_codes=(),
            )
        )
        accepted_ids.append(record.sequence_id)

    sequence_ids = tuple(accepted_ids)
    return QCResult(
        records=tuple(qc_records),
        target_sequence_set=TargetSequenceSet(sequence_ids=sequence_ids),
        evaluation_set=EvaluationSet(sequence_ids=sequence_ids),
    )
