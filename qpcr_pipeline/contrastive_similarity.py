"""Focused and deterministic local-similarity interface for contrastive regions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Protocol

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from qpcr_pipeline.local_input import LocalSequenceRecord


@dataclass(frozen=True, slots=True)
class RegionSimilarity:
    sequence_id: str
    similarity: float
    orientation: Literal["forward", "reverse"]


class RegionSimilarityEngine(Protocol):
    def best_match(
        self,
        query: str,
        records: tuple[LocalSequenceRecord, ...],
    ) -> RegionSimilarity | None: ...


class BiopythonLocalSimilarityEngine:
    """Measure best local sequence representation in either query orientation."""

    def __init__(self) -> None:
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -0.5
        self._aligner = aligner

    def best_match(
        self,
        query: str,
        records: tuple[LocalSequenceRecord, ...],
    ) -> RegionSimilarity | None:
        if not isinstance(query, str) or not query.strip():
            raise ValueError("Contrastive similarity query must be non-blank.")
        if not isinstance(records, tuple):
            raise ValueError("Contrastive similarity records must be a tuple.")
        if not records:
            return None

        normalized_query = query.strip().upper()
        reverse_query = str(Seq(normalized_query).reverse_complement())
        candidates: list[RegionSimilarity] = []
        for record in records:
            if not isinstance(record, LocalSequenceRecord):
                raise ValueError(
                    "Contrastive similarity records must contain LocalSequenceRecord values."
                )
            subject = record.sequence.strip().upper()
            if not subject:
                continue
            for orientation, oriented_query in (
                ("forward", normalized_query),
                ("reverse", reverse_query),
            ):
                score = self._aligner.score(oriented_query, subject)
                similarity = max(
                    0.0,
                    min(1.0, score / len(normalized_query)),
                )
                candidates.append(
                    RegionSimilarity(
                        sequence_id=record.sequence_id,
                        similarity=similarity,
                        orientation=orientation,
                    )
                )

        if not candidates:
            return None
        return min(
            candidates,
            key=lambda item: (
                -item.similarity,
                item.sequence_id,
                0 if item.orientation == "forward" else 1,
            ),
        )
