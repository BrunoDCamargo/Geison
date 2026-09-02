from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ScientificCompleteness:
    complete: bool
    missing_evidence: tuple[str, ...]


def assess_pre_ranking_completeness(
    *,
    evaluation_sequence_count: int,
    assay_count: int,
    inclusivity_status: str,
    specificity_status: str,
) -> ScientificCompleteness:
    missing: list[str] = []
    if evaluation_sequence_count == 0:
        missing.append("EMPTY_EVALUATION_SET")
    if assay_count == 0:
        missing.append("NO_ASSAYS")
    if inclusivity_status != "COMPLETE":
        missing.append("INCLUSIVITY_NOT_COMPLETE")
    if specificity_status != "COMPLETE":
        missing.append("SPECIFICITY_NOT_COMPLETE")
    return ScientificCompleteness(not missing, tuple(missing))


def assess_final_completeness(
    pre_ranking: ScientificCompleteness,
    *,
    ranking_status: str,
) -> ScientificCompleteness:
    missing = list(pre_ranking.missing_evidence)
    if ranking_status != "COMPLETE":
        missing.append("RANKING_NOT_COMPLETE")
    return ScientificCompleteness(not missing, tuple(missing))
