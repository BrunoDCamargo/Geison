from qpcr_pipeline.config import RankingConfig
from qpcr_pipeline.ranking_guard import (
    classify_assays_with_execution_guard,
    rank_assays_with_execution_guard,
)
from ranking_fixtures import make_inclusivity_result, make_primer_result, make_specificity_result


def test_execution_missing_evidence_prevents_in_silico_pass():
    primer = make_primer_result()

    result = classify_assays_with_execution_guard(
        primer,
        make_inclusivity_result(primer),
        make_specificity_result(primer),
        RankingConfig(enabled=True),
        execution_missing_evidence=("EMPTY_EVALUATION_SET",),
    )[0]

    assert result.classification == "REVIEW"
    assert "execution" in result.missing_components
    assert "RUN_EVIDENCE_INCOMPLETE" in {reason.code for reason in result.reasons}


def test_execution_missing_evidence_forces_incomplete_score():
    primer = make_primer_result()

    result = rank_assays_with_execution_guard(
        primer,
        make_inclusivity_result(primer),
        make_specificity_result(primer),
        RankingConfig(enabled=True),
        execution_missing_evidence=("NO_ASSAYS",),
    )[0]

    assert result.score_status == "INCOMPLETE"
    assert result.final_score is None
