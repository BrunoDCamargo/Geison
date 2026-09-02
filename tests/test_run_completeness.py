from qpcr_pipeline.run_recording import assess_final_completeness, assess_pre_ranking_completeness


def test_complete_evidence_is_completed():
    pre = assess_pre_ranking_completeness(
        evaluation_sequence_count=3,
        assay_count=2,
        inclusivity_status="COMPLETE",
        specificity_status="COMPLETE",
    )
    final = assess_final_completeness(pre, ranking_status="COMPLETE")

    assert final.complete is True
    assert final.missing_evidence == ()


def test_missing_evidence_codes_are_stable_and_ordered():
    pre = assess_pre_ranking_completeness(
        evaluation_sequence_count=0,
        assay_count=0,
        inclusivity_status="SKIPPED",
        specificity_status="SKIPPED",
    )
    final = assess_final_completeness(pre, ranking_status="SKIPPED")

    assert final.missing_evidence == (
        "EMPTY_EVALUATION_SET",
        "NO_ASSAYS",
        "INCLUSIVITY_NOT_COMPLETE",
        "SPECIFICITY_NOT_COMPLETE",
        "RANKING_NOT_COMPLETE",
    )
