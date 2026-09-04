import pytest

from qpcr_pipeline.contrastive_similarity import BiopythonLocalSimilarityEngine
from qpcr_pipeline.local_input import LocalSequenceRecord


def _record(sequence_id: str, sequence: str) -> LocalSequenceRecord:
    return LocalSequenceRecord(sequence_id=sequence_id, sequence=sequence)


def test_identical_query_scores_one_and_divergent_scores_lower():
    engine = BiopythonLocalSimilarityEngine()
    identical = engine.best_match(
        "ACGTACGT",
        (_record("same", "TTACGTACGTGG"),),
    )
    divergent = engine.best_match(
        "ACGTACGT",
        (_record("different", "TTTTTTTTTTTT"),),
    )
    assert identical is not None and identical.similarity == 1.0
    assert divergent is not None and divergent.similarity < identical.similarity


def test_reverse_complement_can_win():
    engine = BiopythonLocalSimilarityEngine()
    result = engine.best_match(
        "AGTC",
        (_record("reverse-hit", "TTGACTAA"),),
    )
    assert result is not None
    assert result.similarity == 1.0
    assert result.orientation == "reverse"


def test_blank_query_is_rejected_and_empty_dataset_returns_none():
    engine = BiopythonLocalSimilarityEngine()
    with pytest.raises(ValueError, match="non-blank"):
        engine.best_match("   ", (_record("x", "ACGT"),))
    assert engine.best_match("ACGT", ()) is None


def test_ties_choose_smaller_sequence_id():
    engine = BiopythonLocalSimilarityEngine()
    result = engine.best_match(
        "AGTC",
        (_record("z", "AGTC"), _record("a", "AGTC")),
    )
    assert result is not None
    assert result.sequence_id == "a"
    assert result.orientation == "forward"


def test_orientation_tie_prefers_forward():
    engine = BiopythonLocalSimilarityEngine()
    result = engine.best_match(
        "ACGT",
        (_record("palindrome", "TTACGTAA"),),
    )
    assert result is not None
    assert result.similarity == 1.0
    assert result.orientation == "forward"
