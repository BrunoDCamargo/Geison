import json

import pytest

from qpcr_pipeline.run_recording import StructuredEventLogger, sanitize_diagnostic


def test_sanitizer_redacts_secrets_email_and_long_sequences():
    payload = sanitize_diagnostic({
        "api_key": "private-key",
        "email": "researcher@example.test",
        "message": "failed for " + "ACGT" * 40,
        "stage": "input",
    })

    serialized = json.dumps(payload)
    assert "private-key" not in serialized
    assert "researcher@example.test" not in serialized
    assert "ACGTACGTACGT" not in serialized
    assert payload["stage"] == "input"


def test_event_logger_appends_one_sanitized_json_object_per_line(tmp_path):
    path = tmp_path / "run.log.jsonl"
    logger = StructuredEventLogger(
        path,
        run_id="run-1",
        attempt_id="attempt-1",
        clock=lambda: "2026-08-31T00:00:00Z",
    )
    logger.emit("stage_failed", stage="input", message="boom")
    logger.emit("run_failed", status="FAILED")

    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    assert [row["event"] for row in rows] == ["stage_failed", "run_failed"]
    assert all(row["schema_version"] == 1 for row in rows)
    assert all(row["run_id"] == "run-1" for row in rows)


def test_event_logger_rejects_unknown_event_and_drops_unknown_fields(tmp_path):
    path = tmp_path / "run.log.jsonl"
    logger = StructuredEventLogger(path, run_id="r", attempt_id="a")

    with pytest.raises(ValueError):
        logger.emit("unknown_event")

    logger.emit("stage_started", stage="input", action="RUN", api_key="secret")
    row = json.loads(path.read_text(encoding="utf-8"))
    assert "api_key" not in row
