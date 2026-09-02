import json

from qpcr_pipeline.run_recording import RunRecorder


def test_begin_attempt_logs_environment_and_plan_creation(tmp_path):
    ids = iter(("run-1", "attempt-1"))
    recorder = RunRecorder(
        tmp_path,
        clock=lambda: "2026-08-31T00:00:00Z",
        id_factory=lambda: next(ids),
    )

    recorder.begin_attempt(
        "target",
        {},
        {"resume": False},
        {"python": {"version": "3.12"}},
        [{"stage": "input", "action": "RUN"}],
    )

    events = [
        json.loads(line)["event"]
        for line in (tmp_path / "run.log.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    assert events == ["run_started", "environment_inspected", "plan_created"]


def test_attempt_history_preserves_execution_policy_and_plan(tmp_path):
    ids = iter(("run-1", "attempt-1", "attempt-2"))
    recorder = RunRecorder(
        tmp_path,
        clock=lambda: "2026-08-31T00:00:00Z",
        id_factory=lambda: next(ids),
    )

    recorder.begin_attempt(
        "target",
        {},
        {"resume": False, "from_step": None},
        {},
        [{"stage": "input", "action": "RUN"}],
    )
    recorder.fail(RuntimeError("first failure"), stage="input")
    recorder.begin_attempt(
        "target",
        {},
        {"resume": True, "from_step": None},
        {},
        [{"stage": "input", "action": "REUSE"}],
    )

    payload = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    first, second = payload["attempts"]
    assert first["execution_policy"] == {"resume": False, "from_step": None}
    assert first["plan"] == [{"stage": "input", "action": "RUN"}]
    assert second["execution_policy"] == {"resume": True, "from_step": None}
    assert second["plan"] == [{"stage": "input", "action": "REUSE"}]


def test_resume_retains_run_identity_and_appends_attempt(tmp_path):
    ids = iter(("run-1", "attempt-1", "attempt-2"))
    recorder = RunRecorder(
        tmp_path,
        clock=lambda: "2026-08-31T00:00:00Z",
        id_factory=lambda: next(ids),
    )

    recorder.begin_attempt("target", {}, {"resume": False}, {}, [])
    recorder.fail(RuntimeError("first failure"), stage="alignment")
    recorder.begin_attempt("target", {}, {"resume": True}, {}, [])
    recorder.complete(
        "PARTIAL",
        {"complete": False, "missing_evidence": ["NO_ASSAYS"]},
        {},
        {"id": None, "mode": None},
    )

    payload = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert payload["run_id"] == "run-1"
    assert [attempt["attempt_id"] for attempt in payload["attempts"]] == [
        "attempt-1",
        "attempt-2",
    ]
    assert [attempt["status"] for attempt in payload["attempts"]] == [
        "FAILED",
        "PARTIAL",
    ]


def test_failure_before_first_stage_records_null_stage(tmp_path):
    ids = iter(("run-1", "attempt-1"))
    recorder = RunRecorder(tmp_path, id_factory=lambda: next(ids))
    recorder.begin_attempt("target", {}, {"resume": False}, {}, [])

    recorder.fail(ValueError("invalid environment"), stage=None)

    payload = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert payload["status"] == "FAILED"
    assert payload["failure"]["stage"] is None
    assert payload["failure"]["type"] == "ValueError"


def test_new_attempt_marks_stale_running_attempt_interrupted(tmp_path):
    ids = iter(("run-1", "attempt-1", "attempt-2"))
    recorder = RunRecorder(tmp_path, id_factory=lambda: next(ids))
    recorder.begin_attempt("target", {}, {"resume": False}, {}, [])

    recorder.begin_attempt("target", {}, {"resume": True}, {}, [])

    payload = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    first = payload["attempts"][0]
    assert first["status"] == "FAILED"
    assert first["failure"]["code"] == "INTERRUPTED"
