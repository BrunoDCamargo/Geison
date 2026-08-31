# Issue #12 Reproducibility and Execution Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add auditable run manifests, structured sanitized logs, explicit scientific completion states, environment diagnostics, and a side-effect-free dry-run to Geison.

**Architecture:** A new `diagnostics.py` inspects Python, Geison, Git and external tools through an injectable command boundary. A new `run_recording.py` owns scientific completeness, event sanitization, attempt history and atomic run-manifest publication. Existing checkpoint policy is exposed as a shared read-only plan so `run_pipeline()` and `--dry-run` use the same `RUN`/`REUSE`/`FORCED` decisions without duplicating issue #11.

**Tech Stack:** Python 3.10+, standard library dataclasses/JSON/subprocess/pathlib/uuid, pytest, unittest, existing Biopython and PyYAML dependencies, CircleCI.

**Spec:** `docs/superpowers/specs/2026-08-31-issue-12-reproducibility-diagnostics-design.md`

## Global Constraints

- Work on `feature/issue-12-reproducibility-diagnostics`; never implement directly on `main` or `develop`.
- Preserve the stage graph, checkpoint schemas and selective invalidation semantics from issue #11.
- Do not introduce BLAST+ as a runtime dependency or specificity backend; report it as `NOT_USED` and `required: false`.
- Add no package dependency. Use only the existing project dependencies and Python standard library.
- Never write API keys, authorization values, NCBI e-mail values, raw NCBI payloads or complete biological sequences to manifests, JSONL logs or new console errors.
- Timestamps, run IDs, attempt IDs and run-level status never participate in checkpoint fingerprints.
- `doctor` performs no network calls and requires no configuration.
- `--dry-run` performs no NCBI/scientific execution and creates or modifies no files or directories.
- A final `PARTIAL` result can never contain `IN SILICO PASS`.
- Use the normal CircleCI filter for `develop`/`main`. If local execution remains unavailable, request explicit approval before temporarily adding the feature branch, and restore the original filter before review or merge.
- Use `python -m pytest`, not a bare `pytest`, in all verification commands.

---

### Task 1: Typed environment diagnostics

**Files:**
- Create: `qpcr_pipeline/diagnostics.py`
- Create: `tests/test_diagnostics.py`

**Interfaces:**
- Produces: `CommandResult(returncode: int, stdout: str, stderr: str)`.
- Produces: `CommandRunner.run(argv: tuple[str, ...]) -> CommandResult` protocol.
- Produces: `ComponentReport(name, status, required, installed, version, commit=None, dirty=None)`.
- Produces: `EnvironmentReport(python, geison, git, tools)` with `missing_required_tools`.
- Produces: `EnvironmentInspector.inspect(config: PipelineConfig | None = None) -> EnvironmentReport`.

- [ ] **Step 1: Write failing tests for tool requirement and BLAST semantics**

```python
# tests/test_diagnostics.py
from qpcr_pipeline.config import AlignmentConfig, ClusteringConfig, PipelineConfig, PrimerDesignConfig
from qpcr_pipeline.diagnostics import CommandResult, EnvironmentInspector


class FakeRunner:
    def __init__(self, responses):
        self.responses = responses
        self.calls = []

    def run(self, argv):
        self.calls.append(argv)
        return self.responses.get(argv, CommandResult(127, "", "not found"))


def test_inspector_marks_enabled_tools_required_and_blast_not_used(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    config = PipelineConfig(
        target_name="target",
        input_fasta=fasta,
        clustering=ClusteringConfig(enabled=True),
        alignment=AlignmentConfig(enabled=True),
        primer_design=PrimerDesignConfig(enabled=True),
    )
    runner = FakeRunner({
        ("cd-hit-est", "-h"): CommandResult(0, "CD-HIT version 4.8.1", ""),
        ("mafft", "--version"): CommandResult(0, "v7.526", ""),
        ("primer3_core", "--about"): CommandResult(0, "primer3 release 2.6.1", ""),
        ("blastn", "-version"): CommandResult(127, "", "not found"),
    })

    report = EnvironmentInspector(runner=runner).inspect(config)

    assert report.tools["cd-hit-est"].required is True
    assert report.tools["mafft"].required is True
    assert report.tools["primer3_core"].required is True
    assert report.tools["blast+"].status == "NOT_USED"
    assert report.tools["blast+"].required is False
    assert report.tools["blast+"].installed is False
```

- [ ] **Step 2: Run the focused tests and verify RED**

Run: `python -m pytest tests/test_diagnostics.py -q`

Expected: FAIL during import because `qpcr_pipeline.diagnostics` does not exist.

- [ ] **Step 3: Implement the typed reports and injectable inspector**

```python
# qpcr_pipeline/diagnostics.py
from __future__ import annotations

import platform
import subprocess
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from typing import Mapping, Protocol

from qpcr_pipeline.config import PipelineConfig


@dataclass(frozen=True, slots=True)
class CommandResult:
    returncode: int
    stdout: str
    stderr: str


class CommandRunner(Protocol):
    def run(self, argv: tuple[str, ...]) -> CommandResult: ...


class SubprocessCommandRunner:
    def run(self, argv: tuple[str, ...]) -> CommandResult:
        try:
            completed = subprocess.run(argv, capture_output=True, text=True, check=False)
        except OSError as error:
            return CommandResult(127, "", str(error))
        return CommandResult(completed.returncode, completed.stdout, completed.stderr)


@dataclass(frozen=True, slots=True)
class ComponentReport:
    name: str
    status: str
    required: bool
    installed: bool
    version: str | None
    commit: str | None = None
    dirty: bool | None = None


@dataclass(frozen=True, slots=True)
class EnvironmentReport:
    python: ComponentReport
    geison: ComponentReport
    git: ComponentReport
    tools: Mapping[str, ComponentReport]

    @property
    def missing_required_tools(self) -> tuple[str, ...]:
        return tuple(name for name, report in self.tools.items() if report.required and not report.installed)
```

Implement `EnvironmentInspector` with explicit probes:

```python
TOOL_PROBES = {
    "cd-hit-est": ("cd-hit-est", "-h"),
    "mafft": ("mafft", "--version"),
    "primer3_core": ("primer3_core", "--about"),
    "blast+": ("blastn", "-version"),
}

def inspect(self, config=None):
    required = {
        "cd-hit-est": bool(config and config.clustering.enabled),
        "mafft": bool(config and config.alignment.enabled),
        "primer3_core": bool(config and config.primer_design.enabled),
        "blast+": False,
    }
```

Normalize version output by joining stripped non-empty stdout/stderr lines and bounding it to 500 characters. Determine Git commit with `git rev-parse HEAD` and dirty state with `git status --porcelain`; return `UNAVAILABLE` with `null` metadata outside a checkout. Obtain the Geison version with `importlib.metadata.version("geison-qpcr")`; catch `PackageNotFoundError` and report the component as `UNAVAILABLE` with `version=None`.

- [ ] **Step 4: Add absence, Git-unavailable and disabled-tool tests**

```python
def test_inspector_reports_missing_required_tool(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    config = PipelineConfig(
        target_name="target",
        input_fasta=fasta,
        alignment=AlignmentConfig(enabled=True),
    )
    report = EnvironmentInspector(runner=FakeRunner({})).inspect(config)
    assert report.tools["mafft"].required is True
    assert report.tools["mafft"].installed is False
    assert report.missing_required_tools == ("mafft",)


def test_doctor_context_treats_external_tools_as_optional():
    report = EnvironmentInspector(runner=FakeRunner({})).inspect(None)
    assert report.missing_required_tools == ()
    assert all(not item.required for item in report.tools.values())
```

- [ ] **Step 5: Run tests and commit**

Run: `python -m pytest tests/test_diagnostics.py -q`

Expected: PASS.

```bash
git add qpcr_pipeline/diagnostics.py tests/test_diagnostics.py
git commit -m "feat: inspect Geison execution environment"
```

---

### Task 2: `doctor` CLI command

**Files:**
- Modify: `qpcr_pipeline/cli.py:1-65`
- Modify: `tests/test_cli.py:1-143`
- Test: `tests/test_diagnostics.py`

**Interfaces:**
- Consumes: `EnvironmentInspector.inspect(None)` from Task 1.
- Produces: `render_environment_report(report: EnvironmentReport) -> str`.
- Produces: `doctor_exit_code(report: EnvironmentReport) -> int`.

- [ ] **Step 1: Write failing CLI and rendering tests**

```python
# tests/test_diagnostics.py
from qpcr_pipeline.diagnostics import ComponentReport, EnvironmentReport, doctor_exit_code, render_environment_report


def test_doctor_rendering_names_missing_optional_tools_without_failure():
    missing = ComponentReport("mafft", "UNAVAILABLE", False, False, None)
    used = ComponentReport("Python", "USED", True, True, "3.12")
    report = EnvironmentReport(used, used, ComponentReport("Git", "UNAVAILABLE", False, False, None), {"mafft": missing})
    rendered = render_environment_report(report)
    assert "mafft" in rendered
    assert "UNAVAILABLE" in rendered
    assert doctor_exit_code(report) == 0
```

```python
# tests/test_cli.py, inside PipelineCliTests
def test_doctor_command_runs_without_configuration(self):
    result = subprocess.run(
        [self._executable(), "doctor"],
        capture_output=True,
        text=True,
        check=False,
    )
    self.assertEqual(result.returncode, 0, result.stderr)
    self.assertIn("Python", result.stdout)
    self.assertIn("BLAST+", result.stdout)
    self.assertIn("NOT_USED", result.stdout)
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_diagnostics.py tests/test_cli.py -q`

Expected: FAIL because rendering functions and the `doctor` subcommand do not exist.

- [ ] **Step 3: Implement rendering and CLI routing**

Add a `doctor` subparser in `build_parser()` and handle it before the existing run branch:

```python
subparsers.add_parser("doctor", help="Inspect the Geison execution environment")

if args.command == "doctor":
    report = EnvironmentInspector().inspect()
    print(render_environment_report(report))
    return doctor_exit_code(report)
```

Render a fixed four-column table: component, status, required, version. `doctor_exit_code()` returns non-zero only when Python or Geison is not installed/usable, never for absent external tools or Git.

- [ ] **Step 4: Verify GREEN and commit**

Run: `python -m pytest tests/test_diagnostics.py tests/test_cli.py -q`

Expected: PASS.

```bash
git add qpcr_pipeline/diagnostics.py qpcr_pipeline/cli.py tests/test_diagnostics.py tests/test_cli.py
git commit -m "feat: add environment doctor command"
```

---

### Task 3: Scientific completeness and ranking safety

**Files:**
- Create: `qpcr_pipeline/run_recording.py`
- Create: `tests/test_run_completeness.py`
- Modify: `qpcr_pipeline/ranking.py:115-220,362-430`
- Modify: `tests/test_ranking_classification.py:1-150`

**Interfaces:**
- Produces: `ScientificCompleteness(complete: bool, missing_evidence: tuple[str, ...])`.
- Produces: `assess_pre_ranking_completeness(evaluation_sequence_count, assay_count, inclusivity_status, specificity_status) -> ScientificCompleteness`.
- Produces: `assess_final_completeness(pre_ranking, ranking_status) -> ScientificCompleteness`.
- Extends: `classify_assays(..., execution_missing_evidence: tuple[str, ...] = ())`.
- Extends: `evaluate_ranking(..., execution_missing_evidence: tuple[str, ...] = ())`.

- [ ] **Step 1: Write failing completeness tests**

```python
# tests/test_run_completeness.py
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
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_run_completeness.py -q`

Expected: FAIL because completeness functions do not exist.

- [ ] **Step 3: Implement the minimal pure completeness functions**

```python
@dataclass(frozen=True, slots=True)
class ScientificCompleteness:
    complete: bool
    missing_evidence: tuple[str, ...]


def assess_pre_ranking_completeness(*, evaluation_sequence_count, assay_count, inclusivity_status, specificity_status):
    missing = []
    if evaluation_sequence_count == 0:
        missing.append("EMPTY_EVALUATION_SET")
    if assay_count == 0:
        missing.append("NO_ASSAYS")
    if inclusivity_status != "COMPLETE":
        missing.append("INCLUSIVITY_NOT_COMPLETE")
    if specificity_status != "COMPLETE":
        missing.append("SPECIFICITY_NOT_COMPLETE")
    return ScientificCompleteness(not missing, tuple(missing))


def assess_final_completeness(pre_ranking, *, ranking_status):
    missing = list(pre_ranking.missing_evidence)
    if ranking_status != "COMPLETE":
        missing.append("RANKING_NOT_COMPLETE")
    return ScientificCompleteness(not missing, tuple(missing))
```

- [ ] **Step 4: Write the failing ranking guard test**

```python
# tests/test_ranking_classification.py
def test_execution_missing_evidence_prevents_in_silico_pass(self):
    primer = make_primer_result()
    result = classify_assays(
        primer,
        make_inclusivity_result(primer),
        make_specificity_result(primer),
        RankingConfig(enabled=True),
        execution_missing_evidence=("EMPTY_EVALUATION_SET",),
    )[0]
    self.assertEqual(result.classification, "REVIEW")
    self.assertIn("execution", result.missing_components)
    self.assertIn("RUN_EVIDENCE_INCOMPLETE", {reason.code for reason in result.reasons})
```

In `tests/test_ranking_scoring.py`, add the same `execution_missing_evidence` input to a ranked-result fixture and assert its `score_status == "INCOMPLETE"`.

- [ ] **Step 5: Verify RED, then implement the guard before artifact publication**

Run: `python -m pytest tests/test_ranking_classification.py::RankingClassificationTests::test_execution_missing_evidence_prevents_in_silico_pass -q`

Expected: FAIL because `execution_missing_evidence` is not accepted.

Add an optional keyword-only tuple to `classify_assays()` and `evaluate_ranking()`. When non-empty, append one `RUN_EVIDENCE_INCOMPLETE` reason carrying the stable codes, add `execution` to `missing_components`, force `REVIEW` unless an existing rule is already `HIGH_RISK`, and force score status `INCOMPLETE`. Perform the invariant check before `_atomic_write_text()` promotes TSV/JSON/HTML artifacts:

```python
if execution_missing_evidence and any(
    assay.classification == "IN SILICO PASS" for assay in ranked
):
    raise RankingError("Incomplete run evidence cannot produce IN SILICO PASS.")
```

- [ ] **Step 6: Run focused and regression tests, then commit**

Run: `python -m pytest tests/test_run_completeness.py tests/test_ranking_classification.py tests/test_ranking_scoring.py -q`

Expected: PASS.

```bash
git add qpcr_pipeline/run_recording.py qpcr_pipeline/ranking.py tests/test_run_completeness.py tests/test_ranking_classification.py
git commit -m "feat: enforce scientific run completeness"
```

---

### Task 4: Sanitized structured event log

**Files:**
- Modify: `qpcr_pipeline/run_recording.py`
- Create: `tests/test_run_logging.py`

**Interfaces:**
- Produces: `sanitize_diagnostic(value: object, *, field_name: str | None = None) -> object`.
- Produces: `StructuredEventLogger(path: Path, run_id: str, attempt_id: str)`.
- Produces: `StructuredEventLogger.emit(event, level="INFO", **fields) -> None`.

- [ ] **Step 1: Write failing sanitizer tests**

```python
# tests/test_run_logging.py
import json
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
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_run_logging.py -q`

Expected: FAIL because sanitizer/logger do not exist.

- [ ] **Step 3: Implement recursive sanitization with bounded strings**

Use normalized field-name segments and a DNA-like full-string expression:

```python
SENSITIVE_SEGMENTS = frozenset({"token", "api_key", "secret", "password", "authorization", "email"})
DNA_LIKE = re.compile(r"^[ACGTURYSWKMBDHVNacgturyswkmbdhvn\s-]{40,}$")
MAX_DIAGNOSTIC_STRING = 1000


def sanitize_diagnostic(value, *, field_name=None):
    normalized = (field_name or "").lower().replace("-", "_")
    if normalized in SENSITIVE_SEGMENTS or any(normalized.endswith("_" + item) for item in SENSITIVE_SEGMENTS):
        return "[REDACTED]"
    if isinstance(value, str):
        if DNA_LIKE.fullmatch(value):
            return "[SEQUENCE_REDACTED]"
        return value[:MAX_DIAGNOSTIC_STRING]
    if isinstance(value, dict):
        return {str(key): sanitize_diagnostic(item, field_name=str(key)) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize_diagnostic(item) for item in value]
    if value is None or isinstance(value, (bool, int, float)):
        return value
    return f"<{type(value).__name__}>"
```

- [ ] **Step 4: Write and run a failing event append test**

```python
def test_event_logger_appends_one_sanitized_json_object_per_line(tmp_path):
    path = tmp_path / "run.log.jsonl"
    logger = StructuredEventLogger(path, run_id="run-1", attempt_id="attempt-1", clock=lambda: "2026-08-31T00:00:00Z")
    logger.emit("stage_failed", stage="input", api_key="secret", message="boom")
    logger.emit("run_failed", status="FAILED")
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    assert [row["event"] for row in rows] == ["stage_failed", "run_failed"]
    assert rows[0]["api_key"] == "[REDACTED]"
    assert all(row["schema_version"] == 1 for row in rows)
```

Run: `python -m pytest tests/test_run_logging.py -q`

Expected: FAIL because `StructuredEventLogger` is missing.

- [ ] **Step 5: Implement event allowlists and append/flush behavior**

Define `EVENT_FIELDS` for the nine events in the spec. Reject unknown event names and drop/reject fields outside the event allowlist. Serialize compact sorted JSON, append exactly one newline and flush before closing. Include `schema_version`, UTC timestamp, level, event, run ID and attempt ID automatically.

- [ ] **Step 6: Verify and commit**

Run: `python -m pytest tests/test_run_logging.py -q`

Expected: PASS.

```bash
git add qpcr_pipeline/run_recording.py tests/test_run_logging.py
git commit -m "feat: add sanitized structured run events"
```

---

### Task 5: Atomic run manifest and attempt lifecycle

**Files:**
- Modify: `qpcr_pipeline/run_recording.py`
- Create: `tests/test_run_recording.py`

**Interfaces:**
- Produces: `RunRecorder(outdir: Path, *, clock=utc_now, id_factory=uuid4)`.
- Produces: `begin_attempt(target_name, effective_config, execution_policy, environment, plan) -> str`.
- Produces: `stage_started(stage, action)`, `stage_completed(stage, action, checkpoint_path)`, `stage_reused(...)`.
- Produces: `complete(status, scientific_completeness, input_provenance, reference) -> None`.
- Produces: `fail(error, stage: str | None) -> None`.

- [ ] **Step 1: Write failing run identity and attempt-history tests**

```python
# tests/test_run_recording.py
import json
from qpcr_pipeline.run_recording import RunRecorder


def test_resume_retains_run_identity_and_appends_attempt(tmp_path):
    ids = iter(("run-1", "attempt-1", "attempt-2"))
    recorder = RunRecorder(tmp_path, clock=lambda: "2026-08-31T00:00:00Z", id_factory=lambda: next(ids))
    recorder.begin_attempt("target", {}, {"resume": False}, {}, [])
    recorder.fail(RuntimeError("first failure"), stage="alignment")
    recorder.begin_attempt("target", {}, {"resume": True}, {}, [])
    recorder.complete("PARTIAL", {"complete": False, "missing_evidence": ["NO_ASSAYS"]}, {}, {"id": None, "mode": None})
    payload = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert payload["run_id"] == "run-1"
    assert [attempt["attempt_id"] for attempt in payload["attempts"]] == ["attempt-1", "attempt-2"]
    assert [attempt["status"] for attempt in payload["attempts"]] == ["FAILED", "PARTIAL"]
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_run_recording.py -q`

Expected: FAIL because `RunRecorder` does not exist.

- [ ] **Step 3: Implement atomic manifest publication**

Use a JSON dictionary internally with schema validation on load. Write through `tempfile.mkstemp(dir=outdir)` followed by `Path.replace()`. The initial payload must contain all top-level fields from the spec. Store effective configuration and environment only after `sanitize_diagnostic()`.

Use one active attempt at a time. `begin_attempt()` must convert a previously persisted `RUNNING` attempt to:

```json
{"status": "FAILED", "failure": {"code": "INTERRUPTED", "type": "InterruptedRun", "message": "Previous attempt did not finish."}}
```

before appending the new attempt.

- [ ] **Step 4: Add failure-before-stage and publication-failure tests**

```python
def test_failure_before_first_stage_records_null_stage(tmp_path):
    recorder = RunRecorder(tmp_path, id_factory=lambda: "fixed-id")
    recorder.begin_attempt("target", {}, {"resume": False}, {}, [])
    recorder.fail(ValueError("invalid environment"), stage=None)
    payload = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert payload["status"] == "FAILED"
    assert payload["failure"]["stage"] is None
    assert payload["failure"]["type"] == "ValueError"
```

Patch the private atomic writer to raise during `fail()` and assert the original exception remains the one re-raised by the caller-facing helper; diagnostic failure may be attached as a note but never replaces it.

- [ ] **Step 5: Verify and commit**

Run: `python -m pytest tests/test_run_recording.py tests/test_run_logging.py -q`

Expected: PASS.

```bash
git add qpcr_pipeline/run_recording.py tests/test_run_recording.py
git commit -m "feat: persist auditable run attempts"
```

---

### Task 6: Shared checkpoint-aware execution plan

**Files:**
- Modify: `qpcr_pipeline/execution.py:1-174`
- Modify: `qpcr_pipeline/pipeline.py:54-210,280-360`
- Create: `tests/test_execution_plan.py`
- Modify: `tests/test_pipeline_resume.py:1-220`

**Interfaces:**
- Produces: `PipelineExecutionPlan(decisions, reused_results, reused_manifests)`.
- Produces: `plan_pipeline(config, outdir: Path | None, *, execution=None, tool_identity_provider=None) -> PipelineExecutionPlan`.
- Consumes: existing `plan_from_validity()`, `stage_request()`, `CheckpointManager.validate()` and checkpoint codecs.
- `run_pipeline()` consumes the same plan and does not independently recompute stage action policy.

- [ ] **Step 1: Write failing plan parity tests**

```python
# tests/test_execution_plan.py
from qpcr_pipeline.execution import ExecutionPolicy, STAGE_ORDER
from qpcr_pipeline.pipeline import plan_pipeline, run_pipeline


def test_plan_without_outdir_runs_every_stage(minimal_config):
    plan = plan_pipeline(minimal_config, None)
    assert [(item.stage, item.action) for item in plan.decisions] == [
        (stage, "RUN") for stage in STAGE_ORDER
    ]


def test_resume_plan_matches_real_resume(minimal_config, tmp_path):
    outdir = tmp_path / "run"
    run_pipeline(minimal_config, outdir)
    plan = plan_pipeline(minimal_config, outdir, execution=ExecutionPolicy(resume=True))
    assert [item.action for item in plan.decisions] == ["REUSE"] * len(STAGE_ORDER)
```

Create `minimal_config` in the test file using the same two-record FASTA helper as `tests/test_pipeline_resume.py`; do not add a repository-wide fixture for two tests.

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_execution_plan.py -q`

Expected: FAIL because `plan_pipeline` and `PipelineExecutionPlan` do not exist.

- [ ] **Step 3: Implement the read-only planner**

Add:

```python
@dataclass(frozen=True, slots=True)
class PipelineExecutionPlan:
    decisions: tuple[StageDecision, ...]
    reused_results: Mapping[str, object]
    reused_manifests: Mapping[str, CheckpointManifest]
```

For a normal run or `outdir is None`, use `plan_from_validity(policy, {})`. For resume/from-step, validate checkpoints sequentially. Load state/manifests only for valid `REUSE` decisions so downstream `stage_request()` can be calculated. Record stable invalidity reasons from `CheckpointValidation.invalidity`; once a stage is invalid, mark its transitive descendants non-reusable without probing obsolete descendants.

Do not call `mkdir`, `begin`, `complete`, `fail`, `_run_stage`, NCBI or any scientific runner.

- [ ] **Step 4: Refactor `run_pipeline()` to execute the shared decisions**

Replace the three duplicated policy branches with one loop over `plan.decisions`:

```python
for decision in plan.decisions:
    stage = decision.stage
    if decision.action == "REUSE":
        results[stage] = plan.reused_results[stage]
        manifests[stage] = plan.reused_manifests[stage]
        actions.append(StageActionSummary(stage, "REUSE"))
        continue
    request = stage_request(stage, config, manifests, results, effective_tool_provider)
    result, manifest = _run_and_checkpoint_stage(...)
    results[stage] = result
    manifests[stage] = manifest
    actions.append(StageActionSummary(stage, decision.action))
```

Preserve `refresh_online_input` exactly: normal input run and `--from-step input` refresh; resume refreshes online input only for `--force-step input` or an invalid input checkpoint.

- [ ] **Step 5: Verify planning and all issue #11 regression tests**

Run: `python -m pytest tests/test_execution_plan.py tests/test_pipeline_resume.py tests/test_checkpointing.py tests/test_checkpoint_stages.py -q`

Expected: PASS with unchanged action maps and invalidation behavior.

- [ ] **Step 6: Commit**

```bash
git add qpcr_pipeline/execution.py qpcr_pipeline/pipeline.py tests/test_execution_plan.py tests/test_pipeline_resume.py
git commit -m "refactor: share checkpoint execution planning"
```

---

### Task 7: Side-effect-free `--dry-run`

**Files:**
- Modify: `qpcr_pipeline/cli.py:8-63`
- Modify: `qpcr_pipeline/pipeline.py`
- Create: `tests/test_dry_run.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Produces: `DryRunReport(target_name, decisions, environment)`.
- Produces: `dry_run_pipeline(config, outdir: Path | None, *, execution=None, inspector=None, tool_identity_provider=None) -> DryRunReport`.
- Consumes: `plan_pipeline()` from Task 6 and `EnvironmentInspector` from Task 1.

- [ ] **Step 1: Write failing no-side-effect tests**

```python
# tests/test_dry_run.py
from unittest.mock import Mock
from qpcr_pipeline.pipeline import dry_run_pipeline


def test_dry_run_does_not_create_outdir_or_call_scientific_services(minimal_config, tmp_path):
    outdir = tmp_path / "absent"
    inspector = Mock()
    inspector.inspect.return_value = make_environment_report()
    report = dry_run_pipeline(minimal_config, outdir, inspector=inspector)
    assert not outdir.exists()
    assert report.target_name == "target"
    assert all(item.action == "RUN" for item in report.decisions)
    inspector.inspect.assert_called_once_with(minimal_config)
```

Define `make_environment_report()` locally with installed optional components. Add a second test that snapshots every file under an existing output directory before and after dry-run and asserts byte-for-byte equality.

```python
def test_dry_run_does_not_modify_existing_output(minimal_config, tmp_path):
    outdir = tmp_path / "existing"
    outdir.mkdir()
    checkpoint = outdir / "sentinel.bin"
    checkpoint.write_bytes(b"unchanged")
    before = {path.relative_to(outdir): path.read_bytes() for path in outdir.rglob("*") if path.is_file()}

    dry_run_pipeline(minimal_config, outdir, inspector=make_inspector())

    after = {path.relative_to(outdir): path.read_bytes() for path in outdir.rglob("*") if path.is_file()}
    assert after == before
```

Define `make_inspector()` locally so it returns the same deterministic `EnvironmentReport` as the first test.

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_dry_run.py -q`

Expected: FAIL because `dry_run_pipeline` does not exist.

- [ ] **Step 3: Implement the dry-run API and missing-tool failure**

```python
@dataclass(frozen=True, slots=True)
class DryRunReport:
    target_name: str
    decisions: tuple[StageDecision, ...]
    environment: EnvironmentReport


def dry_run_pipeline(config, outdir=None, *, execution=None, inspector=None, tool_identity_provider=None):
    environment = (inspector or EnvironmentInspector()).inspect(config)
    plan = plan_pipeline(config, outdir, execution=execution, tool_identity_provider=tool_identity_provider)
    return DryRunReport(config.target_name, plan.decisions, environment)
```

Do not call `run_pipeline()`. Provide `render_dry_run_report()` in `cli.py`; if `missing_required_tools` is non-empty, print them and return exit code `2`.

- [ ] **Step 4: Write failing CLI contract tests**

```python
# tests/test_cli.py
def test_dry_run_validates_without_creating_output(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = Path(tmpdir) / "planned-run"
        result = self._run(tmpdir, "--dry-run", "--outdir", str(outdir))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("input", result.stdout)
        self.assertIn("RUN", result.stdout)
        self.assertFalse(outdir.exists())
```

Add `--dry-run` to the run parser. Keep existing resume-control validation: resume/from/force require `--outdir`; plain dry-run without outdir reports all stages as `RUN`.

- [ ] **Step 5: Verify and commit**

Run: `python -m pytest tests/test_dry_run.py tests/test_cli.py tests/test_execution_plan.py -q`

Expected: PASS.

```bash
git add qpcr_pipeline/cli.py qpcr_pipeline/pipeline.py tests/test_dry_run.py tests/test_cli.py
git commit -m "feat: preview pipeline runs without side effects"
```

---

### Task 8: Record real pipeline lifecycle and final status

**Files:**
- Modify: `qpcr_pipeline/pipeline.py:54-280,432-525`
- Modify: `qpcr_pipeline/run_recording.py`
- Create: `tests/test_run_manifest.py`
- Modify: `tests/test_minimal_run.py`
- Modify: `tests/test_pipeline_ranking.py`

**Interfaces:**
- Consumes: `RunRecorder`, `assess_pre_ranking_completeness()` and `assess_final_completeness()`.
- Extends: `run_pipeline(..., environment_inspector: EnvironmentInspector | None = None, recorder_factory: Callable[[Path], RunRecorder] | None = None)` for deterministic tests.
- Produces: `RunSummary.status` of `COMPLETED`, `PARTIAL` or `FAILED` consistent with manifest and `run_summary.json`.

- [ ] **Step 1: Write failing completed/partial manifest tests**

```python
# tests/test_run_manifest.py
import json
from dataclasses import replace
from qpcr_pipeline.config import RankingConfig
from qpcr_pipeline.pipeline import run_pipeline


def test_successful_incomplete_fixture_is_partial_and_has_manifest(minimal_config, tmp_path):
    outdir = tmp_path / "run"
    summary = run_pipeline(minimal_config, outdir)
    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    assert summary.status == "PARTIAL"
    assert manifest["status"] == "PARTIAL"
    assert "RANKING_NOT_COMPLETE" in manifest["scientific_completeness"]["missing_evidence"]
    assert (outdir / "run.log.jsonl").is_file()
```

Use the current minimal configuration, whose scientific stages are disabled by default, to prove `PARTIAL`. In the existing `test_enabled_pipeline_publishes_pass_ranking_and_top_recommendation` in `tests/test_pipeline_ranking.py`, capture the return value from `run_pipeline()`, load `run_manifest.json`, and add these exact assertions:

```python
assert summary.status == "COMPLETED"
assert manifest["status"] == "COMPLETED"
assert manifest["scientific_completeness"]["missing_evidence"] == []
```

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/test_run_manifest.py -q`

Expected: FAIL because the manifest is absent and current status is always `COMPLETED`.

- [ ] **Step 3: Start the recorder after output safety validation**

In `run_pipeline()`:

1. keep `_reject_output_inside_frozen_dataset()` before any write;
2. create the output directory;
3. construct environment and the shared plan;
4. call `recorder.begin_attempt()`;
5. emit start/reuse/run/completion events around each decision;
6. catch `BaseException`, call `recorder.fail(error, current_stage)`, then re-raise.

Serialize effective config from the validated dataclass with a recursive converter that handles dataclasses, `Path`, tuples and mappings without reading raw YAML or environment variables.

- [ ] **Step 4: Calculate pre-ranking completeness before `_run_stage("ranking")`**

Immediately before ranking, derive:

```python
pre_ranking = assess_pre_ranking_completeness(
    evaluation_sequence_count=len(qc_result.evaluation_set.sequence_ids),
    assay_count=len(primer_design.assays),
    inclusivity_status=inclusivity.status,
    specificity_status=specificity.status,
)
```

Pass `pre_ranking.missing_evidence` to `evaluate_ranking()`. After ranking, call `assess_final_completeness()`, choose `COMPLETED` when complete and `PARTIAL` otherwise, and use that same value for `RunSummary`, `run_summary.json` and `RunRecorder.complete()`.

- [ ] **Step 5: Verify focused and ranking integration tests**

Run: `python -m pytest tests/test_run_manifest.py tests/test_minimal_run.py tests/test_pipeline_ranking.py tests/test_ranking_artifacts.py -q`

Expected: PASS. Update old assertions that expected `COMPLETED` from deliberately disabled scientific configurations to expect `PARTIAL`; do not change assertions for fully enabled evidence-complete fixtures.

- [ ] **Step 6: Commit**

```bash
git add qpcr_pipeline/pipeline.py qpcr_pipeline/run_recording.py tests/test_run_manifest.py tests/test_minimal_run.py tests/test_pipeline_ranking.py
git commit -m "feat: record pipeline lifecycle and final status"
```

---

### Task 9: Local, NCBI and reference provenance

**Files:**
- Modify: `qpcr_pipeline/run_recording.py`
- Modify: `qpcr_pipeline/pipeline.py:220-280,526-550`
- Modify: `tests/test_run_manifest.py`
- Modify: `tests/test_ncbi_acquisition.py`

**Interfaces:**
- Produces: `effective_config_payload(config: PipelineConfig) -> dict[str, object]`.
- Produces: `build_input_provenance(config, outdir, qc_result, input_manifest) -> dict[str, object]`.
- Produces: `build_reference_provenance(alignment_result) -> dict[str, str | None]`.
- Consumes: existing `ncbi_dataset_manifest.json`, checkpoint input manifest and alignment result.

- [ ] **Step 1: Write failing local provenance test**

```python
def test_local_run_records_effective_config_hash_counts_and_reference(minimal_config, tmp_path):
    outdir = tmp_path / "run"
    run_pipeline(minimal_config, outdir)
    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["effective_config"]["target_name"] == "target"
    assert manifest["input_provenance"]["kind"] == "fasta"
    assert manifest["input_provenance"]["source_sha256"].startswith("sha256:")
    assert manifest["input_provenance"]["accepted_count"] == 2
    serialized = json.dumps(manifest)
    assert "ACGTACGTACGT" not in serialized
    assert "records" not in manifest["input_provenance"]
    assert set(manifest["reference"]) == {"id", "mode"}
```

- [ ] **Step 2: Verify RED, then implement local/reference builders**

Run: `python -m pytest tests/test_run_manifest.py::test_local_run_records_effective_config_hash_counts_and_reference -q`

Expected: FAIL because provenance fields are empty.

Use the input checkpoint request/manifest identity already calculated by issue #11 for the local source SHA-256; do not hash a second different representation. Record only configured path, format, source identity and QC counts. Record `alignment.reference_id` and `alignment.reference_mode`, or `null` values if unavailable.

- [ ] **Step 3: Write failing frozen/online NCBI provenance tests**

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import NcbiInputConfig, PipelineConfig
from qpcr_pipeline.ncbi_acquisition import NcbiFetchedRecord


class ProvenanceNcbiClient:
    def resolve_query(self, query):
        raise AssertionError("accession mode must not resolve a query")

    def fetch_records(self, identifiers, *, identifier_kind):
        assert identifier_kind == "accession"
        return tuple(
            NcbiFetchedRecord(
                request_id=identifier,
                record=SeqRecord(Seq("ACGTACGTACGT"), id=identifier, name=identifier),
            )
            for identifier in identifiers
        )


def test_ncbi_run_records_request_and_resolved_accession_versions(tmp_path):
    input_config = NcbiInputConfig(accessions=("NC_000001.1",))
    config = PipelineConfig(target_name="target", input_ncbi=input_config)
    outdir = tmp_path / "run"
    run_pipeline(config, outdir, ncbi_client=ProvenanceNcbiClient())
    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    provenance = manifest["input_provenance"]
    assert provenance["kind"] == "ncbi"
    assert provenance["mode"] == "accessions"
    assert provenance["requested_accessions"] == ["NC_000001.1"]
    assert provenance["resolved_accession_versions"]
    serialized = json.dumps(provenance)
    assert "NCBI_API_KEY" not in serialized
    assert "private@example" not in serialized
```

Keep this self-contained fake in `tests/test_run_manifest.py`; do not import one test module from another.

- [ ] **Step 4: Verify RED, then implement NCBI manifest projection**

Read the already copied `<outdir>/ncbi_dataset_manifest.json` after the input stage. Validate only the fields needed for provenance, then project:

- mode (`query`, `accessions`, `frozen_dataset`);
- configured query or requested accessions;
- resolved accession versions from `resolved_entries`;
- source/dataset identity from the checkpoint or existing consolidated manifest hashes.

Never copy completed batches, raw records, request headers, e-mail or API key values.

- [ ] **Step 5: Run focused tests and commit**

Run: `python -m pytest tests/test_run_manifest.py tests/test_ncbi_acquisition.py -q`

Expected: PASS.

```bash
git add qpcr_pipeline/run_recording.py qpcr_pipeline/pipeline.py tests/test_run_manifest.py tests/test_ncbi_acquisition.py
git commit -m "feat: preserve run input and reference provenance"
```

---

### Task 10: Failure, interruption, resume and full verification

**Files:**
- Modify: `tests/test_run_manifest.py`
- Modify: `tests/test_pipeline_resume.py`
- Modify: `README.md`

**Interfaces:**
- Verifies: `RunRecorder.fail()` around real pipeline failures.
- Verifies: one run ID with append-only attempts and log events across resume.
- Documents: `doctor`, `--dry-run`, final statuses and diagnostic artifacts.

- [ ] **Step 1: Write failing mid-stage failure and resume test**

```python
def test_failed_attempt_is_preserved_when_resume_completes(minimal_config, tmp_path, monkeypatch):
    outdir = tmp_path / "run"
    original = pipeline_module._run_stage
    failed_once = {"value": False}

    def interrupt(stage, *args, **kwargs):
        if stage == "alignment" and not failed_once["value"]:
            failed_once["value"] = True
            raise RuntimeError("alignment failed for " + "ACGT" * 50)
        return original(stage, *args, **kwargs)

    monkeypatch.setattr(pipeline_module, "_run_stage", interrupt)
    with pytest.raises(RuntimeError):
        run_pipeline(minimal_config, outdir)
    failed = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    run_id = failed["run_id"]
    assert failed["status"] == "FAILED"
    assert failed["failure"]["stage"] == "alignment"
    assert "ACGTACGTACGT" not in json.dumps(failed)

    resumed = run_pipeline(minimal_config, outdir, execution=ExecutionPolicy(resume=True))
    final = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    assert final["run_id"] == run_id
    assert len(final["attempts"]) == 2
    assert final["attempts"][0]["status"] == "FAILED"
    assert final["attempts"][1]["status"] == resumed.status
```

- [ ] **Step 2: Verify RED, then correct lifecycle gaps only**

Run: `python -m pytest tests/test_run_manifest.py::test_failed_attempt_is_preserved_when_resume_completes -q`

Expected: FAIL at the first lifecycle assertion that Tasks 5/8 did not yet satisfy. Fix recorder/pipeline orchestration without altering checkpoint reuse policy.

- [ ] **Step 3: Add documentation**

Add README sections with these exact user contracts:

```text
qpcr-pipeline doctor
qpcr-pipeline run config.yaml --dry-run [--outdir run]
```

Document:

- `run_manifest.json` and `run.log.jsonl`;
- meanings of `COMPLETED`, `PARTIAL`, `FAILED`;
- that `PARTIAL` cannot contain `IN SILICO PASS`;
- that BLAST+ is reported as `NOT_USED` and is not required;
- that dry-run does not query NCBI or write artifacts;
- that outputs are in silico evidence and do not replace experimental validation.

- [ ] **Step 4: Run the complete unit suite**

Run: `python -m pytest -q`

Expected: PASS with zero failures and no unexpected warnings.

- [ ] **Step 5: Run integration tests when external tools are available**

Run: `python -m unittest discover -s integration_tests -q`

Expected: PASS when CD-HIT, MAFFT and Primer3 are installed. If they are unavailable locally, record that limitation and require CircleCI on `main` before claiming integration completion.

- [ ] **Step 6: Check secrets, sequences and unintended CircleCI changes**

Run:

```bash
git diff --check
git diff -- .circleci/config.yml
rg -n "NCBI_API_KEY|private-key|researcher@example|ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" qpcr_pipeline tests README.md
```

Expected: `git diff --check` succeeds; `.circleci/config.yml` has no permanent feature-branch filter; `rg` finds only explicit redaction test fixtures and no leaked runtime artifacts.

- [ ] **Step 7: Commit final tests and documentation**

```bash
git add README.md tests/test_run_manifest.py tests/test_pipeline_resume.py qpcr_pipeline
git commit -m "docs: explain reproducible run diagnostics"
```

- [ ] **Step 8: Review the complete branch before integration**

Run:

```bash
git status --short --branch
git log --oneline main..HEAD
git diff --stat main...HEAD
git diff main...HEAD -- qpcr_pipeline tests README.md .circleci/config.yml
```

Expected: only issue #12 implementation, tests, spec, plan and documentation are present. No unrelated refactor, dependency addition, BLAST integration or permanent CI-filter change exists.

- [ ] **Step 9: Publish and verify CI**

Push the feature branch and create a PR targeting `develop`. If the normal filter does not run the feature branch and local verification was complete, review may proceed without spending CircleCI credits; integration into `develop` must produce a green normal suite. Merge `develop` to `main` only after review, then require the `main` CircleCI job to pass the normal suite plus real CD-HIT, MAFFT and Primer3 integration tests.

- [ ] **Step 10: Close issue #12 only after main verification**

Mark every acceptance criterion in issue #12 complete only when the verified `main` commit contains the feature and its CircleCI status is success. Link the final `main` SHA and CircleCI build in the issue completion note.

