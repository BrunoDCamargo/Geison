from __future__ import annotations

import json
import re
import tempfile
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable
from uuid import uuid4


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


SENSITIVE_SEGMENTS = frozenset(
    {"token", "api_key", "secret", "password", "authorization", "email"}
)
SEQUENCE_RUN = re.compile(r"[ACGTURYSWKMBDHVNacgturyswkmbdhvn-]{40,}")
MAX_DIAGNOSTIC_STRING = 1000

EVENT_FIELDS: dict[str, frozenset[str]] = {
    "run_started": frozenset({"target_name", "status"}),
    "environment_inspected": frozenset(),
    "plan_created": frozenset(),
    "stage_started": frozenset({"stage", "action"}),
    "stage_completed": frozenset({"stage", "action", "checkpoint_path"}),
    "stage_reused": frozenset({"stage", "action", "checkpoint_path"}),
    "stage_failed": frozenset({"stage", "message", "error_type"}),
    "run_completed": frozenset({"status", "missing_evidence"}),
    "run_action_required": frozenset({"status", "code", "artifact"}),
    "run_failed": frozenset({"status", "stage", "message", "error_type"}),
}


def _normalized_field_name(field_name: str | None) -> str:
    return (field_name or "").lower().replace("-", "_")


def _is_sensitive_field(field_name: str | None) -> bool:
    normalized = _normalized_field_name(field_name)
    return normalized in SENSITIVE_SEGMENTS or any(
        normalized.endswith("_" + item) for item in SENSITIVE_SEGMENTS
    )


def sanitize_diagnostic(value: object, *, field_name: str | None = None) -> object:
    if _is_sensitive_field(field_name):
        return "[REDACTED]"
    if is_dataclass(value) and not isinstance(value, type):
        return sanitize_diagnostic(asdict(value), field_name=field_name)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, str):
        sanitized = SEQUENCE_RUN.sub("[SEQUENCE_REDACTED]", value)
        return sanitized[:MAX_DIAGNOSTIC_STRING]
    if isinstance(value, dict):
        return {
            str(key): sanitize_diagnostic(item, field_name=str(key))
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [sanitize_diagnostic(item) for item in value]
    if value is None or isinstance(value, (bool, int, float)):
        return value
    return f"<{type(value).__name__}>"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


class StructuredEventLogger:
    def __init__(
        self,
        path: Path,
        *,
        run_id: str,
        attempt_id: str,
        clock: Callable[[], str] = utc_now,
    ) -> None:
        self.path = Path(path)
        self.run_id = run_id
        self.attempt_id = attempt_id
        self.clock = clock

    def emit(self, event: str, level: str = "INFO", **fields: object) -> None:
        if event not in EVENT_FIELDS:
            raise ValueError(f"Unknown structured run event: {event}")
        allowed = EVENT_FIELDS[event]
        payload: dict[str, object] = {
            "schema_version": 1,
            "timestamp": self.clock(),
            "level": str(level)[:20],
            "event": event,
            "run_id": self.run_id,
            "attempt_id": self.attempt_id,
        }
        for name, value in fields.items():
            if name not in allowed:
                continue
            payload[name] = sanitize_diagnostic(value, field_name=name)

        self.path.parent.mkdir(parents=True, exist_ok=True)
        line = json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
        with self.path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(line)
            handle.flush()


def _atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        text=True,
    )
    temporary = Path(temporary_name)
    try:
        with open(descriptor, "w", encoding="utf-8", newline="\n", closefd=True) as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _completeness_payload(value: ScientificCompleteness | dict[str, object]) -> dict[str, object]:
    if isinstance(value, ScientificCompleteness):
        return {
            "complete": value.complete,
            "missing_evidence": list(value.missing_evidence),
        }
    sanitized = sanitize_diagnostic(value)
    if not isinstance(sanitized, dict):
        raise ValueError("Scientific completeness must be a mapping.")
    return sanitized


class RunRecorder:
    def __init__(
        self,
        outdir: Path,
        *,
        clock: Callable[[], str] = utc_now,
        id_factory: Callable[[], object] = uuid4,
    ) -> None:
        self.outdir = Path(outdir)
        self.manifest_path = self.outdir / "run_manifest.json"
        self.log_path = self.outdir / "run.log.jsonl"
        self.clock = clock
        self.id_factory = id_factory
        self._payload: dict[str, object] | None = None
        self._active_attempt_id: str | None = None
        self._logger: StructuredEventLogger | None = None

    def _new_id(self) -> str:
        return str(self.id_factory())

    def _load_existing(self) -> dict[str, object] | None:
        if not self.manifest_path.exists():
            return None
        payload = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict) or payload.get("schema_version") != 1:
            raise ValueError("Unsupported or invalid run manifest.")
        if not isinstance(payload.get("attempts"), list):
            raise ValueError("Run manifest attempts must be a list.")
        return payload

    def _write(self) -> None:
        if self._payload is None:
            raise RuntimeError("No active run manifest.")
        _atomic_write_json(self.manifest_path, self._payload)

    def _active_attempt(self) -> dict[str, object]:
        if self._payload is None or self._active_attempt_id is None:
            raise RuntimeError("No active run attempt.")
        attempts = self._payload["attempts"]
        assert isinstance(attempts, list)
        for attempt in reversed(attempts):
            if isinstance(attempt, dict) and attempt.get("attempt_id") == self._active_attempt_id:
                return attempt
        raise RuntimeError("Active run attempt is missing from the manifest.")

    def begin_attempt(
        self,
        target_name: str,
        effective_config: object,
        execution_policy: object,
        environment: object,
        plan: object,
    ) -> str:
        payload = self._load_existing()
        now = self.clock()
        if payload is None:
            payload = {
                "schema_version": 1,
                "run_id": self._new_id(),
                "target_name": target_name,
                "status": "RUNNING",
                "created_at": now,
                "updated_at": now,
                "effective_config": sanitize_diagnostic(effective_config),
                "execution_policy": sanitize_diagnostic(execution_policy),
                "environment": sanitize_diagnostic(environment),
                "plan": sanitize_diagnostic(plan),
                "attempts": [],
                "scientific_completeness": None,
                "input_provenance": {},
                "reference": {"id": None, "mode": None},
                "action_required": None,
                "panel_provenance": {},
                "failure": None,
            }
        else:
            payload.setdefault("action_required", None)
            payload.setdefault("panel_provenance", {})
            payload["action_required"] = None
            attempts = payload["attempts"]
            assert isinstance(attempts, list)
            for attempt in attempts:
                if isinstance(attempt, dict) and attempt.get("status") == "RUNNING":
                    attempt["status"] = "FAILED"
                    attempt["finished_at"] = now
                    attempt["failure"] = {
                        "code": "INTERRUPTED",
                        "type": "InterruptedRun",
                        "message": "Previous attempt did not finish.",
                        "stage": None,
                    }
            payload["target_name"] = target_name
            payload["effective_config"] = sanitize_diagnostic(effective_config)
            payload["execution_policy"] = sanitize_diagnostic(execution_policy)
            payload["environment"] = sanitize_diagnostic(environment)
            payload["plan"] = sanitize_diagnostic(plan)
            payload["failure"] = None

        attempt_id = self._new_id()
        attempts = payload["attempts"]
        assert isinstance(attempts, list)
        attempts.append(
            {
                "attempt_id": attempt_id,
                "execution_policy": sanitize_diagnostic(execution_policy),
                "plan": sanitize_diagnostic(plan),
                "status": "RUNNING",
                "started_at": now,
                "finished_at": None,
                "stages": [],
                "failure": None,
            }
        )
        payload["status"] = "RUNNING"
        payload["updated_at"] = now
        self._payload = payload
        self._active_attempt_id = attempt_id
        self._logger = StructuredEventLogger(
            self.log_path,
            run_id=str(payload["run_id"]),
            attempt_id=attempt_id,
            clock=self.clock,
        )
        self._write()
        self._logger.emit("run_started", target_name=target_name, status="RUNNING")
        self._logger.emit("environment_inspected")
        self._logger.emit("plan_created")
        return attempt_id

    def stage_started(self, stage: str, action: str) -> None:
        attempt = self._active_attempt()
        stages = attempt["stages"]
        assert isinstance(stages, list)
        stages.append(
            {
                "stage": stage,
                "action": action,
                "status": "RUNNING",
                "started_at": self.clock(),
                "finished_at": None,
                "checkpoint_path": None,
            }
        )
        self._write()
        assert self._logger is not None
        self._logger.emit("stage_started", stage=stage, action=action)

    def _finish_stage(self, stage: str, action: str, checkpoint_path: object, event: str) -> None:
        attempt = self._active_attempt()
        stages = attempt["stages"]
        assert isinstance(stages, list)
        matching = [item for item in stages if isinstance(item, dict) and item.get("stage") == stage]
        if matching:
            row = matching[-1]
        else:
            row = {
                "stage": stage,
                "action": action,
                "started_at": self.clock(),
            }
            stages.append(row)
        row["status"] = "COMPLETED"
        row["finished_at"] = self.clock()
        row["checkpoint_path"] = sanitize_diagnostic(checkpoint_path)
        self._write()
        assert self._logger is not None
        self._logger.emit(event, stage=stage, action=action, checkpoint_path=checkpoint_path)

    def stage_completed(self, stage: str, action: str, checkpoint_path: object = None) -> None:
        self._finish_stage(stage, action, checkpoint_path, "stage_completed")

    def stage_reused(self, stage: str, action: str = "REUSE", checkpoint_path: object = None) -> None:
        self._finish_stage(stage, action, checkpoint_path, "stage_reused")

    def complete(
        self,
        status: str,
        scientific_completeness: ScientificCompleteness | dict[str, object],
        *,
        input_provenance: object,
        reference: object,
        panel_provenance: object,
    ) -> None:
        if status not in {"COMPLETED", "PARTIAL"}:
            raise ValueError("Successful run status must be COMPLETED or PARTIAL.")
        if self._payload is None:
            raise RuntimeError("No active run manifest.")
        now = self.clock()
        attempt = self._active_attempt()
        attempt["status"] = status
        attempt["finished_at"] = now
        attempt["failure"] = None
        self._payload["status"] = status
        self._payload["updated_at"] = now
        self._payload["scientific_completeness"] = _completeness_payload(
            scientific_completeness
        )
        self._payload["input_provenance"] = sanitize_diagnostic(input_provenance)
        self._payload["reference"] = sanitize_diagnostic(reference)
        self._payload["panel_provenance"] = sanitize_diagnostic(panel_provenance)
        self._payload["action_required"] = None
        self._payload["failure"] = None
        self._write()
        assert self._logger is not None
        completeness = self._payload["scientific_completeness"]
        assert isinstance(completeness, dict)
        self._logger.emit(
            "run_completed",
            status=status,
            missing_evidence=completeness.get("missing_evidence", []),
        )

    def action_required(self, code: str, artifact: Path) -> None:
        if not isinstance(code, str) or not code:
            raise ValueError("Action-required code must be non-empty.")
        if self._payload is None:
            raise RuntimeError("No active run manifest.")
        now = self.clock()
        action = {
            "code": code,
            "artifact": sanitize_diagnostic(artifact),
        }
        attempt = self._active_attempt()
        attempt["status"] = "ACTION_REQUIRED"
        attempt["finished_at"] = now
        attempt["failure"] = None
        self._payload["status"] = "ACTION_REQUIRED"
        self._payload["updated_at"] = now
        self._payload["action_required"] = action
        self._payload["failure"] = None
        self._write()
        assert self._logger is not None
        self._logger.emit(
            "run_action_required",
            status="ACTION_REQUIRED",
            code=code,
            artifact=artifact,
        )

    def fail(self, error: BaseException, *, stage: str | None) -> None:
        if self._payload is None:
            raise RuntimeError("No active run manifest.")
        now = self.clock()
        failure = {
            "code": "RUN_FAILED",
            "type": type(error).__name__,
            "message": sanitize_diagnostic(str(error), field_name="message"),
            "stage": stage,
        }
        attempt = self._active_attempt()
        attempt["status"] = "FAILED"
        attempt["finished_at"] = now
        attempt["failure"] = failure
        self._payload["status"] = "FAILED"
        self._payload["updated_at"] = now
        self._payload["failure"] = failure
        self._write()
        assert self._logger is not None
        if stage is not None:
            self._logger.emit(
                "stage_failed",
                level="ERROR",
                stage=stage,
                message=str(error),
                error_type=type(error).__name__,
            )
        self._logger.emit(
            "run_failed",
            level="ERROR",
            status="FAILED",
            stage=stage,
            message=str(error),
            error_type=type(error).__name__,
        )
