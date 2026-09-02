from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable


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
    "stage_started": frozenset({"stage", "action"}),
    "stage_completed": frozenset({"stage", "action", "checkpoint_path"}),
    "stage_reused": frozenset({"stage", "action", "checkpoint_path"}),
    "stage_failed": frozenset({"stage", "message", "error_type"}),
    "run_completed": frozenset({"status", "missing_evidence"}),
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
