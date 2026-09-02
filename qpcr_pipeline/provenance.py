from __future__ import annotations

from pathlib import Path
from typing import Mapping

from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.run_recording import sanitize_diagnostic


def effective_config_payload(config: PipelineConfig) -> dict[str, object]:
    """Return the effective configuration in JSON-safe sanitized form."""
    payload = sanitize_diagnostic(config)
    if not isinstance(payload, dict):
        raise ValueError("Effective pipeline configuration must serialize as a mapping.")
    return payload


def _qc_counts(qc_result: object) -> dict[str, int]:
    records = getattr(qc_result, "records", ())
    evaluation_set = getattr(qc_result, "evaluation_set", None)
    sequence_ids = getattr(evaluation_set, "sequence_ids", ())
    accepted = len(sequence_ids)
    return {
        "accepted_count": accepted,
        "rejected_count": max(0, len(records) - accepted),
    }


def build_input_provenance(
    config: PipelineConfig,
    outdir: Path,
    qc_result: object,
    input_manifest: object,
) -> dict[str, object]:
    """Project safe input provenance from existing checkpoint identities."""
    del outdir
    counts = _qc_counts(qc_result)
    selected = config.selected_input
    inputs = getattr(input_manifest, "inputs", {})
    if not isinstance(inputs, Mapping):
        inputs = {}
    source = inputs.get("target_source", {})
    if not isinstance(source, Mapping):
        source = {}

    if isinstance(selected, tuple):
        path, input_format = selected
        source_sha256 = source.get("sha256")
        if not isinstance(source_sha256, str) or not source_sha256.startswith("sha256:"):
            raise ValueError("Local input checkpoint is missing its source SHA-256 identity.")
        return {
            "kind": input_format,
            "configured_path": str(path),
            "source_sha256": source_sha256,
            **counts,
        }

    return {"kind": "ncbi", **counts}


def build_reference_provenance(alignment_result: object) -> dict[str, str | None]:
    return {
        "id": getattr(alignment_result, "reference_id", None),
        "mode": getattr(alignment_result, "reference_mode", None),
    }
