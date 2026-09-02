from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Mapping

from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.run_recording import sanitize_diagnostic


_SHA256_HEX = re.compile(r"^[0-9a-fA-F]{64}$")


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


def _sha256_identity(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    if value.startswith("sha256:") and _SHA256_HEX.fullmatch(value[7:]):
        return "sha256:" + value[7:].lower()
    if _SHA256_HEX.fullmatch(value):
        return "sha256:" + value.lower()
    return None


def _project_source_request(source: Mapping[str, object], source_mode: str) -> dict[str, object]:
    if source_mode == "accessions":
        accessions = source.get("requested_accessions")
        if (
            not isinstance(accessions, list)
            or any(not isinstance(item, str) or not item for item in accessions)
        ):
            raise ValueError("NCBI dataset requested accessions are invalid.")
        return {"requested_accessions": list(accessions)}

    query = source.get("query")
    if not isinstance(query, str) or not query:
        raise ValueError("NCBI dataset query provenance is invalid.")
    return {"query": query}


def _ncbi_manifest_projection(
    config: PipelineConfig,
    outdir: Path,
    checkpoint_source: Mapping[str, object],
    counts: dict[str, int],
) -> dict[str, object]:
    selected = config.input_ncbi
    if selected is None:
        raise ValueError("NCBI provenance requires an NCBI input configuration.")

    manifest_path = Path(outdir) / "ncbi_dataset_manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("NCBI dataset manifest must be a mapping.")
    source = payload.get("source")
    resolved_entries = payload.get("resolved_entries")
    consolidated = payload.get("consolidated")
    if not isinstance(source, dict) or not isinstance(resolved_entries, list):
        raise ValueError("NCBI dataset manifest provenance fields are invalid.")
    if not isinstance(consolidated, dict):
        raise ValueError("NCBI dataset manifest is missing consolidated identity.")

    resolved_versions: list[str] = []
    for entry in resolved_entries:
        if not isinstance(entry, dict):
            raise ValueError("NCBI resolved provenance entry is invalid.")
        accession_version = entry.get("accession_version")
        if not isinstance(accession_version, str) or not accession_version:
            raise ValueError("NCBI resolved provenance accession version is invalid.")
        resolved_versions.append(accession_version)

    dataset_sha256 = _sha256_identity(checkpoint_source.get("records_sha256"))
    if dataset_sha256 is None:
        dataset_sha256 = _sha256_identity(consolidated.get("sha256"))
    if dataset_sha256 is None:
        raise ValueError("NCBI dataset provenance is missing a SHA-256 identity.")

    source_mode = source.get("mode")
    if source_mode not in {"query", "accessions"}:
        raise ValueError("NCBI dataset provenance mode is invalid.")
    mode = "frozen_dataset" if selected.frozen_dataset is not None else source_mode

    result: dict[str, object] = {
        "kind": "ncbi",
        "mode": mode,
        "dataset_sha256": dataset_sha256,
        "resolved_accession_versions": resolved_versions,
        **counts,
    }
    if selected.frozen_dataset is not None:
        result["source_dataset_mode"] = source_mode
        result["configured_path"] = str(selected.frozen_dataset)
        result.update(_project_source_request(source, source_mode))
        manifest_sha256 = _sha256_identity(checkpoint_source.get("manifest_sha256"))
        if manifest_sha256 is not None:
            result["source_manifest_sha256"] = manifest_sha256
    elif selected.query is not None:
        result["query"] = selected.query
    else:
        result["requested_accessions"] = list(selected.accessions)
    return result


def build_input_provenance(
    config: PipelineConfig,
    outdir: Path,
    qc_result: object,
    input_manifest: object,
) -> dict[str, object]:
    """Project safe input provenance from existing checkpoint identities."""
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
        source_sha256 = _sha256_identity(source.get("sha256"))
        if source_sha256 is None:
            raise ValueError("Local input checkpoint is missing its source SHA-256 identity.")
        return {
            "kind": input_format,
            "configured_path": str(path),
            "source_sha256": source_sha256,
            **counts,
        }

    return _ncbi_manifest_projection(config, Path(outdir), source, counts)


def build_reference_provenance(alignment_result: object) -> dict[str, str | None]:
    return {
        "id": getattr(alignment_result, "reference_id", None),
        "mode": getattr(alignment_result, "reference_mode", None),
    }
