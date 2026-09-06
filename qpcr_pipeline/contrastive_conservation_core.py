"""Contrastive target-conservation analysis and auditable artifacts."""

from __future__ import annotations

import csv
import io
import json
import os
import tempfile
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Literal

from qpcr_pipeline.challenge_panel import (
    ChallengeDatasetBinding,
    resolve_challenge_datasets,
)
from qpcr_pipeline.config import (
    ContrastiveConservationConfig,
    OffTargetConfig,
    PrimerDesignConfig,
)
from qpcr_pipeline.conservation import ConservationResult, WindowConservation
from qpcr_pipeline.contrastive_similarity import (
    BiopythonLocalSimilarityEngine,
    RegionSimilarityEngine,
)
from qpcr_pipeline.panel import Criticality
from qpcr_pipeline.panel_manifest import ApprovedPanelManifest
from qpcr_pipeline.region_selection import (
    CandidateRegion,
    candidate_region_from_window,
    is_window_target_eligible,
    overlap_fraction,
)


@dataclass(frozen=True, slots=True)
class ChallengeDatasetSummary:
    name: str
    criticality: Criticality
    source_type: str
    records_sha256: str
    sequence_count: int


@dataclass(frozen=True, slots=True)
class DatasetWindowEvidence:
    reference_start: int
    reference_end: int
    dataset_name: str
    criticality: Criticality
    sequence_count: int
    best_sequence_id: str | None
    best_orientation: str | None
    similarity: float


@dataclass(frozen=True, slots=True)
class ContrastWindowEvidence:
    reference_start: int
    reference_end: int
    target_mean_conservation: float
    target_minimum_conservation: float
    target_mean_coverage: float
    target_mean_gap_frequency: float
    target_mean_entropy_bits: float
    target_eligible: bool
    worst_dataset_name: str | None
    worst_dataset_criticality: Criticality | None
    worst_similarity: float | None
    worst_critical_similarity: float | None
    worst_important_similarity: float | None
    contrast_margin: float | None


@dataclass(frozen=True, slots=True)
class ContrastCandidateRegion:
    region: CandidateRegion
    contributing_windows: tuple[tuple[int, int], ...]
    worst_dataset_name: str | None
    worst_dataset_criticality: Criticality | None
    worst_similarity: float | None
    worst_critical_similarity: float | None
    worst_important_similarity: float | None
    contrast_margin: float | None


@dataclass(frozen=True, slots=True)
class ContrastiveConservationResult:
    status: Literal["SKIPPED", "COMPLETE"]
    reference_id: str | None
    windows: tuple[ContrastWindowEvidence, ...]
    dataset_evidence: tuple[DatasetWindowEvidence, ...]
    candidates: tuple[ContrastCandidateRegion, ...]
    challenge_datasets: tuple[ChallengeDatasetSummary, ...]
    window_metrics_path: Path | None
    dataset_metrics_path: Path | None
    candidate_regions_path: Path | None
    report_path: Path
    html_report_path: Path | None


def analyze_contrastive_conservation(
    conservation: ConservationResult,
    approved_panel: ApprovedPanelManifest | None,
    off_target_configs: tuple[OffTargetConfig, ...],
    config: ContrastiveConservationConfig,
    primer_config: PrimerDesignConfig,
    output_dir: Path,
    *,
    similarity_engine: RegionSimilarityEngine | None = None,
) -> ContrastiveConservationResult:
    """Measure target-window representation across the approved CHALLENGE panel."""
    output_dir = Path(output_dir)
    stage_dir = output_dir / "contrastive_conservation"
    report_path = stage_dir / "contrastive_conservation_report.json"

    if not config.enabled:
        result = ContrastiveConservationResult(
            status="SKIPPED",
            reference_id=None,
            windows=(),
            dataset_evidence=(),
            candidates=(),
            challenge_datasets=(),
            window_metrics_path=None,
            dataset_metrics_path=None,
            candidate_regions_path=None,
            report_path=report_path,
            html_report_path=None,
        )
        _atomic_write_text(
            report_path,
            json.dumps(_report_payload(result, config, approved_panel), indent=2, sort_keys=True)
            + "\n",
        )
        return result

    if conservation.status != "COMPLETE":
        raise ValueError(
            "Enabled contrastive conservation requires COMPLETE conservation."
        )
    if approved_panel is None:
        raise ValueError(
            "Enabled contrastive conservation requires an approved frozen panel."
        )

    bindings = resolve_challenge_datasets(approved_panel, off_target_configs)
    engine = similarity_engine or BiopythonLocalSimilarityEngine()
    challenge_summaries = tuple(_challenge_summary(binding) for binding in bindings)

    window_rows: list[ContrastWindowEvidence] = []
    dataset_rows: list[DatasetWindowEvidence] = []
    raw_candidates: list[
        tuple[WindowConservation, ContrastWindowEvidence, CandidateRegion]
    ] = []

    for window in conservation.windows:
        query = conservation.major_consensus[
            window.reference_start - 1 : window.reference_end
        ]
        region = candidate_region_from_window(conservation, window, primer_config)
        eligible = is_window_target_eligible(conservation, window, primer_config)
        evidence_for_window: list[DatasetWindowEvidence] = []
        for binding in bindings:
            match = engine.best_match(query, binding.dataset.records)
            row = DatasetWindowEvidence(
                reference_start=window.reference_start,
                reference_end=window.reference_end,
                dataset_name=binding.name,
                criticality=binding.criticality,
                sequence_count=len(binding.dataset.records),
                best_sequence_id=None if match is None else match.sequence_id,
                best_orientation=None if match is None else match.orientation,
                similarity=0.0 if match is None else match.similarity,
            )
            dataset_rows.append(row)
            evidence_for_window.append(row)

        summary = _window_evidence(window, eligible, evidence_for_window)
        window_rows.append(summary)
        if eligible:
            raw_candidates.append((window, summary, region))

    candidates = _consolidate_candidates(raw_candidates, primer_config)
    window_metrics_path = stage_dir / "window_metrics.tsv"
    dataset_metrics_path = stage_dir / "dataset_metrics.tsv"
    candidate_regions_path = stage_dir / "candidate_regions.tsv"

    _atomic_write_text(window_metrics_path, _window_tsv(tuple(window_rows)))
    _atomic_write_text(dataset_metrics_path, _dataset_tsv(tuple(dataset_rows)))
    _atomic_write_text(candidate_regions_path, _candidate_tsv(candidates))

    result = ContrastiveConservationResult(
        status="COMPLETE",
        reference_id=conservation.reference_id,
        windows=tuple(window_rows),
        dataset_evidence=tuple(dataset_rows),
        candidates=candidates,
        challenge_datasets=challenge_summaries,
        window_metrics_path=window_metrics_path,
        dataset_metrics_path=dataset_metrics_path,
        candidate_regions_path=candidate_regions_path,
        report_path=report_path,
        html_report_path=None,
    )
    _atomic_write_text(
        report_path,
        json.dumps(_report_payload(result, config, approved_panel), indent=2, sort_keys=True)
        + "\n",
    )
    return result


def _challenge_summary(binding: ChallengeDatasetBinding) -> ChallengeDatasetSummary:
    return ChallengeDatasetSummary(
        name=binding.name,
        criticality=binding.criticality,
        source_type=binding.dataset.source_type,
        records_sha256=binding.dataset.sha256,
        sequence_count=len(binding.dataset.records),
    )


def _window_evidence(
    window: WindowConservation,
    eligible: bool,
    dataset_rows: list[DatasetWindowEvidence],
) -> ContrastWindowEvidence:
    worst = max(dataset_rows, key=lambda row: row.similarity)
    critical = [row.similarity for row in dataset_rows if row.criticality == "CRITICAL"]
    important = [row.similarity for row in dataset_rows if row.criticality == "IMPORTANT"]
    return ContrastWindowEvidence(
        reference_start=window.reference_start,
        reference_end=window.reference_end,
        target_mean_conservation=window.mean_conservation,
        target_minimum_conservation=window.minimum_conservation,
        target_mean_coverage=window.mean_coverage,
        target_mean_gap_frequency=window.mean_gap_frequency,
        target_mean_entropy_bits=window.mean_entropy_bits,
        target_eligible=eligible,
        worst_dataset_name=worst.dataset_name,
        worst_dataset_criticality=worst.criticality,
        worst_similarity=worst.similarity,
        worst_critical_similarity=max(critical) if critical else None,
        worst_important_similarity=max(important) if important else None,
        contrast_margin=window.mean_conservation - worst.similarity,
    )


def _evidence_order_key(row: ContrastWindowEvidence) -> tuple[float | int, ...]:
    return (
        row.worst_critical_similarity
        if row.worst_critical_similarity is not None
        else -1.0,
        row.worst_important_similarity
        if row.worst_important_similarity is not None
        else -1.0,
        row.worst_similarity if row.worst_similarity is not None else -1.0,
        -row.target_mean_conservation,
        row.target_mean_entropy_bits,
        row.reference_start,
        row.reference_end,
    )


def _consolidate_candidates(
    rows: list[tuple[WindowConservation, ContrastWindowEvidence, CandidateRegion]],
    config: PrimerDesignConfig,
) -> tuple[ContrastCandidateRegion, ...]:
    grouped: dict[
        tuple[int, int],
        list[tuple[WindowConservation, ContrastWindowEvidence, CandidateRegion]],
    ] = {}
    for item in rows:
        region = item[2]
        grouped.setdefault((region.reference_start, region.reference_end), []).append(item)

    consolidated: list[ContrastCandidateRegion] = []
    for group in grouped.values():
        representative = min(group, key=lambda item: _evidence_order_key(item[1]))
        _, evidence, region = representative
        contributing = tuple(
            sorted((item[0].reference_start, item[0].reference_end) for item in group)
        )
        consolidated.append(
            ContrastCandidateRegion(
                region=region,
                contributing_windows=contributing,
                worst_dataset_name=evidence.worst_dataset_name,
                worst_dataset_criticality=evidence.worst_dataset_criticality,
                worst_similarity=evidence.worst_similarity,
                worst_critical_similarity=evidence.worst_critical_similarity,
                worst_important_similarity=evidence.worst_important_similarity,
                contrast_margin=evidence.contrast_margin,
            )
        )

    consolidated.sort(
        key=lambda item: (
            item.worst_critical_similarity
            if item.worst_critical_similarity is not None
            else -1.0,
            item.worst_important_similarity
            if item.worst_important_similarity is not None
            else -1.0,
            item.worst_similarity if item.worst_similarity is not None else -1.0,
            -item.region.mean_conservation,
            item.region.mean_entropy_bits,
            item.region.reference_start,
            item.region.reference_end,
        )
    )

    selected: list[ContrastCandidateRegion] = []
    for candidate in consolidated:
        if all(
            overlap_fraction(candidate.region, existing.region)
            <= config.max_region_overlap_fraction
            for existing in selected
        ):
            selected.append(candidate)
            if len(selected) == config.max_candidate_regions:
                break

    return tuple(
        replace(
            item,
            region=replace(
                item.region,
                region_id=f"contrast-region-{rank:03d}",
                rank=rank,
            ),
        )
        for rank, item in enumerate(selected, 1)
    )


def _tsv_text(header: tuple[str, ...], rows: list[tuple[object, ...]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return stream.getvalue()


def _window_tsv(rows: tuple[ContrastWindowEvidence, ...]) -> str:
    header = (
        "reference_start",
        "reference_end",
        "target_mean_conservation",
        "target_minimum_conservation",
        "target_mean_coverage",
        "target_mean_gap_frequency",
        "target_mean_entropy_bits",
        "target_eligible",
        "worst_dataset_name",
        "worst_dataset_criticality",
        "worst_similarity",
        "worst_critical_similarity",
        "worst_important_similarity",
        "contrast_margin",
    )
    return _tsv_text(
        header,
        [tuple(getattr(row, field) for field in header) for row in rows],
    )


def _dataset_tsv(rows: tuple[DatasetWindowEvidence, ...]) -> str:
    header = (
        "reference_start",
        "reference_end",
        "dataset_name",
        "criticality",
        "sequence_count",
        "best_sequence_id",
        "best_orientation",
        "similarity",
    )
    return _tsv_text(
        header,
        [tuple(getattr(row, field) for field in header) for row in rows],
    )


def _candidate_tsv(rows: tuple[ContrastCandidateRegion, ...]) -> str:
    header = (
        "region_id",
        "rank",
        "reference_start",
        "reference_end",
        "peak_start",
        "peak_end",
        "contributing_windows",
        "mean_conservation",
        "minimum_conservation",
        "mean_coverage",
        "mean_gap_frequency",
        "mean_entropy_bits",
        "worst_dataset_name",
        "worst_dataset_criticality",
        "worst_similarity",
        "worst_critical_similarity",
        "worst_important_similarity",
        "contrast_margin",
    )
    values: list[tuple[object, ...]] = []
    for item in rows:
        region = item.region
        values.append(
            (
                region.region_id,
                region.rank,
                region.reference_start,
                region.reference_end,
                region.peak_start,
                region.peak_end,
                json.dumps(item.contributing_windows),
                region.mean_conservation,
                region.minimum_conservation,
                region.mean_coverage,
                region.mean_gap_frequency,
                region.mean_entropy_bits,
                item.worst_dataset_name,
                item.worst_dataset_criticality,
                item.worst_similarity,
                item.worst_critical_similarity,
                item.worst_important_similarity,
                item.contrast_margin,
            )
        )
    return _tsv_text(header, values)


def _report_payload(
    result: ContrastiveConservationResult,
    config: ContrastiveConservationConfig,
    approved_panel: ApprovedPanelManifest | None,
) -> dict[str, object]:
    return {
        "status": result.status,
        "reference_id": result.reference_id,
        "configuration": {"enabled": config.enabled},
        "approved_panel": None
        if approved_panel is None
        else {
            "proposal_sha256": approved_panel.proposal_sha256,
            "target_name": approved_panel.definition.target.name,
        },
        "counts": {
            "windows": len(result.windows),
            "dataset_evidence_rows": len(result.dataset_evidence),
            "candidate_regions": len(result.candidates),
            "challenge_datasets": len(result.challenge_datasets),
        },
        "challenge_datasets": [_jsonable(asdict(item)) for item in result.challenge_datasets],
        "windows": [_jsonable(asdict(item)) for item in result.windows],
        "dataset_evidence": [_jsonable(asdict(item)) for item in result.dataset_evidence],
        "candidates": [_jsonable(asdict(item)) for item in result.candidates],
        "artifacts": {
            "window_metrics": _path_text(result.window_metrics_path),
            "dataset_metrics": _path_text(result.dataset_metrics_path),
            "candidate_regions": _path_text(result.candidate_regions_path),
            "html_report": _path_text(result.html_report_path),
        },
    }


def _path_text(path: Path | None) -> str | None:
    return None if path is None else path.as_posix()


def _jsonable(value: object) -> object:
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        temp_path.replace(path)
    except BaseException:
        temp_path.unlink(missing_ok=True)
        raise
