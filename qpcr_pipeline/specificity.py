"""Public off-target specificity stage and auditable artifacts."""

from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
import json
from pathlib import Path
from typing import Literal
import uuid

from qpcr_pipeline.config import (
    OffTargetConfig,
    SpecificityConfig,
    validate_off_target_config,
    validate_specificity_config,
)
from qpcr_pipeline.off_targets import OffTargetDataset, load_off_target_datasets
from qpcr_pipeline.primer_design import PrimerDesignResult
from qpcr_pipeline.specificity_matching import (
    GeometryAmplicon,
    MatchHit,
    SpecificityMatchingError,
    all_assay_hits,
    find_plausible_amplicons,
)


class SpecificityError(RuntimeError):
    """Raised when specificity cannot produce a trustworthy result."""


@dataclass(frozen=True, slots=True)
class OffTargetHit:
    dataset_name: str
    assay_id: str
    sequence_id: str
    role: str
    orientation: str
    hit_rank: int
    source_start: int
    source_end: int
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    exact_match: bool
    three_prime_mismatch: bool
    compatible: bool


@dataclass(frozen=True, slots=True)
class HitRetentionSummary:
    dataset_name: str
    assay_id: str
    role: str
    total_hit_count: int
    retained_hit_count: int
    truncated: bool


@dataclass(frozen=True, slots=True)
class PlausibleAmplicon:
    dataset_name: str
    assay_id: str
    sequence_id: str
    orientation: str
    source_start: int
    source_end: int
    amplicon_size: int
    forward_source_start: int
    forward_source_end: int
    reverse_source_start: int
    reverse_source_end: int
    probe_source_sites: tuple[tuple[int, int], ...]
    forward_hit_rank: int | None
    reverse_hit_rank: int | None
    probe_hit_ranks: tuple[int, ...]
    primer_amplicon_plausible: bool
    detectable_off_target: bool


@dataclass(frozen=True, slots=True)
class SpecificityResult:
    status: Literal["SKIPPED", "COMPLETE"]
    dataset_names: tuple[str, ...]
    sequence_count: int
    assay_count: int
    hits: tuple[OffTargetHit, ...]
    amplicons: tuple[PlausibleAmplicon, ...]
    retention: tuple[HitRetentionSummary, ...]
    off_target_hits_path: Path | None
    plausible_amplicons_path: Path | None
    report_path: Path


OFF_TARGET_HIT_COLUMNS = (
    "dataset_name",
    "assay_id",
    "sequence_id",
    "role",
    "orientation",
    "hit_rank",
    "source_start",
    "source_end",
    "mismatch_positions",
    "mismatch_count",
    "exact_match",
    "three_prime_mismatch",
    "compatible",
)
PLAUSIBLE_AMPLICON_COLUMNS = (
    "dataset_name",
    "assay_id",
    "sequence_id",
    "orientation",
    "source_start",
    "source_end",
    "amplicon_size",
    "forward_source_start",
    "forward_source_end",
    "reverse_source_start",
    "reverse_source_end",
    "probe_source_sites",
    "forward_hit_rank",
    "reverse_hit_rank",
    "probe_hit_ranks",
    "primer_amplicon_plausible",
    "detectable_off_target",
)
_DATA_ARTIFACT_NAMES = ("off_target_hits.tsv", "plausible_amplicons.tsv")
_ROLE_ORDER = ("FORWARD", "PROBE", "REVERSE")


def _artifact_paths(output_dir: Path) -> dict[str, Path]:
    directory = output_dir / "specificity"
    return {
        "hits": directory / "off_target_hits.tsv",
        "amplicons": directory / "plausible_amplicons.tsv",
        "report": directory / "specificity_report.json",
    }


def _json_value(value: object) -> object:
    if is_dataclass(value) and not isinstance(value, type):
        return {field.name: _json_value(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    return value


def _tsv_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, tuple):
        if value and all(
            isinstance(item, tuple)
            and len(item) == 2
            and all(isinstance(part, int) for part in item)
            for item in value
        ):
            return ",".join(f"{item[0]}-{item[1]}" for item in value)
        return ",".join(_tsv_value(item) for item in value)
    return str(value)


def _tsv_text(columns: tuple[str, ...], rows: tuple[tuple[object, ...], ...]) -> str:
    lines = ["\t".join(columns)]
    lines.extend("\t".join(_tsv_value(value) for value in row) for row in rows)
    return "\n".join(lines) + "\n"


def _hit_text(hits: tuple[OffTargetHit, ...]) -> str:
    return _tsv_text(
        OFF_TARGET_HIT_COLUMNS,
        tuple(
            (
                hit.dataset_name,
                hit.assay_id,
                hit.sequence_id,
                hit.role,
                hit.orientation,
                hit.hit_rank,
                hit.source_start,
                hit.source_end,
                hit.mismatch_positions,
                hit.mismatch_count,
                hit.exact_match,
                hit.three_prime_mismatch,
                hit.compatible,
            )
            for hit in hits
        ),
    )


def _amplicon_text(amplicons: tuple[PlausibleAmplicon, ...]) -> str:
    return _tsv_text(
        PLAUSIBLE_AMPLICON_COLUMNS,
        tuple(
            (
                item.dataset_name,
                item.assay_id,
                item.sequence_id,
                item.orientation,
                item.source_start,
                item.source_end,
                item.amplicon_size,
                item.forward_source_start,
                item.forward_source_end,
                item.reverse_source_start,
                item.reverse_source_end,
                item.probe_source_sites,
                item.forward_hit_rank,
                item.reverse_hit_rank,
                item.probe_hit_ranks,
                item.primer_amplicon_plausible,
                item.detectable_off_target,
            )
            for item in amplicons
        ),
    )


def _relative_path(path: Path, output_dir: Path) -> str:
    try:
        return path.relative_to(output_dir).as_posix()
    except ValueError:
        return path.as_posix()


def _atomic_write_text(destination: Path, content: str) -> None:
    temporary = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.tmp"
    try:
        with temporary.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def _dataset_provenance(dataset: OffTargetDataset) -> dict[str, object]:
    result: dict[str, object] = {
        "name": dataset.name,
        "source_type": dataset.source_type,
        "source_path": dataset.source_path.as_posix(),
        "sha256": dataset.sha256,
        "sequence_ids": list(dataset.sequence_ids),
    }
    if dataset.frozen_manifest is not None:
        result["frozen_manifest"] = {
            "source": _json_value(dataset.frozen_manifest.get("source")),
            "resolved_entries": _json_value(
                dataset.frozen_manifest.get("resolved_entries", [])
            ),
        }
        result["frozen_manifest_path"] = (
            dataset.frozen_manifest_path.as_posix()
            if dataset.frozen_manifest_path is not None
            else None
        )
    return result


def _retained_hits(
    datasets: tuple[OffTargetDataset, ...],
    primer_design: PrimerDesignResult,
    all_hits: tuple[MatchHit, ...],
    config: SpecificityConfig,
) -> tuple[
    tuple[OffTargetHit, ...],
    tuple[HitRetentionSummary, ...],
    dict[MatchHit, int],
]:
    rank_by_hit: dict[MatchHit, int] = {}
    retention: list[HitRetentionSummary] = []
    for dataset in datasets:
        for assay in primer_design.assays:
            for role in _ROLE_ORDER:
                group = tuple(
                    hit
                    for hit in all_hits
                    if hit.dataset_name == dataset.name
                    and hit.assay_id == assay.assay_id
                    and hit.role == role
                )
                retained = group[: config.max_hits_per_oligo_per_dataset]
                for rank, hit in enumerate(retained, 1):
                    rank_by_hit[hit] = rank
                retention.append(
                    HitRetentionSummary(
                        dataset_name=dataset.name,
                        assay_id=assay.assay_id,
                        role=role,
                        total_hit_count=len(group),
                        retained_hit_count=len(retained),
                        truncated=len(group) > len(retained),
                    )
                )

    public_hits = tuple(
        OffTargetHit(
            dataset_name=hit.dataset_name,
            assay_id=hit.assay_id,
            sequence_id=hit.sequence_id,
            role=hit.role,
            orientation=hit.orientation,
            hit_rank=rank_by_hit[hit],
            source_start=hit.source_start,
            source_end=hit.source_end,
            mismatch_positions=hit.mismatch_positions,
            mismatch_count=hit.mismatch_count,
            exact_match=hit.exact_match,
            three_prime_mismatch=hit.three_prime_mismatch,
            compatible=hit.compatible,
        )
        for hit in all_hits
        if hit in rank_by_hit
    )
    return public_hits, tuple(retention), rank_by_hit


def _public_amplicon(
    item: GeometryAmplicon, rank_by_hit: dict[MatchHit, int]
) -> PlausibleAmplicon:
    return PlausibleAmplicon(
        dataset_name=item.dataset_name,
        assay_id=item.assay_id,
        sequence_id=item.sequence_id,
        orientation=item.orientation,
        source_start=item.source_start,
        source_end=item.source_end,
        amplicon_size=item.amplicon_size,
        forward_source_start=item.forward.source_start,
        forward_source_end=item.forward.source_end,
        reverse_source_start=item.reverse.source_start,
        reverse_source_end=item.reverse.source_end,
        probe_source_sites=tuple(
            (probe.source_start, probe.source_end) for probe in item.probes
        ),
        forward_hit_rank=rank_by_hit.get(item.forward),
        reverse_hit_rank=rank_by_hit.get(item.reverse),
        probe_hit_ranks=tuple(
            rank_by_hit[probe] for probe in item.probes if probe in rank_by_hit
        ),
        primer_amplicon_plausible=item.primer_amplicon_plausible,
        detectable_off_target=item.detectable_off_target,
    )


def _report(
    *,
    status: Literal["SKIPPED", "COMPLETE"],
    config: SpecificityConfig,
    datasets: tuple[OffTargetDataset, ...],
    assay_count: int,
    total_hit_count: int,
    hits: tuple[OffTargetHit, ...],
    amplicons: tuple[PlausibleAmplicon, ...],
    retention: tuple[HitRetentionSummary, ...],
    artifacts: dict[str, str | None],
) -> dict[str, object]:
    detectable = tuple(item for item in amplicons if item.detectable_off_target)
    risky_assays = {(item.dataset_name, item.assay_id) for item in detectable}
    return {
        "schema_version": 1,
        "status": status,
        "enabled": config.enabled,
        "configuration": _json_value(config),
        "datasets": [_dataset_provenance(dataset) for dataset in datasets],
        "counts": {
            "datasets": len(datasets),
            "sequences": sum(len(dataset.records) for dataset in datasets),
            "assays": assay_count,
            "total_compatible_hits": total_hit_count,
            "retained_hits": len(hits),
            "plausible_amplicons": len(amplicons),
            "detectable_off_targets": len(detectable),
            "assays_with_detectable_off_target": len(risky_assays),
        },
        "retention": _json_value(retention),
        "artifacts": artifacts,
    }


def evaluate_specificity(
    primer_design: PrimerDesignResult,
    off_target_configs: tuple[OffTargetConfig, ...],
    config: SpecificityConfig,
    output_dir: Path,
) -> SpecificityResult:
    """Evaluate primer assays against offline off-target datasets and publish artifacts."""
    validate_specificity_config(config)
    output_dir = Path(output_dir)
    paths = _artifact_paths(output_dir)
    paths["report"].parent.mkdir(parents=True, exist_ok=True)

    if not config.enabled:
        paths["report"].unlink(missing_ok=True)
        for name in _DATA_ARTIFACT_NAMES:
            (paths["report"].parent / name).unlink(missing_ok=True)
        report = _report(
            status="SKIPPED",
            config=config,
            datasets=(),
            assay_count=0,
            total_hit_count=0,
            hits=(),
            amplicons=(),
            retention=(),
            artifacts={"off_target_hits": None, "plausible_amplicons": None},
        )
        _atomic_write_text(
            paths["report"],
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        )
        return SpecificityResult(
            status="SKIPPED",
            dataset_names=(),
            sequence_count=0,
            assay_count=0,
            hits=(),
            amplicons=(),
            retention=(),
            off_target_hits_path=None,
            plausible_amplicons_path=None,
            report_path=paths["report"],
        )

    if not isinstance(off_target_configs, tuple) or not off_target_configs:
        raise SpecificityError(
            "Enabled specificity requires at least one off-target dataset."
        )
    try:
        for item in off_target_configs:
            validate_off_target_config(item)
    except ValueError as error:
        raise SpecificityError(str(error)) from error
    names = tuple(item.name for item in off_target_configs)
    if len(set(names)) != len(names):
        raise SpecificityError("Off-target dataset names must be unique.")
    if not isinstance(primer_design, PrimerDesignResult) or primer_design.status != "COMPLETE":
        raise SpecificityError(
            "Enabled specificity requires a COMPLETE PrimerDesignResult."
        )

    try:
        datasets = load_off_target_datasets(off_target_configs)
        all_hits = tuple(
            hit
            for dataset in datasets
            for hit in all_assay_hits(
                dataset.name, dataset.records, primer_design.assays, config
            )
        )
        geometry = tuple(
            amplicon
            for dataset in datasets
            for amplicon in find_plausible_amplicons(
                tuple(hit for hit in all_hits if hit.dataset_name == dataset.name),
                config,
            )
        )
    except (ValueError, SpecificityMatchingError) as error:
        raise SpecificityError(str(error)) from error

    hits, retention, rank_by_hit = _retained_hits(
        datasets, primer_design, all_hits, config
    )
    amplicons = tuple(_public_amplicon(item, rank_by_hit) for item in geometry)
    report = _report(
        status="COMPLETE",
        config=config,
        datasets=datasets,
        assay_count=len(primer_design.assays),
        total_hit_count=len(all_hits),
        hits=hits,
        amplicons=amplicons,
        retention=retention,
        artifacts={
            "off_target_hits": _relative_path(paths["hits"], output_dir),
            "plausible_amplicons": _relative_path(paths["amplicons"], output_dir),
        },
    )

    paths["report"].unlink(missing_ok=True)
    _atomic_write_text(paths["hits"], _hit_text(hits))
    _atomic_write_text(paths["amplicons"], _amplicon_text(amplicons))
    _atomic_write_text(
        paths["report"],
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
    )
    return SpecificityResult(
        status="COMPLETE",
        dataset_names=tuple(dataset.name for dataset in datasets),
        sequence_count=sum(len(dataset.records) for dataset in datasets),
        assay_count=len(primer_design.assays),
        hits=hits,
        amplicons=amplicons,
        retention=retention,
        off_target_hits_path=paths["hits"],
        plausible_amplicons_path=paths["amplicons"],
        report_path=paths["report"],
    )
