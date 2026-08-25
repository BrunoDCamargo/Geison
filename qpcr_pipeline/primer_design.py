"""Pure candidate-region selection for Primer3 assay design."""

from __future__ import annotations

import json
import math
import uuid
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Literal

from qpcr_pipeline.config import PrimerDesignConfig, validate_primer_design_config
from qpcr_pipeline.conservation import (
    ConservationResult,
    PositionConservation,
    WindowConservation,
)

if TYPE_CHECKING:
    from qpcr_pipeline.primer3 import Primer3Runner


class PrimerDesignError(RuntimeError):
    """Raised when primer-design inputs cannot produce traceable results."""


@dataclass(frozen=True, slots=True)
class CandidateRegion:
    region_id: str
    rank: int
    reference_start: int
    reference_end: int
    peak_start: int
    peak_end: int
    position_count: int
    usable_length: int
    usable_fraction: float
    mean_conservation: float
    minimum_conservation: float
    mean_coverage: float
    mean_gap_frequency: float
    mean_entropy_bits: float


@dataclass(frozen=True, slots=True)
class DesignedOligo:
    sequence: str
    reference_start: int
    reference_end: int
    length: int
    tm: float
    gc_percent: float
    penalty: float | None
    metrics: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class AssayCandidate:
    assay_id: str
    region_id: str
    primer3_index: int
    forward_primer: DesignedOligo
    probe: DesignedOligo
    reverse_primer: DesignedOligo
    product_size: int
    pair_penalty: float | None
    metrics: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class PrimerDesignResult:
    status: Literal["SKIPPED", "COMPLETE"]
    reference_id: str | None
    candidates: tuple[CandidateRegion, ...]
    assays: tuple[AssayCandidate, ...]
    candidate_regions_path: Path | None
    assays_path: Path | None
    primer3_input_path: Path | None
    primer3_output_path: Path | None
    report_path: Path


_CANDIDATE_HEADER = (
    "region_id\trank\treference_start\treference_end\tpeak_start\tpeak_end\t"
    "position_count\tusable_length\tusable_fraction\tmean_conservation\t"
    "minimum_conservation\tmean_coverage\tmean_gap_frequency\tmean_entropy_bits\n"
)
_ASSAY_HEADER = (
    "assay_id\tregion_id\tprimer3_index\tforward_sequence\t"
    "forward_reference_start\tforward_reference_end\tforward_length\tforward_tm\t"
    "forward_gc_percent\tforward_penalty\tprobe_sequence\tprobe_reference_start\t"
    "probe_reference_end\tprobe_length\tprobe_tm\tprobe_gc_percent\tprobe_penalty\t"
    "reverse_sequence\treverse_reference_start\treverse_reference_end\treverse_length\t"
    "reverse_tm\treverse_gc_percent\treverse_penalty\tproduct_size\tpair_penalty\n"
)


def design_primers(
    conservation: ConservationResult,
    config: PrimerDesignConfig,
    output_dir: Path,
    *,
    runner: Primer3Runner | None = None,
) -> PrimerDesignResult:
    """Select candidate regions and publish auditable Primer3 artifacts."""
    validate_primer_design_config(config)
    output_dir = Path(output_dir)
    paths = _artifact_paths(output_dir)
    if not config.enabled:
        report = _report(
            status="SKIPPED",
            config=config,
            reference_id=None,
            candidates=(),
            assays=(),
            primer3_details={},
            artifacts={
                "candidate_regions": None,
                "assays": None,
                "primer3_input": None,
                "primer3_output": None,
            },
        )
        paths["report"].parent.mkdir(parents=True, exist_ok=True)
        paths["report"].unlink(missing_ok=True)
        for key in ("candidates", "assays", "input", "output"):
            paths[key].unlink(missing_ok=True)
        _atomic_write_text(
            paths["report"], json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        return PrimerDesignResult(
            status="SKIPPED",
            reference_id=None,
            candidates=(),
            assays=(),
            candidate_regions_path=None,
            assays_path=None,
            primer3_input_path=None,
            primer3_output_path=None,
            report_path=paths["report"],
        )

    candidates = _select_candidate_regions(conservation, config)
    if not candidates:
        report = _report(
            status="COMPLETE",
            config=config,
            reference_id=conservation.reference_id,
            candidates=(),
            assays=(),
            primer3_details={},
            artifacts={
                "candidate_regions": _relative_path(paths["candidates"], output_dir),
                "assays": _relative_path(paths["assays"], output_dir),
                "primer3_input": None,
                "primer3_output": None,
            },
        )
        paths["report"].parent.mkdir(parents=True, exist_ok=True)
        paths["report"].unlink(missing_ok=True)
        _atomic_write_text(paths["candidates"], _candidate_text(()))
        _atomic_write_text(paths["assays"], _assay_text(()))
        paths["input"].unlink(missing_ok=True)
        paths["output"].unlink(missing_ok=True)
        _atomic_write_text(
            paths["report"], json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        return PrimerDesignResult(
            status="COMPLETE",
            reference_id=conservation.reference_id,
            candidates=(),
            assays=(),
            candidate_regions_path=paths["candidates"],
            assays_path=paths["assays"],
            primer3_input_path=None,
            primer3_output_path=None,
            report_path=paths["report"],
        )

    from qpcr_pipeline.primer3 import (
        SubprocessPrimer3Runner,
        build_primer3_input,
        parse_primer3_output,
    )

    paths["report"].parent.mkdir(parents=True, exist_ok=True)
    paths["report"].unlink(missing_ok=True)
    input_text = build_primer3_input(
        conservation.major_consensus, candidates, config
    )
    if runner is None:
        runner = SubprocessPrimer3Runner()
    output_text = runner.run(input_text)
    assays, primer3_details = parse_primer3_output(
        output_text, candidates, conservation.major_consensus
    )
    report = _report(
        status="COMPLETE",
        config=config,
        reference_id=conservation.reference_id,
        candidates=candidates,
        assays=assays,
        primer3_details=primer3_details,
        artifacts={
            "candidate_regions": _relative_path(paths["candidates"], output_dir),
            "assays": _relative_path(paths["assays"], output_dir),
            "primer3_input": _relative_path(paths["input"], output_dir),
            "primer3_output": _relative_path(paths["output"], output_dir),
        },
    )
    _atomic_write_text(paths["candidates"], _candidate_text(candidates))
    _atomic_write_text(paths["assays"], _assay_text(assays))
    _atomic_write_text(paths["input"], input_text)
    _atomic_write_text(paths["output"], output_text)
    _atomic_write_text(
        paths["report"], json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    return PrimerDesignResult(
        status="COMPLETE",
        reference_id=conservation.reference_id,
        candidates=candidates,
        assays=assays,
        candidate_regions_path=paths["candidates"],
        assays_path=paths["assays"],
        primer3_input_path=paths["input"],
        primer3_output_path=paths["output"],
        report_path=paths["report"],
    )


def _artifact_paths(output_dir: Path) -> dict[str, Path]:
    directory = output_dir / "primer_design"
    return {
        "candidates": directory / "candidate_regions.tsv",
        "assays": directory / "assays.tsv",
        "input": directory / "primer3_input.txt",
        "output": directory / "primer3_output.txt",
        "report": directory / "primer_design_report.json",
    }


def _candidate_text(candidates: tuple[CandidateRegion, ...]) -> str:
    lines = [_CANDIDATE_HEADER]
    for candidate in candidates:
        lines.append(
            "\t".join(str(value) for value in asdict(candidate).values()) + "\n"
        )
    return "".join(lines)


def _assay_text(assays: tuple[AssayCandidate, ...]) -> str:
    lines = [_ASSAY_HEADER]
    for assay in assays:
        values: list[object] = [
            assay.assay_id,
            assay.region_id,
            assay.primer3_index,
        ]
        for oligo in (
            assay.forward_primer,
            assay.probe,
            assay.reverse_primer,
        ):
            values.extend(
                (
                    oligo.sequence,
                    oligo.reference_start,
                    oligo.reference_end,
                    oligo.length,
                    oligo.tm,
                    oligo.gc_percent,
                    "" if oligo.penalty is None else oligo.penalty,
                )
            )
        values.extend(
            (
                assay.product_size,
                "" if assay.pair_penalty is None else assay.pair_penalty,
            )
        )
        lines.append("\t".join(str(value) for value in values) + "\n")
    return "".join(lines)


def _report(
    *,
    status: Literal["SKIPPED", "COMPLETE"],
    config: PrimerDesignConfig,
    reference_id: str | None,
    candidates: tuple[CandidateRegion, ...],
    assays: tuple[AssayCandidate, ...],
    primer3_details: dict[str, dict[str, str]],
    artifacts: dict[str, str | None],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "status": status,
        "enabled": config.enabled,
        "configuration": asdict(config),
        "reference_id": reference_id,
        "counts": {"candidates": len(candidates), "assays": len(assays)},
        "candidates": [asdict(candidate) for candidate in candidates],
        "assays": [_assay_report(assay) for assay in assays],
        "primer3_details": primer3_details,
        "artifacts": artifacts,
    }


def _assay_report(assay: AssayCandidate) -> dict[str, object]:
    return {
        "assay_id": assay.assay_id,
        "region_id": assay.region_id,
        "primer3_index": assay.primer3_index,
        "forward_primer": _oligo_report(assay.forward_primer),
        "probe": _oligo_report(assay.probe),
        "reverse_primer": _oligo_report(assay.reverse_primer),
        "product_size": assay.product_size,
        "pair_penalty": assay.pair_penalty,
        "metrics": dict(assay.metrics),
    }


def _oligo_report(oligo: DesignedOligo) -> dict[str, object]:
    return {
        "sequence": oligo.sequence,
        "reference_start": oligo.reference_start,
        "reference_end": oligo.reference_end,
        "length": oligo.length,
        "tm": oligo.tm,
        "gc_percent": oligo.gc_percent,
        "penalty": oligo.penalty,
        "metrics": dict(oligo.metrics),
    }


def _relative_path(path: Path, output_dir: Path) -> str:
    return path.relative_to(output_dir).as_posix()


def _atomic_write_text(destination: Path, content: str) -> None:
    temporary = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.tmp"
    try:
        with temporary.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def _select_candidate_regions(
    conservation: ConservationResult, config: PrimerDesignConfig
) -> tuple[CandidateRegion, ...]:
    """Expand conservation windows into fixed-length reference intervals."""
    _validate_conservation_input(conservation)
    if not conservation.positions:
        return ()

    reference_length = len(conservation.positions)
    generated: list[tuple[CandidateRegion, WindowConservation]] = []
    for window in conservation.windows:
        reference_start, reference_end = _expanded_interval(
            window.reference_start,
            window.reference_end,
            reference_length,
            config.candidate_region_length,
        )
        region = _candidate_region(
            region_id="",
            rank=0,
            reference_start=reference_start,
            reference_end=reference_end,
            peak_start=window.reference_start,
            peak_end=window.reference_end,
            positions=conservation.positions[reference_start - 1:reference_end],
            config=config,
        )
        if _is_eligible(region, config):
            generated.append((region, window))

    unique_regions: dict[
        tuple[int, int], tuple[CandidateRegion, WindowConservation]
    ] = {}
    for region, window in generated:
        interval = region.reference_start, region.reference_end
        existing = unique_regions.get(interval)
        if existing is None or _window_ranking_key(window) < _window_ranking_key(
            existing[1]
        ):
            unique_regions[interval] = region, window

    selected: list[CandidateRegion] = []
    for region, _ in sorted(
        unique_regions.values(), key=lambda item: _ranking_key(item[0])
    ):
        if all(
            _overlap_fraction(region, accepted) <= config.max_region_overlap_fraction
            for accepted in selected
        ):
            selected.append(region)
            if len(selected) == config.max_candidate_regions:
                break

    return tuple(
        replace(region, region_id=f"region-{rank:03d}", rank=rank)
        for rank, region in enumerate(selected, 1)
    )


def _expanded_interval(
    peak_start: int,
    peak_end: int,
    reference_length: int,
    requested_length: int,
) -> tuple[int, int]:
    if reference_length <= requested_length:
        return 1, reference_length

    midpoint = (peak_start + peak_end) // 2
    reference_start = midpoint - (requested_length - 1) // 2
    reference_end = reference_start + requested_length - 1
    if reference_start < 1:
        return 1, requested_length
    if reference_end > reference_length:
        return reference_length - requested_length + 1, reference_length
    return reference_start, reference_end


def _candidate_region(
    *,
    region_id: str,
    rank: int,
    reference_start: int,
    reference_end: int,
    peak_start: int,
    peak_end: int,
    positions: tuple[PositionConservation, ...],
    config: PrimerDesignConfig,
) -> CandidateRegion:
    position_count = len(positions)
    usable_length = sum(
        position.major_allele_frequency >= config.min_minimum_conservation
        and position.coverage >= config.min_mean_coverage
        and position.gap_frequency <= config.max_mean_gap_frequency
        and position.entropy_bits <= config.max_mean_entropy_bits
        for position in positions
    )
    return CandidateRegion(
        region_id=region_id,
        rank=rank,
        reference_start=reference_start,
        reference_end=reference_end,
        peak_start=peak_start,
        peak_end=peak_end,
        position_count=position_count,
        usable_length=usable_length,
        usable_fraction=usable_length / position_count,
        mean_conservation=math.fsum(
            position.major_allele_frequency for position in positions
        ) / position_count,
        minimum_conservation=min(
            position.major_allele_frequency for position in positions
        ),
        mean_coverage=(
            math.fsum(position.coverage for position in positions) / position_count
        ),
        mean_gap_frequency=math.fsum(
            position.gap_frequency for position in positions
        ) / position_count,
        mean_entropy_bits=math.fsum(
            position.entropy_bits for position in positions
        ) / position_count,
    )


def _is_eligible(region: CandidateRegion, config: PrimerDesignConfig) -> bool:
    return (
        region.mean_conservation >= config.min_mean_conservation
        and region.minimum_conservation >= config.min_minimum_conservation
        and region.mean_coverage >= config.min_mean_coverage
        and region.mean_gap_frequency <= config.max_mean_gap_frequency
        and region.mean_entropy_bits <= config.max_mean_entropy_bits
        and region.usable_fraction >= config.min_usable_fraction
    )


def _ranking_key(region: CandidateRegion) -> tuple[float | int, ...]:
    return (
        -region.mean_conservation,
        -region.minimum_conservation,
        -region.mean_coverage,
        region.mean_entropy_bits,
        region.mean_gap_frequency,
        -region.usable_length,
        region.reference_start,
        region.reference_end,
    )


def _window_ranking_key(window: WindowConservation) -> tuple[float | int, ...]:
    return (
        -window.mean_conservation,
        -window.minimum_conservation,
        -window.mean_coverage,
        window.mean_entropy_bits,
        window.mean_gap_frequency,
        window.reference_start,
        window.reference_end,
    )


def _overlap_fraction(first: CandidateRegion, second: CandidateRegion) -> float:
    overlap_length = max(
        0,
        min(first.reference_end, second.reference_end)
        - max(first.reference_start, second.reference_start)
        + 1,
    )
    shorter_length = min(
        first.reference_end - first.reference_start + 1,
        second.reference_end - second.reference_start + 1,
    )
    return overlap_length / shorter_length


def _validate_conservation_input(conservation: ConservationResult) -> None:
    if conservation.status != "COMPLETE":
        raise PrimerDesignError("Candidate selection requires COMPLETE conservation.")

    if (
        not isinstance(conservation.reference_id, str)
        or not conservation.reference_id.strip()
    ):
        raise PrimerDesignError("Complete conservation requires a reference ID.")

    if not conservation.positions:
        if conservation.windows or conservation.major_consensus or conservation.iupac_consensus:
            raise PrimerDesignError("Empty conservation cannot contain positions or windows.")
        return

    reference_length = len(conservation.positions)
    if (
        len(conservation.major_consensus) != reference_length
        or len(conservation.iupac_consensus) != reference_length
    ):
        raise PrimerDesignError(
            "Conservation consensus length must match reference positions."
        )

    for expected_position, position in enumerate(conservation.positions, 1):
        if position.reference_position != expected_position:
            raise PrimerDesignError(
                "Conservation reference positions must be contiguous from 1."
            )

    for window in conservation.windows:
        if not (
            1 <= window.reference_start <= window.reference_end <= reference_length
        ):
            raise PrimerDesignError(
                "Conservation windows must lie within reference bounds."
            )
