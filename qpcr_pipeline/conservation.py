"""Deterministic conservation analysis and atomic scientific artifacts."""

from __future__ import annotations

import json
import math
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Mapping

from Bio.SeqFeature import SeqFeature

from qpcr_pipeline.alignment import (
    AlignedSequence,
    AlignmentCoordinate,
    AlignmentResult,
)
from qpcr_pipeline.config import ConservationConfig, validate_conservation_config
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.report_html import render_conservation_html


class ConservationError(RuntimeError):
    """Raised when conservation inputs cannot produce traceable results."""


@dataclass(frozen=True, slots=True)
class PositionConservation:
    alignment_position: int
    reference_position: int | None
    reference_base: str | None
    depth: int
    coverage: float
    frequency_a: float
    frequency_c: float
    frequency_g: float
    frequency_t: float
    gap_frequency: float
    major_allele_frequency: float
    entropy_bits: float
    major_consensus: str
    iupac_consensus: str


@dataclass(frozen=True, slots=True)
class WindowConservation:
    reference_start: int
    reference_end: int
    position_count: int
    mean_conservation: float
    minimum_conservation: float
    mean_coverage: float
    mean_gap_frequency: float
    mean_entropy_bits: float


@dataclass(frozen=True, slots=True)
class ReferenceAnnotation:
    feature_type: str
    start: int
    end: int
    strand: int | None
    label: str


@dataclass(frozen=True, slots=True)
class ConservationResult:
    status: Literal["SKIPPED", "COMPLETE"]
    reference_id: str | None
    positions: tuple[PositionConservation, ...]
    windows: tuple[WindowConservation, ...]
    annotations: tuple[ReferenceAnnotation, ...]
    major_consensus: str
    iupac_consensus: str
    position_metrics_path: Path | None
    window_metrics_path: Path | None
    major_consensus_path: Path | None
    iupac_consensus_path: Path | None
    html_report_path: Path | None
    report_path: Path


_IUPAC_BASES = {
    "A": frozenset("A"), "C": frozenset("C"),
    "G": frozenset("G"), "T": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"),
    "S": frozenset("CG"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"),
    "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
_BASES_TO_IUPAC = {bases: code for code, bases in _IUPAC_BASES.items()}
_CANONICAL_BASES = "ACGT"
_IUPAC_VOTE_SCALE = 12
_VALID_ALIGNMENT_SYMBOLS = frozenset(_IUPAC_BASES) | {"-"}
_POSITION_HEADER = (
    "alignment_position\treference_position\treference_base\tdepth\tcoverage\t"
    "frequency_a\tfrequency_c\tfrequency_g\tfrequency_t\tgap_frequency\t"
    "major_allele_frequency\tentropy_bits\tmajor_consensus\tiupac_consensus\n"
)
_WINDOW_HEADER = (
    "reference_start\treference_end\tposition_count\tmean_conservation\t"
    "minimum_conservation\tmean_coverage\tmean_gap_frequency\tmean_entropy_bits\n"
)
_METRIC_DEFINITIONS = {
    "base_frequencies": "fractional IUPAC support normalized by non-gap depth",
    "conservation": "major allele frequency",
    "coverage": "non-gap depth divided by sequence count",
    "entropy_bits": "Shannon entropy of A/C/G/T frequencies in bits",
    "gap_frequency": "gap count divided by sequence count",
}


def analyze_conservation(
    records: tuple[LocalSequenceRecord, ...],
    alignment: AlignmentResult,
    config: ConservationConfig,
    output_dir: Path,
    *,
    target_name: str,
) -> ConservationResult:
    """Calculate conservation and atomically publish its complete artifact set."""
    validate_conservation_config(config)
    output_dir = Path(output_dir)
    paths = _artifact_paths(output_dir)
    discovery_ids = tuple(alignment.discovery_set.sequence_ids)

    if not config.enabled:
        report = _report(
            status="SKIPPED",
            enabled=False,
            config=config,
            discovery_ids=discovery_ids,
            reference_id=None,
            positions=(),
            windows=(),
            annotations=(),
            annotation_skipped=0,
            major_consensus="",
            iupac_consensus="",
            artifacts={key: None for key in _artifact_report_keys()},
        )
        _publish_skipped(paths, report)
        return ConservationResult(
            "SKIPPED", None, (), (), (), "", "", None, None, None, None,
            None, paths["report"],
        )

    records_by_id, sequences, coordinates = _validate_enabled_inputs(
        records, alignment, discovery_ids
    )
    positions = _calculate_positions(sequences, coordinates)
    windows = _calculate_windows(positions, config)
    reference_id = alignment.reference_id
    annotations, annotation_skipped = _extract_annotations(
        records_by_id.get(reference_id) if reference_id is not None else None,
        sum(position.reference_position is not None for position in positions),
    )
    major_consensus = "".join(
        position.major_consensus
        for position in positions
        if position.reference_position is not None
    )
    iupac_consensus = "".join(
        position.iupac_consensus
        for position in positions
        if position.reference_position is not None
    )

    html_text = render_conservation_html(
        target_name=target_name,
        reference_id=reference_id,
        sequence_count=len(discovery_ids),
        config=config,
        windows=windows,
        annotations=annotations,
    )
    artifacts = {
        "position_metrics": _relative_path(paths["positions"], output_dir),
        "window_metrics": _relative_path(paths["windows"], output_dir),
        "major_consensus": _relative_path(paths["major"], output_dir),
        "iupac_consensus": _relative_path(paths["iupac"], output_dir),
        "html_report": _relative_path(paths["html"], output_dir),
    }
    report = _report(
        status="COMPLETE",
        enabled=True,
        config=config,
        discovery_ids=discovery_ids,
        reference_id=reference_id,
        positions=positions,
        windows=windows,
        annotations=annotations,
        annotation_skipped=annotation_skipped,
        major_consensus=major_consensus,
        iupac_consensus=iupac_consensus,
        artifacts=artifacts,
    )
    _publish_complete(
        paths,
        _position_text(positions),
        _window_text(windows),
        _fasta_text("geison-major-consensus", major_consensus),
        _fasta_text("geison-iupac-consensus", iupac_consensus),
        html_text,
        report,
    )
    return ConservationResult(
        "COMPLETE",
        reference_id,
        positions,
        windows,
        annotations,
        major_consensus,
        iupac_consensus,
        paths["positions"],
        paths["windows"],
        paths["major"],
        paths["iupac"],
        paths["html"],
        paths["report"],
    )


def _artifact_paths(output_dir: Path) -> dict[str, Path]:
    directory = output_dir / "conservation"
    return {
        "positions": directory / "position_metrics.tsv",
        "windows": directory / "window_metrics.tsv",
        "major": directory / "consensus_major.fasta",
        "iupac": directory / "consensus_iupac.fasta",
        "html": output_dir / "report.html",
        "report": directory / "conservation_report.json",
    }


def _artifact_report_keys() -> tuple[str, ...]:
    return (
        "position_metrics", "window_metrics", "major_consensus",
        "iupac_consensus", "html_report",
    )


def _validate_enabled_inputs(
    records: tuple[LocalSequenceRecord, ...],
    alignment: AlignmentResult,
    discovery_ids: tuple[str, ...],
) -> tuple[
    dict[str, LocalSequenceRecord],
    tuple[AlignedSequence, ...],
    tuple[AlignmentCoordinate, ...],
]:
    if alignment.status != "COMPLETE":
        raise ConservationError("Enabled conservation requires a COMPLETE alignment result.")
    if len(set(discovery_ids)) != len(discovery_ids):
        raise ConservationError("Discovery Set IDs must be unique.")

    sequences = tuple(alignment.sequences)
    aligned_ids = tuple(item.sequence_id for item in sequences)
    if aligned_ids != discovery_ids or len(set(aligned_ids)) != len(aligned_ids):
        raise ConservationError(
            "Alignment sequence IDs must equal the Discovery Set IDs exactly once and in order."
        )

    record_ids = tuple(item.sequence_id for item in records)
    if len(set(record_ids)) != len(record_ids) or set(record_ids) != set(discovery_ids):
        raise ConservationError("Supplied records must match the Discovery Set records exactly.")
    records_by_id = {item.sequence_id: item for item in records}

    coordinates = tuple(alignment.coordinates)
    if not discovery_ids:
        if alignment.reference_id is not None or sequences or coordinates:
            raise ConservationError("An empty Discovery Set must have an empty alignment and no reference.")
        return records_by_id, sequences, coordinates

    if alignment.reference_id is None or alignment.reference_id not in discovery_ids:
        raise ConservationError("A non-empty alignment requires a selected reference.")
    reference = sequences[discovery_ids.index(alignment.reference_id)]
    if reference.orientation != "forward":
        raise ConservationError("The selected reference must not be reversed.")

    lengths = {len(item.aligned_sequence) for item in sequences}
    if len(lengths) != 1:
        raise ConservationError("Aligned sequences must have equal length.")
    for item in sequences:
        unknown = set(item.aligned_sequence) - _VALID_ALIGNMENT_SYMBOLS
        if unknown:
            rendered = "".join(sorted(unknown))
            raise ConservationError(f"Alignment sequence {item.sequence_id!r} contains unknown symbol(s): {rendered}.")

    alignment_length = lengths.pop()
    if len(coordinates) != alignment_length:
        raise ConservationError("Alignment coordinate count does not match alignment length.")
    expected_coordinates = _coordinates_for_reference(reference.aligned_sequence)
    if coordinates != expected_coordinates:
        raise ConservationError("Alignment coordinates do not match the selected reference.")
    reference_record = records_by_id[alignment.reference_id]
    if reference.aligned_sequence.replace("-", "") != reference_record.sequence:
        raise ConservationError("Aligned reference does not match the supplied reference record.")
    return records_by_id, sequences, coordinates


def _coordinates_for_reference(sequence: str) -> tuple[AlignmentCoordinate, ...]:
    reference_position = 0
    result = []
    for alignment_position, base in enumerate(sequence, 1):
        if base == "-":
            result.append(AlignmentCoordinate(alignment_position, None, None))
        else:
            reference_position += 1
            result.append(AlignmentCoordinate(alignment_position, reference_position, base))
    return tuple(result)


def _calculate_positions(
    sequences: tuple[AlignedSequence, ...],
    coordinates: tuple[AlignmentCoordinate, ...],
) -> tuple[PositionConservation, ...]:
    if not sequences:
        return ()
    sequence_count = len(sequences)
    positions = []
    for index, coordinate in enumerate(coordinates):
        symbols = tuple(item.aligned_sequence[index] for item in sequences)
        gap_count = symbols.count("-")
        depth = sequence_count - gap_count
        if depth == 0:
            raise ConservationError(
                f"Alignment column {coordinate.alignment_position} is all-gap."
            )
        vote_units = {base: 0 for base in _CANONICAL_BASES}
        for symbol in symbols:
            if symbol == "-":
                continue
            support = _IUPAC_BASES[symbol]
            contribution = _IUPAC_VOTE_SCALE // len(support)
            for base in support:
                vote_units[base] += contribution
        vote_unit_depth = depth * _IUPAC_VOTE_SCALE
        frequencies = {
            base: vote_units[base] / vote_unit_depth for base in _CANONICAL_BASES
        }
        maximum_units = max(vote_units.values())
        maximum = maximum_units / vote_unit_depth
        tied = tuple(
            base for base in _CANONICAL_BASES if vote_units[base] == maximum_units
        )
        reference_base = coordinate.reference_base
        major = reference_base if reference_base in tied else tied[0]
        positive_bases = frozenset(
            base for base in _CANONICAL_BASES if vote_units[base] > 0
        )
        entropy = -math.fsum(
            frequency * math.log2(frequency)
            for frequency in frequencies.values()
            if frequency > 0
        )
        if entropy == 0:
            entropy = 0.0
        positions.append(
            PositionConservation(
                coordinate.alignment_position,
                coordinate.reference_position,
                coordinate.reference_base,
                depth,
                depth / sequence_count,
                frequencies["A"],
                frequencies["C"],
                frequencies["G"],
                frequencies["T"],
                gap_count / sequence_count,
                maximum,
                entropy,
                major,
                _BASES_TO_IUPAC[positive_bases],
            )
        )
    return tuple(positions)


def _calculate_windows(
    positions: tuple[PositionConservation, ...],
    config: ConservationConfig,
) -> tuple[WindowConservation, ...]:
    reference_positions = tuple(
        position for position in positions if position.reference_position is not None
    )
    length = len(reference_positions)
    if length == 0:
        return ()
    if length <= config.window_size:
        starts = (1,)
        window_size = length
    else:
        window_size = config.window_size
        starts_list = list(range(1, length - window_size + 2, config.step_size))
        anchored = length - window_size + 1
        if starts_list[-1] != anchored:
            starts_list.append(anchored)
        starts = tuple(starts_list)

    windows = []
    for start in starts:
        selected = reference_positions[start - 1 : start - 1 + window_size]
        count = len(selected)
        windows.append(
            WindowConservation(
                start,
                start + count - 1,
                count,
                math.fsum(item.major_allele_frequency for item in selected) / count,
                min(item.major_allele_frequency for item in selected),
                math.fsum(item.coverage for item in selected) / count,
                math.fsum(item.gap_frequency for item in selected) / count,
                math.fsum(item.entropy_bits for item in selected) / count,
            )
        )
    return tuple(windows)


def _extract_annotations(
    reference_record: LocalSequenceRecord | None,
    reference_length: int,
) -> tuple[tuple[ReferenceAnnotation, ...], int]:
    if reference_record is None:
        return (), 0
    features = reference_record.metadata.get("features", ())
    if not isinstance(features, (tuple, list)):
        return (), 1
    annotations = []
    skipped = 0
    for feature in features:
        if not isinstance(feature, SeqFeature):
            skipped += 1
            continue
        location = feature.location
        if location is None:
            skipped += 1
            continue
        parts = tuple(location.parts)
        if feature.type == "source":
            skipped += max(1, len(parts))
            continue
        label = _feature_label(feature)
        for part in parts:
            try:
                if part.ref is not None:
                    skipped += 1
                    continue
                start = max(1, int(part.start) + 1)
                end = min(reference_length, int(part.end))
                if start > end:
                    skipped += 1
                    continue
                strand = part.strand if part.strand in (-1, 0, 1, None) else None
                annotations.append(
                    ReferenceAnnotation(str(feature.type), start, end, strand, label)
                )
            except (TypeError, ValueError, AttributeError):
                skipped += 1
    annotations.sort(key=lambda item: (item.start, item.end, item.feature_type, item.label))
    return tuple(annotations), skipped


def _feature_label(feature: SeqFeature) -> str:
    qualifiers = feature.qualifiers if isinstance(feature.qualifiers, Mapping) else {}
    for key in ("gene", "locus_tag", "product"):
        value = qualifiers.get(key)
        if isinstance(value, str) and value:
            return value
        if isinstance(value, (tuple, list)) and value and isinstance(value[0], str) and value[0]:
            return value[0]
    return str(feature.type)


def _position_text(positions: tuple[PositionConservation, ...]) -> str:
    lines = [_POSITION_HEADER]
    for item in positions:
        values = (
            item.alignment_position,
            "" if item.reference_position is None else item.reference_position,
            "" if item.reference_base is None else item.reference_base,
            item.depth,
            item.coverage,
            item.frequency_a,
            item.frequency_c,
            item.frequency_g,
            item.frequency_t,
            item.gap_frequency,
            item.major_allele_frequency,
            item.entropy_bits,
            item.major_consensus,
            item.iupac_consensus,
        )
        lines.append("\t".join(str(value) for value in values) + "\n")
    return "".join(lines)


def _window_text(windows: tuple[WindowConservation, ...]) -> str:
    lines = [_WINDOW_HEADER]
    for item in windows:
        values = (
            item.reference_start,
            item.reference_end,
            item.position_count,
            item.mean_conservation,
            item.minimum_conservation,
            item.mean_coverage,
            item.mean_gap_frequency,
            item.mean_entropy_bits,
        )
        lines.append("\t".join(str(value) for value in values) + "\n")
    return "".join(lines)


def _fasta_text(sequence_id: str, sequence: str) -> str:
    return "" if not sequence else f">{sequence_id}\n{sequence}\n"


def _report(
    *,
    status: Literal["SKIPPED", "COMPLETE"],
    enabled: bool,
    config: ConservationConfig,
    discovery_ids: tuple[str, ...],
    reference_id: str | None,
    positions: tuple[PositionConservation, ...],
    windows: tuple[WindowConservation, ...],
    annotations: tuple[ReferenceAnnotation, ...],
    annotation_skipped: int,
    major_consensus: str,
    iupac_consensus: str,
    artifacts: dict[str, str | None],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "status": status,
        "enabled": enabled,
        "window_parameters": {
            "window_size": config.window_size,
            "step_size": config.step_size,
        },
        "metric_definitions": dict(_METRIC_DEFINITIONS),
        "reference_id": reference_id,
        "discovery_set_ids": list(discovery_ids),
        "counts": {
            "sequences": len(discovery_ids),
            "alignment_columns": len(positions),
            "reference_positions": sum(
                item.reference_position is not None for item in positions
            ),
            "windows": len(windows),
            "annotations": len(annotations),
        },
        "annotation_counts": {
            "published": len(annotations),
            "skipped": annotation_skipped,
        },
        "consensus_lengths": {
            "major": len(major_consensus),
            "iupac": len(iupac_consensus),
        },
        "artifacts": artifacts,
    }


def _publish_complete(
    paths: dict[str, Path],
    position_text: str,
    window_text: str,
    major_text: str,
    iupac_text: str,
    html_text: str,
    report: dict[str, object],
) -> None:
    paths["report"].parent.mkdir(parents=True, exist_ok=True)
    paths["report"].unlink(missing_ok=True)
    _atomic_write_text(paths["positions"], position_text)
    _atomic_write_text(paths["windows"], window_text)
    _atomic_write_text(paths["major"], major_text)
    _atomic_write_text(paths["iupac"], iupac_text)
    _atomic_write_text(paths["html"], html_text)
    _atomic_write_text(
        paths["report"], json.dumps(report, indent=2, sort_keys=True) + "\n"
    )


def _publish_skipped(paths: dict[str, Path], report: dict[str, object]) -> None:
    paths["report"].parent.mkdir(parents=True, exist_ok=True)
    paths["report"].unlink(missing_ok=True)
    for key in ("positions", "windows", "major", "iupac", "html"):
        paths[key].unlink(missing_ok=True)
    _atomic_write_text(
        paths["report"], json.dumps(report, indent=2, sort_keys=True) + "\n"
    )


def _atomic_write_text(destination: Path, content: str) -> None:
    temporary = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.tmp"
    try:
        with temporary.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def _relative_path(path: Path, output_dir: Path) -> str:
    return path.relative_to(output_dir).as_posix()
