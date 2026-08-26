"""Bounded local oligo matching for qPCR assay inclusivity evaluation."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

from qpcr_pipeline.config import InclusivityConfig
from qpcr_pipeline.iupac import (
    IupacError,
    mismatch_positions,
    normalize_iupac,
    reverse_complement_iupac,
)
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import EvaluationSet
from qpcr_pipeline.primer_design import AssayCandidate, DesignedOligo


class InclusivityError(RuntimeError):
    """Raised when an inclusivity input cannot be evaluated safely."""


OligoRole = Literal["FORWARD", "PROBE", "REVERSE"]
Orientation = Literal["FORWARD", "REVERSE_COMPLEMENT"]


@dataclass(frozen=True, slots=True)
class OligoMatch:
    assay_id: str
    sequence_id: str
    role: OligoRole
    orientation: Orientation
    hit_rank: int
    source_start: int
    source_end: int
    expected_start: int
    expected_end: int
    displacement: int
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    exact_match: bool
    three_prime_mismatch: bool
    probe_mismatch: bool
    compatible: bool
    selected: bool


@dataclass(frozen=True, slots=True)
class ProposedOligoCompatibility:
    role: OligoRole
    effective_sequence: str
    exact_match: bool
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    three_prime_mismatch: bool
    probe_mismatch: bool
    compatible: bool


@dataclass(frozen=True, slots=True)
class AssayInclusivity:
    assay_id: str
    sequence_id: str
    orientation: Orientation | None
    geometry_found: bool
    source_amplicon_start: int | None
    source_amplicon_end: int | None
    amplicon_size: int | None
    forward_match: OligoMatch | None
    probe_match: OligoMatch | None
    reverse_match: OligoMatch | None
    original_compatible: bool
    proposed_forward: ProposedOligoCompatibility | None
    proposed_probe: ProposedOligoCompatibility | None
    proposed_reverse: ProposedOligoCompatibility | None
    proposed_compatible: bool


@dataclass(frozen=True, slots=True)
class OligoVariation:
    assay_id: str
    role: OligoRole
    oligo_position: int
    original_symbol: str
    original_support: tuple[str, ...]
    observed_symbol: str
    observed_support: tuple[str, ...]
    affected_sequence_ids: tuple[str, ...]
    affected_sequence_count: int
    affected_fraction: float | None
    primer_3_prime_position: bool


@dataclass(frozen=True, slots=True)
class DegeneracyProposal:
    assay_id: str
    role: OligoRole
    original_sequence: str
    proposed_sequence: str
    status: Literal["ACCEPTED", "UNCHANGED", "REJECTED"]
    reason: str
    original_degeneracy: int
    proposed_degeneracy: int
    changed_positions: tuple[int, ...]
    binding_site_count: int
    original_exact_count: int
    original_exact_fraction: float | None
    proposed_exact_count: int
    proposed_exact_fraction: float | None


@dataclass(frozen=True, slots=True)
class InclusivityResult:
    status: Literal["SKIPPED", "COMPLETE"]
    evaluation_sequence_ids: tuple[str, ...]
    oligo_matches: tuple[OligoMatch, ...]
    assay_results: tuple[AssayInclusivity, ...]
    variations: tuple[OligoVariation, ...]
    proposals: tuple[DegeneracyProposal, ...]
    oligo_matches_path: Path | None
    assay_inclusivity_path: Path | None
    oligo_variations_path: Path | None
    degeneracy_proposals_path: Path | None
    report_path: Path


@dataclass(frozen=True, slots=True)
class _Hit:
    public: OligoMatch
    oriented_start: int
    oriented_end: int
    target_in_synthesis_orientation: str


@dataclass(frozen=True, slots=True)
class _SelectedBinding:
    assay: AssayCandidate
    sequence_id: str
    orientation: Orientation | None
    forward: _Hit | None
    probe: _Hit | None
    reverse: _Hit | None
    geometry_found: bool
    source_amplicon_start: int | None
    source_amplicon_end: int | None
    amplicon_size: int | None
    original_compatible: bool
    retained_hits: tuple[_Hit, ...]


def _oligo_for_role(assay: AssayCandidate, role: OligoRole) -> DesignedOligo:
    return {
        "FORWARD": assay.forward_primer,
        "PROBE": assay.probe,
        "REVERSE": assay.reverse_primer,
    }[role]


def _enumerate_hits(
    *,
    assay_id: str,
    sequence_id: str,
    oriented_sequence: str,
    orientation: Orientation,
    role: OligoRole,
    oligo: DesignedOligo,
    config: InclusivityConfig,
) -> tuple[_Hit, ...]:
    """Return the bounded, deterministically ordered local sites for one oligo."""
    try:
        sequence = normalize_iupac(
            oriented_sequence,
            context=(
                f"record {sequence_id!r} for assay {assay_id!r} role {role!r}"
            ),
        )
        oligo_sequence = normalize_iupac(
            oligo.sequence,
            context=f"assay {assay_id!r} role {role!r} oligo",
        )
    except IupacError as error:
        raise InclusivityError(
            f"Invalid IUPAC input for assay {assay_id!r}, record {sequence_id!r}, "
            f"role {role!r}: {error}"
        ) from error

    oligo_length = len(oligo_sequence)
    record_length = len(sequence)
    if record_length < oligo_length:
        return ()

    window_start = max(1, oligo.reference_start - config.search_flank)
    window_end = min(record_length, oligo.reference_end + config.search_flank)
    first_start = window_start
    last_start = min(record_length - oligo_length + 1, window_end - oligo_length + 1)
    if first_start > last_start:
        return ()

    hits: list[tuple[_Hit, int, str]] = []
    for oriented_start in range(first_start, last_start + 1):
        oriented_end = oriented_start + oligo_length - 1
        target_segment = sequence[oriented_start - 1 : oriented_end]
        try:
            target_in_synthesis_orientation = (
                reverse_complement_iupac(target_segment)
                if role == "REVERSE"
                else target_segment
            )
            positions = mismatch_positions(oligo_sequence, target_in_synthesis_orientation)
        except IupacError as error:
            raise InclusivityError(
                f"Invalid IUPAC input for assay {assay_id!r}, record {sequence_id!r}, "
                f"role {role!r}: {error}"
            ) from error

        three_prime_mismatch = role != "PROBE" and any(
            position > oligo_length - config.primer_3_prime_bases
            for position in positions
        )
        compatible = (
            len(positions) <= (
                config.max_probe_mismatches
                if role == "PROBE"
                else config.max_primer_mismatches
            )
            and not (
                role != "PROBE"
                and config.reject_primer_3_prime_mismatch
                and three_prime_mismatch
            )
        )
        if orientation == "REVERSE_COMPLEMENT":
            source_start = record_length - oriented_end + 1
            source_end = record_length - oriented_start + 1
        else:
            source_start = oriented_start
            source_end = oriented_end

        public = OligoMatch(
            assay_id=assay_id,
            sequence_id=sequence_id,
            role=role,
            orientation=orientation,
            hit_rank=0,
            source_start=source_start,
            source_end=source_end,
            expected_start=oligo.reference_start,
            expected_end=oligo.reference_end,
            displacement=abs(oriented_start - oligo.reference_start),
            mismatch_positions=positions,
            mismatch_count=len(positions),
            exact_match=not positions,
            three_prime_mismatch=three_prime_mismatch,
            probe_mismatch=role == "PROBE" and bool(positions),
            compatible=compatible,
            selected=False,
        )
        hits.append(
            (
                _Hit(
                    public=public,
                    oriented_start=oriented_start,
                    oriented_end=oriented_end,
                    target_in_synthesis_orientation=target_in_synthesis_orientation,
                ),
                sum(
                    position > oligo_length - config.primer_3_prime_bases
                    for position in positions
                )
                if role != "PROBE"
                else 0,
                target_segment,
            )
        )

    hits.sort(
        key=lambda item: (
            item[0].public.mismatch_count,
            item[1],
            item[0].public.displacement,
            item[0].oriented_start,
            item[2],
        )
    )
    return tuple(
        replace(hit, public=replace(hit.public, hit_rank=rank))
        for rank, (hit, _, _) in enumerate(hits[: config.max_hits_per_oligo], 1)
    )


def _validate_enabled_inputs(
    records: tuple[LocalSequenceRecord, ...],
    evaluation_set: EvaluationSet,
    assays: tuple[AssayCandidate, ...],
) -> tuple[LocalSequenceRecord, ...]:
    """Validate the exact inputs required for local inclusivity evaluation."""
    if not isinstance(evaluation_set, EvaluationSet):
        raise InclusivityError("Inclusivity evaluation requires an EvaluationSet.")
    evaluation_ids = evaluation_set.sequence_ids
    if not isinstance(evaluation_ids, tuple):
        raise InclusivityError("Evaluation Set IDs must be a tuple.")
    if any(not isinstance(sequence_id, str) or not sequence_id.strip() for sequence_id in evaluation_ids):
        raise InclusivityError("Evaluation Set IDs must be non-blank strings.")
    if len(set(evaluation_ids)) != len(evaluation_ids):
        raise InclusivityError("Evaluation Set IDs must be unique.")

    if not isinstance(records, tuple):
        raise InclusivityError("Evaluation records must be a tuple.")
    for index, record in enumerate(records, 1):
        if not isinstance(record, LocalSequenceRecord):
            raise InclusivityError(
                f"Evaluation record {index} must be a LocalSequenceRecord."
            )
        if not isinstance(record.sequence_id, str) or not record.sequence_id.strip():
            raise InclusivityError(
                f"Evaluation record {index} must have a non-blank string ID."
            )
    record_ids = tuple(record.sequence_id for record in records)
    if record_ids != evaluation_ids:
        if len(record_ids) != len(evaluation_ids) or set(record_ids) != set(evaluation_ids):
            raise InclusivityError("Evaluation records must match Evaluation Set IDs exactly once.")
        raise InclusivityError("Evaluation records must follow Evaluation Set ID order.")

    if not isinstance(assays, tuple):
        raise InclusivityError("Inclusivity assays must be a tuple.")
    assay_ids: list[str] = []
    for index, assay in enumerate(assays, 1):
        if not isinstance(assay, AssayCandidate):
            raise InclusivityError(f"Inclusivity assay {index} must be an AssayCandidate.")
        if not isinstance(assay.assay_id, str) or not assay.assay_id.strip():
            raise InclusivityError(f"Inclusivity assay {index} has a blank assay ID.")
        assay_ids.append(assay.assay_id)
        _validate_assay(assay)
    if len(set(assay_ids)) != len(assay_ids):
        raise InclusivityError("Inclusivity assay IDs must be unique.")
    return records


def _validate_assay(assay: AssayCandidate) -> None:
    oligos = (
        ("FORWARD", assay.forward_primer),
        ("PROBE", assay.probe),
        ("REVERSE", assay.reverse_primer),
    )
    for role, oligo in oligos:
        if not isinstance(oligo, DesignedOligo):
            raise InclusivityError(
                f"Assay {assay.assay_id!r} role {role!r} must be a DesignedOligo."
            )
        if not isinstance(oligo.sequence, str):
            raise InclusivityError(
                f"Assay {assay.assay_id!r} role {role!r} sequence must be a string."
            )
        try:
            normalized = normalize_iupac(
                oligo.sequence, context=f"assay {assay.assay_id!r} role {role!r} oligo"
            )
        except IupacError as error:
            raise InclusivityError(
                f"Invalid IUPAC input for assay {assay.assay_id!r}, role {role!r}: {error}"
            ) from error
        for name, value in (
            ("reference_start", oligo.reference_start),
            ("reference_end", oligo.reference_end),
            ("length", oligo.length),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise InclusivityError(
                    f"Assay {assay.assay_id!r} role {role!r} {name} must be positive."
                )
        if oligo.length != len(normalized):
            raise InclusivityError(
                f"Assay {assay.assay_id!r} role {role!r} length must match its sequence length."
            )
        if oligo.reference_end - oligo.reference_start + 1 != oligo.length:
            raise InclusivityError(
                f"Assay {assay.assay_id!r} role {role!r} coordinates must match its length."
            )

    forward = assay.forward_primer
    probe = assay.probe
    reverse = assay.reverse_primer
    if not (
        forward.reference_end < probe.reference_start <= probe.reference_end < reverse.reference_start
    ):
        raise InclusivityError(
            f"Assay {assay.assay_id!r} has invalid forward/probe/reverse geometry."
        )
    if isinstance(assay.product_size, bool) or not isinstance(assay.product_size, int) or assay.product_size < 1:
        raise InclusivityError(f"Assay {assay.assay_id!r} product_size must be a positive integer.")
    expected_size = reverse.reference_end - forward.reference_start + 1
    if assay.product_size != expected_size:
        raise InclusivityError(
            f"Assay {assay.assay_id!r} product_size must match its primer coordinates."
        )


def _select_binding(
    record: LocalSequenceRecord,
    assay: AssayCandidate,
    config: InclusivityConfig,
) -> _SelectedBinding:
    """Select one deterministic complete assay binding or per-role fallback hits."""
    if not isinstance(record, LocalSequenceRecord):
        raise InclusivityError("Inclusivity selection requires a LocalSequenceRecord.")
    if not isinstance(assay, AssayCandidate):
        raise InclusivityError("Inclusivity selection requires an AssayCandidate.")
    if not isinstance(config, InclusivityConfig):
        raise InclusivityError("Inclusivity selection requires an InclusivityConfig.")
    _validate_assay(assay)
    if not isinstance(record.sequence_id, str) or not record.sequence_id.strip():
        raise InclusivityError(f"Assay {assay.assay_id!r} has a record with a blank ID.")
    if not isinstance(record.sequence, str):
        raise InclusivityError(
            f"Record {record.sequence_id!r} for assay {assay.assay_id!r} sequence must be a string."
        )
    try:
        sequence = normalize_iupac(
            record.sequence, context=f"record {record.sequence_id!r} for assay {assay.assay_id!r}"
        )
    except IupacError as error:
        raise InclusivityError(
            f"Invalid IUPAC input for assay {assay.assay_id!r}, record {record.sequence_id!r}: {error}"
        ) from error

    oriented_sequences: tuple[tuple[Orientation, str], ...] = (
        ("FORWARD", sequence),
        ("REVERSE_COMPLEMENT", reverse_complement_iupac(sequence)),
    )
    roles: tuple[OligoRole, ...] = ("FORWARD", "PROBE", "REVERSE")
    hits_by_role: dict[OligoRole, tuple[_Hit, ...]] = {}
    retained: list[_Hit] = []
    for role in roles:
        role_hits: list[_Hit] = []
        for orientation, oriented_sequence in oriented_sequences:
            hits = _enumerate_hits(
                assay_id=assay.assay_id,
                sequence_id=record.sequence_id,
                oriented_sequence=oriented_sequence,
                orientation=orientation,
                role=role,
                oligo=_oligo_for_role(assay, role),
                config=config,
            )
            role_hits.extend(hits)
            retained.extend(hits)
        hits_by_role[role] = tuple(role_hits)

    candidates: list[tuple[tuple[object, ...], _Hit, _Hit, _Hit]] = []
    for orientation_rank, orientation in enumerate(("FORWARD", "REVERSE_COMPLEMENT")):
        forward_hits = tuple(hit for hit in hits_by_role["FORWARD"] if hit.public.orientation == orientation)
        probe_hits = tuple(hit for hit in hits_by_role["PROBE"] if hit.public.orientation == orientation)
        reverse_hits = tuple(hit for hit in hits_by_role["REVERSE"] if hit.public.orientation == orientation)
        for forward in forward_hits:
            for probe in probe_hits:
                for reverse in reverse_hits:
                    amplicon_size = reverse.oriented_end - forward.oriented_start + 1
                    if not (
                        forward.oriented_end < probe.oriented_start <= probe.oriented_end < reverse.oriented_start
                        and assay.product_size - config.max_amplicon_size_delta <= amplicon_size <= assay.product_size + config.max_amplicon_size_delta
                    ):
                        continue
                    candidates.append((
                        (
                            forward.public.mismatch_count + probe.public.mismatch_count + reverse.public.mismatch_count,
                            _primer_three_prime_mismatch_count(forward, config) + _primer_three_prime_mismatch_count(reverse, config),
                            abs(amplicon_size - assay.product_size),
                            forward.public.displacement + probe.public.displacement + reverse.public.displacement,
                            orientation_rank,
                            forward.oriented_start,
                            forward.oriented_end,
                            probe.oriented_start,
                            probe.oriented_end,
                            reverse.oriented_start,
                            reverse.oriented_end,
                        ),
                        forward,
                        probe,
                        reverse,
                    ))

    if candidates:
        _, forward, probe, reverse = min(candidates, key=lambda item: item[0])
        selected_ids = {id(forward), id(probe), id(reverse)}
        selected_retained = _mark_selected(retained, selected_ids)
        selected_by_id = {id(original): marked for original, marked in zip(retained, selected_retained)}
        selected_forward = selected_by_id[id(forward)]
        selected_probe = selected_by_id[id(probe)]
        selected_reverse = selected_by_id[id(reverse)]
        return _SelectedBinding(
            assay=assay,
            sequence_id=record.sequence_id,
            orientation=selected_forward.public.orientation,
            forward=selected_forward,
            probe=selected_probe,
            reverse=selected_reverse,
            geometry_found=True,
            source_amplicon_start=min(selected_forward.public.source_start, selected_reverse.public.source_start),
            source_amplicon_end=max(selected_forward.public.source_end, selected_reverse.public.source_end),
            amplicon_size=selected_reverse.oriented_end - selected_forward.oriented_start + 1,
            original_compatible=(
                selected_forward.public.compatible
                and selected_probe.public.compatible
                and selected_reverse.public.compatible
            ),
            retained_hits=selected_retained,
        )

    fallback = tuple(
        _best_fallback(hits_by_role[role], role, config)
        for role in roles
    )
    selected_ids = {id(hit) for hit in fallback if hit is not None}
    selected_retained = _mark_selected(retained, selected_ids)
    selected_by_id = {id(original): marked for original, marked in zip(retained, selected_retained)}
    forward, probe, reverse = tuple(
        selected_by_id[id(hit)] if hit is not None else None
        for hit in fallback
    )
    return _SelectedBinding(
        assay=assay,
        sequence_id=record.sequence_id,
        orientation=None,
        forward=forward,
        probe=probe,
        reverse=reverse,
        geometry_found=False,
        source_amplicon_start=None,
        source_amplicon_end=None,
        amplicon_size=None,
        original_compatible=False,
        retained_hits=selected_retained,
    )


def _primer_three_prime_mismatch_count(hit: _Hit, config: InclusivityConfig) -> int:
    if hit.public.role == "PROBE":
        return 0
    length = len(hit.target_in_synthesis_orientation)
    return sum(
        position > length - config.primer_3_prime_bases
        for position in hit.public.mismatch_positions
    )


def _best_fallback(
    hits: tuple[_Hit, ...], role: OligoRole, config: InclusivityConfig
) -> _Hit | None:
    if not hits:
        return None
    return min(
        hits,
        key=lambda hit: (
            hit.public.mismatch_count,
            _primer_three_prime_mismatch_count(hit, config) if role != "PROBE" else 0,
            hit.public.displacement,
            0 if hit.public.orientation == "FORWARD" else 1,
            hit.oriented_start,
            hit.oriented_end,
            hit.target_in_synthesis_orientation,
        ),
    )


def _mark_selected(hits: list[_Hit], selected_ids: set[int]) -> tuple[_Hit, ...]:
    return tuple(
        replace(hit, public=replace(hit.public, selected=id(hit) in selected_ids))
        for hit in hits
    )


def _evaluate_original(
    records: tuple[LocalSequenceRecord, ...],
    evaluation_set: EvaluationSet,
    assays: tuple[AssayCandidate, ...],
    config: InclusivityConfig,
) -> tuple[_SelectedBinding, ...]:
    validated = _validate_enabled_inputs(records, evaluation_set, assays)
    return tuple(
        _select_binding(record, assay, config)
        for assay in assays
        for record in validated
    )
