"""Bounded local oligo matching for qPCR assay inclusivity evaluation."""

from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass, replace
import json
from pathlib import Path
from typing import Literal
import uuid

from qpcr_pipeline.config import InclusivityConfig, validate_inclusivity_config
from qpcr_pipeline.iupac import (
    IupacError,
    iupac_support,
    minimal_iupac_symbol,
    mismatch_positions,
    normalize_iupac,
    reverse_complement_iupac,
    sequence_degeneracy,
)
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import EvaluationSet
from qpcr_pipeline.primer_design import AssayCandidate, DesignedOligo, PrimerDesignResult


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


OLIGO_MATCH_COLUMNS = (
    "assay_id", "sequence_id", "role", "orientation", "hit_rank",
    "source_start", "source_end", "expected_start", "expected_end",
    "displacement", "mismatch_positions", "mismatch_count", "exact_match",
    "three_prime_mismatch", "probe_mismatch", "compatible", "selected",
)
ASSAY_COLUMNS = (
    "assay_id", "sequence_id", "orientation", "geometry_found",
    "source_amplicon_start", "source_amplicon_end", "amplicon_size",
    "forward_hit_rank", "probe_hit_rank", "reverse_hit_rank",
    "original_forward_compatible", "original_probe_compatible",
    "original_reverse_compatible", "original_compatible",
    "proposed_forward_compatible", "proposed_probe_compatible",
    "proposed_reverse_compatible", "proposed_compatible",
)
VARIATION_COLUMNS = (
    "assay_id", "role", "oligo_position", "original_symbol",
    "original_support", "observed_symbol", "observed_support",
    "affected_sequence_ids", "affected_sequence_count", "affected_fraction",
    "primer_3_prime_position",
)
PROPOSAL_COLUMNS = (
    "assay_id", "role", "original_sequence", "proposed_sequence", "status",
    "reason", "original_degeneracy", "proposed_degeneracy", "changed_positions",
    "binding_site_count", "original_exact_count", "original_exact_fraction",
    "proposed_exact_count", "proposed_exact_fraction",
)

_DATA_ARTIFACT_NAMES = (
    "oligo_matches.tsv",
    "assay_inclusivity.tsv",
    "oligo_variations.tsv",
    "degeneracy_proposals.tsv",
)


@dataclass(frozen=True, slots=True)
class _Compatibility:
    exact_match: bool
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    three_prime_mismatch: bool
    probe_mismatch: bool
    compatible: bool


@dataclass(frozen=True, slots=True)
class _Hit:
    public: OligoMatch
    oriented_start: int
    oriented_end: int
    oriented_target_segment: str
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


def _compatibility(
    effective_sequence: str,
    target_in_synthesis_orientation: str,
    role: OligoRole,
    config: InclusivityConfig,
) -> _Compatibility:
    positions = mismatch_positions(
        effective_sequence, target_in_synthesis_orientation
    )
    three_prime_mismatch = role != "PROBE" and any(
        position > len(effective_sequence) - config.primer_3_prime_bases
        for position in positions
    )
    compatible = (
        len(positions)
        <= (
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
    return _Compatibility(
        exact_match=not positions,
        mismatch_positions=positions,
        mismatch_count=len(positions),
        three_prime_mismatch=three_prime_mismatch,
        probe_mismatch=role == "PROBE" and bool(positions),
        compatible=compatible,
    )


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
            compatibility = _compatibility(
                oligo_sequence,
                target_in_synthesis_orientation,
                role,
                config,
            )
        except IupacError as error:
            raise InclusivityError(
                f"Invalid IUPAC input for assay {assay_id!r}, record {sequence_id!r}, "
                f"role {role!r}: {error}"
            ) from error

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
            mismatch_positions=compatibility.mismatch_positions,
            mismatch_count=compatibility.mismatch_count,
            exact_match=compatibility.exact_match,
            three_prime_mismatch=compatibility.three_prime_mismatch,
            probe_mismatch=compatibility.probe_mismatch,
            compatible=compatibility.compatible,
            selected=False,
        )
        hits.append(
            (
                _Hit(
                    public=public,
                    oriented_start=oriented_start,
                    oriented_end=oriented_end,
                    oriented_target_segment=target_segment,
                    target_in_synthesis_orientation=target_in_synthesis_orientation,
                ),
                sum(
                    position > oligo_length - config.primer_3_prime_bases
                    for position in compatibility.mismatch_positions
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
            hit.oriented_target_segment,
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


def _variation_rows(
    selected: tuple[_SelectedBinding, ...],
    evaluation_set: EvaluationSet,
    config: InclusivityConfig,
) -> tuple[OligoVariation, ...]:
    """Aggregate positional support changes from selected geometric sites only."""
    assay_order: list[AssayCandidate] = []
    seen_assays: set[str] = set()
    for binding in selected:
        if binding.assay.assay_id not in seen_assays:
            assay_order.append(binding.assay)
            seen_assays.add(binding.assay.assay_id)

    evaluation_order = {
        sequence_id: index
        for index, sequence_id in enumerate(evaluation_set.sequence_ids)
    }
    denominator = len(evaluation_set.sequence_ids)
    rows: list[OligoVariation] = []
    for assay in assay_order:
        assay_bindings = tuple(
            binding
            for binding in selected
            if binding.assay.assay_id == assay.assay_id and binding.geometry_found
        )
        for role in ("FORWARD", "PROBE", "REVERSE"):
            oligo_sequence = normalize_iupac(
                _oligo_for_role(assay, role).sequence,
                context=f"assay {assay.assay_id!r} role {role!r} oligo",
            )
            sites = tuple(
                (binding.sequence_id, getattr(binding, role.lower()))
                for binding in assay_bindings
                if getattr(binding, role.lower()) is not None
            )
            for position, original_symbol in enumerate(oligo_sequence, 1):
                original_support = iupac_support(original_symbol)
                affected = tuple(
                    sequence_id
                    for sequence_id, hit in sites
                    if iupac_support(
                        hit.target_in_synthesis_orientation[position - 1]
                    ) != original_support
                )
                if not affected:
                    continue
                affected = tuple(
                    sorted(
                        affected,
                        key=lambda sequence_id: evaluation_order[sequence_id],
                    )
                )
                observed_support = frozenset().union(
                    *(
                        iupac_support(
                            hit.target_in_synthesis_orientation[position - 1]
                        )
                        for _, hit in sites
                    )
                )
                rows.append(
                    OligoVariation(
                        assay_id=assay.assay_id,
                        role=role,
                        oligo_position=position,
                        original_symbol=original_symbol,
                        original_support=tuple(sorted(original_support)),
                        observed_symbol=minimal_iupac_symbol(observed_support),
                        observed_support=tuple(sorted(observed_support)),
                        affected_sequence_ids=affected,
                        affected_sequence_count=len(affected),
                        affected_fraction=(len(affected) / denominator if denominator else None),
                        primer_3_prime_position=(
                            role != "PROBE"
                            and position > len(oligo_sequence) - config.primer_3_prime_bases
                        ),
                    )
                )
    return tuple(rows)


def _proposal_for_role(
    assay: AssayCandidate,
    role: OligoRole,
    selected: tuple[_SelectedBinding, ...],
    evaluation_set: EvaluationSet,
    config: InclusivityConfig,
) -> DegeneracyProposal:
    """Choose the best cap-bounded expansion for one immutable assay oligo."""
    del evaluation_set
    original = normalize_iupac(
        _oligo_for_role(assay, role).sequence,
        context=f"assay {assay.assay_id!r} role {role!r} oligo",
    )
    original_degeneracy = sequence_degeneracy(original)
    sites = tuple(
        getattr(binding, role.lower())
        for binding in selected
        if binding.assay.assay_id == assay.assay_id
        and binding.geometry_found
        and getattr(binding, role.lower()) is not None
    )
    binding_site_count = len(sites)

    def exact_count(candidate: str) -> int:
        return sum(
            not mismatch_positions(candidate, hit.target_in_synthesis_orientation)
            for hit in sites
        )

    original_exact_count = exact_count(original)
    original_fraction = (
        original_exact_count / binding_site_count if binding_site_count else None
    )

    def unchanged(status: Literal["UNCHANGED", "REJECTED"], reason: str) -> DegeneracyProposal:
        return DegeneracyProposal(
            assay_id=assay.assay_id,
            role=role,
            original_sequence=original,
            proposed_sequence=original,
            status=status,
            reason=reason,
            original_degeneracy=original_degeneracy,
            proposed_degeneracy=original_degeneracy,
            changed_positions=(),
            binding_site_count=binding_site_count,
            original_exact_count=original_exact_count,
            original_exact_fraction=original_fraction,
            proposed_exact_count=original_exact_count,
            proposed_exact_fraction=original_fraction,
        )

    if not sites:
        return unchanged("UNCHANGED", "NO_GEOMETRIC_SITES")

    expansions: list[tuple[int, str]] = []
    variation_found = False
    for position, original_symbol in enumerate(original, 1):
        original_support = iupac_support(original_symbol)
        observed_support = frozenset().union(
            *(
                iupac_support(hit.target_in_synthesis_orientation[position - 1])
                for hit in sites
            )
        )
        if any(
            iupac_support(hit.target_in_synthesis_orientation[position - 1])
            != original_support
            for hit in sites
        ):
            variation_found = True
        expanded_symbol = minimal_iupac_symbol(original_support | observed_support)
        if expanded_symbol != original_symbol:
            expansions.append((position, expanded_symbol))

    if not variation_found:
        return unchanged("UNCHANGED", "NO_VARIATION")

    def expanded_sequence(changes: tuple[tuple[int, str], ...]) -> str:
        symbols = list(original)
        for position, symbol in changes:
            symbols[position - 1] = symbol
        return "".join(symbols)

    unrestricted = expanded_sequence(tuple(expansions))
    if exact_count(unrestricted) <= original_exact_count:
        return unchanged("UNCHANGED", "NO_IMPROVEMENT")

    protected_start = len(original) - config.primer_3_prime_bases + 1
    allowed_expansions = tuple(
        expansion
        for expansion in expansions
        if not (
            role != "PROBE"
            and not config.allow_primer_3_prime_degeneracy
            and expansion[0] >= protected_start
        )
    )
    allowed_full = expanded_sequence(allowed_expansions)
    if exact_count(allowed_full) <= original_exact_count:
        return unchanged("REJECTED", "REJECTED_3_PRIME")

    limit = (
        config.max_probe_degeneracy
        if role == "PROBE"
        else config.max_primer_degeneracy
    )
    best: tuple[tuple[object, ...], str, tuple[int, ...], int] | None = None

    def search(index: int, candidate: str, changed: tuple[int, ...]) -> None:
        nonlocal best
        if sequence_degeneracy(candidate) > limit:
            return
        if index == len(allowed_expansions):
            candidate_exact_count = exact_count(candidate)
            if candidate_exact_count <= original_exact_count:
                return
            key = (
                -candidate_exact_count,
                sequence_degeneracy(candidate),
                len(changed),
                candidate,
            )
            if best is None or key < best[0]:
                best = (key, candidate, changed, candidate_exact_count)
            return

        position, expanded_symbol = allowed_expansions[index]
        search(index + 1, candidate, changed)
        symbols = list(candidate)
        symbols[position - 1] = expanded_symbol
        search(index + 1, "".join(symbols), changed + (position,))

    search(0, original, ())
    if best is None:
        return unchanged("REJECTED", "REJECTED_LIMIT")

    _, proposed, changed_positions, proposed_exact_count = best
    return DegeneracyProposal(
        assay_id=assay.assay_id,
        role=role,
        original_sequence=original,
        proposed_sequence=proposed,
        status="ACCEPTED",
        reason="ACCEPTED_IMPROVEMENT",
        original_degeneracy=original_degeneracy,
        proposed_degeneracy=sequence_degeneracy(proposed),
        changed_positions=changed_positions,
        binding_site_count=binding_site_count,
        original_exact_count=original_exact_count,
        original_exact_fraction=original_fraction,
        proposed_exact_count=proposed_exact_count,
        proposed_exact_fraction=proposed_exact_count / binding_site_count,
    )


def _proposals(
    assays: tuple[AssayCandidate, ...],
    selected: tuple[_SelectedBinding, ...],
    evaluation_set: EvaluationSet,
    config: InclusivityConfig,
) -> tuple[DegeneracyProposal, ...]:
    """Return exactly three role proposals for every assay in input order."""
    return tuple(
        _proposal_for_role(assay, role, selected, evaluation_set, config)
        for assay in assays
        for role in ("FORWARD", "PROBE", "REVERSE")
    )


def _assay_results_with_proposals(
    selected: tuple[_SelectedBinding, ...],
    proposals: tuple[DegeneracyProposal, ...],
    config: InclusivityConfig,
) -> tuple[AssayInclusivity, ...]:
    """Recompute proposed compatibility at each existing selected geometry."""
    proposal_by_role = {
        (proposal.assay_id, proposal.role): proposal for proposal in proposals
    }
    results: list[AssayInclusivity] = []
    for binding in selected:
        proposed_by_role: dict[OligoRole, ProposedOligoCompatibility | None] = {}
        for role in ("FORWARD", "PROBE", "REVERSE"):
            hit = getattr(binding, role.lower())
            if not binding.geometry_found or hit is None:
                proposed_by_role[role] = None
                continue
            proposal = proposal_by_role[(binding.assay.assay_id, role)]
            effective_sequence = (
                proposal.proposed_sequence
                if proposal.status == "ACCEPTED"
                else proposal.original_sequence
            )
            detail = _compatibility(
                effective_sequence,
                hit.target_in_synthesis_orientation,
                role,
                config,
            )
            proposed_by_role[role] = ProposedOligoCompatibility(
                role=role,
                effective_sequence=effective_sequence,
                exact_match=detail.exact_match,
                mismatch_positions=detail.mismatch_positions,
                mismatch_count=detail.mismatch_count,
                three_prime_mismatch=detail.three_prime_mismatch,
                probe_mismatch=detail.probe_mismatch,
                compatible=detail.compatible,
            )

        proposed_roles = tuple(proposed_by_role[role] for role in ("FORWARD", "PROBE", "REVERSE"))
        results.append(
            AssayInclusivity(
                assay_id=binding.assay.assay_id,
                sequence_id=binding.sequence_id,
                orientation=binding.orientation,
                geometry_found=binding.geometry_found,
                source_amplicon_start=binding.source_amplicon_start,
                source_amplicon_end=binding.source_amplicon_end,
                amplicon_size=binding.amplicon_size,
                forward_match=binding.forward.public if binding.forward is not None else None,
                probe_match=binding.probe.public if binding.probe is not None else None,
                reverse_match=binding.reverse.public if binding.reverse is not None else None,
                original_compatible=binding.original_compatible,
                proposed_forward=proposed_by_role["FORWARD"],
                proposed_probe=proposed_by_role["PROBE"],
                proposed_reverse=proposed_by_role["REVERSE"],
                proposed_compatible=(
                    binding.geometry_found
                    and all(
                        proposed is not None and proposed.compatible
                        for proposed in proposed_roles
                    )
                ),
            )
        )
    return tuple(results)


def _artifact_paths(output_dir: Path) -> dict[str, Path]:
    directory = output_dir / "inclusivity"
    return {
        "oligo_matches": directory / "oligo_matches.tsv",
        "assay_inclusivity": directory / "assay_inclusivity.tsv",
        "oligo_variations": directory / "oligo_variations.tsv",
        "degeneracy_proposals": directory / "degeneracy_proposals.tsv",
        "report": directory / "inclusivity_report.json",
    }


def _relative_path(path: Path, output_dir: Path) -> str:
    return path.relative_to(output_dir).as_posix()


def _tsv_scalar(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, tuple):
        return ",".join(_tsv_scalar(item) for item in value)
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _tsv_text(columns: tuple[str, ...], rows: tuple[tuple[object, ...], ...]) -> str:
    return "".join(
        ("\t".join(columns) + "\n",)
        + tuple("\t".join(_tsv_scalar(value) for value in row) + "\n" for row in rows)
    )


def _oligo_match_text(matches: tuple[OligoMatch, ...]) -> str:
    return _tsv_text(
        OLIGO_MATCH_COLUMNS,
        tuple(
            (
                match.assay_id, match.sequence_id, match.role, match.orientation,
                match.hit_rank, match.source_start, match.source_end,
                match.expected_start, match.expected_end, match.displacement,
                match.mismatch_positions, match.mismatch_count, match.exact_match,
                match.three_prime_mismatch, match.probe_mismatch, match.compatible,
                match.selected,
            )
            for match in matches
        ),
    )


def _assay_text(results: tuple[AssayInclusivity, ...]) -> str:
    return _tsv_text(
        ASSAY_COLUMNS,
        tuple(
            (
                result.assay_id, result.sequence_id, result.orientation,
                result.geometry_found, result.source_amplicon_start,
                result.source_amplicon_end, result.amplicon_size,
                result.forward_match.hit_rank if result.forward_match else None,
                result.probe_match.hit_rank if result.probe_match else None,
                result.reverse_match.hit_rank if result.reverse_match else None,
                result.forward_match.compatible if result.forward_match else None,
                result.probe_match.compatible if result.probe_match else None,
                result.reverse_match.compatible if result.reverse_match else None,
                result.original_compatible,
                result.proposed_forward.compatible if result.proposed_forward else None,
                result.proposed_probe.compatible if result.proposed_probe else None,
                result.proposed_reverse.compatible if result.proposed_reverse else None,
                result.proposed_compatible,
            )
            for result in results
        ),
    )


def _variation_text(variations: tuple[OligoVariation, ...]) -> str:
    return _tsv_text(
        VARIATION_COLUMNS,
        tuple(
            (
                variation.assay_id, variation.role, variation.oligo_position,
                variation.original_symbol, variation.original_support,
                variation.observed_symbol, variation.observed_support,
                variation.affected_sequence_ids, variation.affected_sequence_count,
                variation.affected_fraction, variation.primer_3_prime_position,
            )
            for variation in variations
        ),
    )


def _proposal_text(proposals: tuple[DegeneracyProposal, ...]) -> str:
    return _tsv_text(
        PROPOSAL_COLUMNS,
        tuple(
            (
                proposal.assay_id, proposal.role, proposal.original_sequence,
                proposal.proposed_sequence, proposal.status, proposal.reason,
                proposal.original_degeneracy, proposal.proposed_degeneracy,
                proposal.changed_positions, proposal.binding_site_count,
                proposal.original_exact_count, proposal.original_exact_fraction,
                proposal.proposed_exact_count, proposal.proposed_exact_fraction,
            )
            for proposal in proposals
        ),
    )


def _json_value(value: object) -> object:
    if is_dataclass(value) and not isinstance(value, type):
        return {field.name: _json_value(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    return value


def _counts(
    evaluation_sequence_ids: tuple[str, ...],
    assays: tuple[AssayCandidate, ...],
    assay_results: tuple[AssayInclusivity, ...],
    oligo_matches: tuple[OligoMatch, ...],
    variations: tuple[OligoVariation, ...],
    proposals: tuple[DegeneracyProposal, ...],
) -> dict[str, int]:
    return {
        "evaluation_sequences": len(evaluation_sequence_ids),
        "assays": len(assays),
        "assay_evaluations": len(assay_results),
        "retained_oligo_hits": len(oligo_matches),
        "variations": len(variations),
        "proposals": len(proposals),
        "accepted_proposals": sum(proposal.status == "ACCEPTED" for proposal in proposals),
        "original_compatible": sum(result.original_compatible for result in assay_results),
        "proposed_compatible": sum(result.proposed_compatible for result in assay_results),
    }


def _report(
    *,
    status: Literal["SKIPPED", "COMPLETE"],
    config: InclusivityConfig,
    evaluation_sequence_ids: tuple[str, ...],
    assays: tuple[AssayCandidate, ...],
    oligo_matches: tuple[OligoMatch, ...],
    assay_results: tuple[AssayInclusivity, ...],
    variations: tuple[OligoVariation, ...],
    proposals: tuple[DegeneracyProposal, ...],
    artifacts: dict[str, str | None],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "status": status,
        "enabled": config.enabled,
        "configuration": _json_value(config),
        "evaluation_sequence_ids": _json_value(evaluation_sequence_ids),
        "counts": _counts(
            evaluation_sequence_ids, assays, assay_results, oligo_matches, variations, proposals
        ),
        "oligo_matches": _json_value(oligo_matches),
        "assay_results": _json_value(assay_results),
        "variations": _json_value(variations),
        "proposals": _json_value(proposals),
        "artifacts": artifacts,
    }


def _atomic_write_text(destination: Path, content: str) -> None:
    temporary = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.tmp"
    try:
        with temporary.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def evaluate_inclusivity(
    records: tuple[LocalSequenceRecord, ...],
    evaluation_set: EvaluationSet,
    primer_design: PrimerDesignResult,
    config: InclusivityConfig,
    output_dir: Path,
) -> InclusivityResult:
    """Evaluate local assay inclusivity and atomically publish public artifacts."""
    validate_inclusivity_config(config)
    output_dir = Path(output_dir)
    paths = _artifact_paths(output_dir)

    if not config.enabled:
        report = _report(
            status="SKIPPED", config=config, evaluation_sequence_ids=(), assays=(),
            oligo_matches=(), assay_results=(), variations=(), proposals=(),
            artifacts={
                "oligo_matches": None,
                "assay_inclusivity": None,
                "oligo_variations": None,
                "degeneracy_proposals": None,
            },
        )
        paths["report"].parent.mkdir(parents=True, exist_ok=True)
        paths["report"].unlink(missing_ok=True)
        for name in _DATA_ARTIFACT_NAMES:
            (paths["report"].parent / name).unlink(missing_ok=True)
        _atomic_write_text(
            paths["report"], json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        return InclusivityResult(
            status="SKIPPED", evaluation_sequence_ids=(), oligo_matches=(),
            assay_results=(), variations=(), proposals=(), oligo_matches_path=None,
            assay_inclusivity_path=None, oligo_variations_path=None,
            degeneracy_proposals_path=None, report_path=paths["report"],
        )

    if not isinstance(primer_design, PrimerDesignResult) or primer_design.status != "COMPLETE":
        raise InclusivityError("Inclusivity evaluation requires a COMPLETE PrimerDesignResult.")

    selected = _evaluate_original(records, evaluation_set, primer_design.assays, config)
    variations = _variation_rows(selected, evaluation_set, config)
    proposals = _proposals(primer_design.assays, selected, evaluation_set, config)
    assay_results = _assay_results_with_proposals(selected, proposals, config)
    oligo_matches = tuple(hit.public for binding in selected for hit in binding.retained_hits)
    report = _report(
        status="COMPLETE", config=config,
        evaluation_sequence_ids=evaluation_set.sequence_ids, assays=primer_design.assays,
        oligo_matches=oligo_matches, assay_results=assay_results,
        variations=variations, proposals=proposals,
        artifacts={
            key: _relative_path(paths[key], output_dir)
            for key in ("oligo_matches", "assay_inclusivity", "oligo_variations", "degeneracy_proposals")
        },
    )
    paths["report"].parent.mkdir(parents=True, exist_ok=True)
    paths["report"].unlink(missing_ok=True)
    _atomic_write_text(paths["oligo_matches"], _oligo_match_text(oligo_matches))
    _atomic_write_text(paths["assay_inclusivity"], _assay_text(assay_results))
    _atomic_write_text(paths["oligo_variations"], _variation_text(variations))
    _atomic_write_text(paths["degeneracy_proposals"], _proposal_text(proposals))
    _atomic_write_text(
        paths["report"], json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return InclusivityResult(
        status="COMPLETE", evaluation_sequence_ids=evaluation_set.sequence_ids,
        oligo_matches=oligo_matches, assay_results=assay_results,
        variations=variations, proposals=proposals,
        oligo_matches_path=paths["oligo_matches"],
        assay_inclusivity_path=paths["assay_inclusivity"],
        oligo_variations_path=paths["oligo_variations"],
        degeneracy_proposals_path=paths["degeneracy_proposals"], report_path=paths["report"],
    )
