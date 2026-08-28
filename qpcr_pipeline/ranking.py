"""Classify final qPCR assay evidence before any quantitative ranking."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from qpcr_pipeline.config import RankingConfig, validate_ranking_config
from qpcr_pipeline.inclusivity import DegeneracyProposal, InclusivityResult
from qpcr_pipeline.primer_design import AssayCandidate, CandidateRegion, PrimerDesignResult
from qpcr_pipeline.specificity import SpecificityResult


ReasonSeverity = Literal["HIGH_RISK", "REVIEW", "ADVISORY"]
AssayClassification = Literal["IN SILICO PASS", "REVIEW", "HIGH_RISK"]

_ROLE_ORDER = {"FORWARD": 0, "PROBE": 1, "REVERSE": 2}
_SEVERITY_ORDER = {"HIGH_RISK": 0, "REVIEW": 1, "ADVISORY": 2}


class RankingError(RuntimeError):
    """Raised when ranking evidence is structurally untrustworthy."""


@dataclass(frozen=True, slots=True)
class RankingReason:
    code: str
    severity: ReasonSeverity
    source: str
    message: str
    evidence: tuple[tuple[str, object], ...]


@dataclass(frozen=True, slots=True)
class ClassifiedAssay:
    assay: AssayCandidate
    region: CandidateRegion
    classification: AssayClassification
    reasons: tuple[RankingReason, ...]
    original_compatible_count: int | None
    evaluation_sequence_count: int | None
    compatible_off_target_hit_count: int | None
    plausible_off_target_count: int | None
    detectable_off_target_count: int | None
    proposals: tuple[DegeneracyProposal, ...]
    inclusivity_available: bool
    specificity_available: bool
    missing_components: tuple[str, ...]


def classify_assays(
    primer_design: PrimerDesignResult,
    inclusivity: InclusivityResult,
    specificity: SpecificityResult,
    config: RankingConfig,
) -> tuple[ClassifiedAssay, ...]:
    """Validate evidence and classify every original Primer3 assay before scoring."""
    validate_ranking_config(config)
    assays_by_id, regions_by_id = _validate_primer_design(primer_design)
    inclusivity_rows, proposals_by_assay = _validate_inclusivity(
        inclusivity, assays_by_id
    )
    specificity_by_assay = _validate_specificity(specificity, assays_by_id)

    classified: list[ClassifiedAssay] = []
    for assay in primer_design.assays:
        region = regions_by_id[assay.region_id]
        assay_rows = inclusivity_rows.get(assay.assay_id, ())
        assay_proposals = proposals_by_assay.get(assay.assay_id, ())
        specificity_summary = specificity_by_assay.get(assay.assay_id)

        inclusivity_available = (
            inclusivity.status == "COMPLETE"
            and bool(inclusivity.evaluation_sequence_ids)
        )
        specificity_available = specificity.status == "COMPLETE"

        missing_components: set[str] = set()
        reasons: list[RankingReason] = []

        if inclusivity.status == "SKIPPED":
            missing_components.update(("inclusivity", "robustness"))
            reasons.append(
                _reason(
                    "INCLUSIVITY_EVIDENCE_MISSING",
                    "REVIEW",
                    "inclusivity",
                    "Inclusivity evidence was not evaluated.",
                )
            )
        elif not inclusivity.evaluation_sequence_ids:
            missing_components.update(("inclusivity", "robustness"))
            reasons.append(
                _reason(
                    "INCLUSIVITY_EVIDENCE_MISSING",
                    "REVIEW",
                    "inclusivity",
                    "Inclusivity evidence has an empty Evaluation Set.",
                    (("evaluation_sequence_count", 0),),
                )
            )

        if specificity.status == "SKIPPED":
            missing_components.add("specificity")
            reasons.append(
                _reason(
                    "SPECIFICITY_EVIDENCE_MISSING",
                    "REVIEW",
                    "specificity",
                    "Specificity evidence was not evaluated.",
                )
            )

        if assay.pair_penalty is None:
            missing_components.add("primer3_quality")

        original_compatible_count: int | None = None
        evaluation_sequence_count: int | None = None
        if inclusivity.status == "COMPLETE":
            evaluation_sequence_count = len(inclusivity.evaluation_sequence_ids)
            if evaluation_sequence_count:
                original_compatible_count = sum(
                    row.original_compatible for row in assay_rows
                )
                fraction = original_compatible_count / evaluation_sequence_count
                if fraction < config.min_inclusivity_before_high_risk:
                    reasons.append(
                        _reason(
                            "INCLUSIVITY_BELOW_MINIMUM",
                            "HIGH_RISK",
                            "inclusivity",
                            "Original assay inclusivity is below the configured high-risk floor.",
                            (
                                ("compatible_count", original_compatible_count),
                                ("evaluation_sequence_count", evaluation_sequence_count),
                                ("fraction", fraction),
                                (
                                    "threshold",
                                    config.min_inclusivity_before_high_risk,
                                ),
                            ),
                        )
                    )
                elif fraction < config.min_inclusivity_for_pass:
                    reasons.append(
                        _reason(
                            "INCLUSIVITY_BELOW_PASS",
                            "REVIEW",
                            "inclusivity",
                            "Original assay inclusivity is below the configured PASS threshold.",
                            (
                                ("compatible_count", original_compatible_count),
                                ("evaluation_sequence_count", evaluation_sequence_count),
                                ("fraction", fraction),
                                ("threshold", config.min_inclusivity_for_pass),
                            ),
                        )
                    )

        compatible_hit_count: int | None = None
        plausible_count: int | None = None
        detectable_count: int | None = None
        if specificity_summary is not None:
            compatible_hit_count, plausible_count, detectable_count = specificity_summary
            if detectable_count:
                reasons.append(
                    _reason(
                        "DETECTABLE_OFF_TARGET",
                        "HIGH_RISK",
                        "specificity",
                        "At least one plausible off-target amplicon contains a compatible probe.",
                        (
                            ("detectable_off_target_count", detectable_count),
                            ("plausible_off_target_count", plausible_count),
                        ),
                    )
                )
            elif plausible_count:
                reasons.append(
                    _reason(
                        "PLAUSIBLE_OFF_TARGET_AMPLICON",
                        "REVIEW",
                        "specificity",
                        "At least one plausible forward/reverse off-target amplicon was found without a compatible probe.",
                        (("plausible_off_target_count", plausible_count),),
                    )
                )
            elif compatible_hit_count:
                reasons.append(
                    _reason(
                        "ISOLATED_OFF_TARGET_HITS",
                        "ADVISORY",
                        "specificity",
                        "Compatible isolated off-target oligo hits were found without a plausible amplicon.",
                        (("compatible_hit_count", compatible_hit_count),),
                    )
                )

        reasons.extend(_proposal_reasons(assay_proposals))

        if missing_components:
            reasons.append(
                _reason(
                    "EVIDENCE_INCOMPLETE",
                    "REVIEW",
                    "ranking",
                    "One or more score components cannot be computed from the available evidence.",
                    (("components", tuple(sorted(missing_components))),),
                )
            )

        final_reasons = _deduplicate_and_sort_reasons(reasons)
        classification = _classification_from_reasons(final_reasons)
        classified.append(
            ClassifiedAssay(
                assay=assay,
                region=region,
                classification=classification,
                reasons=final_reasons,
                original_compatible_count=original_compatible_count,
                evaluation_sequence_count=evaluation_sequence_count,
                compatible_off_target_hit_count=compatible_hit_count,
                plausible_off_target_count=plausible_count,
                detectable_off_target_count=detectable_count,
                proposals=assay_proposals,
                inclusivity_available=inclusivity_available,
                specificity_available=specificity_available,
                missing_components=tuple(sorted(missing_components)),
            )
        )
    return tuple(classified)


def _validate_primer_design(
    primer_design: PrimerDesignResult,
) -> tuple[dict[str, AssayCandidate], dict[str, CandidateRegion]]:
    if primer_design.status != "COMPLETE":
        raise RankingError("Enabled ranking requires a COMPLETE primer design result.")

    assays_by_id: dict[str, AssayCandidate] = {}
    for assay in primer_design.assays:
        if assay.assay_id in assays_by_id:
            raise RankingError(f"Duplicate Primer3 assay_id: {assay.assay_id}")
        assays_by_id[assay.assay_id] = assay

    regions_by_id: dict[str, CandidateRegion] = {}
    for region in primer_design.candidates:
        if region.region_id in regions_by_id:
            raise RankingError(f"Duplicate candidate region_id: {region.region_id}")
        regions_by_id[region.region_id] = region

    for assay in primer_design.assays:
        if assay.region_id not in regions_by_id:
            raise RankingError(
                f"Assay {assay.assay_id} references missing candidate region {assay.region_id}."
            )
    return assays_by_id, regions_by_id


def _validate_inclusivity(
    result: InclusivityResult,
    assays_by_id: dict[str, AssayCandidate],
) -> tuple[
    dict[str, tuple[object, ...]],
    dict[str, tuple[DegeneracyProposal, ...]],
]:
    if result.status not in ("SKIPPED", "COMPLETE"):
        raise RankingError(f"Unknown inclusivity status: {result.status}")
    if result.status == "SKIPPED":
        return {}, {}

    if len(set(result.evaluation_sequence_ids)) != len(result.evaluation_sequence_ids):
        raise RankingError("Inclusivity Evaluation Set IDs must be unique.")

    assay_ids = set(assays_by_id)
    sequence_ids = set(result.evaluation_sequence_ids)
    expected_pairs = {
        (assay_id, sequence_id)
        for assay_id in assay_ids
        for sequence_id in result.evaluation_sequence_ids
    }
    observed_pairs: list[tuple[str, str]] = []
    rows_by_assay: dict[str, list[object]] = {assay_id: [] for assay_id in assay_ids}
    for row in result.assay_results:
        if row.assay_id not in assay_ids:
            raise RankingError(
                f"Inclusivity assay result references unknown assay {row.assay_id}."
            )
        if row.sequence_id not in sequence_ids:
            raise RankingError(
                f"Inclusivity assay result references unknown Evaluation Set sequence {row.sequence_id}."
            )
        observed_pairs.append((row.assay_id, row.sequence_id))
        rows_by_assay[row.assay_id].append(row)
    if len(set(observed_pairs)) != len(observed_pairs):
        raise RankingError("Inclusivity assay matrix contains duplicate assay/sequence rows.")
    if set(observed_pairs) != expected_pairs:
        raise RankingError("Inclusivity assay matrix is missing required assay/sequence rows.")

    for hit in result.oligo_matches:
        if hit.assay_id not in assay_ids:
            raise RankingError(
                f"Inclusivity oligo match references unknown assay {hit.assay_id}."
            )
    for variation in result.variations:
        if variation.assay_id not in assay_ids:
            raise RankingError(
                f"Inclusivity variation references unknown assay {variation.assay_id}."
            )

    proposals_by_assay: dict[str, list[DegeneracyProposal]] = {
        assay_id: [] for assay_id in assay_ids
    }
    proposal_keys: set[tuple[str, str]] = set()
    for proposal in result.proposals:
        if proposal.assay_id not in assay_ids:
            raise RankingError(
                f"Degeneracy proposal references unknown assay {proposal.assay_id}."
            )
        if proposal.role not in _ROLE_ORDER:
            raise RankingError(f"Invalid degeneracy proposal role: {proposal.role}")
        key = (proposal.assay_id, proposal.role)
        if key in proposal_keys:
            raise RankingError(
                f"Duplicate degeneracy proposal for assay/role {proposal.assay_id}/{proposal.role}."
            )
        proposal_keys.add(key)
        proposals_by_assay[proposal.assay_id].append(proposal)

    return (
        {
            assay_id: tuple(rows_by_assay[assay_id])
            for assay_id in assays_by_id
        },
        {
            assay_id: tuple(
                sorted(
                    proposals_by_assay[assay_id],
                    key=lambda item: (_ROLE_ORDER[item.role], item.status, item.proposed_sequence),
                )
            )
            for assay_id in assays_by_id
        },
    )


def _validate_specificity(
    result: SpecificityResult,
    assays_by_id: dict[str, AssayCandidate],
) -> dict[str, tuple[int, int, int]]:
    if result.status not in ("SKIPPED", "COMPLETE"):
        raise RankingError(f"Unknown specificity status: {result.status}")
    if result.status == "SKIPPED":
        return {}

    assay_ids = set(assays_by_id)
    if result.assay_count != len(assay_ids):
        raise RankingError(
            "Specificity assay_count does not match the Primer3 assay count."
        )
    if len(set(result.dataset_names)) != len(result.dataset_names):
        raise RankingError("Specificity dataset names must be unique.")
    dataset_names = set(result.dataset_names)

    expected_keys = {
        (dataset_name, assay_id, role)
        for dataset_name in result.dataset_names
        for assay_id in assay_ids
        for role in _ROLE_ORDER
    }
    observed_keys: list[tuple[str, str, str]] = []
    hit_totals_by_assay = {assay_id: 0 for assay_id in assay_ids}
    for row in result.retention:
        if row.assay_id not in assay_ids:
            raise RankingError(
                f"Specificity retention references unknown assay {row.assay_id}."
            )
        if row.dataset_name not in dataset_names:
            raise RankingError(
                f"Specificity retention references unknown dataset {row.dataset_name}."
            )
        if row.role not in _ROLE_ORDER:
            raise RankingError(f"Invalid specificity retention role: {row.role}")
        for name, value in (
            ("total_hit_count", row.total_hit_count),
            ("retained_hit_count", row.retained_hit_count),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise RankingError(
                    f"Specificity retention {name} must be a non-negative integer."
                )
        if row.retained_hit_count > row.total_hit_count:
            raise RankingError(
                "Specificity retained_hit_count cannot exceed total_hit_count."
            )
        key = (row.dataset_name, row.assay_id, row.role)
        observed_keys.append(key)
        hit_totals_by_assay[row.assay_id] += row.total_hit_count
    if len(set(observed_keys)) != len(observed_keys):
        raise RankingError("Specificity retention matrix contains duplicate rows.")
    if set(observed_keys) != expected_keys:
        raise RankingError("Specificity retention matrix is missing required rows.")

    for hit in result.hits:
        if hit.assay_id not in assay_ids:
            raise RankingError(
                f"Specificity hit references unknown assay {hit.assay_id}."
            )
        if hit.dataset_name not in dataset_names:
            raise RankingError(
                f"Specificity hit references unknown dataset {hit.dataset_name}."
            )
    for amplicon in result.amplicons:
        if amplicon.assay_id not in assay_ids:
            raise RankingError(
                f"Specificity amplicon references unknown assay {amplicon.assay_id}."
            )
        if amplicon.dataset_name not in dataset_names:
            raise RankingError(
                f"Specificity amplicon references unknown dataset {amplicon.dataset_name}."
            )

    plausible_by_assay = {assay_id: 0 for assay_id in assay_ids}
    detectable_by_assay = {assay_id: 0 for assay_id in assay_ids}
    for amplicon in result.amplicons:
        if amplicon.primer_amplicon_plausible:
            plausible_by_assay[amplicon.assay_id] += 1
        if amplicon.detectable_off_target:
            detectable_by_assay[amplicon.assay_id] += 1

    return {
        assay_id: (
            hit_totals_by_assay[assay_id],
            plausible_by_assay[assay_id],
            detectable_by_assay[assay_id],
        )
        for assay_id in assays_by_id
    }


def _proposal_reasons(
    proposals: tuple[DegeneracyProposal, ...],
) -> tuple[RankingReason, ...]:
    reasons: list[RankingReason] = []
    for status, code, message in (
        (
            "ACCEPTED",
            "IUPAC_PROPOSAL_ACCEPTED",
            "An accepted IUPAC degeneracy proposal is available as contextual evidence.",
        ),
        (
            "REJECTED",
            "IUPAC_PROPOSAL_REJECTED",
            "A rejected IUPAC degeneracy proposal is preserved as contextual evidence.",
        ),
    ):
        selected = tuple(item for item in proposals if item.status == status)
        if not selected:
            continue
        roles = tuple(
            item.role
            for item in sorted(selected, key=lambda item: _ROLE_ORDER[item.role])
        )
        reasons.append(
            _reason(
                code,
                "ADVISORY",
                "degeneracy",
                message,
                (("roles", roles),),
            )
        )
    return tuple(reasons)


def _reason(
    code: str,
    severity: ReasonSeverity,
    source: str,
    message: str,
    evidence: tuple[tuple[str, object], ...] = (),
) -> RankingReason:
    return RankingReason(
        code=code,
        severity=severity,
        source=source,
        message=message,
        evidence=evidence,
    )


def _deduplicate_and_sort_reasons(
    reasons: list[RankingReason],
) -> tuple[RankingReason, ...]:
    by_identity: dict[tuple[str, str], RankingReason] = {}
    for reason in reasons:
        by_identity.setdefault((reason.code, reason.source), reason)
    return tuple(
        sorted(
            by_identity.values(),
            key=lambda item: (
                _SEVERITY_ORDER[item.severity],
                item.source,
                item.code,
            ),
        )
    )


def _classification_from_reasons(
    reasons: tuple[RankingReason, ...],
) -> AssayClassification:
    if any(reason.severity == "HIGH_RISK" for reason in reasons):
        return "HIGH_RISK"
    if any(reason.severity == "REVIEW" for reason in reasons):
        return "REVIEW"
    return "IN SILICO PASS"
