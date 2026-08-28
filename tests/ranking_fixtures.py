from pathlib import Path

from qpcr_pipeline.inclusivity import (
    AssayInclusivity,
    DegeneracyProposal,
    InclusivityResult,
)
from qpcr_pipeline.primer_design import (
    AssayCandidate,
    CandidateRegion,
    DesignedOligo,
    PrimerDesignResult,
)
from qpcr_pipeline.specificity import (
    HitRetentionSummary,
    OffTargetHit,
    PlausibleAmplicon,
    SpecificityResult,
)


ROLES = ("FORWARD", "PROBE", "REVERSE")


def make_region(region_id: str = "r1") -> CandidateRegion:
    return CandidateRegion(
        region_id=region_id,
        rank=1,
        reference_start=1,
        reference_end=100,
        peak_start=1,
        peak_end=100,
        position_count=100,
        usable_length=100,
        usable_fraction=1.0,
        mean_conservation=1.0,
        minimum_conservation=1.0,
        mean_coverage=1.0,
        mean_gap_frequency=0.0,
        mean_entropy_bits=0.0,
    )


def make_oligo(sequence: str, start: int) -> DesignedOligo:
    return DesignedOligo(
        sequence=sequence,
        reference_start=start,
        reference_end=start + len(sequence) - 1,
        length=len(sequence),
        tm=60.0,
        gc_percent=50.0,
        penalty=0.0,
        metrics=(),
    )


def make_assay(
    assay_id: str = "a1",
    region_id: str = "r1",
    primer3_index: int = 0,
    pair_penalty: float | None = 1.0,
) -> AssayCandidate:
    return AssayCandidate(
        assay_id=assay_id,
        region_id=region_id,
        primer3_index=primer3_index,
        forward_primer=make_oligo("ACGT", 1),
        probe=make_oligo("TTAA", 9),
        reverse_primer=make_oligo("AGTC", 17),
        product_size=20,
        pair_penalty=pair_penalty,
        metrics=(),
    )


def make_primer_result(
    assays: tuple[AssayCandidate, ...] | None = None,
    candidates: tuple[CandidateRegion, ...] | None = None,
    status: str = "COMPLETE",
) -> PrimerDesignResult:
    actual_assays = assays if assays is not None else (make_assay(),)
    actual_candidates = candidates if candidates is not None else (make_region(),)
    return PrimerDesignResult(
        status=status,
        reference_id="ref" if status == "COMPLETE" else None,
        candidates=actual_candidates,
        assays=actual_assays,
        candidate_regions_path=None,
        assays_path=None,
        primer3_input_path=None,
        primer3_output_path=None,
        report_path=Path("primer_design_report.json"),
    )


def make_inclusivity_result(
    primer: PrimerDesignResult,
    compatibility: dict[tuple[str, str], bool] | None = None,
    sequence_ids: tuple[str, ...] = ("s1",),
    proposals: tuple[DegeneracyProposal, ...] = (),
    status: str = "COMPLETE",
) -> InclusivityResult:
    if status == "SKIPPED":
        return InclusivityResult(
            status="SKIPPED",
            evaluation_sequence_ids=(),
            oligo_matches=(),
            assay_results=(),
            variations=(),
            proposals=(),
            oligo_matches_path=None,
            assay_inclusivity_path=None,
            oligo_variations_path=None,
            degeneracy_proposals_path=None,
            report_path=Path("inclusivity_report.json"),
        )
    compatibility = compatibility or {}
    rows = []
    for assay in primer.assays:
        for sequence_id in sequence_ids:
            compatible = compatibility.get((assay.assay_id, sequence_id), True)
            rows.append(
                AssayInclusivity(
                    assay_id=assay.assay_id,
                    sequence_id=sequence_id,
                    orientation="FORWARD" if compatible else None,
                    geometry_found=compatible,
                    source_amplicon_start=1 if compatible else None,
                    source_amplicon_end=assay.product_size if compatible else None,
                    amplicon_size=assay.product_size if compatible else None,
                    forward_match=None,
                    probe_match=None,
                    reverse_match=None,
                    original_compatible=compatible,
                    proposed_forward=None,
                    proposed_probe=None,
                    proposed_reverse=None,
                    proposed_compatible=compatible,
                )
            )
    return InclusivityResult(
        status="COMPLETE",
        evaluation_sequence_ids=sequence_ids,
        oligo_matches=(),
        assay_results=tuple(rows),
        variations=(),
        proposals=proposals,
        oligo_matches_path=None,
        assay_inclusivity_path=None,
        oligo_variations_path=None,
        degeneracy_proposals_path=None,
        report_path=Path("inclusivity_report.json"),
    )


def make_specificity_result(
    primer: PrimerDesignResult,
    dataset_names: tuple[str, ...] = ("off",),
    hit_totals: dict[tuple[str, str, str], int] | None = None,
    amplicons: tuple[PlausibleAmplicon, ...] = (),
    hits: tuple[OffTargetHit, ...] = (),
    status: str = "COMPLETE",
) -> SpecificityResult:
    if status == "SKIPPED":
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
            report_path=Path("specificity_report.json"),
        )
    hit_totals = hit_totals or {}
    retention = []
    for dataset_name in dataset_names:
        for assay in primer.assays:
            for role in ROLES:
                total = hit_totals.get((dataset_name, assay.assay_id, role), 0)
                retained = min(total, 20)
                retention.append(
                    HitRetentionSummary(
                        dataset_name=dataset_name,
                        assay_id=assay.assay_id,
                        role=role,
                        total_hit_count=total,
                        retained_hit_count=retained,
                        truncated=retained < total,
                    )
                )
    return SpecificityResult(
        status="COMPLETE",
        dataset_names=dataset_names,
        sequence_count=len(dataset_names),
        assay_count=len(primer.assays),
        hits=hits,
        amplicons=amplicons,
        retention=tuple(retention),
        off_target_hits_path=None,
        plausible_amplicons_path=None,
        report_path=Path("specificity_report.json"),
    )


def make_amplicon(
    assay_id: str = "a1",
    dataset_name: str = "off",
    detectable: bool = False,
) -> PlausibleAmplicon:
    return PlausibleAmplicon(
        dataset_name=dataset_name,
        assay_id=assay_id,
        sequence_id="off-seq",
        orientation="FORWARD",
        source_start=1,
        source_end=20,
        amplicon_size=20,
        forward_source_start=1,
        forward_source_end=4,
        reverse_source_start=17,
        reverse_source_end=20,
        probe_source_sites=((9, 12),) if detectable else (),
        forward_hit_rank=1,
        reverse_hit_rank=1,
        probe_hit_ranks=(1,) if detectable else (),
        primer_amplicon_plausible=True,
        detectable_off_target=detectable,
    )


def make_proposal(
    assay_id: str = "a1",
    role: str = "FORWARD",
    status: str = "ACCEPTED",
    original_degeneracy: int = 1,
    proposed_degeneracy: int = 2,
) -> DegeneracyProposal:
    return DegeneracyProposal(
        assay_id=assay_id,
        role=role,
        original_sequence="ACGT",
        proposed_sequence="ARGT",
        status=status,
        reason="fixture",
        original_degeneracy=original_degeneracy,
        proposed_degeneracy=proposed_degeneracy,
        changed_positions=(2,),
        binding_site_count=10,
        original_exact_count=9,
        original_exact_fraction=0.9,
        proposed_exact_count=10,
        proposed_exact_fraction=1.0,
    )


def make_hit(assay_id: str = "a1", dataset_name: str = "off") -> OffTargetHit:
    return OffTargetHit(
        dataset_name=dataset_name,
        assay_id=assay_id,
        sequence_id="off-seq",
        role="FORWARD",
        orientation="FORWARD",
        hit_rank=1,
        source_start=1,
        source_end=4,
        mismatch_positions=(),
        mismatch_count=0,
        exact_match=True,
        three_prime_mismatch=False,
        compatible=True,
    )
