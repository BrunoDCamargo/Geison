from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from qpcr_pipeline.assay_report_html import render_assay_report_html
from qpcr_pipeline.config import RankingConfig, validate_ranking_config
from qpcr_pipeline.inclusivity import InclusivityResult
from qpcr_pipeline.primer_design import PrimerDesignResult
from qpcr_pipeline.ranking import (
    ClassifiedAssay,
    RankedAssay,
    RankingError,
    RankingReason,
    RankingResult,
    _artifact_paths,
    _atomic_write_text,
    _classification_from_reasons,
    _deduplicate_and_sort_reasons,
    _ranking_report,
    _ranking_sort_key,
    _ranking_tsv_text,
    _score_components,
    classify_assays,
    evaluate_ranking,
)
from qpcr_pipeline.specificity import SpecificityResult


def classify_assays_with_execution_guard(
    primer_design: PrimerDesignResult,
    inclusivity: InclusivityResult,
    specificity: SpecificityResult,
    config: RankingConfig,
    *,
    execution_missing_evidence: tuple[str, ...] = (),
) -> tuple[ClassifiedAssay, ...]:
    classified = classify_assays(primer_design, inclusivity, specificity, config)
    if not execution_missing_evidence:
        return classified

    guarded: list[ClassifiedAssay] = []
    for item in classified:
        missing_components = tuple(sorted(set(item.missing_components) | {"execution"}))
        reason = RankingReason(
            code="RUN_EVIDENCE_INCOMPLETE",
            severity="REVIEW",
            source="execution",
            message="The run is missing scientific evidence required for IN SILICO PASS.",
            evidence=(("missing_evidence", execution_missing_evidence),),
        )
        reasons = _deduplicate_and_sort_reasons([*item.reasons, reason])
        guarded.append(
            replace(
                item,
                classification=_classification_from_reasons(reasons),
                reasons=reasons,
                missing_components=missing_components,
            )
        )
    return tuple(guarded)


def rank_assays_with_execution_guard(
    primer_design: PrimerDesignResult,
    inclusivity: InclusivityResult,
    specificity: SpecificityResult,
    config: RankingConfig,
    *,
    execution_missing_evidence: tuple[str, ...] = (),
) -> tuple[RankedAssay, ...]:
    classified = classify_assays_with_execution_guard(
        primer_design,
        inclusivity,
        specificity,
        config,
        execution_missing_evidence=execution_missing_evidence,
    )
    scored: list[RankedAssay] = []
    for item in classified:
        components = _score_components(item)
        component_values = (
            components.inclusivity,
            components.specificity,
            components.conservation,
            components.primer3_quality,
            components.robustness,
        )
        if execution_missing_evidence or any(value is None for value in component_values):
            score_status = "INCOMPLETE"
            final_score = None
        else:
            inclusivity_score = components.inclusivity
            specificity_score = components.specificity
            conservation_score = components.conservation
            primer3_score = components.primer3_quality
            robustness_score = components.robustness
            assert inclusivity_score is not None
            assert specificity_score is not None
            assert conservation_score is not None
            assert primer3_score is not None
            assert robustness_score is not None
            score_status = "COMPLETE"
            final_score = 100.0 * (
                config.weights.inclusivity * inclusivity_score
                + config.weights.specificity * specificity_score
                + config.weights.conservation * conservation_score
                + config.weights.primer3_quality * primer3_score
                + config.weights.robustness * robustness_score
            )
        scored.append(
            RankedAssay(
                rank=0,
                assay_id=item.assay.assay_id,
                region_id=item.assay.region_id,
                primer3_index=item.assay.primer3_index,
                classification=item.classification,
                final_score=final_score,
                score_status=score_status,
                components=components,
                reasons=item.reasons,
                original_compatible_count=item.original_compatible_count,
                evaluation_sequence_count=item.evaluation_sequence_count,
                compatible_off_target_hit_count=item.compatible_off_target_hit_count,
                plausible_off_target_count=item.plausible_off_target_count,
                detectable_off_target_count=item.detectable_off_target_count,
                pair_penalty=item.assay.pair_penalty,
            )
        )

    ordered = sorted(scored, key=_ranking_sort_key)
    return tuple(replace(item, rank=index) for index, item in enumerate(ordered, 1))


def evaluate_ranking_with_execution_guard(
    primer_design: PrimerDesignResult,
    inclusivity: InclusivityResult,
    specificity: SpecificityResult,
    config: RankingConfig,
    output_dir: Path,
    *,
    target_name: str,
    execution_missing_evidence: tuple[str, ...] = (),
) -> RankingResult:
    """Publish ranking artifacts while enforcing run-level evidence completeness."""
    validate_ranking_config(config)
    output_dir = Path(output_dir)

    if not config.enabled:
        return evaluate_ranking(
            primer_design,
            inclusivity,
            specificity,
            config,
            output_dir,
            target_name=target_name,
        )

    paths = _artifact_paths(output_dir)
    paths["report"].parent.mkdir(parents=True, exist_ok=True)
    for key in ("tsv", "report", "html"):
        paths[key].unlink(missing_ok=True)

    assays = rank_assays_with_execution_guard(
        primer_design,
        inclusivity,
        specificity,
        config,
        execution_missing_evidence=execution_missing_evidence,
    )
    if execution_missing_evidence and any(
        assay.classification == "IN SILICO PASS" for assay in assays
    ):
        raise RankingError("Incomplete run evidence cannot produce IN SILICO PASS.")

    tsv_text = _ranking_tsv_text(assays)
    html_text = render_assay_report_html(
        target_name=target_name,
        primer_design=primer_design,
        inclusivity=inclusivity,
        specificity=specificity,
        assays=assays,
    )
    report = _ranking_report(
        status="COMPLETE",
        config=config,
        assays=assays,
        artifacts={
            "ranking_tsv": "ranking/assay_ranking.tsv",
            "ranking_report": "ranking/ranking_report.json",
            "html_report": "report.html",
        },
    )

    _atomic_write_text(paths["tsv"], tsv_text)
    _atomic_write_text(paths["html"], html_text)
    _atomic_write_text(
        paths["report"], json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    return RankingResult(
        status="COMPLETE",
        assays=assays,
        ranking_tsv_path=paths["tsv"],
        ranking_report_path=paths["report"],
        html_report_path=paths["html"],
    )
