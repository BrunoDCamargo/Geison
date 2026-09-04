"""Artifact-driven model for the final researcher-facing run report."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path


class ResearcherReportError(RuntimeError):
    """Raised when published run evidence cannot be loaded safely."""


@dataclass(frozen=True, slots=True)
class ResearcherReportData:
    output_dir: Path
    run_manifest: dict[str, object]
    run_summary: dict[str, object] | None
    qc_report: dict[str, object] | None
    panel: dict[str, object] | None
    conservation: dict[str, object] | None
    contrastive: dict[str, object] | None
    primer_design: dict[str, object] | None
    inclusivity: dict[str, object] | None
    specificity: dict[str, object] | None
    ranking: dict[str, object] | None


def _load_json_object(
    output_dir: Path,
    relative_path: str,
    *,
    required: bool = False,
) -> dict[str, object] | None:
    path = output_dir / relative_path
    if not path.is_file():
        if required:
            raise ResearcherReportError(
                f"Required researcher report artifact is missing: {relative_path}"
            )
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ResearcherReportError(
            f"Could not read researcher report artifact {relative_path}: "
            f"{type(error).__name__}: {error}"
        ) from error
    if not isinstance(payload, dict):
        raise ResearcherReportError(
            f"Researcher report artifact {relative_path} must contain a JSON object."
        )
    return payload


def load_researcher_report_data(output_dir: Path) -> ResearcherReportData:
    """Load only published Geison artifacts used by the final report."""
    root = Path(output_dir)
    run_manifest = _load_json_object(root, "run_manifest.json", required=True)
    assert run_manifest is not None
    return ResearcherReportData(
        output_dir=root,
        run_manifest=run_manifest,
        run_summary=_load_json_object(root, "run_summary.json"),
        qc_report=_load_json_object(root, "qc_report.json"),
        panel=_load_json_object(root, "panel/approved_panel.json"),
        conservation=_load_json_object(root, "conservation/conservation_report.json"),
        contrastive=_load_json_object(
            root,
            "contrastive_conservation/contrastive_conservation_report.json",
        ),
        primer_design=_load_json_object(
            root,
            "primer_design/primer_design_report.json",
        ),
        inclusivity=_load_json_object(root, "inclusivity/inclusivity_report.json"),
        specificity=_load_json_object(root, "specificity/specificity_report.json"),
        ranking=_load_json_object(root, "ranking/ranking_report.json"),
    )


def _ranking_assays(
    ranking_report: dict[str, object] | None,
) -> tuple[dict[str, object], ...]:
    if ranking_report is None:
        return ()
    assays = ranking_report.get("assays")
    if not isinstance(assays, list):
        return ()
    return tuple(item for item in assays if isinstance(item, dict))


def scientific_outcome(
    run_status: str,
    ranking_report: dict[str, object] | None,
) -> tuple[str, str]:
    """Map recorded execution/classification evidence to report wording."""
    if run_status == "FAILED":
        return (
            "Execution failed - no conclusive scientific outcome",
            "The computational run failed before a conclusive in-silico assay outcome "
            "could be established. Review the recorded failure evidence before reuse.",
        )
    if run_status == "ACTION_REQUIRED":
        return (
            "Inconclusive - scientific review required before execution",
            "The run is waiting for required scientific approval and has not completed "
            "the assay-discovery workflow.",
        )
    if run_status != "COMPLETED":
        return (
            "Inconclusive - insufficient evidence",
            "The run does not contain all scientific evidence required for a conclusive "
            "in-silico assay classification.",
        )

    assays = _ranking_assays(ranking_report)
    if not assays:
        return (
            "Inconclusive - insufficient evidence",
            "The computational run recorded no ranked assay candidates, so no conclusive "
            "assay acceptance statement can be made.",
        )

    classifications = [item.get("classification") for item in assays]
    if "IN SILICO PASS" in classifications:
        return (
            "In-silico candidate(s) identified",
            "At least one designed assay is classified IN SILICO PASS based on the "
            "recorded computational evidence. Experimental validation is still required.",
        )
    if "REVIEW" in classifications:
        return (
            "No in-silico pass; candidate(s) require review",
            "No assay reached IN SILICO PASS. At least one candidate is classified REVIEW "
            "and requires scientific assessment before experimental follow-up.",
        )

    if all(item == "HIGH_RISK" for item in classifications):
        has_detectable_off_target = any(
            reason.get("code") == "DETECTABLE_OFF_TARGET"
            for assay in assays
            for reason in (
                assay.get("reasons")
                if isinstance(assay.get("reasons"), list)
                else []
            )
            if isinstance(reason, dict)
        )
        if has_detectable_off_target:
            detail = (
                "The computational workflow completed successfully, but all designed "
                "assays were classified HIGH_RISK because detectable off-target "
                "amplification was identified."
            )
        else:
            detail = (
                "The computational workflow completed successfully, but all designed "
                "assays were classified HIGH_RISK by the recorded scientific evidence."
            )
        return "No in-silico acceptable assay candidates identified", detail

    return (
        "Inconclusive - insufficient evidence",
        "The recorded ranking classifications are insufficient for a conclusive in-silico "
        "assay outcome.",
    )
