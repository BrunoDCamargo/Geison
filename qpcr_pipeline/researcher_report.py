"""Artifact-driven model for the final researcher-facing run report."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import os
from pathlib import Path
import tempfile


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
    conservation_windows: tuple[dict[str, object], ...] = ()
    specificity_hits: tuple[dict[str, object], ...] = ()
    specificity_amplicons: tuple[dict[str, object], ...] = ()


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


def _coerce_tsv_scalar(value: str | None) -> object:
    if value is None or value == "":
        return None
    lowered = value.casefold()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _load_tsv_rows(
    output_dir: Path,
    relative_path: str,
) -> tuple[dict[str, object], ...]:
    path = output_dir / relative_path
    if not path.is_file():
        return ()
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ResearcherReportError(
                    f"Researcher report TSV artifact has no header: {relative_path}"
                )
            return tuple(
                {
                    str(key): _coerce_tsv_scalar(value)
                    for key, value in row.items()
                    if key is not None
                }
                for row in reader
            )
    except ResearcherReportError:
        raise
    except (OSError, UnicodeError, csv.Error) as error:
        raise ResearcherReportError(
            f"Could not read researcher report artifact {relative_path}: "
            f"{type(error).__name__}: {error}"
        ) from error


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
        conservation_windows=_load_tsv_rows(
            root,
            "conservation/window_metrics.tsv",
        ),
        specificity_hits=_load_tsv_rows(
            root,
            "specificity/off_target_hits.tsv",
        ),
        specificity_amplicons=_load_tsv_rows(
            root,
            "specificity/plausible_amplicons.tsv",
        ),
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


def generate_researcher_report(output_dir: Path) -> Path:
    """Render and atomically publish the final researcher-facing HTML report."""
    root = Path(output_dir)
    data = load_researcher_report_data(root)
    from qpcr_pipeline.researcher_report_html import render_researcher_report_html

    html = render_researcher_report_html(data)
    destination = root / "report.html"
    descriptor, temporary_name = tempfile.mkstemp(
        dir=root,
        prefix=".report.html.",
        suffix=".tmp",
        text=True,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(
            descriptor,
            "w",
            encoding="utf-8",
            newline="\n",
        ) as handle:
            handle.write(html)
            if not html.endswith("\n"):
                handle.write("\n")
            handle.flush()
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination
