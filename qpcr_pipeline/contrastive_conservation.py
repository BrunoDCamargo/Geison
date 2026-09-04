"""Contrastive conservation public API with presentation artifact publishing."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from qpcr_pipeline.contrastive_conservation_core import *  # noqa: F401,F403
from qpcr_pipeline.contrastive_conservation_core import (
    _atomic_write_text,
    _report_payload,
    analyze_contrastive_conservation as _analyze_contrastive_conservation,
)


def analyze_contrastive_conservation(
    conservation,
    approved_panel,
    off_target_configs,
    config,
    primer_config,
    output_dir,
    *,
    similarity_engine=None,
):
    """Run the core analysis and publish its read-only HTML artifact."""
    result = _analyze_contrastive_conservation(
        conservation,
        approved_panel,
        off_target_configs,
        config,
        primer_config,
        output_dir,
        similarity_engine=similarity_engine,
    )
    html_path = Path(output_dir) / "contrastive_conservation" / "report.html"
    if result.status != "COMPLETE":
        html_path.unlink(missing_ok=True)
        return result

    from qpcr_pipeline.contrastive_report_html import render_contrastive_html

    updated = replace(result, html_report_path=html_path)
    _atomic_write_text(
        html_path,
        render_contrastive_html(
            target_name=approved_panel.definition.target.name,
            reference_id=updated.reference_id,
            windows=updated.windows,
            dataset_evidence=updated.dataset_evidence,
            candidates=updated.candidates,
        ),
    )
    _atomic_write_text(
        updated.report_path,
        json.dumps(
            _report_payload(updated, config, approved_panel),
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    return updated
