"""Self-contained HTML boundary for the conservation report."""

from __future__ import annotations

import json

from qpcr_pipeline.config import ConservationConfig


def render_conservation_html(
    *,
    target_name: str,
    reference_id: str | None,
    sequence_count: int,
    config: ConservationConfig,
    windows: tuple[object, ...],
    annotations: tuple[object, ...],
) -> str:
    """Render a deterministic basic shell; Task 3 supplies the interactive plot."""
    payload = {
        "identity": {
            "targetName": target_name,
            "referenceId": reference_id,
            "sequenceCount": sequence_count,
            "windowSize": config.window_size,
            "stepSize": config.step_size,
        },
        "windows": [
            [
                item.reference_start,
                item.reference_end,
                item.position_count,
                item.mean_conservation,
                item.minimum_conservation,
                item.mean_coverage,
                item.mean_gap_frequency,
                item.mean_entropy_bits,
            ]
            for item in windows
        ],
        "annotations": [
            [item.feature_type, item.start, item.end, item.strand, item.label]
            for item in annotations
        ],
    }
    payload_json = json.dumps(
        payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    empty_state = (
        '<p id="empty-state">No conservation windows are available.</p>'
        if not windows
        else '<p id="empty-state" hidden>No conservation windows are available.</p>'
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Geison conservation report</title>
<style>body{{font-family:system-ui,sans-serif;max-width:72rem;margin:2rem auto;padding:0 1rem;color:#172033}}dt{{font-weight:700}}dd{{margin:0 0 .5rem}}code{{font-family:ui-monospace,monospace}}</style>
</head>
<body>
<main>
<h1>Conservation report</h1>
<dl><dt>Target</dt><dd id="target"></dd><dt>Reference</dt><dd id="reference"></dd><dt>Discovery sequences</dt><dd id="sequence-count"></dd><dt>Windows</dt><dd id="window-parameters"></dd></dl>
{empty_state}
<p>Per-position metrics and reference-coordinate consensus sequences are available in the conservation artifact directory.</p>
</main>
<script>
"use strict";
const reportData={payload_json};
document.getElementById("target").textContent=reportData.identity.targetName;
document.getElementById("reference").textContent=reportData.identity.referenceId ?? "None";
document.getElementById("sequence-count").textContent=String(reportData.identity.sequenceCount);
document.getElementById("window-parameters").textContent=`${{reportData.identity.windowSize}} bases, step ${{reportData.identity.stepSize}}`;
</script>
</body>
</html>
"""
