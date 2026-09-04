"""Self-contained HTML renderer for the final Geison researcher report."""

from __future__ import annotations

from html import escape
import json
from typing import Iterable

from qpcr_pipeline.researcher_report import ResearcherReportData, scientific_outcome


def _text(value: object) -> str:
    if value is None or value == "":
        return "Not available"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, float):
        return f"{value:.4g}"
    return escape(str(value), quote=True)


def _range(start: object, end: object) -> str:
    if start is None or end is None:
        return "Not available"
    return f"{_text(start)}–{_text(end)}"


def _dict(value: object) -> dict[str, object]:
    return value if isinstance(value, dict) else {}


def _list(value: object) -> list[object]:
    return value if isinstance(value, list) else []


def _evidence_unavailable(path: str) -> str:
    return (
        '<div class="unavailable"><strong>Evidence unavailable</strong>'
        f"<p>The published artifact <code>{escape(path)}</code> is not available for "
        "this run.</p></div>"
    )


def _summary_card(label: str, value: object) -> str:
    return (
        '<div class="summary-card">'
        f'<span class="summary-label">{escape(label)}</span>'
        f'<strong>{_text(value)}</strong>'
        "</div>"
    )


def _metric_table(mapping: dict[str, object], fields: Iterable[tuple[str, str]]) -> str:
    rows = "".join(
        f"<tr><th>{escape(label)}</th><td>{_text(mapping.get(key))}</td></tr>"
        for key, label in fields
    )
    return f'<table class="compact"><tbody>{rows}</tbody></table>'


def _polyline_chart(
    rows: list[dict[str, object]],
    series: tuple[tuple[str, str, str], ...],
    *,
    title: str,
) -> str:
    numeric_rows: list[tuple[float, dict[str, object]]] = []
    for row in rows:
        start = row.get("reference_start")
        end = row.get("reference_end")
        if isinstance(start, (int, float)) and isinstance(end, (int, float)):
            numeric_rows.append(((float(start) + float(end)) / 2.0, row))
    if not numeric_rows:
        return '<p class="muted">No recorded window series available for this view.</p>'
    min_x = min(item[0] for item in numeric_rows)
    max_x = max(item[0] for item in numeric_rows)
    span_x = max(max_x - min_x, 1.0)
    left, top, width, height = 54.0, 18.0, 630.0, 150.0

    polylines: list[str] = []
    legend: list[str] = []
    dash_patterns = ("", "6 5", "2 4")
    for index, (key, label, stroke_class) in enumerate(series):
        points: list[str] = []
        for center, row in numeric_rows:
            value = row.get(key)
            if not isinstance(value, (int, float)):
                continue
            x = left + (center - min_x) / span_x * width
            bounded = max(0.0, min(1.0, float(value)))
            y = top + (1.0 - bounded) * height
            points.append(f"{x:.1f},{y:.1f}")
        if points:
            dash = dash_patterns[index % len(dash_patterns)]
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            polylines.append(
                f'<polyline class="{stroke_class}" points="{" ".join(points)}"'
                f'{dash_attr} fill="none" stroke-width="3" />'
            )
            legend.append(
                f'<span><i class="legend-line {stroke_class}"></i>{escape(label)}</span>'
            )

    return (
        '<div class="chart-wrap">'
        f'<h4>{escape(title)}</h4>'
        '<svg viewBox="0 0 720 205" role="img" aria-label="Recorded window metrics">'
        '<line x1="54" y1="18" x2="54" y2="168" class="axis" />'
        '<line x1="54" y1="168" x2="684" y2="168" class="axis" />'
        '<text x="12" y="24" class="axis-label">1.0</text>'
        '<text x="18" y="171" class="axis-label">0.0</text>'
        f'<text x="54" y="192" class="axis-label">{_text(int(min_x))}</text>'
        f'<text x="650" y="192" class="axis-label">{_text(int(max_x))}</text>'
        + "".join(polylines)
        + "</svg>"
        + '<div class="legend">'
        + "".join(legend)
        + "</div></div>"
    )


def _target_name(data: ResearcherReportData) -> str:
    panel = _dict(data.panel)
    definition = _dict(panel.get("definition"))
    target = _dict(definition.get("target"))
    value = target.get("name") or data.run_manifest.get("target_name")
    return str(value) if value is not None else "Not available"


def _run_summary_section(data: ResearcherReportData) -> str:
    ranking_counts = _dict(_dict(data.ranking).get("counts"))
    contrast_counts = _dict(_dict(data.contrastive).get("counts"))
    primer_counts = _dict(_dict(data.primer_design).get("counts"))
    run_summary = _dict(data.run_summary)
    panel_definition = _dict(_dict(data.panel).get("definition"))
    non_targets = [
        item
        for item in _list(panel_definition.get("non_targets"))
        if isinstance(item, dict) and "CHALLENGE" in _list(item.get("dataset_roles"))
    ]
    challenge_count = contrast_counts.get("challenge_datasets", len(non_targets))
    assay_count = ranking_counts.get("assays", primer_counts.get("assays"))
    cards = "".join(
        (
            _summary_card("Technical status", data.run_manifest.get("status")),
            _summary_card("Target sequences", run_summary.get("sequence_count")),
            _summary_card("Challenge datasets", challenge_count),
            _summary_card("Contrastive regions", contrast_counts.get("candidate_regions")),
            _summary_card("Designed assays", assay_count),
            _summary_card("IN SILICO PASS", ranking_counts.get("in_silico_pass", 0)),
            _summary_card("REVIEW", ranking_counts.get("review", 0)),
            _summary_card("HIGH_RISK", ranking_counts.get("high_risk", 0)),
        )
    )
    completeness = _dict(data.run_manifest.get("scientific_completeness"))
    missing = _list(completeness.get("missing_evidence"))
    missing_html = ""
    if missing:
        missing_html = (
            '<p class="warning-line"><strong>Missing evidence:</strong> '
            + ", ".join(escape(str(item)) for item in missing)
            + "</p>"
        )
    return f'<section><h2>Run summary</h2><div class="summary-grid">{cards}</div>{missing_html}</section>'


def _outcome_section(data: ResearcherReportData) -> str:
    status = str(data.run_manifest.get("status", "UNKNOWN"))
    title, body = scientific_outcome(status, data.ranking)
    if title == "In-silico candidate(s) identified":
        css = "outcome-pass"
    elif "acceptable" in title or "failed" in title.lower():
        css = "outcome-risk"
    elif "review" in title.lower():
        css = "outcome-review"
    else:
        css = "outcome-inconclusive"
    return (
        f'<section class="outcome {css}"><p class="eyebrow">Scientific outcome</p>'
        f'<h2>{escape(title)}</h2><p>{escape(body)}</p></section>'
    )


def _panel_section(data: ResearcherReportData) -> str:
    if data.panel is None:
        return (
            '<section><h2>Approved panel and study context</h2>'
            + _evidence_unavailable("panel/approved_panel.json")
            + "</section>"
        )
    panel = data.panel
    definition = _dict(panel.get("definition"))
    target = _dict(definition.get("target"))
    groups = [item for item in _list(target.get("groups")) if isinstance(item, dict)]
    non_targets = [
        item for item in _list(definition.get("non_targets")) if isinstance(item, dict)
    ]
    context = _dict(definition.get("diagnostic_context"))

    group_rows = "".join(
        "<tr>"
        f"<td>{_text(item.get('name'))}</td>"
        f"<td>{_text(', '.join(str(role) for role in _list(item.get('dataset_roles'))))}</td>"
        f"<td>{_text(item.get('required'))}</td>"
        f"<td>{_text('; '.join(str(reason) for reason in _list(item.get('reasons'))))}</td>"
        "</tr>"
        for item in groups
    ) or '<tr><td colspan="4">No DESIGN groups recorded.</td></tr>'
    challenge_rows = "".join(
        "<tr>"
        f"<td>{_text(item.get('name'))}</td>"
        f"<td>{_text(item.get('criticality'))}</td>"
        f"<td>{_text(', '.join(str(role) for role in _list(item.get('dataset_roles'))))}</td>"
        f"<td>{_text('; '.join(str(reason) for reason in _list(item.get('reasons'))))}</td>"
        "</tr>"
        for item in non_targets
    ) or '<tr><td colspan="4">No challenge datasets recorded.</td></tr>'

    context_table = _metric_table(
        context,
        (
            ("syndrome", "Syndrome"),
            ("geography", "Geography"),
            ("sample_type", "Sample type"),
            ("vector", "Vector / context"),
        ),
    )
    return (
        '<section><h2>Approved panel and study context</h2>'
        '<div class="two-col">'
        '<div><h3>Panel provenance</h3><dl>'
        f'<dt>Status</dt><dd>{_text(panel.get("status"))}</dd>'
        f'<dt>Target</dt><dd>{_text(target.get("name"))}</dd>'
        f'<dt>Target mode</dt><dd>{_text(target.get("mode"))}</dd>'
        f'<dt>Proposal hash</dt><dd class="mono">{_text(panel.get("proposal_sha256"))}</dd>'
        '</dl></div><div><h3>Diagnostic context</h3>'
        + context_table
        + "</div></div>"
        '<h3>Target groups</h3><div class="table-scroll"><table><thead><tr>'
        '<th>Group</th><th>Roles</th><th>Required</th><th>Rationale</th>'
        f'</tr></thead><tbody>{group_rows}</tbody></table></div>'
        '<h3>Challenge panel</h3><div class="table-scroll"><table><thead><tr>'
        '<th>Dataset</th><th>Criticality</th><th>Roles</th><th>Rationale</th>'
        f'</tr></thead><tbody>{challenge_rows}</tbody></table></div></section>'
    )


def _conservation_section(data: ResearcherReportData) -> str:
    if data.conservation is None:
        return (
            '<section><h2>Target conservation</h2>'
            + _evidence_unavailable("conservation/conservation_report.json")
            + "</section>"
        )
    report = data.conservation
    counts = _dict(report.get("counts"))
    windows = [item for item in _list(report.get("windows")) if isinstance(item, dict)]
    chart = _polyline_chart(
        windows,
        (("mean_conservation", "Target conservation", "series-a"),),
        title="Recorded target conservation by reference position",
    )
    return (
        '<section><h2>Target conservation</h2>'
        '<p class="section-intro">Conservation describes agreement within the target '
        'population. It does not establish specificity against non-target organisms.</p>'
        '<div class="summary-grid small">'
        + _summary_card("Status", report.get("status"))
        + _summary_card("Discovery sequences", counts.get("sequences"))
        + _summary_card("Reference positions", counts.get("reference_positions", counts.get("positions")))
        + _summary_card("Windows", counts.get("windows"))
        + "</div>"
        + chart
        + "</section>"
    )


def _contrast_section(data: ResearcherReportData) -> str:
    if data.contrastive is None:
        return (
            '<section><h2>Target vs non-target contrast</h2>'
            + _evidence_unavailable(
                "contrastive_conservation/contrastive_conservation_report.json"
            )
            + "</section>"
        )
    report = data.contrastive
    windows = [item for item in _list(report.get("windows")) if isinstance(item, dict)]
    candidates = [
        item for item in _list(report.get("candidates")) if isinstance(item, dict)
    ]
    chart = _polyline_chart(
        windows,
        (
            ("target_mean_conservation", "Target conservation", "series-a"),
            ("worst_similarity", "Worst challenge similarity", "series-b"),
        ),
        title="Recorded target conservation and challenge similarity",
    )
    rows: list[str] = []
    for item in candidates:
        region = _dict(item.get("region"))
        rows.append(
            "<tr>"
            f"<td>{_text(region.get('rank'))}</td>"
            f"<td>{_text(region.get('region_id'))}</td>"
            f"<td>{_range(region.get('reference_start'), region.get('reference_end'))}</td>"
            f"<td><strong>{_range(region.get('peak_start'), region.get('peak_end'))}</strong></td>"
            f"<td>{_text(item.get('worst_dataset_name'))}</td>"
            f"<td>{_text(item.get('worst_dataset_criticality'))}</td>"
            f"<td>{_text(item.get('worst_similarity'))}</td>"
            f"<td>{_text(item.get('contrast_margin'))}</td>"
            "</tr>"
        )
    body = "".join(rows) or '<tr><td colspan="8">No contrastive candidate regions recorded.</td></tr>'
    return (
        '<section><h2>Target vs non-target contrast</h2>'
        '<p class="section-intro">The candidate region is the broad interval available to '
        'Primer3. The <strong>contrast anchor</strong> is the discriminant window that '
        'must be contained by a contrastive assay.</p>'
        + chart
        + '<div class="table-scroll"><table><thead><tr><th>Rank</th><th>Region</th>'
        '<th>Candidate region</th><th>Contrast anchor</th><th>Worst challenge</th>'
        '<th>Criticality</th><th>Similarity</th><th>Contrast margin</th></tr></thead><tbody>'
        + body
        + "</tbody></table></div></section>"
    )


def _candidate_map(data: ResearcherReportData) -> dict[str, dict[str, object]]:
    report = _dict(data.primer_design)
    return {
        str(item.get("region_id")): item
        for item in _list(report.get("candidates"))
        if isinstance(item, dict) and item.get("region_id") is not None
    }


def _oligo_row(role: str, oligo: dict[str, object]) -> str:
    return (
        "<tr>"
        f"<th>{escape(role)}</th>"
        f'<td class="mono sequence">{_text(oligo.get("sequence"))}</td>'
        f"<td>{_range(oligo.get('reference_start'), oligo.get('reference_end'))}</td>"
        f"<td>{_text(oligo.get('tm'))}</td>"
        f"<td>{_text(oligo.get('gc_percent'))}</td>"
        f"<td>{_text(oligo.get('penalty'))}</td>"
        "</tr>"
    )


def _assay_section(data: ResearcherReportData) -> str:
    if data.primer_design is None:
        return (
            '<section><h2>Assay design</h2>'
            + _evidence_unavailable("primer_design/primer_design_report.json")
            + "</section>"
        )
    report = data.primer_design
    source = report.get("candidate_source")
    candidates = _candidate_map(data)
    assays = [item for item in _list(report.get("assays")) if isinstance(item, dict)]
    blocks: list[str] = []
    for assay in assays:
        region = candidates.get(str(assay.get("region_id")), {})
        forward = _dict(assay.get("forward_primer"))
        probe = _dict(assay.get("probe"))
        reverse = _dict(assay.get("reverse_primer"))
        peak_start = region.get("peak_start")
        peak_end = region.get("peak_end")
        anchor_known = isinstance(peak_start, (int, float)) and isinstance(peak_end, (int, float))
        anchor_contained = (
            anchor_known
            and isinstance(forward.get("reference_start"), (int, float))
            and isinstance(reverse.get("reference_end"), (int, float))
            and float(forward["reference_start"]) <= float(peak_start)
            and float(peak_end) <= float(reverse["reference_end"])
        )
        anchor_text = _range(peak_start, peak_end) if source == "CONTRASTIVE_CONSERVATION" else "Not applicable"
        contained_text = "Yes" if anchor_contained else ("No" if anchor_known else "Not available")
        oligos = (
            _oligo_row("Forward primer", forward)
            + _oligo_row("Probe", probe)
            + _oligo_row("Reverse primer", reverse)
        )
        blocks.append(
            '<article class="assay-card">'
            f'<h3>{_text(assay.get("assay_id"))}</h3><dl>'
            f'<dt>Candidate source</dt><dd>{_text(source)}</dd>'
            f'<dt>Candidate region</dt><dd>{_range(region.get("reference_start"), region.get("reference_end"))}</dd>'
            f'<dt>Contrast anchor</dt><dd>{anchor_text}</dd>'
            f'<dt>Anchor contained</dt><dd><strong>{contained_text}</strong></dd>'
            f'<dt>Product size</dt><dd>{_text(assay.get("product_size"))}</dd>'
            f'<dt>Pair penalty</dt><dd>{_text(assay.get("pair_penalty"))}</dd></dl>'
            f'<p class="anchor-check">Anchor contained: {contained_text}</p>'
            '<div class="table-scroll"><table><thead><tr><th>Role</th><th>Sequence</th>'
            '<th>Reference coordinates</th><th>Tm</th><th>GC %</th><th>Penalty</th>'
            f'</tr></thead><tbody>{oligos}</tbody></table></div></article>'
        )
    content = "".join(blocks) or '<p class="muted">No designed assays were recorded.</p>'
    return '<section><h2>Assay design</h2>' + content + "</section>"


def _inclusivity_section(data: ResearcherReportData) -> str:
    if data.inclusivity is None:
        return (
            '<section><h2>Target coverage / inclusivity</h2>'
            + _evidence_unavailable("inclusivity/inclusivity_report.json")
            + "</section>"
        )
    report = data.inclusivity
    counts = _dict(report.get("counts"))
    proposals = [item for item in _list(report.get("proposals")) if isinstance(item, dict)]
    proposal_rows = "".join(
        "<tr>"
        f"<td>{_text(item.get('assay_id'))}</td>"
        f"<td>{_text(item.get('role'))}</td>"
        f'<td class="mono">{_text(item.get("original_sequence"))}</td>'
        f'<td class="mono">{_text(item.get("proposed_sequence"))}</td>'
        f"<td>{_text(item.get('status'))}</td>"
        f"<td>{_text(item.get('reason'))}</td>"
        "</tr>"
        for item in proposals
    ) or '<tr><td colspan="6">No degeneracy proposal changes recorded.</td></tr>'
    return (
        '<section><h2>Target coverage / inclusivity</h2><div class="summary-grid small">'
        + _summary_card("Status", report.get("status"))
        + _summary_card("Evaluation sequences", counts.get("evaluation_sequences"))
        + _summary_card("Assays", counts.get("assays"))
        + _summary_card("Original compatible", counts.get("original_compatible"))
        + _summary_card("Accepted proposals", counts.get("accepted_proposals"))
        + "</div>"
        '<h3>Degeneracy proposals</h3><div class="table-scroll"><table><thead><tr>'
        '<th>Assay</th><th>Role</th><th>Original</th><th>Proposed</th><th>Status</th><th>Reason</th>'
        f'</tr></thead><tbody>{proposal_rows}</tbody></table></div></section>'
    )


def _specificity_section(data: ResearcherReportData) -> str:
    if data.specificity is None:
        return (
            '<section><h2>Specificity</h2>'
            + _evidence_unavailable("specificity/specificity_report.json")
            + "</section>"
        )
    report = data.specificity
    counts = _dict(report.get("counts"))
    retention = [item for item in _list(report.get("retention")) if isinstance(item, dict)]
    rows = "".join(
        "<tr>"
        f"<td>{_text(item.get('dataset_name'))}</td>"
        f"<td>{_text(item.get('assay_id'))}</td>"
        f"<td>{_text(item.get('role'))}</td>"
        f"<td>{_text(item.get('total_hit_count'))}</td>"
        f"<td>{_text(item.get('retained_hit_count'))}</td>"
        f"<td>{_text(item.get('truncated'))}</td>"
        "</tr>"
        for item in retention
    ) or '<tr><td colspan="6">No retained off-target hit summary recorded.</td></tr>'
    return (
        '<section><h2>Specificity</h2>'
        '<p class="section-intro">Final specificity is assay-level evidence and remains '
        'independent from region-level contrast.</p><div class="summary-grid small">'
        + _summary_card("Status", report.get("status"))
        + _summary_card("Datasets", counts.get("datasets"))
        + _summary_card("Sequences", counts.get("sequences"))
        + _summary_card("Plausible amplicons", counts.get("plausible_amplicons"))
        + _summary_card("Detectable off-targets", counts.get("detectable_off_targets"))
        + "</div>"
        '<h3>Compatible hit retention</h3><div class="table-scroll"><table><thead><tr>'
        '<th>Dataset</th><th>Assay</th><th>Role</th><th>Total hits</th><th>Retained</th><th>Truncated</th>'
        f'</tr></thead><tbody>{rows}</tbody></table></div></section>'
    )


def _ranking_section(data: ResearcherReportData) -> str:
    if data.ranking is None:
        return (
            '<section><h2>Final candidates</h2>'
            + _evidence_unavailable("ranking/ranking_report.json")
            + "</section>"
        )
    assays = [item for item in _list(data.ranking.get("assays")) if isinstance(item, dict)]
    recommended = next(
        (item for item in assays if item.get("classification") == "IN SILICO PASS"),
        None,
    )
    recommendation = ""
    if recommended is not None:
        recommendation = (
            '<div class="recommendation"><strong>Recommended in-silico candidate</strong>'
            f'<p>{_text(recommended.get("assay_id"))} · rank {_text(recommended.get("rank"))}</p>'
            '<small>This recommendation is computational only and requires independent review and experimental validation.</small></div>'
        )
    rows: list[str] = []
    details: list[str] = []
    for assay in assays:
        reasons = [item for item in _list(assay.get("reasons")) if isinstance(item, dict)]
        reason_codes = ", ".join(str(item.get("code")) for item in reasons) or "None"
        rows.append(
            "<tr>"
            f"<td>{_text(assay.get('rank'))}</td>"
            f"<td>{_text(assay.get('assay_id'))}</td>"
            f"<td><span class=\"badge\">{_text(assay.get('classification'))}</span></td>"
            f"<td>{_text(assay.get('score_status'))}</td>"
            f"<td>{_text(assay.get('final_score'))}</td>"
            f"<td>{_text(reason_codes)}</td>"
            "</tr>"
        )
        components = _dict(assay.get("components"))
        component_table = _metric_table(
            components,
            (
                ("inclusivity", "Inclusivity"),
                ("specificity", "Specificity"),
                ("conservation", "Conservation"),
                ("primer3_quality", "Primer3 quality"),
                ("robustness", "Robustness"),
            ),
        )
        reason_items = "".join(
            "<li>"
            f"<strong>{_text(item.get('code'))}</strong> "
            f"[{_text(item.get('severity'))} / {_text(item.get('source'))}] — "
            f"{_text(item.get('message'))}</li>"
            for item in reasons
        ) or "<li>None</li>"
        details.append(
            '<details><summary>'
            f"Rank {_text(assay.get('rank'))} · {_text(assay.get('assay_id'))} · {_text(assay.get('classification'))}"
            '</summary><div class="two-col"><div><h4>Score components</h4>'
            + component_table
            + '</div><div><h4>Reason codes</h4><ul class="reasons">'
            + reason_items
            + "</ul></div></div></details>"
        )
    body = "".join(rows) or '<tr><td colspan="6">No ranked assays recorded.</td></tr>'
    return (
        '<section><h2>Final candidates</h2>'
        + recommendation
        + '<div class="table-scroll"><table><thead><tr><th>Rank</th><th>Assay</th>'
        '<th>Classification</th><th>Score status</th><th>Final score</th><th>Reason codes</th>'
        + f'</tr></thead><tbody>{body}</tbody></table></div>'
        + "".join(details)
        + "</section>"
    )


def _limitations_section() -> str:
    return (
        '<section><h2>Interpretation and limitations</h2><div class="limitations"><ul>'
        '<li>This report summarizes <strong>in-silico</strong> assay-discovery evidence.</li>'
        '<li>Region-level target/non-target contrast is not equivalent to final oligo specificity.</li>'
        '<li><code>IN SILICO PASS</code> does not establish experimental, clinical, or diagnostic validity.</li>'
        '<li>Candidate assays require independent scientific review, wet-lab optimization, and experimental validation before use.</li>'
        '</ul></div></section>'
    )


def _json_preview(value: object) -> str:
    try:
        text = json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)
    except (TypeError, ValueError):
        text = str(value)
    if len(text) > 12000:
        text = text[:12000] + "\n[truncated for report readability]"
    return escape(text)


def _reproducibility_section(data: ResearcherReportData) -> str:
    manifest = data.run_manifest
    reference = _dict(manifest.get("reference"))
    panel_provenance = _dict(manifest.get("panel_provenance"))
    return (
        '<section><h2>Reproducibility</h2><div class="two-col"><div><dl>'
        f'<dt>Run ID</dt><dd class="mono">{_text(manifest.get("run_id"))}</dd>'
        f'<dt>Created</dt><dd>{_text(manifest.get("created_at"))}</dd>'
        f'<dt>Updated</dt><dd>{_text(manifest.get("updated_at"))}</dd>'
        f'<dt>Reference</dt><dd>{_text(reference.get("id"))}</dd>'
        f'<dt>Reference mode</dt><dd>{_text(reference.get("mode"))}</dd>'
        f'<dt>Panel manifest</dt><dd class="mono">{_text(panel_provenance.get("manifest_sha256"))}</dd>'
        '</dl></div><div><p class="muted">The following values are taken from the sanitized run manifest.</p></div></div>'
        '<details><summary>Effective configuration</summary>'
        f'<pre>{_json_preview(manifest.get("effective_config", {}))}</pre></details>'
        '<details><summary>Environment and tool information</summary>'
        f'<pre>{_json_preview(manifest.get("environment", {}))}</pre></details>'
        '</section>'
    )


def render_researcher_report_html(data: ResearcherReportData) -> str:
    """Render the complete static researcher report from published artifacts."""
    target = _target_name(data)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<link rel="icon" href="data:,">
<title>Geison Researcher Report — {escape(target, quote=True)}</title>
<style>
:root{{--ink:#17202d;--muted:#647083;--line:#d9e0e9;--panel:#f7f9fc;--soft:#eef3f8;--pass:#e7f6ed;--pass-line:#9bc9ab;--review:#fff7dc;--review-line:#d9bd5b;--risk:#fdecec;--risk-line:#df9b9b;--accent:#245e72;--accent2:#8c5a75}}
*{{box-sizing:border-box}}body{{margin:0;font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--ink);background:#fff;line-height:1.5}}main{{max-width:1120px;margin:0 auto;padding:42px 28px 72px}}header{{padding-bottom:24px;border-bottom:1px solid var(--line);margin-bottom:28px}}h1{{font-size:2.25rem;letter-spacing:-.03em;margin:.15rem 0}}h2{{font-size:1.4rem;margin:0 0 1rem}}h3{{font-size:1.02rem;margin:1.35rem 0 .65rem}}h4{{margin:.6rem 0}}section{{margin:34px 0}}.eyebrow{{text-transform:uppercase;letter-spacing:.11em;font-size:.72rem;font-weight:750;color:var(--muted);margin:0 0 .4rem}}.muted,.section-intro{{color:var(--muted)}}.outcome{{border:1px solid var(--line);border-radius:14px;padding:22px 24px}}.outcome h2{{font-size:1.65rem;margin:.2rem 0 .45rem}}.outcome p:last-child{{margin-bottom:0}}.outcome-pass{{background:var(--pass);border-color:var(--pass-line)}}.outcome-review{{background:var(--review);border-color:var(--review-line)}}.outcome-risk{{background:var(--risk);border-color:var(--risk-line)}}.outcome-inconclusive{{background:var(--soft)}}.summary-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px}}.summary-grid.small{{grid-template-columns:repeat(auto-fit,minmax(135px,1fr))}}.summary-card{{border:1px solid var(--line);border-radius:10px;padding:12px;background:#fff;min-height:76px}}.summary-card strong{{font-size:1.25rem;display:block;margin-top:7px}}.summary-label{{font-size:.75rem;color:var(--muted);text-transform:uppercase;letter-spacing:.05em}}.warning-line,.unavailable{{border:1px solid var(--review-line);background:var(--review);padding:12px 14px;border-radius:9px}}.unavailable p{{margin:.25rem 0 0}}.two-col{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:22px}}dl{{display:grid;grid-template-columns:minmax(125px,max-content) 1fr;gap:7px 14px;margin:.4rem 0}}dt{{font-weight:700}}dd{{margin:0}}table{{width:100%;border-collapse:collapse;font-size:.88rem}}th,td{{padding:9px 10px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}}thead th{{background:var(--panel);font-size:.78rem;text-transform:uppercase;letter-spacing:.04em;color:#465268}}.compact th{{width:48%}}.table-scroll{{overflow-x:auto;border:1px solid var(--line);border-radius:10px}}.table-scroll table tr:last-child td{{border-bottom:0}}.mono,code,pre,.sequence{{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}}.sequence{{overflow-wrap:anywhere}}.assay-card{{border:1px solid var(--line);border-radius:12px;padding:16px 18px;margin:14px 0;background:#fff}}.anchor-check{{display:inline-block;background:var(--soft);border-radius:999px;padding:5px 10px;font-weight:700;font-size:.82rem}}.chart-wrap{{border:1px solid var(--line);border-radius:12px;padding:12px 16px;margin:16px 0;background:#fff}}svg{{display:block;width:100%;height:auto}}.axis{{stroke:#9aa6b4;stroke-width:1}}.axis-label{{font-size:12px;fill:#6d7987}}.series-a{{stroke:var(--accent)}}.series-b{{stroke:var(--accent2)}}.series-c{{stroke:#65734c}}.legend{{display:flex;gap:18px;flex-wrap:wrap;font-size:.8rem;color:var(--muted)}}.legend span{{display:flex;align-items:center;gap:6px}}.legend-line{{display:inline-block;width:24px;border-top:3px solid}}.recommendation{{background:var(--pass);border:1px solid var(--pass-line);border-radius:10px;padding:13px 15px;margin-bottom:14px}}.recommendation p{{margin:.3rem 0}}.badge{{font-weight:750}}details{{border:1px solid var(--line);border-radius:9px;margin:10px 0;padding:10px 12px}}summary{{cursor:pointer;font-weight:700}}.reasons{{padding-left:18px}}.limitations{{border-left:4px solid var(--accent);background:var(--panel);padding:12px 18px}}pre{{background:#111827;color:#eef2f7;padding:14px;border-radius:9px;overflow:auto;font-size:.78rem;white-space:pre-wrap}}footer{{margin-top:44px;padding-top:18px;border-top:1px solid var(--line);color:var(--muted);font-size:.8rem}}
@media(max-width:720px){{main{{padding:24px 14px 48px}}h1{{font-size:1.8rem}}.two-col{{grid-template-columns:1fr}}dl{{grid-template-columns:1fr}}}}
@media print{{body{{font-size:10pt}}main{{max-width:none;padding:0}}details{{break-inside:avoid}}.assay-card,.outcome,.table-scroll{{break-inside:avoid}}}}
</style>
</head>
<body>
<main>
<header><p class="eyebrow">Geison · qPCR assay-discovery evidence</p><h1>Geison Researcher Report</h1><p><strong>Target:</strong> {escape(target)}</p><p class="muted">A readable view of the scientific evidence published by this Geison run.</p></header>
{_outcome_section(data)}
{_run_summary_section(data)}
{_panel_section(data)}
{_conservation_section(data)}
{_contrast_section(data)}
{_assay_section(data)}
{_inclusivity_section(data)}
{_specificity_section(data)}
{_ranking_section(data)}
{_limitations_section()}
{_reproducibility_section(data)}
<footer>Generated from published Geison artifacts. This report does not replace the underlying JSON/TSV evidence or experimental validation.</footer>
</main>
</body>
</html>
"""
