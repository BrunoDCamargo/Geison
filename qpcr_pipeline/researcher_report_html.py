"""Self-contained HTML renderer for the final Geison researcher report."""

from __future__ import annotations

from html import escape
import json

from qpcr_pipeline.researcher_report import ResearcherReportData, scientific_outcome


def _obj(value: object) -> dict[str, object]:
    return value if isinstance(value, dict) else {}


def _arr(value: object) -> list[object]:
    return value if isinstance(value, list) else []


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


def _unavailable(path: str) -> str:
    return (
        '<div class="unavailable"><strong>Evidence unavailable</strong>'
        f'<p>The published artifact <code>{escape(path)}</code> is not available for this run.</p>'
        '</div>'
    )


def _card(label: str, value: object) -> str:
    return (
        '<div class="card">'
        f'<span>{escape(label)}</span><strong>{_text(value)}</strong>'
        '</div>'
    )


def _target_name(data: ResearcherReportData) -> str:
    target = _obj(_obj(_obj(data.panel).get("definition")).get("target")).get("name")
    value = target or data.run_manifest.get("target_name")
    return str(value) if value is not None else "Not available"


def _run_summary(data: ResearcherReportData) -> str:
    ranking_counts = _obj(_obj(data.ranking).get("counts"))
    contrast_counts = _obj(_obj(data.contrastive).get("counts"))
    primer_counts = _obj(_obj(data.primer_design).get("counts"))
    run_summary = _obj(data.run_summary)
    panel_definition = _obj(_obj(data.panel).get("definition"))
    challenge_count = sum(
        1
        for item in _arr(panel_definition.get("non_targets"))
        if isinstance(item, dict) and "CHALLENGE" in _arr(item.get("dataset_roles"))
    )
    cards = "".join(
        (
            _card("Technical status", data.run_manifest.get("status")),
            _card("Target sequences", run_summary.get("sequence_count")),
            _card(
                "Challenge datasets",
                contrast_counts.get("challenge_datasets", challenge_count),
            ),
            _card("Contrastive regions", contrast_counts.get("candidate_regions")),
            _card("Designed assays", ranking_counts.get("assays", primer_counts.get("assays"))),
            _card("IN SILICO PASS", ranking_counts.get("in_silico_pass", 0)),
            _card("REVIEW", ranking_counts.get("review", 0)),
            _card("HIGH_RISK", ranking_counts.get("high_risk", 0)),
        )
    )
    completeness = _obj(data.run_manifest.get("scientific_completeness"))
    missing = _arr(completeness.get("missing_evidence"))
    missing_html = ""
    if missing:
        missing_html = (
            '<p class="warning"><strong>Missing evidence:</strong> '
            + ", ".join(escape(str(item)) for item in missing)
            + "</p>"
        )
    return f'<section><h2>Run summary</h2><div class="cards">{cards}</div>{missing_html}</section>'


def _outcome(data: ResearcherReportData) -> str:
    status = str(data.run_manifest.get("status", "UNKNOWN"))
    title, body = scientific_outcome(status, data.ranking)
    if title == "In-silico candidate(s) identified":
        css = "pass"
    elif "review" in title.lower():
        css = "review"
    elif "acceptable" in title.lower() or "failed" in title.lower():
        css = "risk"
    else:
        css = "neutral"
    return (
        f'<section class="outcome {css}"><p class="eyebrow">Scientific outcome</p>'
        f'<h2>{escape(title)}</h2><p>{escape(body)}</p></section>'
    )


def _panel(data: ResearcherReportData) -> str:
    if data.panel is None:
        return '<section><h2>Approved panel and study context</h2>' + _unavailable(
            "panel/approved_panel.json"
        ) + '</section>'
    definition = _obj(data.panel.get("definition"))
    target = _obj(definition.get("target"))
    groups = [x for x in _arr(target.get("groups")) if isinstance(x, dict)]
    non_targets = [x for x in _arr(definition.get("non_targets")) if isinstance(x, dict)]
    context = _obj(definition.get("diagnostic_context"))
    group_rows = "".join(
        '<tr>'
        f'<td>{_text(x.get("name"))}</td>'
        f'<td>{_text(", ".join(str(v) for v in _arr(x.get("dataset_roles"))))}</td>'
        f'<td>{_text(x.get("required"))}</td>'
        f'<td>{_text("; ".join(str(v) for v in _arr(x.get("reasons"))))}</td>'
        '</tr>'
        for x in groups
    ) or '<tr><td colspan="4">No target groups recorded.</td></tr>'
    challenge_rows = "".join(
        '<tr>'
        f'<td>{_text(x.get("name"))}</td>'
        f'<td>{_text(x.get("criticality"))}</td>'
        f'<td>{_text(", ".join(str(v) for v in _arr(x.get("dataset_roles"))))}</td>'
        f'<td>{_text("; ".join(str(v) for v in _arr(x.get("reasons"))))}</td>'
        '</tr>'
        for x in non_targets
    ) or '<tr><td colspan="4">No challenge datasets recorded.</td></tr>'
    return (
        '<section><h2>Approved panel and study context</h2>'
        '<div class="two"><div><h3>Panel</h3><dl>'
        f'<dt>Status</dt><dd>{_text(data.panel.get("status"))}</dd>'
        f'<dt>Target</dt><dd>{_text(target.get("name"))}</dd>'
        f'<dt>Mode</dt><dd>{_text(target.get("mode"))}</dd>'
        f'<dt>Proposal hash</dt><dd class="mono">{_text(data.panel.get("proposal_sha256"))}</dd>'
        '</dl></div><div><h3>Diagnostic context</h3><dl>'
        f'<dt>Syndrome</dt><dd>{_text(context.get("syndrome"))}</dd>'
        f'<dt>Geography</dt><dd>{_text(context.get("geography"))}</dd>'
        f'<dt>Sample type</dt><dd>{_text(context.get("sample_type"))}</dd>'
        f'<dt>Vector / context</dt><dd>{_text(context.get("vector"))}</dd>'
        '</dl></div></div>'
        '<h3>Target groups</h3><div class="scroll"><table><thead><tr>'
        '<th>Group</th><th>Roles</th><th>Required</th><th>Rationale</th>'
        f'</tr></thead><tbody>{group_rows}</tbody></table></div>'
        '<h3>Challenge panel</h3><div class="scroll"><table><thead><tr>'
        '<th>Dataset</th><th>Criticality</th><th>Roles</th><th>Rationale</th>'
        f'</tr></thead><tbody>{challenge_rows}</tbody></table></div></section>'
    )


def _line_chart(
    rows: list[dict[str, object]],
    series: tuple[tuple[str, str, str], ...],
    title: str,
) -> str:
    points_by_row: list[tuple[float, dict[str, object]]] = []
    for row in rows:
        start, end = row.get("reference_start"), row.get("reference_end")
        if isinstance(start, (int, float)) and isinstance(end, (int, float)):
            points_by_row.append(((float(start) + float(end)) / 2.0, row))
    if not points_by_row:
        return '<p class="muted">No recorded window series available for this view.</p>'
    lo = min(x for x, _ in points_by_row)
    hi = max(x for x, _ in points_by_row)
    span = max(hi - lo, 1.0)
    colors = ("#245e72", "#8c5a75", "#65734c")
    lines: list[str] = []
    legend: list[str] = []
    for idx, (key, label, _) in enumerate(series):
        pts: list[str] = []
        for center, row in points_by_row:
            value = row.get(key)
            if not isinstance(value, (int, float)):
                continue
            x = 52 + (center - lo) / span * 626
            y = 18 + (1 - max(0.0, min(1.0, float(value)))) * 145
            pts.append(f"{x:.1f},{y:.1f}")
        if pts:
            color = colors[idx % len(colors)]
            dash = ' stroke-dasharray="7 5"' if idx else ""
            lines.append(
                f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" '
                f'stroke-width="3"{dash}/>'
            )
            legend.append(
                f'<span><i style="border-color:{color}"></i>{escape(label)}</span>'
            )
    return (
        '<div class="chart"><h4>' + escape(title) + '</h4>'
        '<svg viewBox="0 0 720 200" role="img" aria-label="Recorded metric series">'
        '<line x1="52" y1="18" x2="52" y2="163" class="axis"/>'
        '<line x1="52" y1="163" x2="678" y2="163" class="axis"/>'
        '<text x="12" y="24">1.0</text><text x="18" y="166">0.0</text>'
        f'<text x="52" y="188">{_text(int(lo))}</text>'
        f'<text x="644" y="188">{_text(int(hi))}</text>'
        + "".join(lines)
        + '</svg><div class="legend">' + "".join(legend) + '</div></div>'
    )


def _conservation(data: ResearcherReportData) -> str:
    if data.conservation is None:
        return '<section><h2>Target conservation</h2>' + _unavailable(
            "conservation/conservation_report.json"
        ) + '</section>'
    counts = _obj(data.conservation.get("counts"))
    json_windows = [x for x in _arr(data.conservation.get("windows")) if isinstance(x, dict)]
    windows = json_windows or list(data.conservation_windows)
    chart = _line_chart(
        windows,
        (("mean_conservation", "Target conservation", "a"),),
        "Recorded target conservation by reference position",
    )
    return (
        '<section><h2>Target conservation</h2>'
        '<p class="muted">Conservation describes agreement within the target population; '
        'it does not establish non-target specificity.</p><div class="cards small">'
        + _card("Status", data.conservation.get("status"))
        + _card("Discovery sequences", counts.get("sequences"))
        + _card("Reference positions", counts.get("reference_positions"))
        + _card("Windows", counts.get("windows"))
        + '</div>' + chart + '</section>'
    )


def _contrast(data: ResearcherReportData) -> str:
    if data.contrastive is None:
        return '<section><h2>Target vs non-target contrast</h2>' + _unavailable(
            "contrastive_conservation/contrastive_conservation_report.json"
        ) + '</section>'
    windows = [x for x in _arr(data.contrastive.get("windows")) if isinstance(x, dict)]
    chart = _line_chart(
        windows,
        (
            ("target_mean_conservation", "Target conservation", "a"),
            ("worst_similarity", "Worst challenge similarity", "b"),
        ),
        "Recorded target conservation and challenge similarity",
    )
    rows: list[str] = []
    for item in _arr(data.contrastive.get("candidates")):
        if not isinstance(item, dict):
            continue
        region = _obj(item.get("region"))
        rows.append(
            '<tr>'
            f'<td>{_text(region.get("rank"))}</td><td>{_text(region.get("region_id"))}</td>'
            f'<td>{_range(region.get("reference_start"), region.get("reference_end"))}</td>'
            f'<td><strong>{_range(region.get("peak_start"), region.get("peak_end"))}</strong></td>'
            f'<td>{_text(item.get("worst_dataset_name"))}</td>'
            f'<td>{_text(item.get("worst_dataset_criticality"))}</td>'
            f'<td>{_text(item.get("worst_similarity"))}</td>'
            f'<td>{_text(item.get("contrast_margin"))}</td></tr>'
        )
    body = "".join(rows) or '<tr><td colspan="8">No contrastive candidate regions recorded.</td></tr>'
    return (
        '<section><h2>Target vs non-target contrast</h2>'
        '<p class="muted">The <strong>Candidate region</strong> is the broad Primer3 design '
        'interval. The <strong>Contrast anchor</strong> is the discriminant window the '
        'amplicon must contain.</p>' + chart
        + '<div class="scroll"><table><thead><tr><th>Rank</th><th>Region</th>'
        '<th>Candidate region</th><th>Contrast anchor</th><th>Worst challenge</th>'
        '<th>Criticality</th><th>Similarity</th><th>Margin</th></tr></thead><tbody>'
        + body + '</tbody></table></div></section>'
    )


def _candidate_map(data: ResearcherReportData) -> dict[str, dict[str, object]]:
    return {
        str(x.get("region_id")): x
        for x in _arr(_obj(data.primer_design).get("candidates"))
        if isinstance(x, dict) and x.get("region_id") is not None
    }


def _oligo_row(label: str, oligo: dict[str, object]) -> str:
    return (
        f'<tr><th>{escape(label)}</th><td class="mono">{_text(oligo.get("sequence"))}</td>'
        f'<td>{_range(oligo.get("reference_start"), oligo.get("reference_end"))}</td>'
        f'<td>{_text(oligo.get("tm"))}</td><td>{_text(oligo.get("gc_percent"))}</td>'
        f'<td>{_text(oligo.get("penalty"))}</td></tr>'
    )


def _assays(data: ResearcherReportData) -> str:
    if data.primer_design is None:
        return '<section><h2>Assay design</h2>' + _unavailable(
            "primer_design/primer_design_report.json"
        ) + '</section>'
    source = data.primer_design.get("candidate_source")
    candidates = _candidate_map(data)
    blocks: list[str] = []
    for assay in _arr(data.primer_design.get("assays")):
        if not isinstance(assay, dict):
            continue
        region = candidates.get(str(assay.get("region_id")), {})
        forward = _obj(assay.get("forward_primer"))
        probe = _obj(assay.get("probe"))
        reverse = _obj(assay.get("reverse_primer"))
        ps, pe = region.get("peak_start"), region.get("peak_end")
        known = isinstance(ps, (int, float)) and isinstance(pe, (int, float))
        contained = (
            known
            and isinstance(forward.get("reference_start"), (int, float))
            and isinstance(reverse.get("reference_end"), (int, float))
            and float(forward["reference_start"]) <= float(ps)
            and float(pe) <= float(reverse["reference_end"])
        )
        contained_text = "Yes" if contained else ("No" if known else "Not available")
        anchor = _range(ps, pe) if source == "CONTRASTIVE_CONSERVATION" else "Not applicable"
        blocks.append(
            '<article class="assay"><h3>' + _text(assay.get("assay_id")) + '</h3><dl>'
            f'<dt>Candidate source</dt><dd>{_text(source)}</dd>'
            f'<dt>Candidate region</dt><dd>{_range(region.get("reference_start"), region.get("reference_end"))}</dd>'
            f'<dt>Contrast anchor</dt><dd>{anchor}</dd>'
            f'<dt>Product size</dt><dd>{_text(assay.get("product_size"))}</dd>'
            f'<dt>Pair penalty</dt><dd>{_text(assay.get("pair_penalty"))}</dd></dl>'
            f'<p class="anchor">Anchor contained: {contained_text}</p>'
            '<div class="scroll"><table><thead><tr><th>Role</th><th>Sequence</th>'
            '<th>Coordinates</th><th>Tm</th><th>GC %</th><th>Penalty</th></tr></thead><tbody>'
            + _oligo_row("Forward primer", forward)
            + _oligo_row("Probe", probe)
            + _oligo_row("Reverse primer", reverse)
            + '</tbody></table></div></article>'
        )
    return '<section><h2>Assay design</h2>' + (
        "".join(blocks) or '<p class="muted">No designed assays were recorded.</p>'
    ) + '</section>'


def _inclusivity(data: ResearcherReportData) -> str:
    if data.inclusivity is None:
        return '<section><h2>Target coverage / inclusivity</h2>' + _unavailable(
            "inclusivity/inclusivity_report.json"
        ) + '</section>'
    counts = _obj(data.inclusivity.get("counts"))
    proposal_rows = "".join(
        '<tr>'
        f'<td>{_text(x.get("assay_id"))}</td><td>{_text(x.get("role"))}</td>'
        f'<td class="mono">{_text(x.get("original_sequence"))}</td>'
        f'<td class="mono">{_text(x.get("proposed_sequence"))}</td>'
        f'<td>{_text(x.get("status"))}</td><td>{_text(x.get("reason"))}</td></tr>'
        for x in _arr(data.inclusivity.get("proposals")) if isinstance(x, dict)
    ) or '<tr><td colspan="6">No degeneracy proposal changes recorded.</td></tr>'
    return (
        '<section><h2>Target coverage / inclusivity</h2><div class="cards small">'
        + _card("Status", data.inclusivity.get("status"))
        + _card("Evaluation sequences", counts.get("evaluation_sequences"))
        + _card("Assays", counts.get("assays"))
        + _card("Original compatible", counts.get("original_compatible"))
        + _card("Accepted proposals", counts.get("accepted_proposals"))
        + '</div><h3>Degeneracy proposals</h3><div class="scroll"><table><thead><tr>'
        '<th>Assay</th><th>Role</th><th>Original</th><th>Proposed</th><th>Status</th><th>Reason</th>'
        f'</tr></thead><tbody>{proposal_rows}</tbody></table></div></section>'
    )


def _specificity(data: ResearcherReportData) -> str:
    if data.specificity is None:
        return '<section><h2>Specificity</h2>' + _unavailable(
            "specificity/specificity_report.json"
        ) + '</section>'
    counts = _obj(data.specificity.get("counts"))
    retention_rows = "".join(
        '<tr>'
        f'<td>{_text(x.get("dataset_name"))}</td><td>{_text(x.get("assay_id"))}</td>'
        f'<td>{_text(x.get("role"))}</td><td>{_text(x.get("total_hit_count"))}</td>'
        f'<td>{_text(x.get("retained_hit_count"))}</td><td>{_text(x.get("truncated"))}</td></tr>'
        for x in _arr(data.specificity.get("retention")) if isinstance(x, dict)
    ) or '<tr><td colspan="6">No retained off-target hit summary recorded.</td></tr>'
    amplicon_rows = "".join(
        '<tr>'
        f'<td>{_text(x.get("dataset_name"))}</td><td>{_text(x.get("assay_id"))}</td>'
        f'<td>{_text(x.get("sequence_id"))}</td>'
        f'<td>{_range(x.get("source_start"), x.get("source_end"))}</td>'
        f'<td>{_text(x.get("amplicon_size"))}</td>'
        f'<td>{_text(x.get("primer_amplicon_plausible"))}</td>'
        f'<td>{_text(x.get("detectable_off_target"))}</td></tr>'
        for x in data.specificity_amplicons
    ) or '<tr><td colspan="7">No plausible off-target amplicons recorded.</td></tr>'
    return (
        '<section><h2>Specificity</h2><p class="muted">Specificity is evaluated at the '
        'complete assay level and remains independent from region-level contrast.</p>'
        '<div class="cards small">'
        + _card("Status", data.specificity.get("status"))
        + _card("Datasets", counts.get("datasets"))
        + _card("Sequences", counts.get("sequences"))
        + _card("Plausible amplicons", counts.get("plausible_amplicons"))
        + _card("Detectable off-targets", counts.get("detectable_off_targets"))
        + '</div><h3>Compatible hit retention</h3><div class="scroll"><table><thead><tr>'
        '<th>Dataset</th><th>Assay</th><th>Role</th><th>Total hits</th><th>Retained</th><th>Truncated</th>'
        f'</tr></thead><tbody>{retention_rows}</tbody></table></div>'
        '<h3>Plausible off-target amplicons</h3><div class="scroll"><table><thead><tr>'
        '<th>Dataset</th><th>Assay</th><th>Sequence</th><th>Source range</th><th>Size</th>'
        '<th>Plausible</th><th>Detectable off-target</th></tr></thead><tbody>'
        f'{amplicon_rows}</tbody></table></div></section>'
    )


def _ranking(data: ResearcherReportData) -> str:
    if data.ranking is None:
        return '<section><h2>Final candidates</h2>' + _unavailable(
            "ranking/ranking_report.json"
        ) + '</section>'
    assays = [x for x in _arr(data.ranking.get("assays")) if isinstance(x, dict)]
    recommended = next(
        (x for x in assays if x.get("classification") == "IN SILICO PASS"), None
    )
    recommendation = ""
    if recommended is not None:
        recommendation = (
            '<div class="recommend"><strong>Recommended in-silico candidate</strong>'
            f'<p>{_text(recommended.get("assay_id"))} · rank {_text(recommended.get("rank"))}</p>'
            '<small>Computational recommendation only; independent review and experimental validation are required.</small></div>'
        )
    rows: list[str] = []
    details: list[str] = []
    for assay in assays:
        reasons = [x for x in _arr(assay.get("reasons")) if isinstance(x, dict)]
        codes = ", ".join(str(x.get("code")) for x in reasons) or "None"
        rows.append(
            '<tr>'
            f'<td>{_text(assay.get("rank"))}</td><td>{_text(assay.get("assay_id"))}</td>'
            f'<td><strong>{_text(assay.get("classification"))}</strong></td>'
            f'<td>{_text(assay.get("score_status"))}</td><td>{_text(assay.get("final_score"))}</td>'
            f'<td>{_text(codes)}</td></tr>'
        )
        components = _obj(assay.get("components"))
        reason_items = "".join(
            '<li>'
            f'<strong>{_text(x.get("code"))}</strong> [{_text(x.get("severity"))} / {_text(x.get("source"))}] — '
            f'{_text(x.get("message"))}</li>'
            for x in reasons
        ) or '<li>None</li>'
        details.append(
            '<details><summary>'
            f'Rank {_text(assay.get("rank"))} · {_text(assay.get("assay_id"))} · {_text(assay.get("classification"))}'
            '</summary><div class="two"><div><h4>Score components</h4><dl>'
            f'<dt>Inclusivity</dt><dd>{_text(components.get("inclusivity"))}</dd>'
            f'<dt>Specificity</dt><dd>{_text(components.get("specificity"))}</dd>'
            f'<dt>Conservation</dt><dd>{_text(components.get("conservation"))}</dd>'
            f'<dt>Primer3 quality</dt><dd>{_text(components.get("primer3_quality"))}</dd>'
            f'<dt>Robustness</dt><dd>{_text(components.get("robustness"))}</dd>'
            '</dl></div><div><h4>Reason codes</h4><ul>' + reason_items + '</ul></div></div></details>'
        )
    body = "".join(rows) or '<tr><td colspan="6">No ranked assays recorded.</td></tr>'
    return (
        '<section><h2>Final candidates</h2>' + recommendation
        + '<div class="scroll"><table><thead><tr><th>Rank</th><th>Assay</th><th>Classification</th>'
        '<th>Score status</th><th>Final score</th><th>Reason codes</th></tr></thead><tbody>'
        + body + '</tbody></table></div>' + "".join(details) + '</section>'
    )


def _limitations() -> str:
    return (
        '<section><h2>Interpretation and limitations</h2><div class="limits"><ul>'
        '<li>Results are <strong>in silico</strong>.</li>'
        '<li>Region-level contrast is not equivalent to oligo specificity.</li>'
        '<li><code>IN SILICO PASS</code> does not establish experimental, clinical, or diagnostic validity.</li>'
        '<li>Candidate assays require independent scientific review and experimental validation.</li>'
        '</ul></div></section>'
    )


def _preview(value: object) -> str:
    try:
        raw = json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)
    except (TypeError, ValueError):
        raw = str(value)
    if len(raw) > 12000:
        raw = raw[:12000] + "\n[truncated for report readability]"
    return escape(raw)


def _reproducibility(data: ResearcherReportData) -> str:
    manifest = data.run_manifest
    reference = _obj(manifest.get("reference"))
    panel_provenance = _obj(manifest.get("panel_provenance"))
    return (
        '<section><h2>Reproducibility</h2><dl>'
        f'<dt>Run ID</dt><dd class="mono">{_text(manifest.get("run_id"))}</dd>'
        f'<dt>Created</dt><dd>{_text(manifest.get("created_at"))}</dd>'
        f'<dt>Updated</dt><dd>{_text(manifest.get("updated_at"))}</dd>'
        f'<dt>Reference</dt><dd>{_text(reference.get("id"))}</dd>'
        f'<dt>Reference mode</dt><dd>{_text(reference.get("mode"))}</dd>'
        f'<dt>Panel manifest</dt><dd class="mono">{_text(panel_provenance.get("manifest_sha256"))}</dd>'
        '</dl><details><summary>Effective configuration</summary><pre>'
        + _preview(manifest.get("effective_config", {}))
        + '</pre></details><details><summary>Environment and tool information</summary><pre>'
        + _preview(manifest.get("environment", {}))
        + '</pre></details></section>'
    )


def render_researcher_report_html(data: ResearcherReportData) -> str:
    """Render a static researcher report from already-published Geison evidence."""
    target = _target_name(data)
    return f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<link rel="icon" href="data:,">
<title>Geison Researcher Report — {escape(target, quote=True)}</title>
<style>
:root{{--ink:#17202d;--muted:#647083;--line:#d9e0e9;--panel:#f7f9fc;--pass:#e7f6ed;--review:#fff7dc;--risk:#fdecec;--accent:#245e72}}
*{{box-sizing:border-box}}body{{margin:0;font-family:system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--ink);line-height:1.5}}main{{max-width:1120px;margin:auto;padding:40px 28px 70px}}header{{border-bottom:1px solid var(--line);padding-bottom:22px;margin-bottom:28px}}h1{{font-size:2.2rem;margin:.2rem 0}}h2{{font-size:1.4rem;margin:0 0 1rem}}h3{{font-size:1rem;margin:1.3rem 0 .6rem}}h4{{margin:.5rem 0}}section{{margin:32px 0}}.eyebrow{{font-size:.72rem;text-transform:uppercase;letter-spacing:.1em;font-weight:750;color:var(--muted)}}.muted{{color:var(--muted)}}.outcome{{padding:20px 22px;border:1px solid var(--line);border-radius:13px}}.outcome.pass{{background:var(--pass)}}.outcome.review{{background:var(--review)}}.outcome.risk{{background:var(--risk)}}.outcome.neutral{{background:var(--panel)}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:9px}}.card{{border:1px solid var(--line);border-radius:9px;padding:11px;background:#fff}}.card span{{display:block;color:var(--muted);font-size:.72rem;text-transform:uppercase;letter-spacing:.04em}}.card strong{{display:block;font-size:1.2rem;margin-top:6px}}.warning,.unavailable{{background:var(--review);border:1px solid #d8bd61;border-radius:9px;padding:11px 13px}}.two{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}}dl{{display:grid;grid-template-columns:max-content 1fr;gap:6px 12px}}dt{{font-weight:700}}dd{{margin:0}}table{{width:100%;border-collapse:collapse;font-size:.87rem}}th,td{{text-align:left;padding:8px 9px;border-bottom:1px solid var(--line);vertical-align:top}}thead th{{background:var(--panel);font-size:.75rem;text-transform:uppercase;letter-spacing:.03em;color:#4b5668}}.scroll{{overflow-x:auto;border:1px solid var(--line);border-radius:9px}}.scroll table tr:last-child td{{border-bottom:0}}.mono,code,pre{{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}}.assay{{border:1px solid var(--line);border-radius:11px;padding:15px 17px;margin:12px 0}}.anchor{{display:inline-block;background:var(--panel);padding:5px 10px;border-radius:999px;font-weight:700}}.chart{{border:1px solid var(--line);border-radius:11px;padding:12px 15px;margin:15px 0}}svg{{width:100%;height:auto}}svg text{{font-size:12px;fill:#657181}}.axis{{stroke:#9aa5b2;stroke-width:1}}.legend{{display:flex;gap:18px;flex-wrap:wrap;font-size:.8rem;color:var(--muted)}}.legend span{{display:flex;gap:6px;align-items:center}}.legend i{{display:inline-block;width:22px;border-top:3px solid}}.recommend{{background:var(--pass);border:1px solid #9cc8ab;border-radius:9px;padding:12px 14px;margin-bottom:12px}}.recommend p{{margin:.25rem 0}}details{{border:1px solid var(--line);border-radius:8px;padding:9px 11px;margin:9px 0}}summary{{font-weight:700;cursor:pointer}}.limits{{border-left:4px solid var(--accent);background:var(--panel);padding:10px 16px}}pre{{white-space:pre-wrap;background:#111827;color:#eef2f7;padding:13px;border-radius:8px;overflow:auto;font-size:.78rem}}footer{{border-top:1px solid var(--line);padding-top:16px;margin-top:40px;color:var(--muted);font-size:.8rem}}
@media(max-width:720px){{main{{padding:24px 14px 45px}}.two{{grid-template-columns:1fr}}dl{{grid-template-columns:1fr}}}}
@media print{{main{{max-width:none;padding:0}}.outcome,.assay,.scroll{{break-inside:avoid}}}}
</style>
</head>
<body><main>
<header><p class="eyebrow">Geison · qPCR assay-discovery evidence</p><h1>Geison Researcher Report</h1><p><strong>Target:</strong> {escape(target)}</p><p class="muted">Readable view of the scientific evidence published by this Geison run.</p></header>
{_outcome(data)}
{_run_summary(data)}
{_panel(data)}
{_conservation(data)}
{_contrast(data)}
{_assays(data)}
{_inclusivity(data)}
{_specificity(data)}
{_ranking(data)}
{_limitations()}
{_reproducibility(data)}
<footer>Generated from published Geison artifacts. This report does not replace the underlying evidence or experimental validation.</footer>
</main></body></html>'''
