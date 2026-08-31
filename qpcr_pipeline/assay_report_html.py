"""Deterministic, self-contained HTML for final assay ranking results."""

from __future__ import annotations

from html import escape

from qpcr_pipeline.inclusivity import InclusivityResult
from qpcr_pipeline.primer_design import PrimerDesignResult
from qpcr_pipeline.ranking import RankedAssay
from qpcr_pipeline.specificity import SpecificityResult


_ROLE_ORDER = {"FORWARD": 0, "PROBE": 1, "REVERSE": 2}


def _text(value: object) -> str:
    if value is None:
        return "Not available"
    if isinstance(value, float):
        return f"{value:.6g}"
    return escape(str(value), quote=True)


def _score(value: float | None) -> str:
    return "Not available" if value is None else f"{value:.2f}"


def render_assay_report_html(
    *,
    target_name: str,
    primer_design: PrimerDesignResult,
    inclusivity: InclusivityResult,
    specificity: SpecificityResult,
    assays: tuple[RankedAssay, ...],
) -> str:
    """Render a static, escaped, offline final assay report."""
    primer_by_id = {item.assay_id: item for item in primer_design.assays}
    region_by_id = {item.region_id: item for item in primer_design.candidates}
    proposals_by_assay: dict[str, list[object]] = {}
    for proposal in inclusivity.proposals:
        proposals_by_assay.setdefault(proposal.assay_id, []).append(proposal)
    for proposals in proposals_by_assay.values():
        proposals.sort(key=lambda item: (_ROLE_ORDER.get(item.role, 99), item.status))

    specificity_counts: dict[tuple[str, str], list[int]] = {}
    for dataset_name in specificity.dataset_names:
        for ranked in assays:
            specificity_counts[(ranked.assay_id, dataset_name)] = [0, 0, 0]
    for row in specificity.retention:
        key = (row.assay_id, row.dataset_name)
        if key in specificity_counts:
            specificity_counts[key][0] += row.total_hit_count
    for amplicon in specificity.amplicons:
        key = (amplicon.assay_id, amplicon.dataset_name)
        if key not in specificity_counts:
            continue
        if amplicon.primer_amplicon_plausible:
            specificity_counts[key][1] += 1
        if amplicon.detectable_off_target:
            specificity_counts[key][2] += 1

    summary_rows: list[str] = []
    detail_blocks: list[str] = []
    for ranked in assays:
        assay = primer_by_id[ranked.assay_id]
        region = region_by_id[ranked.region_id]
        reason_codes = "; ".join(reason.code for reason in ranked.reasons) or "None"
        summary_rows.append(
            "<tr>"
            f"<td>{ranked.rank}</td>"
            f"<td>{escape(ranked.assay_id)}</td>"
            f"<td>{escape(ranked.classification)}</td>"
            f"<td>{_score(ranked.final_score)}</td>"
            f"<td>{escape(reason_codes)}</td>"
            "</tr>"
        )

        oligo_rows = []
        for role, oligo in (
            ("FORWARD", assay.forward_primer),
            ("PROBE", assay.probe),
            ("REVERSE", assay.reverse_primer),
        ):
            oligo_rows.append(
                "<tr>"
                f"<th>{role}</th>"
                f"<td class=sequence>{escape(oligo.sequence)}</td>"
                f"<td>{oligo.reference_start}-{oligo.reference_end}</td>"
                f"<td>{_text(oligo.tm)}</td>"
                f"<td>{_text(oligo.gc_percent)}</td>"
                f"<td>{_text(oligo.penalty)}</td>"
                "</tr>"
            )

        proposal_rows = []
        for proposal in proposals_by_assay.get(ranked.assay_id, []):
            proposal_rows.append(
                "<tr>"
                f"<td>{escape(proposal.role)}</td>"
                f"<td class=sequence>{escape(proposal.original_sequence)}</td>"
                f"<td class=sequence>{escape(proposal.proposed_sequence)}</td>"
                f"<td>{escape(proposal.status)}</td>"
                f"<td>{proposal.original_degeneracy} → {proposal.proposed_degeneracy}</td>"
                f"<td>{escape(proposal.reason)}</td>"
                "</tr>"
            )
        if not proposal_rows:
            proposal_rows.append('<tr><td colspan="6">No IUPAC proposal</td></tr>')

        specificity_rows = []
        for dataset_name in specificity.dataset_names:
            compatible_hits, plausible, detectable = specificity_counts.get(
                (ranked.assay_id, dataset_name), [0, 0, 0]
            )
            specificity_rows.append(
                "<tr>"
                f"<td>{escape(dataset_name)}</td>"
                f"<td>{compatible_hits}</td>"
                f"<td>{plausible}</td>"
                f"<td>{detectable}</td>"
                "</tr>"
            )
        if not specificity_rows:
            specificity_rows.append(
                '<tr><td colspan="4">Specificity evidence not available</td></tr>'
            )

        components = ranked.components
        component_rows = "".join(
            "<tr>"
            f"<td>{name}</td>"
            f"<td>{_text(value)}</td>"
            "</tr>"
            for name, value in (
                ("inclusivity", components.inclusivity),
                ("specificity", components.specificity),
                ("conservation", components.conservation),
                ("primer3_quality", components.primer3_quality),
                ("robustness", components.robustness),
            )
        )

        inclusivity_value = (
            "Not available"
            if ranked.original_compatible_count is None
            or ranked.evaluation_sequence_count is None
            else f"{ranked.original_compatible_count} / {ranked.evaluation_sequence_count}"
        )
        inclusivity_fraction = (
            None
            if ranked.components.inclusivity is None
            else ranked.components.inclusivity
        )

        reason_items = "".join(
            "<li>"
            f"<strong>{escape(reason.code)}</strong> "
            f"[{escape(reason.severity)} / {escape(reason.source)}]: "
            f"{escape(reason.message)}"
            "</li>"
            for reason in ranked.reasons
        ) or "<li>None</li>"

        detail_blocks.append(
            f'<details class="assay"><summary>Rank {ranked.rank} · '
            f'{escape(ranked.assay_id)} · {escape(ranked.classification)} · '
            f'Final score {_score(ranked.final_score)}</summary>'
            '<div class="grid">'
            '<section><h3>Assay design</h3>'
            '<table><thead><tr><th>Role</th><th>Sequence</th><th>Coordinates</th>'
            '<th>Tm</th><th>GC %</th><th>Penalty</th></tr></thead><tbody>'
            + "".join(oligo_rows)
            + "</tbody></table>"
            '<dl>'
            f'<dt>Product size</dt><dd>{assay.product_size}</dd>'
            f'<dt>Pair penalty</dt><dd>{_text(assay.pair_penalty)}</dd>'
            '</dl></section>'
            '<section><h3>Conservation</h3><dl>'
            f'<dt>Mean conservation</dt><dd>{_text(region.mean_conservation)}</dd>'
            f'<dt>Minimum conservation</dt><dd>{_text(region.minimum_conservation)}</dd>'
            f'<dt>Mean coverage</dt><dd>{_text(region.mean_coverage)}</dd>'
            f'<dt>Mean gap frequency</dt><dd>{_text(region.mean_gap_frequency)}</dd>'
            f'<dt>Mean entropy</dt><dd>{_text(region.mean_entropy_bits)}</dd>'
            '</dl></section>'
            '<section><h3>Inclusivity and degeneracy</h3><dl>'
            f'<dt>Original compatible</dt><dd>{inclusivity_value}</dd>'
            f'<dt>Original inclusivity fraction</dt><dd>{_text(inclusivity_fraction)}</dd>'
            '</dl><table><thead><tr><th>Role</th><th>Original</th><th>Proposal</th>'
            '<th>Status</th><th>Degeneracy</th><th>Reason</th></tr></thead><tbody>'
            + "".join(proposal_rows)
            + "</tbody></table></section>"
            '<section><h3>Specificity and score</h3>'
            '<table><thead><tr><th>Dataset</th><th>Compatible hits</th>'
            '<th>Plausible amplicons</th><th>Detectable off-targets</th></tr></thead><tbody>'
            + "".join(specificity_rows)
            + "</tbody></table>"
            f'<p><strong>Score status:</strong> {escape(ranked.score_status)} · '
            f'<strong>Final score:</strong> {_score(ranked.final_score)}</p>'
            '<table><thead><tr><th>Component</th><th>Value</th></tr></thead><tbody>'
            + component_rows
            + "</tbody></table>"
            '<h4>Reason codes</h4><ul>'
            + reason_items
            + "</ul></section></div></details>"
        )

    empty_state = "" if assays else '<p class="empty">No assay candidates</p>'
    summary_body = "".join(summary_rows) or '<tr><td colspan="5">No assay candidates</td></tr>'
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<link rel="icon" href="data:,">
<title>Geison final assay report</title>
<style>
:root{{--ink:#182132;--muted:#5d6878;--line:#d8dee8;--panel:#f6f8fb;--pass:#e8f6ed;--review:#fff7db;--risk:#fdeceb}}
*{{box-sizing:border-box}}body{{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;max-width:86rem;margin:0 auto;padding:1.5rem;color:var(--ink)}}
h1,h2,h3,h4{{line-height:1.2}}.muted{{color:var(--muted)}}.panel,details{{border:1px solid var(--line);border-radius:.6rem;margin:1rem 0;padding:1rem;background:var(--panel)}}
table{{width:100%;border-collapse:collapse;font-size:.9rem}}th,td{{padding:.45rem;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}}.sequence{{font-family:ui-monospace,monospace;overflow-wrap:anywhere}}
summary{{font-weight:700;cursor:pointer}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(27rem,1fr));gap:1rem}}section{{overflow-x:auto}}dl{{display:grid;grid-template-columns:max-content 1fr;gap:.35rem .8rem}}dt{{font-weight:650}}dd{{margin:0}}.empty{{padding:1rem;border:1px solid var(--line)}}
@media(max-width:40rem){{body{{padding:.8rem}}.grid{{grid-template-columns:1fr}}dl{{grid-template-columns:1fr}}}}
</style>
</head>
<body>
<main>
<h1>Final assay ranking</h1>
<p><strong>Target:</strong> {escape(target_name, quote=True)}</p>
<p class="muted">In silico classification and ranking of the original Primer3 assays. IUPAC proposals are contextual evidence and do not replace the original oligos.</p>
{empty_state}
<section class="panel"><h2>Ordered assays</h2><table><thead><tr><th>Rank</th><th>Assay</th><th>Class</th><th>Final score</th><th>Reason codes</th></tr></thead><tbody>{summary_body}</tbody></table></section>
{''.join(detail_blocks)}
</main>
</body>
</html>
"""
