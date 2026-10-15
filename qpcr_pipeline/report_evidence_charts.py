"""Offline views of recorded window metrics, without scientific classification."""
from __future__ import annotations

from dataclasses import dataclass
from html import escape
import math

from qpcr_pipeline.researcher_report import ResearcherReportError


@dataclass(frozen=True)
class WindowMatrix:
    windows: tuple[tuple[int, int], ...]
    datasets: tuple[str, ...]
    values: tuple[tuple[float | None, ...], ...]
    states: tuple[tuple[str, ...], ...]
    evidence: tuple[tuple[dict | None, ...], ...]


def fraction(value: object, context: str) -> float | None:
    """Accept a recorded fraction or absence; never silently clamp invalid data."""
    if value is None:
        return None
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or not 0 <= value <= 1):
        raise ResearcherReportError(f'{context} must be a finite number in [0, 1].')
    return float(value)


def _rows(value: object) -> list[dict]:
    if value is None:
        return []
    if not isinstance(value, (list, tuple)) or any(not isinstance(r, dict) for r in value):
        raise ResearcherReportError("Metric rows must be a sequence of objects.")
    return list(value)


def interval(row: dict) -> tuple[int, int]:
    start, end = row.get('reference_start'), row.get('reference_end')
    if (type(start) is not int or type(end) is not int or start < 1 or end < start):
        raise ResearcherReportError('Window coordinates must be positive inclusive integers.')
    return start, end


def window_matrix(report: dict | None) -> WindowMatrix:
    """Keep per-dataset values and unavailable states distinct from the max curve."""
    empty = WindowMatrix((), (), (), (), ())
    if not report or report.get('status') != 'COMPLETE':
        return empty
    windows = tuple(sorted(interval(row) for row in _rows(report.get('windows'))))
    if len(set(windows)) != len(windows):
        raise ResearcherReportError('Duplicate windows in contrastive evidence.')
    declarations = _rows(report.get('challenge_datasets'))
    names = [row.get('name') for row in declarations]
    rows = _rows(report.get('dataset_evidence'))
    if not names:
        names = sorted({row.get('dataset_name') for row in rows}, key=str)
    if any(not isinstance(name, str) or not name.strip() for name in names):
        raise ResearcherReportError('Each dataset needs a nonblank name.')
    if len(set(names)) != len(names):
        raise ResearcherReportError('Duplicate dataset declarations.')
    datasets = tuple(names)
    lookup = {}
    for row in rows:
        name = row.get('dataset_name')
        span = interval(row)
        key = (name, span)
        if name not in datasets or span not in windows:
            raise ResearcherReportError('Dataset evidence is outside the declared matrix.')
        if key in lookup:
            raise ResearcherReportError('Duplicate dataset/window evidence.')
        value = fraction(row.get('similarity'), f'{name} {span} similarity')
        count = row.get('sequence_count')
        if count is not None and (type(count) is not int or count < 0):
            raise ResearcherReportError('Sequence count must be a nonnegative integer.')
        state = 'observed'
        if count == 0:
            state, value = 'unevaluated', None
        elif 'best_sequence_id' in row and row['best_sequence_id'] is None:
            state, value = 'no-match', None
        elif value is None:
            state = 'missing'
        lookup[key] = (value, state, row)
    cells = [[lookup.get((name, span), (None, 'missing', None)) for span in windows]
             for name in datasets]
    return WindowMatrix(windows, datasets,
                        tuple(tuple(c[0] for c in row) for row in cells),
                        tuple(tuple(c[1] for c in row) for row in cells),
                        tuple(tuple(c[2] for c in row) for row in cells))


def _shade(value: float) -> str:
    # One fixed numerical scale across datasets; no stacked transparency.
    low, high = (239, 246, 255), (35, 94, 154)
    channels = [round(a + value * (b - a)) for a, b in zip(low, high)]
    return 'rgb(' + ','.join(map(str, channels)) + ')'


def render_window_heatmap(report: dict | None) -> str:
    matrix = window_matrix(report)
    if not matrix.windows or not matrix.datasets:
        return '<p class="muted">Per-dataset window matrix: Not available.</p>'
    left, width, row_height = 215.0, 740.0, 30
    cell_width = width / len(matrix.windows)
    height = 85 + row_height * len(matrix.datasets)
    svg, table, evidence_table = [], [], []
    labels = {'missing': 'Not available', 'no-match': 'No recorded best match',
              'unevaluated': 'Not evaluated'}
    for i, dataset in enumerate(matrix.datasets):
        y = 20 + i * row_height
        svg.append(f'<text x="5" y="{y + 20}">{escape(dataset)}</text>')
        cells = []
        for j, span in enumerate(matrix.windows):
            value, state = matrix.values[i][j], matrix.states[i][j]
            row = matrix.evidence[i][j] or {}
            text = format(value, '.4g') if value is not None else labels[state]
            details = '; '.join((dataset, f'window {span[0]}-{span[1]}',
                                f'normalized local alignment score: {text}',
                                f'sequences: {row.get("sequence_count", "Not available")}',
                                f'best sequence: {row.get("best_sequence_id", "Not available")}',
                                f'orientation: {row.get("best_orientation", "Not available")}'))
            raw_values = (dataset, f'{span[0]}-{span[1]}',
                          value if value is not None else labels[state],
                          row.get('sequence_count'), row.get('best_sequence_id'),
                          row.get('best_orientation'))
            evidence_table.append('<tr>' + ''.join(
                '<td>' + escape('Not available' if v is None else str(v)) + '</td>'
                for v in raw_values) + '</tr>')
            fill = _shade(value) if value is not None else 'url(#missing-window)'
            attrs = f'data-state="{state}"'
            if value is not None:
                attrs += f' data-value="{value}"'
            svg.append(f'<rect x="{left + j * cell_width:.3f}" y="{y}" '
                       f'width="{cell_width:.3f}" height="26" fill="{fill}" {attrs}>'
                       f'<title>{escape(details)}</title></rect>')
            cells.append(f'<td title="{escape(details, quote=True)}">{escape(text)}</td>')
        table.append('<tr><th scope="row">' + escape(dataset) + '</th>' + ''.join(cells) + '</tr>')
    headings = ''.join(f'<th scope="col">{a}-{b}</th>' for a, b in matrix.windows)
    return (
        '<figure><h3>Non-target window heatmap</h3>'
        '<p>One column per window, not per nucleotide. Columns follow reference order; '
        'equal column widths do not represent genomic distances. Overlapping windows remain separate.</p>'
        f'<svg viewBox="0 0 1000 {height}" role="img" aria-label="Recorded per-dataset window scores">'
        '<defs><pattern id="missing-window" width="5" height="5" patternUnits="userSpaceOnUse">'
        '<rect width="5" height="5" fill="#eee"/>'
        '<path d="M0 5L5 0" stroke="#aaa" stroke-width="1"/></pattern>'
        '<linearGradient id="window-scale"><stop stop-color="rgb(239,246,255)"/>'
        '<stop offset="1" stop-color="rgb(35,94,154)"/></linearGradient></defs>'
        + ''.join(svg)
        + f'<text x="215" y="{height-39}">First window: {matrix.windows[0][0]}-{matrix.windows[0][1]}</text>'
        + f'<text x="705" y="{height-39}">Last: {matrix.windows[-1][0]}-{matrix.windows[-1][1]}</text>'
        + f'<rect x="215" y="{height-26}" width="160" height="10" fill="url(#window-scale)"/>'
        + f'<text x="215" y="{height-3}">0</text><text x="367" y="{height-3}">1</text>'
        + f'<text x="395" y="{height-15}">Hatched: missing, no match or not evaluated; see table.</text></svg>'
        '<p class="muted">Metric: local alignment score divided by query length. '
        'It is not percent identity or a probability. Source: '
        '<code>contrastive_conservation/contrastive_conservation_report.json: dataset_evidence</code>.</p>'
        '<details><summary>Window values and evidence details</summary>'
        '<div class="scroll" tabindex="0" role="region" aria-label="Per-window numeric evidence">'
        '<table><caption>Recorded scores and distinct unavailable states</caption>'
        '<thead><tr><th scope="col">Dataset</th>' + headings + '</tr></thead><tbody>'
        + ''.join(table) + '</tbody></table></div>'
        '<h4>Exact recorded values and identifiers</h4><div class="scroll" tabindex="0" role="region" aria-label="Cell evidence">'
        '<table><thead><tr><th scope="col">Dataset</th><th scope="col">Window</th>'
        '<th scope="col">Recorded value / state</th><th scope="col">Sequence count</th>'
        '<th scope="col">Best sequence</th><th scope="col">Orientation</th></tr></thead><tbody>'
        + ''.join(evidence_table) + '</tbody></table></div></details></figure>'
    )


def render_curve(rows: list[dict], series: tuple[tuple[str, str], ...]) -> str:
    """Plot recorded fractions at window centres, breaking at missing values."""
    if not rows:
        return '<p class="muted">Window series: Not available.</p>'
    ordered = sorted(((interval(row), row) for row in _rows(rows)), key=lambda item: item[0])
    spans = [span for span, _ in ordered]
    if len(set(spans)) != len(spans):
        raise ResearcherReportError('Duplicate window coordinates in curve.')
    lo, hi = min(a for a, _ in spans), max(b for _, b in spans)
    span = max(hi - lo, 1)
    elements, legends = [], []
    colors = ('#246fa0', '#c76627')
    for index, (key, label) in enumerate(series):
        segments, points = [], []
        for (a, b), row in ordered:
            value = fraction(row.get(key), key)
            if value is None:
                if points:
                    segments.append(points)
                points = []
                continue
            x, y = 48 + (((a + b) / 2 - lo) / span) * 904, 15 + (1 - value) * 180
            points.append((x, y))
        if points:
            segments.append(points)
        color = colors[index % len(colors)]
        for points in segments:
            if len(points) == 1:
                x, y = points[0]
                elements.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3" fill="{color}"/>')
            else:
                coords = ' '.join(f'{x:.3f},{y:.3f}' for x, y in points)
                elements.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2"/>')
        legends.append(f'<span style="border-left:4px solid {color};padding-left:8px">{escape(label)}</span>')
    return (
        '<figure><svg viewBox="0 0 1000 230" role="img" aria-label="Recorded window metrics along the reference">'
        '<path d="M48 15V195H952" fill="none" stroke="#999"/>'
        '<text x="7" y="20">1.0</text><text x="7" y="198">0.0</text>'
        f'<text x="48" y="219">{lo}</text><text x="910" y="219">{hi}</text>'
        + ''.join(elements) + '</svg><figcaption class="legend">' + ' '.join(legends)
        + '</figcaption><p class="muted">Reference coordinates, 1-based inclusive. '
        'Points represent window centres. Missing values break the curves.</p></figure>'
    )
