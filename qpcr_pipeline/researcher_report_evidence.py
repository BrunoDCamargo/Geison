"""Read-only evidence report. No assay selection or scientific classification."""
from __future__ import annotations

from html import escape
import json
from typing import Iterable

from qpcr_pipeline.researcher_report import ResearcherReportData, ResearcherReportError
from qpcr_pipeline.report_evidence_charts import render_curve, render_window_heatmap


def _obj(value: object) -> dict:
    return value if isinstance(value, dict) else {}


def _rows(value: object) -> list[dict]:
    if value is None:
        return []
    if not isinstance(value, (list, tuple)) or any(not isinstance(x, dict) for x in value):
        raise ResearcherReportError('Evidence rows must be a sequence of objects.')
    return list(value)


def _identifier(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ResearcherReportError(f'{name} must be a nonblank string.')
    return value


def _count(value: object, name: str) -> int | None:
    if value is None:
        return None
    if type(value) is not int or value < 0:
        raise ResearcherReportError(f'{name} must be a nonnegative integer.')
    return value


def _indexed(rows: list[dict], key: str) -> dict[str, dict]:
    indexed = {}
    for row in rows:
        identifier = _identifier(row.get(key), key)
        if identifier in indexed:
            raise ResearcherReportError(f'Duplicate {key}: {identifier}')
        indexed[identifier] = row
    return indexed


def _complete(report: dict | None) -> bool:
    return _obj(report).get('status') == 'COMPLETE'


def assay_observations(data: ResearcherReportData) -> tuple[dict, ...]:
    """Use recorded matrix counts, or explicit legacy counts, never rounded scores."""
    if data.primer_design is not None and not _complete(data.primer_design):
        return ()
    designs = _indexed(_rows(_obj(data.primer_design).get('assays')), 'assay_id')
    legacy = _indexed(_rows(_obj(data.ranking).get('assays')), 'assay_id') if _complete(data.ranking) else {}
    report = _obj(data.inclusivity)
    active = _complete(report)
    matrix = _rows(report.get('assay_results')) if active else []
    groups: dict[str, list[dict]] = {}
    seen = set()
    for row in matrix:
        assay = _identifier(row.get('assay_id'), 'assay_id')
        if data.primer_design is not None and assay not in designs:
            raise ResearcherReportError(f'Unknown assay in compatibility matrix: {assay}')
        seq = _identifier(row.get('sequence_id'), 'sequence_id')
        if (assay, seq) in seen:
            raise ResearcherReportError('Duplicate assay/sequence observation.')
        seen.add((assay, seq))
        if type(row.get('original_compatible')) is not bool:
            raise ResearcherReportError('original_compatible must be boolean.')
        groups.setdefault(assay, []).append(row)
    expected = _count(_obj(report.get('counts')).get('evaluation_sequences'), 'evaluation_sequences') if active else None
    if active and 'evaluation_sequence_ids' in report:
        declared = report['evaluation_sequence_ids']
        if not isinstance(declared, (list, tuple)):
            raise ResearcherReportError('Declared population must be a sequence of identifiers.')
        declared_ids = [_identifier(v, 'population sequence_id') for v in declared]
        if len(set(declared_ids)) != len(declared_ids):
            raise ResearcherReportError('Duplicate sequence in declared population.')
        if expected is not None and len(declared_ids) != expected:
            raise ResearcherReportError('Declared population and denominator disagree.')
        if any(seq not in declared_ids for _, seq in seen):
            raise ResearcherReportError('Observation outside the declared population.')
    populations = [{r['sequence_id'] for r in rows} for rows in groups.values()
                   if expected is not None and len(rows) == expected]
    if populations and any(ids != populations[0] for ids in populations[1:]):
        raise ResearcherReportError('Complete matrix rows refer to different populations.')
    if data.primer_design is not None and set(legacy) - set(designs):
        raise ResearcherReportError('Unknown assay in legacy evidence summary.')
    result = []
    for assay in sorted(set(designs) | set(legacy) | set(groups)):
        n, denominator, source, note = None, None, 'Not available', 'Not evaluated or evidence unavailable.'
        observations = groups.get(assay, [])
        if observations:
            if expected is not None and len(observations) > expected:
                raise ResearcherReportError('Observed rows exceed the recorded evaluation count.')
            source = 'inclusivity/inclusivity_report.json: assay_results'
            if expected is not None and len(observations) == expected:
                n, denominator = sum(r['original_compatible'] for r in observations), expected
                summary = _obj(legacy.get(assay, {}).get('evidence_summary'))
                legacy_n = _count(summary.get('original_compatible_count'), 'legacy compatible count')
                legacy_d = _count(summary.get('evaluation_sequence_count'), 'legacy evaluation count')
                if ((legacy_n is not None and legacy_n != n)
                        or (legacy_d is not None and legacy_d != denominator)):
                    raise ResearcherReportError('Matrix and legacy counts disagree.')
                note = 'Original assay; compatibility under the recorded rules.'
            else:
                note = f'{len(observations)} recorded rows; full evaluation denominator unavailable or matrix incomplete.'
        elif active and not matrix:
            summary = _obj(legacy.get(assay, {}).get('evidence_summary'))
            n = _count(summary.get('original_compatible_count'), 'original_compatible_count')
            denominator = _count(summary.get('evaluation_sequence_count'), 'evaluation_sequence_count')
            if n is not None and denominator is not None and n > denominator:
                raise ResearcherReportError('Compatible count exceeds evaluation count.')
            if denominator is not None and expected is not None and denominator != expected:
                raise ResearcherReportError('Legacy and inclusivity denominators disagree.')
            if n is not None and denominator is not None:
                source = 'ranking/ranking_report.json: evidence_summary (legacy counts)'
                note = 'Recorded counts only; underlying assay/sequence matrix unavailable.'
            else:
                n, denominator = None, None
        result.append(dict(assay_id=assay, compatible_count=n, evaluation_count=denominator,
                           source=source, note=note))
    return tuple(result)


def _text(value: object) -> str:
    if value is None or value == '':
        return 'Not available'
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    if isinstance(value, (list, dict, tuple)):
        return escape(json.dumps(value, ensure_ascii=False, sort_keys=True), quote=True)
    return escape(str(value), quote=True)


def _table(headers: tuple[str, ...], rows: Iterable[tuple]) -> str:
    body = ''.join('<tr>' + ''.join('<td>' + _text(cell) + '</td>' for cell in row) + '</tr>' for row in rows)
    if not body:
        body = '<tr><td colspan="' + str(len(headers)) + '">No rows recorded.</td></tr>'
    return ('<div class="scroll" tabindex="0" role="region" aria-label="Evidence table"><table><thead><tr>'
            + ''.join('<th scope="col">' + escape(h) + '</th>' for h in headers)
            + '</tr></thead><tbody>' + body + '</tbody></table></div>')


def _preview(value: object) -> str:
    if value is None:
        return '<p>Not available.</p>'
    return '<pre>' + escape(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)) + '</pre>'


def _status(report: dict | None) -> str:
    return '<p><strong>Recorded stage status:</strong> ' + _text(_obj(report).get('status')) + '</p>'


def _span(row: dict, start: str = 'reference_start', end: str = 'reference_end') -> str | None:
    if row.get(start) is None or row.get(end) is None:
        return None
    return f'{row[start]}-{row[end]}'


def _coverage(data: ResearcherReportData, observations: tuple[dict, ...]) -> str:
    aggregate = ''
    if observations and all(r['compatible_count'] is not None and r['evaluation_count'] is not None for r in observations):
        total = sum(r['compatible_count'] for r in observations)
        denominator = sum(r['evaluation_count'] for r in observations)
        aggregate = (f'<p><strong>{total} / {denominator}</strong> compatible assay-sequence pairs. '
                     'This sum spans assays; it is not the number of unique target sequences.</p>')
    table = _table(('Assay', 'Compatible / evaluated sequences', 'Evidence source', 'Scope'), (
        (r['assay_id'], None if r['compatible_count'] is None else
         f"{r['compatible_count']} / {r['evaluation_count']}", r['source'], r['note']) for r in observations))
    proposals = ''
    if _complete(data.inclusivity):
        proposals = '<h3>Recorded proposals</h3><p>Policy outcomes below are recorded software rules, not researcher approval. Original and proposed sequences remain separate.</p>' + _table(
            ('Assay', 'Role', 'Original', 'Proposed', 'Policy state', 'Rule code'), (
                tuple(r.get(k) for k in ('assay_id', 'role', 'original_sequence', 'proposed_sequence', 'status', 'reason'))
                for r in _rows(_obj(data.inclusivity).get('proposals'))))
        for label, key in (('Assay-sequence records', 'assay_results'),
                           ('Oligo observations', 'oligo_matches'),
                           ('Recorded variations', 'variations')):
            proposals += '<details><summary>' + label + '</summary>' + _preview(_obj(data.inclusivity).get(key)) + '</details>'
    return '<section><h2>Target compatibility by assay</h2>' + _status(data.inclusivity) + aggregate + table + proposals + '</section>'


def _design(data: ResearcherReportData) -> str:
    report = _obj(data.primer_design)
    content = _status(data.primer_design)
    if not _complete(report):
        return '<section><h2>Recorded assay candidates</h2>' + content + '<p>Not available for a complete design stage.</p></section>'
    regions = _indexed(_rows(report.get('candidates')), 'region_id')
    for assay in sorted(_rows(report.get('assays')), key=lambda r: r['assay_id']):
        region = regions.get(assay.get('region_id'), {})
        content += '<article><h3>' + _text(assay.get('assay_id')) + '</h3>' + _table(
            ('Candidate region', 'Representative window', 'Product size', 'Primer3 pair penalty'), [
                (_span(region), _span(region, 'peak_start', 'peak_end'), assay.get('product_size'), assay.get('pair_penalty'))])
        content += _table(('Role', 'Recorded sequence', 'Reference coordinates', 'Tm', 'GC %', 'Primer3 penalty'), (
            (role, oligo.get('sequence'), _span(oligo), oligo.get('tm'), oligo.get('gc_percent'), oligo.get('penalty'))
            for role, oligo in ((name, _obj(assay.get(key))) for name, key in (
                ('Forward', 'forward_primer'), ('Probe', 'probe'), ('Reverse', 'reverse_primer'))))) + '</article>'
    return '<section><h2>Recorded assay candidates</h2><p>Identifier order, without recommendation. Primer3 penalties retain their original units.</p>' + content + '</section>'


def _conservation(data: ResearcherReportData) -> str:
    report = _obj(data.conservation)
    chart = ''
    if _complete(report):
        rows = _rows(report.get('windows')) or list(data.conservation_windows)
        chart = render_curve(rows, (('mean_conservation', 'Target conservation'),))
    return '<section><h2>Target conservation</h2>' + _status(data.conservation) + chart + '</section>'


def _contrast(data: ResearcherReportData) -> str:
    report = _obj(data.contrastive)
    content = _status(data.contrastive)
    if _complete(report):
        content += render_curve(_rows(report.get('windows')), (
            ('target_mean_conservation', 'Target conservation'),
            ('worst_similarity', 'Maximum non-target normalized local alignment score')))
        content += render_window_heatmap(report)
        content += '<h3>Recorded regions</h3><p>The window value below does not describe the entire expanded region. No ranking or discrimination claim is made here.</p>'
        content += _table(('Region', 'Expanded interval', 'Representative window', 'Contributing windows', 'Dataset at recorded maximum', 'Window score'), (
            (_obj(r.get('region')).get('region_id'), _span(_obj(r.get('region'))),
             _span(_obj(r.get('region')), 'peak_start', 'peak_end'), r.get('contributing_windows'),
             r.get('worst_dataset_name'), r.get('worst_similarity')) for r in _rows(report.get('candidates'))))
    return '<section><h2>Target and non-target window evidence</h2>' + content + '</section>'


def _specificity(data: ResearcherReportData) -> str:
    report = _obj(data.specificity)
    content = _status(data.specificity)
    if _complete(report):
        counts = _obj(report.get('counts'))
        content += _table(('Recorded datasets', 'Recorded sequences', 'F/R geometries', 'Geometries with probe correspondence'), [
            tuple(counts.get(k) for k in ('datasets', 'sequences', 'plausible_amplicons', 'detectable_off_targets'))])
        content += '<details><summary>Hit retention and recorded geometries</summary>' + _table(
            ('Dataset', 'Assay', 'Role', 'Total hits', 'Retained', 'Truncated'), (
                tuple(r.get(k) for k in ('dataset_name', 'assay_id', 'role', 'total_hit_count', 'retained_hit_count', 'truncated'))
                for r in _rows(report.get('retention'))))
        content += '<h3>Retained oligo correspondences</h3>' + _table(
            ('Dataset', 'Assay', 'Sequence', 'Role', 'Orientation', 'Coordinates',
             'Mismatch positions', 'Mismatch count', '3-prime mismatch', 'Compatible under recorded rules'), (
                (r.get('dataset_name'), r.get('assay_id'), r.get('sequence_id'), r.get('role'),
                 r.get('orientation'), _span(r, 'source_start', 'source_end'),
                 r.get('mismatch_positions'), r.get('mismatch_count'),
                 r.get('three_prime_mismatch'), r.get('compatible'))
                for r in data.specificity_hits))
        content += _table(('Dataset', 'Assay', 'Sequence', 'Coordinates', 'Size'), (
            (r.get('dataset_name'), r.get('assay_id'), r.get('sequence_id'), _span(r, 'source_start', 'source_end'), r.get('amplicon_size'))
            for r in data.specificity_amplicons)) + '</details>'
    else:
        content += '<p>Not available. Absence of evidence is not a zero result.</p>'
    return ('<section><h2>Recorded non-target correspondences</h2><p>Counts describe the recorded search scope and filters. '
            'Geometries are computational observations, not observed amplification. Zero is not proof of absence.</p>' + content + '</section>')


def render_researcher_report_html(data: ResearcherReportData) -> str:
    """Render published data without changing measurements, candidates or files."""
    observations = assay_observations(data)
    manifest = data.run_manifest
    target = manifest.get('target_name') or _obj(_obj(_obj(data.panel).get('definition')).get('target')).get('name')
    design_rows = _obj(data.primer_design).get('assays')
    candidate_count = (len(design_rows) if _complete(data.primer_design)
                       and isinstance(design_rows, (list, tuple))
                       else len(observations) if observations else None)
    summary = '<section class="summary"><h2>Evidence summary</h2>' + _table(
        ('Technical run status', 'Recorded target sequences', 'Candidate identifiers'), [
            (manifest.get('status'), _obj(data.run_summary).get('sequence_count'), candidate_count)])
    summary += '<p>No scientific acceptance or experimental performance is inferred from the execution status.</p>'
    missing = _obj(manifest.get('scientific_completeness')).get('missing_evidence')
    if missing:
        summary += '<p><strong>Missing evidence:</strong> ' + _text(missing) + '</p>'
    summary += '</section>'
    panel = '<section><h2>Recorded panel and context</h2><p>Panel approval records a human input decision. Dataset roles alone do not establish independence from development.</p>'
    panel += '<details><summary>Panel records</summary>' + (_preview(data.panel) if data.panel is not None else '<p>Not available.</p>') + '</details></section>'
    audit = '<section><h2>Evidence audit</h2><p>Legacy ranking artifacts remain unchanged. Their classifications and composite scores are not used to select, order or approve candidates in this report.</p>'
    audit += '<details><summary>Recorded rule codes and supporting values</summary>' + _table(
        ('Assay', 'Code', 'Source', 'Recorded evidence'), (
            (assay.get('assay_id'), r.get('code'), r.get('source'), r.get('evidence'))
            for assay in (_rows(_obj(data.ranking).get('assays')) if _complete(data.ranking) else [])
            for r in _rows(assay.get('reasons')))) + '</details></section>'
    provenance = '<section><h2>Reproducibility</h2>' + _table(('Run ID', 'Created', 'Updated', 'Reference'), [
        (manifest.get('run_id'), manifest.get('created_at'), manifest.get('updated_at'), _obj(manifest.get('reference')).get('id'))])
    for title, key in (('Effective configuration', 'effective_config'), ('Environment and tool information', 'environment')):
        provenance += '<details><summary>' + title + '</summary>' + _preview(manifest.get(key)) + '</details>'
    provenance += '<p>Source: published run JSON/TSV artifacts. Missing metadata are not reconstructed.</p></section>'
    body = summary + panel + _conservation(data) + _contrast(data) + _design(data) + _coverage(data, observations) + _specificity(data) + audit + provenance
    return '''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><link rel="icon" href="data:,">
<title>Geison Evidence Report</title><style>
*{box-sizing:border-box}body{margin:0;color:#203144;font:16px/1.5 system-ui,sans-serif;background:#fff}
main{max-width:1160px;margin:auto;padding:32px 28px 60px}h1{font-size:2rem;margin-bottom:8px}
h2{font-size:1.4rem}h3{font-size:1.05rem}section{margin:32px 0}header{border-bottom:1px solid #cbd5df;padding-bottom:20px}
.summary,figure,article{border:1px solid #d8e2ed;border-radius:8px;padding:16px;margin:16px 0}.summary{background:#f3f7fa}
figure svg{width:100%;height:auto}svg text{font:12px system-ui,sans-serif;fill:#304454}table{border-collapse:collapse;width:100%;font-size:.86rem}
th,td{padding:9px;text-align:left;border-bottom:1px solid #dbe3eb;vertical-align:top}th{background:#f5f8fa}
.scroll{overflow:auto;max-height:480px;border:1px solid #dbe3eb}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f5f8fa;padding:12px;font-size:.8rem}
code{overflow-wrap:anywhere}
.muted{font-size:.86rem;color:#4c6275}.legend{display:flex;flex-wrap:wrap;gap:18px;font-size:.86rem}
details{border:1px solid #dbe3eb;border-radius:6px;padding:10px;margin:12px 0}summary{cursor:pointer;font-weight:600}
:focus-visible{outline:3px solid #267ba5;outline-offset:3px}footer{border-top:1px solid #dbe3eb;padding-top:16px;color:#4c6275}
@media(max-width:700px){main{padding:20px 12px}figure{padding:8px}table{font-size:.8rem}}
@media print{.scroll{max-height:none}details{display:block}section{break-inside:auto}}
</style></head><body><main><header><h1>Geison Evidence Report</h1><p><strong>Target:</strong> ''' + _text(target) + '''</p>
<p>Recorded measurements, scope and parameters for researcher inspection.</p></header>''' + body + '''
<footer>In-silico records only. This presentation change does not validate the scientific algorithms or experimental performance.
Review the underlying artifacts and their limitations before interpreting the results.</footer></main></body></html>'''
