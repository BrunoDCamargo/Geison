"""Synthetic presentation tests; these do not validate biological performance."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import importlib
import json
from pathlib import Path
import tempfile
import unittest

from qpcr_pipeline.researcher_report import ResearcherReportData, ResearcherReportError


def sample_data() -> ResearcherReportData:
    assay_rows = [
        {'assay_id': assay, 'sequence_id': f's{n}', 'original_compatible': n < count}
        for assay, count in (('A', 3), ('B', 2)) for n in range(4)
    ]
    return ResearcherReportData(
        output_dir=Path('/synthetic'),
        run_manifest={'status': 'COMPLETED', 'target_name': 'Synthetic fixture',
                      'run_id': 'fixture-only', 'effective_config': {'synthetic': True},
                      'scientific_completeness': {'missing_evidence': []}},
        run_summary={'sequence_count': 4}, qc_report=None, panel=None,
        conservation={'status': 'COMPLETE', 'counts': {'windows': 2}, 'windows': [
            {'reference_start': 1, 'reference_end': 100, 'mean_conservation': 0.9},
            {'reference_start': 21, 'reference_end': 120, 'mean_conservation': 1.0}]},
        contrastive={'status': 'COMPLETE', 'windows': [
            {'reference_start': 1, 'reference_end': 100, 'target_mean_conservation': 0.9,
             'worst_similarity': 0.75},
            {'reference_start': 21, 'reference_end': 120, 'target_mean_conservation': 1.0,
             'worst_similarity': 0.5}],
            'challenge_datasets': [{'name': 'X'}, {'name': 'Y'}],
            'dataset_evidence': [
                {'dataset_name': 'X', 'reference_start': 1, 'reference_end': 100,
                 'similarity': 0.0, 'sequence_count': 3, 'best_sequence_id': '0001'},
                {'dataset_name': 'X', 'reference_start': 21, 'reference_end': 120,
                 'similarity': 0.5, 'sequence_count': 3, 'best_sequence_id': '0002'},
                {'dataset_name': 'Y', 'reference_start': 1, 'reference_end': 100,
                 'similarity': 0.75, 'sequence_count': 2, 'best_sequence_id': 's'}],
            'candidates': []},
        primer_design={'status': 'COMPLETE', 'assays': [
            {'assay_id': name, 'region_id': 'r', 'product_size': 100, 'pair_penalty': .25,
             'forward_primer': {'sequence': 'ACGT', 'reference_start': 1, 'reference_end': 4},
             'probe': {'sequence': 'TGCA', 'reference_start': 20, 'reference_end': 23},
             'reverse_primer': {'sequence': 'CATG', 'reference_start': 97, 'reference_end': 100}}
            for name in ('B', 'A')],
            'candidates': [{'region_id': 'r', 'reference_start': 1, 'reference_end': 120,
                            'peak_start': 21, 'peak_end': 100}]},
        inclusivity={'status': 'COMPLETE', 'counts': {'assays': 2, 'evaluation_sequences': 4,
                      'original_compatible': 5}, 'assay_results': assay_rows, 'proposals': []},
        specificity={'status': 'COMPLETE', 'counts': {'sequences': 6, 'datasets': 2,
                     'plausible_amplicons': 0, 'detectable_off_targets': 0}, 'retention': []},
        ranking={'status': 'COMPLETE', 'assays': [
            {'assay_id': 'B', 'rank': 1, 'classification': 'HIGH_RISK', 'final_score': 81.03,
             'evidence_summary': {'original_compatible_count': 2, 'evaluation_sequence_count': 4}},
            {'assay_id': 'A', 'rank': 2, 'classification': 'IN SILICO PASS', 'final_score': 79,
             'evidence_summary': {'original_compatible_count': 3, 'evaluation_sequence_count': 4}}]},
    )


class EvidenceReportTests(unittest.TestCase):
    def api(self):
        try:
            return importlib.import_module('qpcr_pipeline.researcher_report_evidence')
        except ModuleNotFoundError:
            self.fail('The neutral evidence renderer is not implemented')

    def test_default_output_does_not_recommend_or_classify(self):
        html = self.api().render_researcher_report_html(sample_data())
        for phrase in ('IN SILICO PASS', 'HIGH_RISK', 'Scientific outcome',
                       'Recommended in-silico', 'Final score', '81.03'):
            self.assertNotIn(phrase, html)
        self.assertIn('COMPLETED', html)
        self.assertIn('5 / 8', html)
        self.assertIn('assay-sequence pairs', html)

    def test_counts_are_per_assay_and_order_is_not_legacy_rank(self):
        rows = self.api().assay_observations(sample_data())
        self.assertEqual([(r['assay_id'], r['compatible_count'], r['evaluation_count'])
                          for r in rows], [('A', 3, 4), ('B', 2, 4)])

    def test_candidates_and_counts_do_not_require_ranking(self):
        data = replace(sample_data(), ranking=None)
        rows = self.api().assay_observations(data)
        self.assertEqual([r['compatible_count'] for r in rows], [3, 2])
        self.assertIn('ACGT', self.api().render_researcher_report_html(data))

    def test_missing_counts_are_not_reconstructed_from_fraction_or_target_size(self):
        data = sample_data()
        data.inclusivity['assay_results'] = []
        for row in data.ranking['assays']:
            row.pop('evidence_summary')
            row['components'] = {'inclusivity': .75}
        rows = self.api().assay_observations(data)
        self.assertTrue(all(r['compatible_count'] is None for r in rows))
        self.assertNotIn('5 / 8', self.api().render_researcher_report_html(data))

    def test_missing_sequence_ids_do_not_inflate_denominator(self):
        data = sample_data()
        data.inclusivity['assay_results'][0].pop('sequence_id')
        with self.assertRaises(ResearcherReportError):
            self.api().assay_observations(data)

    def test_duplicate_matrix_rows_are_rejected(self):
        data = sample_data()
        data.inclusivity['assay_results'].append(dict(data.inclusivity['assay_results'][0]))
        with self.assertRaises(ResearcherReportError):
            self.api().assay_observations(data)

    def test_skipped_stage_does_not_reuse_stale_counts(self):
        data = sample_data()
        data.inclusivity['status'] = 'SKIPPED'
        data.specificity['status'] = 'SKIPPED'
        rows = self.api().assay_observations(data)
        self.assertTrue(all(r['compatible_count'] is None for r in rows))
        html = self.api().render_researcher_report_html(data)
        self.assertIn('SKIPPED', html)
        self.assertNotIn('5 / 8', html)

    def test_partial_and_missing_evidence_are_explicit(self):
        data = sample_data()
        data.run_manifest['status'] = 'PARTIAL'
        data.run_manifest['scientific_completeness']['missing_evidence'] = ['SPECIFICITY']
        data = replace(data, specificity=None)
        html = self.api().render_researcher_report_html(data)
        self.assertIn('PARTIAL', html)
        self.assertIn('SPECIFICITY', html)
        self.assertIn('Not available', html)

    def test_escaping_and_input_immutability(self):
        data = sample_data()
        data.run_manifest['target_name'] = '<img src=x onerror=alert(1)>'
        original = deepcopy(data)
        html = self.api().render_researcher_report_html(data)
        self.assertEqual(data, original)
        self.assertNotIn('<img', html)
        self.assertIn('&lt;img', html)
        self.assertNotIn('<script', html.lower())
        self.assertNotIn('https://', html.lower())

    def test_legacy_summary_counts_have_explicit_source(self):
        data = sample_data()
        data.inclusivity['assay_results'] = []
        rows = self.api().assay_observations(data)
        self.assertEqual(rows[0]['compatible_count'], 3)
        self.assertIn('ranking', rows[0]['source'])

    def test_invalid_legacy_counts_are_not_clamped(self):
        for count in (-1, 8, True, float('nan')):
            with self.subTest(count=count):
                data = sample_data()
                data.inclusivity['assay_results'] = []
                data.ranking['assays'][0]['evidence_summary']['original_compatible_count'] = count
                with self.assertRaises(ResearcherReportError):
                    self.api().assay_observations(data)

    def test_regeneration_uses_neutral_renderer_and_preserves_json(self):
        from qpcr_pipeline.researcher_report import generate_researcher_report
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            paths = {'run_manifest.json': {'status': 'COMPLETED'},
                     'primer_design/primer_design_report.json': sample_data().primer_design,
                     'inclusivity/inclusivity_report.json': sample_data().inclusivity}
            for name, payload in paths.items():
                p = root / name
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(json.dumps(payload), encoding='utf-8')
            before = {name: (root / name).read_bytes() for name in paths}
            try:
                result = generate_researcher_report(root)
            except ModuleNotFoundError:
                self.fail('The default generator has not been wired to the neutral renderer')
            html = result.read_text()
            self.assertIn('Evidence summary', html)
            self.assertNotIn('Scientific outcome', html)
            self.assertEqual(before, {name: (root / name).read_bytes() for name in paths})


class WindowMatrixTests(unittest.TestCase):
    def api(self):
        try:
            return importlib.import_module('qpcr_pipeline.report_evidence_charts')
        except ModuleNotFoundError:
            self.fail('The numeric evidence heatmap is not implemented')

    def test_matrix_uses_each_dataset_not_global_maximum(self):
        matrix = self.api().window_matrix(sample_data().contrastive)
        self.assertEqual(matrix.windows, ((1, 100), (21, 120)))
        self.assertEqual(matrix.datasets, ('X', 'Y'))
        self.assertEqual(matrix.values, ((0.0, 0.5), (0.75, None)))

    def test_missing_and_numeric_zero_have_distinct_cells(self):
        html = self.api().render_window_heatmap(sample_data().contrastive)
        self.assertIn('data-value="0.0"', html)
        self.assertIn('data-state="missing"', html)
        self.assertIn('0001', html)
        self.assertIn('one column per window', html.lower())

    def test_no_best_match_is_not_zero_similarity(self):
        report = sample_data().contrastive
        report['dataset_evidence'][0]['best_sequence_id'] = None
        matrix = self.api().window_matrix(report)
        self.assertIsNone(matrix.values[0][0])
        self.assertIn('no-match', self.api().render_window_heatmap(report))

    def test_duplicate_cells_are_rejected(self):
        report = sample_data().contrastive
        report['dataset_evidence'].append(dict(report['dataset_evidence'][0]))
        with self.assertRaises(ResearcherReportError):
            self.api().window_matrix(report)

    def test_nonfinite_out_of_range_and_boolean_values_are_rejected(self):
        for value in (float('nan'), float('inf'), -0.01, 1.01, True):
            with self.subTest(value=value):
                report = sample_data().contrastive
                report['dataset_evidence'][0]['similarity'] = value
                with self.assertRaises(ResearcherReportError):
                    self.api().window_matrix(report)

    def test_empty_skipped_and_single_window(self):
        self.assertEqual(self.api().window_matrix(None).values, ())
        report = sample_data().contrastive
        report['status'] = 'SKIPPED'
        self.assertEqual(self.api().window_matrix(report).values, ())
        report = sample_data().contrastive
        report['windows'] = report['windows'][:1]
        report['dataset_evidence'] = [r for r in report['dataset_evidence']
                                      if r['reference_start'] == 1]
        self.assertEqual(self.api().window_matrix(report).values, ((0.0,), (0.75,)))

    def test_cells_outside_declared_windows_are_rejected(self):
        report = sample_data().contrastive
        report['dataset_evidence'][0]['reference_end'] = 99
        with self.assertRaises(ResearcherReportError):
            self.api().window_matrix(report)

    def test_missing_curve_values_break_segments(self):
        rows = [dict(reference_start=i, reference_end=i, metric=v)
                for i, v in ((1, .2), (2, .4), (3, None), (4, .5), (5, .6))]
        html = self.api().render_curve(rows, (('metric', 'Recorded metric'),))
        self.assertEqual(html.count('<polyline'), 2)


if __name__ == '__main__':
    unittest.main()
