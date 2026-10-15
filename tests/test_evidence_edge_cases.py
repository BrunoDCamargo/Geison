"""Boundary tests for report trustworthiness, using synthetic records only."""
from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest

from test_evidence_report import sample_data
from qpcr_pipeline.researcher_report import ResearcherReportError, generate_researcher_report
from qpcr_pipeline.researcher_report_evidence import assay_observations, render_researcher_report_html
from qpcr_pipeline.report_evidence_charts import render_curve, render_window_heatmap, window_matrix


class EvidenceBoundaryTests(unittest.TestCase):
    def test_skipped_design_never_counts_stale_candidates(self):
        data = sample_data()
        data.primer_design['status'] = 'SKIPPED'
        data.inclusivity['status'] = 'SKIPPED'
        self.assertEqual(assay_observations(data), ())

    def test_conflicting_matrix_and_legacy_counts_are_reported(self):
        data = sample_data()
        data.ranking['assays'][0]['evidence_summary']['original_compatible_count'] = 4
        with self.assertRaisesRegex(ResearcherReportError, 'disagree'):
            assay_observations(data)

    def test_unknown_assay_is_not_silently_added_to_design(self):
        data = sample_data()
        data.inclusivity['assay_results'][0]['assay_id'] = 'unexpected'
        with self.assertRaisesRegex(ResearcherReportError, 'Unknown assay'):
            assay_observations(data)

    def test_incomplete_matrix_does_not_fall_back_to_complete_legacy_counts(self):
        data = sample_data()
        data.inclusivity['assay_results'] = data.inclusivity['assay_results'][:-1]
        rows = assay_observations(data)
        self.assertIsNone(rows[1]['compatible_count'])
        self.assertIn('incomplete', rows[1]['note'])
        self.assertNotIn('5 / 8', render_researcher_report_html(data))

    def test_nonboolean_compatibility_and_count_overflow_are_errors(self):
        data = sample_data()
        data.inclusivity['assay_results'][0]['original_compatible'] = 'true'
        with self.assertRaises(ResearcherReportError):
            assay_observations(data)
        data = sample_data()
        data.inclusivity['counts']['evaluation_sequences'] = 2
        with self.assertRaises(ResearcherReportError):
            assay_observations(data)

    def test_inconsistent_populations_are_rejected(self):
        data = sample_data()
        data.inclusivity['assay_results'][-1]['sequence_id'] = 'different-population'
        with self.assertRaisesRegex(ResearcherReportError, 'population'):
            assay_observations(data)

    def test_duplicate_curves_raise_domain_error(self):
        rows = [{'reference_start': 1, 'reference_end': 100, 'metric': .5},
                {'reference_start': 1, 'reference_end': 100, 'metric': .7}]
        with self.assertRaisesRegex(ResearcherReportError, 'Duplicate'):
            render_curve(rows, (('metric', 'Metric'),))

    def test_invalid_matrix_shapes_raise_domain_error(self):
        for key in ('windows', 'dataset_evidence', 'challenge_datasets'):
            with self.subTest(key=key):
                report = sample_data().contrastive
                report[key] = [None]
                with self.assertRaises(ResearcherReportError):
                    window_matrix(report)

    def test_no_values_are_mutated_and_dataset_text_is_escaped(self):
        report = sample_data().contrastive
        report['challenge_datasets'][0]['name'] = '<script>bad</script>'
        for row in report['dataset_evidence']:
            if row['dataset_name'] == 'X':
                row['dataset_name'] = '<script>bad</script>'
        original = deepcopy(report)
        rendered = render_window_heatmap(report)
        self.assertEqual(report, original)
        self.assertNotIn('<script>', rendered)
        self.assertIn('&lt;script&gt;', rendered)

    def test_empty_dataset_is_unevaluated_not_observed_zero(self):
        report = sample_data().contrastive
        report['dataset_evidence'][0]['sequence_count'] = 0
        matrix = window_matrix(report)
        self.assertEqual(matrix.states[0][0], 'unevaluated')
        self.assertIsNone(matrix.values[0][0])

    def test_numeric_evidence_is_readable_without_hover(self):
        rendered = render_window_heatmap(sample_data().contrastive)
        self.assertIn('<td>0001</td>', rendered)
        self.assertIn('Best sequence', rendered)

    def test_failed_regeneration_preserves_previous_html(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / 'run_manifest.json').write_text('{"status":"COMPLETED"}')
            (root / 'contrastive_conservation').mkdir()
            report = sample_data().contrastive
            report['dataset_evidence'].append(dict(report['dataset_evidence'][0]))
            (root / 'contrastive_conservation/contrastive_conservation_report.json').write_text(json.dumps(report))
            old = root / 'report.html'
            old.write_text('original report bytes')
            with self.assertRaises(ResearcherReportError):
                generate_researcher_report(root)
            self.assertEqual(old.read_text(), 'original report bytes')
            self.assertFalse(list(root.glob('.report.html.*')))


if __name__ == '__main__':
    unittest.main()
