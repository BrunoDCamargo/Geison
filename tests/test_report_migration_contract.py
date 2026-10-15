"""Regression tests for descriptive output and lossless artifact loading."""
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest

from test_evidence_report import sample_data
from qpcr_pipeline.researcher_report import load_researcher_report_data
from qpcr_pipeline.researcher_report_evidence import render_researcher_report_html


class ReportMigrationContractTests(unittest.TestCase):
    def test_completed_empty_design_reports_zero_not_missing(self):
        data = replace(sample_data(), primer_design={
            'status': 'COMPLETE', 'assays': [], 'candidates': []},
            inclusivity=None, ranking=None)
        summary = render_researcher_report_html(data).split('</section>', 1)[0]
        self.assertIn('<td>0</td>', summary)

    def test_missing_design_does_not_invent_zero_candidates(self):
        data = replace(sample_data(), primer_design=None, inclusivity=None, ranking=None)
        summary = render_researcher_report_html(data).split('</section>', 1)[0]
        self.assertNotIn('<td>0</td>', summary)
        self.assertIn('Not available', summary)

    def test_retained_correspondences_remain_inspectable(self):
        hit = {'dataset_name': 'reference-set', 'assay_id': 'A',
               'sequence_id': 'retained-site-001', 'role': 'FORWARD',
               'orientation': 'FORWARD', 'source_start': 10, 'source_end': 29,
               'mismatch_positions': '2,5', 'mismatch_count': 2,
               'three_prime_mismatch': False, 'compatible': True}
        data = replace(sample_data(), specificity_hits=(hit,))
        html = render_researcher_report_html(data)
        self.assertIn('retained-site-001', html)
        self.assertIn('Mismatch positions', html)
        self.assertIn('2,5', html)
        data.specificity['status'] = 'SKIPPED'
        self.assertNotIn('retained-site-001', render_researcher_report_html(data))

    def test_skipped_ranking_does_not_publish_stale_rule_evidence(self):
        data = sample_data()
        data.ranking['status'] = 'SKIPPED'
        data.ranking['assays'][0]['reasons'] = [
            {'code': 'STALE_RULE_SENTINEL', 'source': 'old-run', 'evidence': {}}]
        self.assertNotIn('STALE_RULE_SENTINEL', render_researcher_report_html(data))

    def test_tsv_identifiers_and_position_lists_keep_their_original_text(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'run_manifest.json').write_text(json.dumps({'status': 'COMPLETED'}))
            directory = root / 'specificity'
            directory.mkdir()
            path = directory / 'off_target_hits.tsv'
            raw = ('dataset_name\tassay_id\tsequence_id\tmismatch_positions\tsource_start\tcompatible\n'
                   'false\t001\t0007\t2\t10\ttrue\n')
            path.write_text(raw)
            row, = load_researcher_report_data(root).specificity_hits
            self.assertEqual(row['dataset_name'], 'false')
            self.assertEqual(row['assay_id'], '001')
            self.assertEqual(row['sequence_id'], '0007')
            self.assertEqual(row['mismatch_positions'], '2')
            self.assertEqual(row['source_start'], 10)
            self.assertIs(row['compatible'], True)
            self.assertEqual(path.read_text(), raw)
