"""Validate declared populations and preserve published inspection details."""
import unittest
from test_evidence_report import sample_data
from qpcr_pipeline.researcher_report import ResearcherReportError
from qpcr_pipeline.researcher_report_evidence import assay_observations, render_researcher_report_html


class PopulationTests(unittest.TestCase):
    def test_declared_population_rejects_unknown_sequence(self):
        data = sample_data()
        data.inclusivity['evaluation_sequence_ids'] = ['s0', 's1', 's2', 's3']
        for row in data.inclusivity['assay_results']:
            if row['sequence_id'] == 's3':
                row['sequence_id'] = 'unknown'
        with self.assertRaisesRegex(ResearcherReportError, 'population'):
            assay_observations(data)

    def test_duplicate_declared_ids_are_rejected(self):
        data = sample_data()
        data.inclusivity['evaluation_sequence_ids'] = ['s0', 's1', 's2', 's2']
        with self.assertRaises(ResearcherReportError):
            assay_observations(data)

    def test_population_length_and_recorded_denominator_must_agree(self):
        data = sample_data()
        data.inclusivity['evaluation_sequence_ids'] = ['s0', 's1', 's2']
        with self.assertRaises(ResearcherReportError):
            assay_observations(data)

    def test_complete_population_preserves_counts(self):
        data = sample_data()
        data.inclusivity['evaluation_sequence_ids'] = ['s0', 's1', 's2', 's3']
        self.assertEqual(assay_observations(data)[0]['compatible_count'], 3)

    def test_published_detail_is_available_without_recomputation(self):
        data = sample_data()
        data.inclusivity['oligo_matches'] = [{'sequence_id': 'synthetic-detail-0001', 'mismatch_positions': [2]}]
        data.inclusivity['variations'] = [{'note': 'synthetic-variation-detail'}]
        html = render_researcher_report_html(data)
        self.assertIn('synthetic-detail-0001', html)
        self.assertIn('mismatch_positions', html)
        self.assertIn('synthetic-variation-detail', html)


if __name__ == '__main__':
    unittest.main()
