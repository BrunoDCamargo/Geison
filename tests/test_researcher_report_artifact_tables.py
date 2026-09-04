from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from qpcr_pipeline.researcher_report import load_researcher_report_data


class ResearcherReportArtifactTableTests(unittest.TestCase):
    def _write(self, root: Path, relative: str, text: str) -> None:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def test_loads_conservation_window_metrics_from_published_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write(root, "run_manifest.json", json.dumps({"status": "COMPLETED"}))
            self._write(
                root,
                "conservation/window_metrics.tsv",
                "reference_start\treference_end\tposition_count\tmean_conservation\t"
                "minimum_conservation\tmean_coverage\tmean_gap_frequency\tmean_entropy_bits\n"
                "601\t700\t100\t1.0\t1.0\t1.0\t0.0\t0.0\n",
            )

            data = load_researcher_report_data(root)

        self.assertEqual(len(data.conservation_windows), 1)
        row = data.conservation_windows[0]
        self.assertEqual(row["reference_start"], 601)
        self.assertEqual(row["reference_end"], 700)
        self.assertEqual(row["position_count"], 100)
        self.assertEqual(row["mean_conservation"], 1.0)
        self.assertEqual(row["mean_gap_frequency"], 0.0)

    def test_loads_specificity_amplicons_and_hits_from_published_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write(root, "run_manifest.json", json.dumps({"status": "COMPLETED"}))
            self._write(
                root,
                "specificity/plausible_amplicons.tsv",
                "dataset_name\tassay_id\tsequence_id\torientation\tsource_start\tsource_end\t"
                "amplicon_size\tforward_source_start\tforward_source_end\treverse_source_start\t"
                "reverse_source_end\tprobe_source_sites\tforward_hit_rank\treverse_hit_rank\t"
                "probe_hit_ranks\tprimer_amplicon_plausible\tdetectable_off_target\n"
                "Related A\ta1\toff-1\tFORWARD\t724\t795\t72\t724\t743\t775\t795\t"
                "744-769\t1\t1\t1\ttrue\ttrue\n",
            )
            self._write(
                root,
                "specificity/off_target_hits.tsv",
                "dataset_name\tassay_id\tsequence_id\trole\torientation\thit_rank\t"
                "source_start\tsource_end\tmismatch_positions\tmismatch_count\texact_match\t"
                "three_prime_mismatch\tcompatible\n"
                "Related A\ta1\toff-1\tFORWARD\tFORWARD\t1\t724\t743\t\t0\ttrue\tfalse\ttrue\n",
            )

            data = load_researcher_report_data(root)

        self.assertEqual(len(data.specificity_amplicons), 1)
        amplicon = data.specificity_amplicons[0]
        self.assertEqual(amplicon["amplicon_size"], 72)
        self.assertIs(amplicon["primer_amplicon_plausible"], True)
        self.assertIs(amplicon["detectable_off_target"], True)
        self.assertEqual(len(data.specificity_hits), 1)
        hit = data.specificity_hits[0]
        self.assertEqual(hit["mismatch_count"], 0)
        self.assertIs(hit["compatible"], True)
        self.assertIs(hit["three_prime_mismatch"], False)


if __name__ == "__main__":
    unittest.main()
