from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from qpcr_pipeline.researcher_report import load_researcher_report_data
from qpcr_pipeline.researcher_report_html import render_researcher_report_html


class ResearcherReportRealTableHtmlTests(unittest.TestCase):
    def _write(self, root: Path, relative: str, text: str) -> None:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def test_report_renders_conservation_windows_and_specificity_amplicons_from_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write(
                root,
                "run_manifest.json",
                json.dumps(
                    {
                        "status": "COMPLETED",
                        "target_name": "Synthetic target",
                        "scientific_completeness": {"complete": True, "missing_evidence": []},
                    }
                ),
            )
            self._write(root, "run_summary.json", json.dumps({"sequence_count": 4}))
            self._write(
                root,
                "conservation/conservation_report.json",
                json.dumps(
                    {
                        "status": "COMPLETE",
                        "counts": {"sequences": 4, "reference_positions": 1200, "windows": 2},
                    }
                ),
            )
            self._write(
                root,
                "conservation/window_metrics.tsv",
                "reference_start\treference_end\tposition_count\tmean_conservation\t"
                "minimum_conservation\tmean_coverage\tmean_gap_frequency\tmean_entropy_bits\n"
                "581\t680\t100\t0.98\t0.95\t1.0\t0.0\t0.01\n"
                "601\t700\t100\t1.0\t1.0\t1.0\t0.0\t0.0\n",
            )
            self._write(
                root,
                "specificity/specificity_report.json",
                json.dumps(
                    {
                        "status": "COMPLETE",
                        "counts": {
                            "datasets": 1,
                            "sequences": 1,
                            "assays": 1,
                            "plausible_amplicons": 1,
                            "detectable_off_targets": 1,
                        },
                        "retention": [],
                    }
                ),
            )
            self._write(
                root,
                "specificity/plausible_amplicons.tsv",
                "dataset_name\tassay_id\tsequence_id\torientation\tsource_start\tsource_end\t"
                "amplicon_size\tforward_source_start\tforward_source_end\treverse_source_start\t"
                "reverse_source_end\tprobe_source_sites\tforward_hit_rank\treverse_hit_rank\t"
                "probe_hit_ranks\tprimer_amplicon_plausible\tdetectable_off_target\n"
                "Related synthetic A\ta1\toff-1\tFORWARD\t724\t795\t72\t724\t743\t775\t795\t"
                "744-769\t1\t1\t1\ttrue\ttrue\n",
            )
            self._write(
                root,
                "ranking/ranking_report.json",
                json.dumps(
                    {
                        "status": "COMPLETE",
                        "counts": {"assays": 1, "in_silico_pass": 0, "review": 0, "high_risk": 1},
                        "assays": [
                            {
                                "rank": 1,
                                "assay_id": "a1",
                                "region_id": "r1",
                                "classification": "HIGH_RISK",
                                "score_status": "COMPLETE",
                                "final_score": 69.4,
                                "components": {},
                                "reasons": [
                                    {
                                        "code": "DETECTABLE_OFF_TARGET",
                                        "severity": "HIGH_RISK",
                                        "source": "specificity",
                                        "message": "Detected.",
                                    }
                                ],
                            }
                        ],
                    }
                ),
            )

            html = render_researcher_report_html(load_researcher_report_data(root))

        self.assertIn("Recorded target conservation by reference position", html)
        self.assertNotIn("No recorded window series available for this view.", html)
        self.assertIn("Plausible off-target amplicons", html)
        self.assertIn("Related synthetic A", html)
        self.assertIn("724–795", html)
        self.assertIn("Detectable off-target", html)
        self.assertIn(">Yes<", html)


if __name__ == "__main__":
    unittest.main()
