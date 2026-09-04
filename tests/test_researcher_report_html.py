from __future__ import annotations

from pathlib import Path
import unittest

from qpcr_pipeline.researcher_report import ResearcherReportData
from qpcr_pipeline.researcher_report_html import render_researcher_report_html


def report_data(
    *,
    run_status: str = "COMPLETED",
    classification: str = "IN SILICO PASS",
    target_name: str = "Synthetic target",
) -> ResearcherReportData:
    reasons = []
    specificity_score = 1.0
    detectable_count = 0
    if classification == "HIGH_RISK":
        reasons = [
            {
                "code": "DETECTABLE_OFF_TARGET",
                "severity": "HIGH_RISK",
                "source": "specificity",
                "message": "A detectable off-target amplicon was found.",
                "evidence": {"detectable_off_target_count": 1},
            }
        ]
        specificity_score = 0.0
        detectable_count = 1
    elif classification == "REVIEW":
        reasons = [
            {
                "code": "PLAUSIBLE_OFF_TARGET_AMPLICON",
                "severity": "REVIEW",
                "source": "specificity",
                "message": "A plausible primer amplicon requires review.",
                "evidence": {"plausible_off_target_count": 1},
            }
        ]
        specificity_score = 0.4

    ranking_assay = {
        "rank": 1,
        "assay_id": "contrast-region-001-assay-001",
        "region_id": "contrast-region-001",
        "classification": classification,
        "score_status": "COMPLETE",
        "final_score": 92.5 if classification == "IN SILICO PASS" else 69.4,
        "components": {
            "inclusivity": 1.0,
            "specificity": specificity_score,
            "conservation": 1.0,
            "primer3_quality": 0.44,
            "robustness": 1.0,
        },
        "evidence_summary": {
            "original_compatible_count": 4,
            "evaluation_sequence_count": 4,
            "compatible_off_target_hit_count": 2,
            "plausible_off_target_count": detectable_count,
            "detectable_off_target_count": detectable_count,
            "pair_penalty": 1.25,
        },
        "reasons": reasons,
    }
    counts = {
        "assays": 1,
        "in_silico_pass": 1 if classification == "IN SILICO PASS" else 0,
        "review": 1 if classification == "REVIEW" else 0,
        "high_risk": 1 if classification == "HIGH_RISK" else 0,
    }

    return ResearcherReportData(
        output_dir=Path("/tmp/output"),
        run_manifest={
            "status": run_status,
            "target_name": target_name,
            "run_id": "run-123",
            "created_at": "2026-09-04T19:00:00Z",
            "updated_at": "2026-09-04T19:05:00Z",
            "effective_config": {"target_name": target_name},
            "environment": {
                "python": {"version": "3.12"},
                "tools": {"primer3_core": {"version": "2.6"}},
            },
            "scientific_completeness": {
                "complete": run_status == "COMPLETED",
                "missing_evidence": [] if run_status == "COMPLETED" else ["NO_ASSAYS"],
            },
            "panel_provenance": {"manifest_sha256": "sha256:abc"},
            "reference": {"id": "synthetic-target-reference", "mode": "explicit"},
        },
        run_summary={"status": run_status, "sequence_count": 4},
        qc_report={
            "evaluation_set": {
                "sequence_ids": ["ref", "variant-1", "variant-2", "variant-3"]
            }
        },
        panel={
            "status": "APPROVED",
            "proposal_sha256": "sha256:panel",
            "definition": {
                "target": {
                    "name": target_name,
                    "mode": "broad_detection",
                    "groups": [
                        {
                            "name": "synthetic target diversity",
                            "required": True,
                            "dataset_roles": ["DESIGN"],
                            "reasons": ["synthetic demonstration"],
                        }
                    ],
                },
                "non_targets": [
                    {
                        "name": "Related synthetic A",
                        "criticality": "CRITICAL",
                        "dataset_roles": ["CHALLENGE"],
                        "reasons": ["close differential"],
                    }
                ],
                "diagnostic_context": {
                    "syndrome": "synthetic febrile scenario",
                    "geography": "synthetic setting",
                    "sample_type": "synthetic serum-like material",
                    "vector": "synthetic mosquito-like context",
                },
            },
        },
        conservation={
            "status": "COMPLETE",
            "counts": {"positions": 1200, "windows": 56},
            "windows": [
                {
                    "reference_start": 601,
                    "reference_end": 700,
                    "mean_conservation": 1.0,
                    "minimum_conservation": 1.0,
                    "mean_coverage": 1.0,
                    "mean_gap_frequency": 0.0,
                    "mean_entropy_bits": 0.0,
                }
            ],
        },
        contrastive={
            "status": "COMPLETE",
            "counts": {
                "windows": 56,
                "candidate_regions": 2,
                "challenge_datasets": 1,
            },
            "challenge_datasets": [
                {"name": "Related synthetic A", "criticality": "CRITICAL", "sequence_count": 1}
            ],
            "windows": [
                {
                    "reference_start": 601,
                    "reference_end": 700,
                    "target_mean_conservation": 1.0,
                    "worst_similarity": 0.53,
                    "contrast_margin": 0.47,
                }
            ],
            "candidates": [
                {
                    "region": {
                        "region_id": "contrast-region-001",
                        "rank": 1,
                        "reference_start": 501,
                        "reference_end": 800,
                        "peak_start": 601,
                        "peak_end": 700,
                        "mean_conservation": 1.0,
                        "minimum_conservation": 1.0,
                        "mean_coverage": 1.0,
                        "mean_gap_frequency": 0.0,
                        "mean_entropy_bits": 0.0,
                    },
                    "contributing_windows": [[601, 700]],
                    "worst_dataset_name": "Related synthetic A",
                    "worst_dataset_criticality": "CRITICAL",
                    "worst_similarity": 0.53,
                    "contrast_margin": 0.47,
                }
            ],
        },
        primer_design={
            "status": "COMPLETE",
            "candidate_source": "CONTRASTIVE_CONSERVATION",
            "counts": {"candidates": 2, "assays": 1},
            "candidates": [
                {
                    "region_id": "contrast-region-001",
                    "reference_start": 501,
                    "reference_end": 800,
                    "peak_start": 601,
                    "peak_end": 700,
                }
            ],
            "assays": [
                {
                    "assay_id": "contrast-region-001-assay-001",
                    "region_id": "contrast-region-001",
                    "forward_primer": {
                        "sequence": "TACAAGGCTCTCATGCACCC",
                        "reference_start": 580,
                        "reference_end": 599,
                        "tm": 59.7,
                        "gc_percent": 55.0,
                        "penalty": 0.25,
                    },
                    "probe": {
                        "sequence": "AACCGGTTAACCGGTTAACCGGTTA",
                        "reference_start": 620,
                        "reference_end": 644,
                        "tm": 70.0,
                        "gc_percent": 48.0,
                        "penalty": 0.4,
                    },
                    "reverse_primer": {
                        "sequence": "CGCTTGAAGTGTCCTACCAGT",
                        "reference_start": 701,
                        "reference_end": 721,
                        "tm": 60.0,
                        "gc_percent": 52.4,
                        "penalty": 1.0,
                    },
                    "product_size": 142,
                    "pair_penalty": 1.25,
                }
            ],
        },
        inclusivity={
            "status": "COMPLETE",
            "counts": {"assays": 1, "evaluation_sequences": 4, "original_compatible": 4},
            "assay_results": [],
            "proposals": [],
        },
        specificity={
            "status": "COMPLETE",
            "counts": {
                "assays": 1,
                "datasets": 1,
                "detectable_off_targets": detectable_count,
                "plausible_amplicons": detectable_count,
            },
            "retention": [],
            "amplicons": [],
        },
        ranking={"status": "COMPLETE", "counts": counts, "assays": [ranking_assay]},
    )


class ResearcherReportHtmlTests(unittest.TestCase):
    def test_report_has_readable_summary_and_full_evidence_hierarchy(self) -> None:
        html = render_researcher_report_html(report_data())

        for expected in (
            "Geison Researcher Report",
            "Run summary",
            "Scientific outcome",
            "In-silico candidate(s) identified",
            "Approved panel and study context",
            "Target conservation",
            "Target vs non-target contrast",
            "Assay design",
            "Target coverage / inclusivity",
            "Specificity",
            "Final candidates",
            "Interpretation and limitations",
            "Reproducibility",
            "Recommended in-silico candidate",
        ):
            with self.subTest(expected=expected):
                self.assertIn(expected, html)

    def test_report_exposes_candidate_region_and_contrast_anchor_separately(self) -> None:
        html = render_researcher_report_html(report_data())

        self.assertIn("Candidate region", html)
        self.assertIn("501–800", html)
        self.assertIn("Contrast anchor", html)
        self.assertIn("601–700", html)
        self.assertIn("CONTRASTIVE_CONSERVATION", html)
        self.assertIn("Anchor contained: Yes", html)

    def test_high_risk_run_is_negative_and_never_recommended(self) -> None:
        html = render_researcher_report_html(
            report_data(classification="HIGH_RISK")
        )

        self.assertIn("No in-silico acceptable assay candidates identified", html)
        self.assertIn("DETECTABLE_OFF_TARGET", html)
        self.assertNotIn("Recommended in-silico candidate", html)

    def test_partial_zero_evidence_is_visible_as_inconclusive(self) -> None:
        data = report_data(run_status="PARTIAL")
        data = ResearcherReportData(
            output_dir=data.output_dir,
            run_manifest=data.run_manifest,
            run_summary=data.run_summary,
            qc_report=data.qc_report,
            panel=data.panel,
            conservation=None,
            contrastive=None,
            primer_design=None,
            inclusivity=None,
            specificity=None,
            ranking=None,
        )

        html = render_researcher_report_html(data)

        self.assertIn("Inconclusive - insufficient evidence", html)
        self.assertIn("Evidence unavailable", html)
        self.assertIn("NO_ASSAYS", html)

    def test_report_escapes_dynamic_text_and_has_no_remote_resources_or_script(self) -> None:
        html = render_researcher_report_html(
            report_data(target_name='<img src=x onerror="boom">')
        )

        self.assertNotIn('<img src=x onerror="boom">', html)
        self.assertIn("&lt;img", html)
        lowered = html.lower()
        self.assertNotIn("http://", lowered)
        self.assertNotIn("https://", lowered)
        self.assertNotIn("<script", lowered)


if __name__ == "__main__":
    unittest.main()
