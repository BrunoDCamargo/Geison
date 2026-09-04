from qpcr_pipeline.contrastive_conservation import (
    ContrastCandidateRegion,
    ContrastWindowEvidence,
    DatasetWindowEvidence,
)
from qpcr_pipeline.contrastive_report_html import render_contrastive_html
from qpcr_pipeline.region_selection import CandidateRegion


def _candidate() -> ContrastCandidateRegion:
    region = CandidateRegion(
        region_id="contrast-region-001",
        rank=1,
        reference_start=101,
        reference_end=300,
        peak_start=150,
        peak_end=200,
        position_count=200,
        usable_length=200,
        usable_fraction=1.0,
        mean_conservation=0.99,
        minimum_conservation=0.95,
        mean_coverage=1.0,
        mean_gap_frequency=0.0,
        mean_entropy_bits=0.02,
    )
    return ContrastCandidateRegion(
        region=region,
        contributing_windows=((150, 200),),
        worst_dataset_name="challenge<script>alert(1)</script>",
        worst_dataset_criticality="CRITICAL",
        worst_similarity=0.12,
        worst_critical_similarity=0.12,
        worst_important_similarity=None,
        contrast_margin=0.87,
    )


def test_report_is_self_contained_safe_and_traceable():
    windows = (
        ContrastWindowEvidence(
            150, 200, 0.99, 0.95, 1.0, 0.0, 0.02, True,
            "challenge<script>alert(1)</script>", "CRITICAL", 0.12, 0.12, None, 0.87,
        ),
    )
    dataset = (
        DatasetWindowEvidence(
            150, 200, "challenge<script>alert(1)</script>", "CRITICAL", 2,
            "seq-1", "forward", 0.12,
        ),
    )
    html = render_contrastive_html(
        target_name="Synthetic target",
        reference_id="ref-1",
        windows=windows,
        dataset_evidence=dataset,
        candidates=(_candidate(),),
    )
    assert "Target vs non-target contrast" in html
    assert "<canvas" in html
    assert "contrast-region-001" in html
    assert "CRITICAL" in html
    assert "http://" not in html
    assert "https://" not in html
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    assert "const targetPoints" in html
    assert "const challengePoints" in html
    assert "pointMap" not in html


def test_report_rendering_is_deterministic():
    html_a = render_contrastive_html(
        target_name="target",
        reference_id=None,
        windows=(),
        dataset_evidence=(),
        candidates=(),
    )
    html_b = render_contrastive_html(
        target_name="target",
        reference_id=None,
        windows=(),
        dataset_evidence=(),
        candidates=(),
    )
    assert html_a == html_b
