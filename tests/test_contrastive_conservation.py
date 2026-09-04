import csv
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from qpcr_pipeline.config import ContrastiveConservationConfig, OffTargetConfig, PrimerDesignConfig
from qpcr_pipeline.conservation import ConservationResult, PositionConservation, WindowConservation
from qpcr_pipeline.contrastive_conservation import analyze_contrastive_conservation
from qpcr_pipeline.contrastive_similarity import RegionSimilarity
from qpcr_pipeline.panel import DiagnosticContext, PanelDefinition, PanelNonTarget, PanelTarget, TargetGroup
from qpcr_pipeline.panel_manifest import ApprovedPanelManifest


def _position(index: int, base: str) -> PositionConservation:
    return PositionConservation(
        alignment_position=index,
        reference_position=index,
        reference_base=base,
        depth=10,
        coverage=1.0,
        frequency_a=1.0 if base == "A" else 0.0,
        frequency_c=1.0 if base == "C" else 0.0,
        frequency_g=0.0,
        frequency_t=0.0,
        gap_frequency=0.0,
        major_allele_frequency=1.0,
        entropy_bits=0.0,
        major_consensus=base,
        iupac_consensus=base,
    )


def _window(start: int, end: int) -> WindowConservation:
    return WindowConservation(
        reference_start=start,
        reference_end=end,
        position_count=end - start + 1,
        mean_conservation=1.0,
        minimum_conservation=1.0,
        mean_coverage=1.0,
        mean_gap_frequency=0.0,
        mean_entropy_bits=0.0,
    )


def _conservation() -> ConservationResult:
    bases = "A" * 200 + "C" * 200
    positions = tuple(_position(index, base) for index, base in enumerate(bases, 1))
    return ConservationResult(
        status="COMPLETE",
        reference_id="synthetic-ref",
        positions=positions,
        windows=(_window(80, 120), _window(280, 320), _window(290, 310)),
        annotations=(),
        major_consensus=bases,
        iupac_consensus=bases,
        position_metrics_path=None,
        window_metrics_path=None,
        major_consensus_path=None,
        iupac_consensus_path=None,
        html_report_path=None,
        report_path=Path("unused.json"),
    )


def _target() -> PanelTarget:
    return PanelTarget(
        name="target",
        taxid=None,
        mode="broad_detection",
        subtype=None,
        groups=(TargetGroup("all", True, ("DESIGN",), ("fixture",), ("test",)),),
    )


def _manifest(*items: PanelNonTarget) -> ApprovedPanelManifest:
    return ApprovedPanelManifest(
        schema_version=1,
        status="APPROVED",
        approved_by_user=True,
        proposal_sha256="sha256:" + "a" * 64,
        definition=PanelDefinition(_target(), tuple(items), DiagnosticContext()),
    )


def _non_target(name: str, criticality: str) -> PanelNonTarget:
    return PanelNonTarget(
        name=name,
        taxid=None,
        criticality=criticality,
        dataset_roles=("CHALLENGE",),
        reasons=("fixture",),
        proposed_by=("test",),
    )


def _primer_config(**changes) -> PrimerDesignConfig:
    return replace(
        PrimerDesignConfig(
            candidate_region_length=100,
            max_candidate_regions=4,
            min_mean_conservation=0.0,
            min_minimum_conservation=0.0,
            min_mean_coverage=0.0,
            max_mean_gap_frequency=1.0,
            max_mean_entropy_bits=2.0,
            min_usable_fraction=0.0,
            max_region_overlap_fraction=0.5,
        ),
        **changes,
    )


def _write_fasta(path: Path, sequence_id: str) -> Path:
    path.write_text(f">{sequence_id}\nACGTACGTACGT\n", encoding="utf-8")
    return path


class _FakeEngine:
    def best_match(self, query, records):
        record = records[0]
        is_shared = query.startswith("A")
        if "critical" in record.sequence_id:
            similarity = 0.90 if is_shared else 0.10
        else:
            similarity = 0.80 if is_shared else 0.20
        return RegionSimilarity(record.sequence_id, similarity, "forward")


def test_disabled_path_is_skipped_without_loading_panel_datasets(tmp_path):
    with patch("qpcr_pipeline.contrastive_conservation.resolve_challenge_datasets") as resolver:
        result = analyze_contrastive_conservation(
            object(),
            None,
            (),
            ContrastiveConservationConfig(enabled=False),
            _primer_config(),
            tmp_path,
        )
    resolver.assert_not_called()
    assert result.status == "SKIPPED"
    assert result.windows == ()
    assert result.dataset_evidence == ()
    assert result.candidates == ()
    assert result.html_report_path is None
    assert result.report_path.exists()


def test_discriminant_region_ranks_first_and_consolidates_windows(tmp_path):
    critical = _write_fasta(tmp_path / "critical.fasta", "critical-seq")
    important = _write_fasta(tmp_path / "important.fasta", "important-seq")
    manifest = _manifest(
        _non_target("critical", "CRITICAL"),
        _non_target("important", "IMPORTANT"),
    )
    result = analyze_contrastive_conservation(
        _conservation(),
        manifest,
        (
            OffTargetConfig("critical", fasta=critical),
            OffTargetConfig("important", fasta=important),
        ),
        ContrastiveConservationConfig(enabled=True),
        _primer_config(),
        tmp_path,
        similarity_engine=_FakeEngine(),
    )
    assert result.status == "COMPLETE"
    assert len(result.windows) == 3
    assert len(result.dataset_evidence) == 6
    assert result.candidates[0].region.peak_start >= 280
    assert result.candidates[0].worst_critical_similarity == 0.10
    assert result.candidates[0].contributing_windows == ((280, 320), (290, 310))
    assert result.candidates[0].region.region_id == "contrast-region-001"
    assert result.window_metrics_path.exists()
    assert result.dataset_metrics_path.exists()
    assert result.candidate_regions_path.exists()
    assert result.report_path.exists()


def test_tsv_serialization_round_trips_control_characters_in_dataset_name(tmp_path):
    dataset_name = "Challenge\tName\nLine"
    fasta = _write_fasta(tmp_path / "special.fasta", "critical-special")
    result = analyze_contrastive_conservation(
        _conservation(),
        _manifest(_non_target(dataset_name, "CRITICAL")),
        (OffTargetConfig(dataset_name, fasta=fasta),),
        ContrastiveConservationConfig(enabled=True),
        _primer_config(max_candidate_regions=1),
        tmp_path,
        similarity_engine=_FakeEngine(),
    )
    with result.dataset_metrics_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert {row["dataset_name"] for row in rows} == {dataset_name}
