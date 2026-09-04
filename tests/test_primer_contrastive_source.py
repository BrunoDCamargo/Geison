import json
from dataclasses import replace
from pathlib import Path

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.conservation import ConservationResult, PositionConservation, WindowConservation
from qpcr_pipeline.contrastive_conservation import (
    ContrastCandidateRegion,
    ContrastiveConservationResult,
)
from qpcr_pipeline.primer_design import design_primers, primer3_required
from qpcr_pipeline.region_selection import CandidateRegion, select_conservation_candidate_regions


def _position(index: int) -> PositionConservation:
    return PositionConservation(
        alignment_position=index,
        reference_position=index,
        reference_base="A",
        depth=10,
        coverage=1.0,
        frequency_a=1.0,
        frequency_c=0.0,
        frequency_g=0.0,
        frequency_t=0.0,
        gap_frequency=0.0,
        major_allele_frequency=1.0,
        entropy_bits=0.0,
        major_consensus="A",
        iupac_consensus="A",
    )


def _conservation() -> ConservationResult:
    positions = tuple(_position(i) for i in range(1, 301))
    window = WindowConservation(150, 150, 1, 1.0, 1.0, 1.0, 0.0, 0.0)
    consensus = "A" * 300
    return ConservationResult(
        status="COMPLETE",
        reference_id="ref",
        positions=positions,
        windows=(window,),
        annotations=(),
        major_consensus=consensus,
        iupac_consensus=consensus,
        position_metrics_path=None,
        window_metrics_path=None,
        major_consensus_path=None,
        iupac_consensus_path=None,
        html_report_path=None,
        report_path=Path("unused.json"),
    )


def _config() -> PrimerDesignConfig:
    return replace(
        PrimerDesignConfig(
            enabled=True,
            candidate_region_length=200,
            min_mean_conservation=0.0,
            min_minimum_conservation=0.0,
            min_mean_coverage=0.0,
            max_mean_gap_frequency=1.0,
            max_mean_entropy_bits=2.0,
            min_usable_fraction=0.0,
        ),
        product_size_min=70,
        product_size_max=200,
    )


def _contrastive(candidates: tuple[ContrastCandidateRegion, ...]) -> ContrastiveConservationResult:
    return ContrastiveConservationResult(
        status="COMPLETE",
        reference_id="ref",
        windows=(),
        dataset_evidence=(),
        candidates=candidates,
        challenge_datasets=(),
        window_metrics_path=Path("contrast/window.tsv"),
        dataset_metrics_path=Path("contrast/dataset.tsv"),
        candidate_regions_path=Path("contrast/candidates.tsv"),
        report_path=Path("contrast/report.json"),
        html_report_path=None,
    )


class _ZeroPairRunner:
    def __init__(self, sequence_id: str):
        self.sequence_id = sequence_id
        self.calls = 0

    def run(self, input_text: str) -> str:
        self.calls += 1
        assert f"SEQUENCE_ID={self.sequence_id}" in input_text
        return (
            f"SEQUENCE_ID={self.sequence_id}\n"
            "PRIMER_WARNING=no viable pair\n"
            "PRIMER_LEFT_EXPLAIN=considered 1, low tm 1\n"
            "PRIMER_INTERNAL_EXPLAIN=considered 1, high tm 1\n"
            "PRIMER_RIGHT_EXPLAIN=considered 1, high tm 1\n"
            "PRIMER_PAIR_EXPLAIN=considered 0\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
        )


def test_legacy_selection_records_conservation_only_source(tmp_path):
    conservation = _conservation()
    config = _config()
    runner = _ZeroPairRunner("region-001")
    result = design_primers(conservation, config, tmp_path, runner=runner)
    assert result.candidate_source == "CONSERVATION_ONLY"
    assert result.candidates[0].region_id == "region-001"
    report = json.loads(result.report_path.read_text(encoding="utf-8"))
    assert report["candidate_source"] == "CONSERVATION_ONLY"


def test_complete_contrastive_result_is_the_only_candidate_source(tmp_path):
    conservation = _conservation()
    config = _config()
    base = select_conservation_candidate_regions(conservation, config)[0]
    contrast_region = replace(base, region_id="contrast-region-001", rank=1)
    contrastive_candidate = ContrastCandidateRegion(
        region=contrast_region,
        contributing_windows=((150, 150),),
        worst_dataset_name="challenge",
        worst_dataset_criticality="CRITICAL",
        worst_similarity=0.1,
        worst_critical_similarity=0.1,
        worst_important_similarity=None,
        contrast_margin=0.9,
    )
    contrastive = _contrastive((contrastive_candidate,))
    runner = _ZeroPairRunner("contrast-region-001")
    result = design_primers(
        conservation,
        config,
        tmp_path,
        contrastive=contrastive,
        runner=runner,
    )
    assert result.candidate_source == "CONTRASTIVE_CONSERVATION"
    assert result.candidates == (contrast_region,)
    assert primer3_required(conservation, config, contrastive=contrastive) is True
    report = json.loads(result.report_path.read_text(encoding="utf-8"))
    assert report["candidate_source"] == "CONTRASTIVE_CONSERVATION"


def test_complete_contrastive_with_zero_candidates_never_falls_back(tmp_path):
    conservation = _conservation()
    config = _config()
    contrastive = _contrastive(())
    result = design_primers(
        conservation,
        config,
        tmp_path,
        contrastive=contrastive,
    )
    assert result.candidate_source == "CONTRASTIVE_CONSERVATION"
    assert result.candidates == ()
    assert result.assays == ()
    assert primer3_required(conservation, config, contrastive=contrastive) is False
