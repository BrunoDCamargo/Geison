from dataclasses import replace
from pathlib import Path

import pytest

from qpcr_pipeline import checkpoint_stages
from qpcr_pipeline.alignment import AlignmentResult
from qpcr_pipeline.clustering import ClusteringResult
from qpcr_pipeline.config import (
    ClusteringConfig,
    NcbiInputConfig,
    OffTargetConfig,
    PipelineConfig,
    PrimerDesignConfig,
    RankingConfig,
    RankingWeights,
    SpecificityConfig,
)
from qpcr_pipeline.conservation import ConservationResult, PositionConservation, WindowConservation
from qpcr_pipeline.models import DiscoverySet, EvaluationSet, TargetSequenceSet
from qpcr_pipeline.ncbi import AcquiredNcbiDataset
from qpcr_pipeline.primer_design import PrimerDesignResult, primer3_required
from qpcr_pipeline.qc import QCResult
from qpcr_pipeline.ranking import RankingResult
from qpcr_pipeline.specificity import SpecificityResult
from qpcr_pipeline.checkpoint_stages import (
    STAGE_DEFINITIONS,
    stage_input_identities,
    stage_outputs,
    stage_parameters,
    stage_tool_identities,
)


class FakeToolIdentityProvider:
    def __init__(self, versions=None):
        self.versions = versions or {
            "cd-hit-est": "4.8.1",
            "mafft": "7.526",
            "primer3_core": "2.6.1",
        }
        self.calls = []

    def identity(self, tool_name):
        self.calls.append(tool_name)
        return {"name": tool_name, "version": self.versions[tool_name]}


def _config(fasta: Path, **changes):
    config = PipelineConfig(target_name="target", input_fasta=fasta)
    return replace(config, **changes)


def test_registry_contains_every_pipeline_stage():
    assert tuple(STAGE_DEFINITIONS) == (
        "input", "qc", "clustering", "alignment", "conservation",
        "primer_design", "inclusivity", "specificity", "ranking",
    )


def test_specificity_parameter_change_is_isolated_to_specificity(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    first = _config(fasta)
    second = replace(first, specificity=SpecificityConfig(max_hits_per_oligo_per_dataset=7))
    unaffected = ("input", "qc", "clustering", "alignment", "conservation", "primer_design", "inclusivity")
    for stage in unaffected:
        assert stage_parameters(stage, first) == stage_parameters(stage, second)
    assert stage_parameters("specificity", first) != stage_parameters("specificity", second)


def test_clustering_parameter_change_changes_only_clustering_projection_directly(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    first = _config(fasta)
    second = replace(first, clustering=ClusteringConfig(identity=0.90))
    assert stage_parameters("clustering", first) != stage_parameters("clustering", second)
    for stage in ("qc", "alignment", "conservation", "primer_design", "inclusivity", "specificity", "ranking"):
        assert stage_parameters(stage, first) == stage_parameters(stage, second)


def test_ranking_weights_change_only_ranking_projection(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    first = _config(fasta)
    second = replace(first, ranking=RankingConfig(weights=RankingWeights(inclusivity=0.30, specificity=0.30, conservation=0.20, primer3_quality=0.10, robustness=0.10)))
    assert stage_parameters("ranking", first) != stage_parameters("ranking", second)
    for stage in ("input", "qc", "clustering", "alignment", "conservation", "primer_design", "inclusivity", "specificity"):
        assert stage_parameters(stage, first) == stage_parameters(stage, second)


def test_local_input_identity_uses_bytes_not_path(tmp_path):
    first_path = tmp_path / "one.fasta"
    second_path = tmp_path / "two.fasta"
    first_path.write_text(">s1\nACGT\n", encoding="utf-8")
    second_path.write_bytes(first_path.read_bytes())
    first_identity = stage_input_identities("input", _config(first_path))
    second_identity = stage_input_identities("input", _config(second_path))
    assert first_identity == second_identity
    first_path.write_text(">s1\nACGA\n", encoding="utf-8")
    assert stage_input_identities("input", _config(first_path)) != second_identity


def test_frozen_ncbi_identity_uses_records_and_manifest_bytes(tmp_path, monkeypatch):
    frozen = tmp_path / "frozen"
    frozen.mkdir()
    records = frozen / "records.gb"
    manifest = frozen / "dataset_manifest.json"
    records.write_text("records-a", encoding="utf-8")
    manifest.write_text('{"a":1}\n', encoding="utf-8")
    monkeypatch.setattr(
        checkpoint_stages,
        "validate_frozen_dataset",
        lambda path: AcquiredNcbiDataset(records_path=records, manifest_path=manifest),
    )
    config = PipelineConfig(target_name="target", input_ncbi=NcbiInputConfig(frozen_dataset=frozen))
    first = stage_input_identities("input", config)
    records.write_text("records-b", encoding="utf-8")
    second = stage_input_identities("input", config)
    assert first != second


def test_specificity_input_identity_changes_with_off_target_bytes(tmp_path):
    target = tmp_path / "target.fasta"
    off_target = tmp_path / "off.fasta"
    target.write_text(">s1\nACGT\n", encoding="utf-8")
    off_target.write_text(">o1\nAAAA\n", encoding="utf-8")
    config = replace(_config(target), off_targets=(OffTargetConfig("db", fasta=off_target),))
    first = stage_input_identities("specificity", config)
    off_target.write_text(">o1\nTTTT\n", encoding="utf-8")
    second = stage_input_identities("specificity", config)
    assert first != second


def test_disabled_clustering_does_not_request_tool_identity(tmp_path):
    provider = FakeToolIdentityProvider()
    config = _config(tmp_path / "target.fasta")
    qc = QCResult((), TargetSequenceSet(("s1",)), EvaluationSet(("s1",)))
    assert stage_tool_identities("clustering", config, {"qc": qc}, provider) == {}
    assert provider.calls == []


def test_enabled_nonempty_clustering_requests_cd_hit_identity(tmp_path):
    provider = FakeToolIdentityProvider()
    config = _config(tmp_path / "target.fasta", clustering=ClusteringConfig(enabled=True))
    qc = QCResult((), TargetSequenceSet(("s1",)), EvaluationSet(("s1",)))
    tools = stage_tool_identities("clustering", config, {"qc": qc}, provider)
    assert tools["cd-hit-est"]["version"] == "4.8.1"
    assert provider.calls == ["cd-hit-est"]


def test_alignment_requests_mafft_only_for_more_than_one_discovery_sequence(tmp_path):
    provider = FakeToolIdentityProvider()
    config = replace(_config(tmp_path / "target.fasta"), alignment=replace(PipelineConfig(target_name="x", input_fasta=tmp_path / "x").alignment, enabled=True))
    one = ClusteringResult(DiscoverySet(("s1",)), (), tmp_path / "d.fa", tmp_path / "r.json", None)
    two = ClusteringResult(DiscoverySet(("s1", "s2")), (), tmp_path / "d.fa", tmp_path / "r.json", None)
    assert stage_tool_identities("alignment", config, {"clustering": one}, provider) == {}
    assert provider.calls == []
    tools = stage_tool_identities("alignment", config, {"clustering": two}, provider)
    assert "mafft" in tools
    assert provider.calls == ["mafft"]


def _perfect_conservation(outdir: Path, position_count: int):
    positions = tuple(
        PositionConservation(i, i, "A", 1, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, "A", "A")
        for i in range(1, position_count + 1)
    )
    windows = (WindowConservation(1, position_count, position_count, 1.0, 1.0, 1.0, 0.0, 0.0),)
    return ConservationResult("COMPLETE", "s1", positions, windows, (), "A" * position_count, "A" * position_count, None, None, None, None, None, outdir / "conservation_report.json")


def test_primer3_required_matches_candidate_selection_and_tool_identity(tmp_path):
    provider = FakeToolIdentityProvider()
    conservation = _perfect_conservation(tmp_path, 100)
    disabled = PrimerDesignConfig(enabled=False, candidate_region_length=100)
    enabled = PrimerDesignConfig(enabled=True, candidate_region_length=100)
    assert primer3_required(conservation, disabled) is False
    assert primer3_required(ConservationResult("COMPLETE", "s1", (), (), (), "", "", None, None, None, None, None, tmp_path / "empty.json"), enabled) is False
    assert primer3_required(conservation, enabled) is True
    config = replace(_config(tmp_path / "target.fasta"), primer_design=enabled)
    assert "primer3_core" in stage_tool_identities("primer_design", config, {"conservation": conservation}, provider)
    assert provider.calls == ["primer3_core"]


def test_conservation_never_declares_root_report_alias(tmp_path):
    outdir = tmp_path / "run"
    result = ConservationResult("COMPLETE", "s1", (), (), (), "", "", None, None, None, None, outdir / "report.html", outdir / "conservation" / "conservation_report.json")
    outputs = stage_outputs("conservation", result, outdir)
    assert outdir / "report.html" not in outputs
    assert result.report_path in outputs


def test_ranking_declares_root_report_only_when_owned(tmp_path):
    outdir = tmp_path / "run"
    skipped = RankingResult("SKIPPED", (), None, outdir / "ranking" / "ranking_report.json", None)
    complete = RankingResult("COMPLETE", (), outdir / "ranking" / "assay_ranking.tsv", outdir / "ranking" / "ranking_report.json", outdir / "report.html")
    assert outdir / "report.html" not in stage_outputs("ranking", skipped, outdir)
    assert outdir / "report.html" in stage_outputs("ranking", complete, outdir)


def test_stage_outputs_omit_none_and_reject_outside_paths(tmp_path):
    outdir = tmp_path / "run"
    result = SpecificityResult("SKIPPED", (), 0, 0, (), (), (), None, None, outdir / "specificity" / "specificity_report.json")
    assert stage_outputs("specificity", result, outdir) == (result.report_path,)
    bad = replace(result, report_path=tmp_path / "outside.json")
    with pytest.raises(ValueError, match="inside"):
        stage_outputs("specificity", bad, outdir)
