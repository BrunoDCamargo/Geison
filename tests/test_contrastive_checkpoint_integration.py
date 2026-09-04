from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

from qpcr_pipeline import checkpoint_stages, pipeline
from qpcr_pipeline.config import (
    ConservationConfig,
    ContrastiveConservationConfig,
    OffTargetConfig,
    PanelConfig,
    PipelineConfig,
    PrimerDesignConfig,
)
from qpcr_pipeline.conservation import (
    ConservationResult,
    PositionConservation,
    WindowConservation,
)
from qpcr_pipeline.contrastive_conservation import ContrastiveConservationResult
from qpcr_pipeline.panel import (
    DiagnosticContext,
    PanelDefinition,
    PanelNonTarget,
    PanelTarget,
    TargetGroup,
)
from qpcr_pipeline.panel_manifest import ApprovedPanelManifest


REGION_SELECTION_FIELDS = {
    "max_candidate_regions",
    "candidate_region_length",
    "max_region_overlap_fraction",
    "min_mean_conservation",
    "min_minimum_conservation",
    "min_mean_coverage",
    "max_mean_gap_frequency",
    "max_mean_entropy_bits",
    "min_usable_fraction",
}


class FakeToolIdentityProvider:
    def __init__(self):
        self.calls = []

    def identity(self, tool_name):
        self.calls.append(tool_name)
        return {"name": tool_name, "version": "test"}


def _manifest(*, criticality="CRITICAL"):
    target = PanelTarget(
        name="target",
        taxid=None,
        mode="broad_detection",
        subtype=None,
        groups=(
            TargetGroup(
                name="all",
                required=True,
                dataset_roles=("DESIGN",),
                reasons=("fixture",),
                proposed_by=("test",),
            ),
        ),
    )
    challenge = PanelNonTarget(
        name="challenge",
        taxid=None,
        criticality=criticality,
        dataset_roles=("CHALLENGE",),
        reasons=("fixture",),
        proposed_by=("test",),
    )
    design_only = PanelNonTarget(
        name="design-only",
        taxid=None,
        criticality="BACKGROUND",
        dataset_roles=("DESIGN",),
        reasons=("fixture",),
        proposed_by=("test",),
    )
    return ApprovedPanelManifest(
        schema_version=1,
        status="APPROVED",
        approved_by_user=True,
        proposal_sha256="sha256:" + "a" * 64,
        definition=PanelDefinition(
            target=target,
            non_targets=(challenge, design_only),
            diagnostic_context=DiagnosticContext(),
        ),
    )


def _config(tmp_path):
    target = tmp_path / "target.fasta"
    target.write_text(">target\n" + "A" * 200 + "\n", encoding="utf-8")
    approved = tmp_path / "approved.json"
    approved.write_text("approved-a\n", encoding="utf-8")
    challenge = tmp_path / "challenge.fasta"
    challenge.write_text(">challenge\n" + "C" * 200 + "\n", encoding="utf-8")
    design_only = tmp_path / "design-only.fasta"
    design_only.write_text(">design\n" + "G" * 200 + "\n", encoding="utf-8")
    return PipelineConfig(
        target_name="target",
        input_fasta=target,
        panel=PanelConfig(frozen_manifest=approved),
        conservation=ConservationConfig(enabled=True),
        contrastive_conservation=ContrastiveConservationConfig(enabled=True),
        primer_design=PrimerDesignConfig(
            enabled=True,
            max_candidate_regions=4,
            candidate_region_length=100,
        ),
        off_targets=(
            OffTargetConfig("challenge", fasta=challenge),
            OffTargetConfig("design-only", fasta=design_only),
        ),
    )


def _perfect_conservation(outdir: Path):
    positions = tuple(
        PositionConservation(
            i, i, "A", 1, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, "A", "A"
        )
        for i in range(1, 201)
    )
    windows = (
        WindowConservation(1, 100, 100, 1.0, 1.0, 1.0, 0.0, 0.0),
    )
    return ConservationResult(
        "COMPLETE",
        "target",
        positions,
        windows,
        (),
        "A" * 200,
        "A" * 200,
        None,
        None,
        None,
        None,
        None,
        outdir / "conservation_report.json",
    )


def _empty_contrastive(outdir: Path):
    return ContrastiveConservationResult(
        status="COMPLETE",
        reference_id="target",
        windows=(),
        dataset_evidence=(),
        candidates=(),
        challenge_datasets=(),
        window_metrics_path=outdir / "contrastive_conservation" / "window_metrics.tsv",
        dataset_metrics_path=outdir / "contrastive_conservation" / "dataset_metrics.tsv",
        candidate_regions_path=outdir / "contrastive_conservation" / "candidate_regions.tsv",
        report_path=outdir / "contrastive_conservation" / "contrastive_conservation_report.json",
        html_report_path=outdir / "contrastive_conservation" / "report.html",
    )


def test_contrastive_parameters_include_only_region_selection_and_panel_semantics(tmp_path, monkeypatch):
    config = _config(tmp_path)
    monkeypatch.setattr(checkpoint_stages, "load_approved_panel_manifest", lambda path: _manifest())

    parameters = checkpoint_stages.stage_parameters("contrastive_conservation", config)

    assert parameters["config"] == {"enabled": True}
    assert set(parameters["region_selection"]) == REGION_SELECTION_FIELDS
    assert parameters["challenge_panel"] == [
        {
            "name": "challenge",
            "criticality": "CRITICAL",
            "dataset_roles": ["CHALLENGE"],
        }
    ]

    tm_only = replace(
        config,
        primer_design=replace(
            config.primer_design,
            primer=replace(config.primer_design.primer, min_tm=57.0),
        ),
    )
    assert checkpoint_stages.stage_parameters("contrastive_conservation", tm_only) == parameters

    region_change = replace(
        config,
        primer_design=replace(config.primer_design, max_candidate_regions=7),
    )
    assert checkpoint_stages.stage_parameters("contrastive_conservation", region_change) != parameters


def test_contrastive_input_identity_hashes_only_approved_challenge_sources(tmp_path, monkeypatch):
    config = _config(tmp_path)
    monkeypatch.setattr(checkpoint_stages, "load_approved_panel_manifest", lambda path: _manifest())

    first = checkpoint_stages.stage_input_identities("contrastive_conservation", config)
    assert [row["name"] for row in first["challenge_datasets"]] == ["challenge"]

    design_only = config.off_targets[1].fasta
    assert design_only is not None
    design_only.write_text(">design\nTTTT\n", encoding="utf-8")
    assert checkpoint_stages.stage_input_identities("contrastive_conservation", config) == first

    challenge = config.off_targets[0].fasta
    assert challenge is not None
    challenge.write_text(">challenge\nAAAA\n", encoding="utf-8")
    assert checkpoint_stages.stage_input_identities("contrastive_conservation", config) != first


def test_primer3_identity_uses_actual_contrastive_candidate_source(tmp_path):
    config = _config(tmp_path)
    provider = FakeToolIdentityProvider()
    conservation = _perfect_conservation(tmp_path)
    contrastive = _empty_contrastive(tmp_path)

    assert checkpoint_stages.stage_tool_identities(
        "primer_design",
        config,
        {"conservation": conservation, "contrastive_conservation": contrastive},
        provider,
    ) == {}
    assert provider.calls == []


def _stage_results(tmp_path):
    return {
        "input": (),
        "qc": SimpleNamespace(evaluation_set=SimpleNamespace(sequence_ids=())),
        "clustering": SimpleNamespace(discovery_set=SimpleNamespace(sequence_ids=())),
        "alignment": object(),
        "conservation": _perfect_conservation(tmp_path),
        "contrastive_conservation": _empty_contrastive(tmp_path),
        "panel": object(),
    }


def test_pipeline_executes_contrastive_stage_with_approved_manifest(tmp_path, monkeypatch):
    config = _config(tmp_path)
    results = _stage_results(tmp_path)
    approved = _manifest()
    sentinel = object()
    captured = {}

    monkeypatch.setattr(pipeline, "load_approved_panel_manifest", lambda path: approved, raising=False)

    def fake_analyze(conservation, panel, off_targets, contrast_config, primer_config, output_dir):
        captured["args"] = (conservation, panel, off_targets, contrast_config, primer_config, output_dir)
        return sentinel

    monkeypatch.setattr(pipeline, "analyze_contrastive_conservation", fake_analyze, raising=False)

    result = pipeline._run_stage(
        "contrastive_conservation",
        config,
        tmp_path,
        results,
        ncbi_client=None,
        cdhit_runner=None,
        mafft_runner=None,
        primer3_runner=None,
        refresh_online_input=False,
    )

    assert result is sentinel
    assert captured["args"][1] is approved
    assert captured["args"][2] == config.off_targets


def test_pipeline_passes_contrastive_result_into_primer_design(tmp_path, monkeypatch):
    config = _config(tmp_path)
    results = _stage_results(tmp_path)
    sentinel = object()
    captured = {}

    def fake_design(conservation, primer_config, output_dir, *, contrastive=None, runner=None):
        captured["contrastive"] = contrastive
        return sentinel

    monkeypatch.setattr(pipeline, "design_primers", fake_design)

    result = pipeline._run_stage(
        "primer_design",
        config,
        tmp_path,
        results,
        ncbi_client=None,
        cdhit_runner=None,
        mafft_runner=None,
        primer3_runner=None,
        refresh_online_input=False,
    )

    assert result is sentinel
    assert captured["contrastive"] is results["contrastive_conservation"]
