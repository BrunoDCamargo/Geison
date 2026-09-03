import json
from dataclasses import replace
from pathlib import Path

import pytest

from qpcr_pipeline.config import PanelConfig, PipelineConfig, PrimerDesignConfig
from qpcr_pipeline.execution import ExecutionPolicy
from qpcr_pipeline.panel_manifest import approve_panel_proposal
from qpcr_pipeline.pipeline import run_pipeline
from panel_fixtures import proposal_panel_config


def panel_proposal_pipeline_config(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGTACGTACGT\n", encoding="utf-8")
    return PipelineConfig(
        target_name="Target virus",
        input_fasta=fasta,
        panel=proposal_panel_config("Target virus"),
    )


def test_panel_proposal_stops_before_input_and_writes_review_artifact(tmp_path):
    config = panel_proposal_pipeline_config(tmp_path)
    outdir = tmp_path / "run"

    summary = run_pipeline(config, outdir)

    proposal_path = outdir / "panel_proposal.yaml"
    assert summary.status == "ACTION_REQUIRED"
    assert summary.action_required_code == "PANEL_APPROVAL_REQUIRED"
    assert summary.action_required_artifact == str(proposal_path)
    assert summary.sequence_count == 0
    assert summary.sequence_ids == []
    assert proposal_path.is_file()
    assert not (outdir / ".checkpoints" / "input" / "manifest.json").exists()
    assert not (outdir / "qc_report.json").exists()

    manifest = json.loads(
        (outdir / "run_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["status"] == "ACTION_REQUIRED"
    assert manifest["action_required"] == {
        "code": "PANEL_APPROVAL_REQUIRED",
        "artifact": str(proposal_path),
    }
    assert manifest["attempts"][0]["status"] == "ACTION_REQUIRED"
    assert manifest["failure"] is None

    rows = [
        json.loads(line)
        for line in (outdir / "run.log.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    assert [row["event"] for row in rows].count("run_action_required") == 1


def test_repeated_proposal_runs_write_identical_yaml(tmp_path):
    config = panel_proposal_pipeline_config(tmp_path)
    outdir = tmp_path / "run"
    run_pipeline(config, outdir)
    first = (outdir / "panel_proposal.yaml").read_bytes()

    run_pipeline(config, outdir)

    assert (outdir / "panel_proposal.yaml").read_bytes() == first


def test_action_required_replaces_stale_summary_and_qc_report(tmp_path):
    config = panel_proposal_pipeline_config(tmp_path)
    outdir = tmp_path / "run"
    outdir.mkdir()
    (outdir / "run_summary.json").write_text(
        '{"status":"COMPLETED"}',
        encoding="utf-8",
    )
    (outdir / "qc_report.json").write_text("stale", encoding="utf-8")

    run_pipeline(config, outdir)

    summary = json.loads(
        (outdir / "run_summary.json").read_text(encoding="utf-8")
    )
    assert summary["status"] == "ACTION_REQUIRED"
    assert summary["action_required_code"] == "PANEL_APPROVAL_REQUIRED"
    assert not (outdir / "qc_report.json").exists()


def test_proposal_approve_rerun_reaches_checkpointed_pipeline(tmp_path):
    proposal_config = panel_proposal_pipeline_config(tmp_path)
    outdir = tmp_path / "run"
    first = run_pipeline(proposal_config, outdir)
    assert first.status == "ACTION_REQUIRED"

    approved = tmp_path / "approved_panel.json"
    approve_panel_proposal(outdir / "panel_proposal.yaml", approved)
    approved_config = replace(
        proposal_config,
        panel=PanelConfig(frozen_manifest=approved),
        primer_design=PrimerDesignConfig(enabled=False),
    )
    second = run_pipeline(
        approved_config,
        outdir,
        execution=ExecutionPolicy(resume=True),
    )

    assert second.status in {"PARTIAL", "COMPLETED"}
    assert (outdir / "panel" / "approved_panel.json").is_file()
    assert (outdir / ".checkpoints" / "panel" / "manifest.json").is_file()
    assert (outdir / ".checkpoints" / "input" / "manifest.json").is_file()


def test_frozen_panel_target_must_match_pipeline_target(tmp_path):
    approved = tmp_path / "approved_panel.json"
    approve_panel_proposal(
        Path("tests/fixtures/panels/west_nile_proposal.yaml"),
        approved,
    )
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGTACGTACGT\n", encoding="utf-8")
    config = PipelineConfig(
        target_name="Zika virus",
        input_fasta=fasta,
        panel=PanelConfig(frozen_manifest=approved),
    )

    with pytest.raises(
        ValueError,
        match="Approved panel target 'West Nile virus'.*pipeline target 'Zika virus'",
    ):
        run_pipeline(config, tmp_path / "run")
