import json

from qpcr_pipeline.config import PipelineConfig
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
