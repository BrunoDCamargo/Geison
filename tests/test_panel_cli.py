import json
import sys
from pathlib import Path
from types import SimpleNamespace

from qpcr_pipeline import cli
from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.pipeline import RunSummary


def test_dry_run_render_explains_panel_approval_destination(tmp_path):
    proposal_path = tmp_path / "run" / "panel_proposal.yaml"
    report = SimpleNamespace(
        target_name="target",
        decisions=(),
        environment=SimpleNamespace(missing_required_tools=()),
        panel_action_required=True,
        panel_proposal_would_be_written=str(proposal_path),
    )

    rendered = cli._render_dry_run(report)

    assert "Panel approval required before scientific execution." in rendered
    assert f"Proposal would be written to: {proposal_path}" in rendered


def test_panel_approve_command_writes_frozen_manifest(tmp_path, monkeypatch, capsys):
    proposal = Path("tests/fixtures/panels/west_nile_proposal.yaml")
    approved = tmp_path / "approved.json"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "qpcr-pipeline",
            "panel",
            "approve",
            str(proposal),
            "--output",
            str(approved),
        ],
    )

    assert cli.main() == 0
    payload = json.loads(approved.read_text(encoding="utf-8"))
    assert payload["status"] == "APPROVED"
    assert payload["approved_by_user"] is True
    assert payload["proposal_sha256"].startswith("sha256:")
    assert "Approved panel manifest:" in capsys.readouterr().out


def test_run_returns_three_when_panel_approval_is_required(
    tmp_path,
    monkeypatch,
    capsys,
):
    config_path = tmp_path / "config.yaml"
    proposal_path = tmp_path / "run" / "panel" / "panel_proposal.yaml"
    monkeypatch.setattr(
        cli,
        "load_config",
        lambda path: PipelineConfig(
            target_name="target",
            input_fasta=tmp_path / "target.fasta",
        ),
    )
    monkeypatch.setattr(
        cli,
        "run_pipeline",
        lambda *args, **kwargs: RunSummary(
            status="ACTION_REQUIRED",
            target_name="target",
            sequence_count=0,
            sequence_ids=[],
            action_required_code="PANEL_APPROVAL_REQUIRED",
            action_required_artifact=str(proposal_path),
        ),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "qpcr-pipeline",
            "run",
            str(config_path),
            "--outdir",
            str(tmp_path / "run"),
        ],
    )

    assert cli.main() == 3
    output = capsys.readouterr().out
    assert "PANEL_APPROVAL_REQUIRED" in output
    assert str(proposal_path) in output
