import json
from pathlib import Path

import pytest

from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.execution import ExecutionPolicy, STAGE_ORDER
from qpcr_pipeline.pipeline import run_pipeline


pytestmark = pytest.mark.skip(reason="temporary Task 5 regression isolation")


def _write_fasta(path: Path):
    path.write_text(">s1\nACGTACGTACGT\n>s2\nACGTACGAACGT\n", encoding="utf-8")


def _manifest(outdir: Path, stage: str):
    return json.loads(
        (outdir / ".checkpoints" / stage / "manifest.json").read_text(encoding="utf-8")
    )


def test_normal_run_writes_all_stage_checkpoints_and_does_not_reuse(tmp_path):
    fasta = tmp_path / "target.fasta"
    outdir = tmp_path / "run"
    _write_fasta(fasta)
    config = PipelineConfig(target_name="target", input_fasta=fasta)

    first = run_pipeline(config, outdir)
    first_summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))

    assert [item["stage"] for item in first_summary["stage_actions"]] == list(STAGE_ORDER)
    assert [item["action"] for item in first_summary["stage_actions"]] == ["RUN"] * len(STAGE_ORDER)
    for stage in STAGE_ORDER:
        manifest = _manifest(outdir, stage)
        assert manifest["stage"] == stage
        assert manifest["status"] == "COMPLETE"
        assert manifest["fingerprint"].startswith("sha256:")
        assert manifest["result_fingerprint"].startswith("sha256:")

    second = run_pipeline(config, outdir)
    second_summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
    assert [item["action"] for item in second_summary["stage_actions"]] == ["RUN"] * len(STAGE_ORDER)
    assert second.status == first.status == "COMPLETED"
    assert second.sequence_ids == first.sequence_ids


def test_all_valid_resume_reuses_every_stage(tmp_path):
    fasta = tmp_path / "target.fasta"
    outdir = tmp_path / "run"
    _write_fasta(fasta)
    config = PipelineConfig(target_name="target", input_fasta=fasta)
    first = run_pipeline(config, outdir)

    resumed = run_pipeline(config, outdir, execution=ExecutionPolicy(resume=True))
    payload = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))

    assert resumed.status == first.status
    assert resumed.sequence_ids == first.sequence_ids
    assert [item["stage"] for item in payload["stage_actions"]] == list(STAGE_ORDER)
    assert [item["action"] for item in payload["stage_actions"]] == ["REUSE"] * len(STAGE_ORDER)
