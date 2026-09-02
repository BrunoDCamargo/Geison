import json
from unittest.mock import patch

import pytest

from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.execution import ExecutionPolicy
from qpcr_pipeline.pipeline import run_pipeline


def test_failed_attempt_then_resume_keeps_run_identity_and_sanitizes_failure(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGTACGTACGT\n>s2\nACGTACGAACGT\n", encoding="utf-8")
    config = PipelineConfig(target_name="target", input_fasta=fasta)
    outdir = tmp_path / "run"
    sensitive_sequence = "ACGT" * 50

    with patch(
        "qpcr_pipeline.pipeline.align_discovery",
        side_effect=RuntimeError("alignment failed for " + sensitive_sequence),
    ):
        with pytest.raises(RuntimeError, match="alignment failed"):
            run_pipeline(config, outdir)

    failed = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    run_id = failed["run_id"]
    assert failed["status"] == "FAILED"
    assert len(failed["attempts"]) == 1
    assert failed["attempts"][0]["status"] == "FAILED"
    assert failed["failure"]["stage"] == "alignment"
    assert not (outdir / "run_summary.json").exists()
    assert not (outdir / "qc_report.json").exists()
    assert sensitive_sequence not in json.dumps(failed)
    assert sensitive_sequence not in (outdir / "run.log.jsonl").read_text(encoding="utf-8")

    summary = run_pipeline(
        config,
        outdir,
        execution=ExecutionPolicy(resume=True),
    )

    resumed = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    assert summary.status == "PARTIAL"
    assert resumed["run_id"] == run_id
    assert resumed["status"] == "PARTIAL"
    assert len(resumed["attempts"]) == 2
    assert [item["status"] for item in resumed["attempts"]] == ["FAILED", "PARTIAL"]
