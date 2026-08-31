import json
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest

from qpcr_pipeline.config import (
    AlignmentConfig,
    ClusteringConfig,
    PipelineConfig,
    SpecificityConfig,
)
from qpcr_pipeline.execution import ExecutionPolicy, STAGE_ORDER
from qpcr_pipeline.pipeline import run_pipeline
from pipeline_checkpoint_fixtures import checkpoint_alignment


def _write_fasta(path: Path):
    path.write_text(">s1\nACGTACGTACGT\n>s2\nACGTACGAACGT\n", encoding="utf-8")


def _manifest(outdir: Path, stage: str):
    return json.loads(
        (outdir / ".checkpoints" / stage / "manifest.json").read_text(encoding="utf-8")
    )


def _actions(outdir: Path):
    payload = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
    return {item["stage"]: item["action"] for item in payload["stage_actions"]}


def _default_run(tmp_path: Path):
    fasta = tmp_path / "target.fasta"
    outdir = tmp_path / "run"
    _write_fasta(fasta)
    config = PipelineConfig(target_name="target", input_fasta=fasta)
    run_pipeline(config, outdir)
    return config, outdir


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
    config, outdir = _default_run(tmp_path)
    resumed = run_pipeline(config, outdir, execution=ExecutionPolicy(resume=True))
    payload = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))

    assert resumed.status == "COMPLETED"
    assert resumed.sequence_ids == ["s1", "s2"]
    assert [item["stage"] for item in payload["stage_actions"]] == list(STAGE_ORDER)
    assert [item["action"] for item in payload["stage_actions"]] == ["REUSE"] * len(STAGE_ORDER)


def test_specificity_parameter_change_does_not_recalculate_alignment_or_conservation(tmp_path):
    config, outdir = _default_run(tmp_path)
    changed = replace(
        config,
        specificity=SpecificityConfig(max_hits_per_oligo_per_dataset=7),
    )

    run_pipeline(changed, outdir, execution=ExecutionPolicy(resume=True))

    assert _actions(outdir) == {
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "REUSE",
        "alignment": "REUSE",
        "conservation": "REUSE",
        "primer_design": "REUSE",
        "inclusivity": "REUSE",
        "specificity": "RUN",
        "ranking": "FORCED",
    }


def test_clustering_parameter_change_invalidates_complete_dependent_chain(tmp_path):
    config, outdir = _default_run(tmp_path)
    changed = replace(config, clustering=ClusteringConfig(identity=0.90))

    run_pipeline(changed, outdir, execution=ExecutionPolicy(resume=True))

    assert _actions(outdir) == {
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "RUN",
        "alignment": "FORCED",
        "conservation": "FORCED",
        "primer_design": "FORCED",
        "inclusivity": "FORCED",
        "specificity": "FORCED",
        "ranking": "FORCED",
    }


def test_force_inclusivity_keeps_independent_specificity_reusable(tmp_path):
    config, outdir = _default_run(tmp_path)

    run_pipeline(
        config,
        outdir,
        execution=ExecutionPolicy(resume=True, force_step="inclusivity"),
    )

    assert _actions(outdir) == {
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "REUSE",
        "alignment": "REUSE",
        "conservation": "REUSE",
        "primer_design": "REUSE",
        "inclusivity": "FORCED",
        "specificity": "REUSE",
        "ranking": "FORCED",
    }


def test_valid_from_specificity_reuses_boundary_and_forces_selected_subgraph(tmp_path):
    config, outdir = _default_run(tmp_path)

    run_pipeline(
        config,
        outdir,
        execution=ExecutionPolicy(from_step="specificity"),
    )

    assert _actions(outdir) == {
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "REUSE",
        "alignment": "REUSE",
        "conservation": "REUSE",
        "primer_design": "REUSE",
        "inclusivity": "REUSE",
        "specificity": "FORCED",
        "ranking": "FORCED",
    }


def test_from_specificity_fails_before_scientific_write_when_boundary_is_invalid(tmp_path):
    config, outdir = _default_run(tmp_path)
    specificity_before = (outdir / "specificity" / "specificity_report.json").read_bytes()
    ranking_before = (outdir / "ranking" / "ranking_report.json").read_bytes()
    state_path = outdir / ".checkpoints" / "inclusivity" / "state.json"
    state_path.write_text('{"corrupt":true}\n', encoding="utf-8")

    with pytest.raises(ValueError) as error:
        run_pipeline(
            config,
            outdir,
            execution=ExecutionPolicy(from_step="specificity"),
        )

    message = str(error.value)
    assert "inclusivity" in message
    assert "STATE_HASH_MISMATCH" in message
    assert "--resume" in message
    assert (outdir / "specificity" / "specificity_report.json").read_bytes() == specificity_before
    assert (outdir / "ranking" / "ranking_report.json").read_bytes() == ranking_before


def test_corrupt_alignment_output_reruns_alignment_and_descendants_only(tmp_path):
    config, outdir = _default_run(tmp_path)
    alignment_report = outdir / "alignment" / "alignment_report.json"
    alignment_report.write_text('{"corrupt":true}\n', encoding="utf-8")

    run_pipeline(config, outdir, execution=ExecutionPolicy(resume=True))

    assert _actions(outdir) == {
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "REUSE",
        "alignment": "RUN",
        "conservation": "FORCED",
        "primer_design": "FORCED",
        "inclusivity": "FORCED",
        "specificity": "FORCED",
        "ranking": "FORCED",
    }


def test_interrupted_run_preserves_completed_upstream_and_resume_finishes(tmp_path):
    fasta = tmp_path / "target.fasta"
    outdir = tmp_path / "run"
    _write_fasta(fasta)
    config = PipelineConfig(target_name="target", input_fasta=fasta)

    with patch(
        "qpcr_pipeline.pipeline.evaluate_specificity",
        side_effect=RuntimeError("simulated interruption"),
    ):
        with pytest.raises(RuntimeError, match="simulated interruption"):
            run_pipeline(config, outdir)

    for stage in STAGE_ORDER[:7]:
        assert _manifest(outdir, stage)["status"] == "COMPLETE"
    assert _manifest(outdir, "specificity")["status"] == "FAILED"
    assert not (outdir / ".checkpoints" / "ranking" / "manifest.json").exists()

    run_pipeline(config, outdir, execution=ExecutionPolicy(resume=True))

    assert _actions(outdir) == {
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "REUSE",
        "alignment": "REUSE",
        "conservation": "REUSE",
        "primer_design": "REUSE",
        "inclusivity": "REUSE",
        "specificity": "RUN",
        "ranking": "FORCED",
    }
    assert _manifest(outdir, "specificity")["status"] == "COMPLETE"
    assert _manifest(outdir, "ranking")["status"] == "COMPLETE"


def test_geison_version_change_invalidates_entire_chain(tmp_path):
    config, outdir = _default_run(tmp_path)

    with patch("qpcr_pipeline.checkpoint_stages.geison_version", return_value="999.0.0"):
        run_pipeline(config, outdir, execution=ExecutionPolicy(resume=True))

    assert _actions(outdir)["input"] == "RUN"
    assert all(
        _actions(outdir)[stage] == "FORCED"
        for stage in STAGE_ORDER[1:]
    )


class _ToolProvider:
    def __init__(self, mafft_version: str):
        self.mafft_version = mafft_version

    def identity(self, tool_name: str):
        assert tool_name == "mafft"
        return {"name": "mafft", "version": self.mafft_version}


def test_mafft_version_change_invalidates_alignment_and_descendants_only(tmp_path):
    fasta = tmp_path / "target.fasta"
    outdir = tmp_path / "run"
    _write_fasta(fasta)
    config = PipelineConfig(
        target_name="target",
        input_fasta=fasta,
        alignment=AlignmentConfig(enabled=True),
    )

    with patch(
        "qpcr_pipeline.pipeline.align_discovery",
        side_effect=lambda *args, **kwargs: checkpoint_alignment(outdir),
    ):
        run_pipeline(
            config,
            outdir,
            tool_identity_provider=_ToolProvider("7.1"),
        )

    with patch(
        "qpcr_pipeline.pipeline.align_discovery",
        side_effect=lambda *args, **kwargs: checkpoint_alignment(outdir),
    ):
        run_pipeline(
            config,
            outdir,
            execution=ExecutionPolicy(resume=True),
            tool_identity_provider=_ToolProvider("7.2"),
        )

    assert _actions(outdir) == {
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "REUSE",
        "alignment": "RUN",
        "conservation": "FORCED",
        "primer_design": "FORCED",
        "inclusivity": "FORCED",
        "specificity": "FORCED",
        "ranking": "FORCED",
    }
