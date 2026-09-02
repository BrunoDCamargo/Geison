from pathlib import Path

from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.execution import ExecutionPolicy, STAGE_ORDER
from qpcr_pipeline.pipeline import run_pipeline
from qpcr_pipeline.planning import plan_pipeline


def _config(tmp_path: Path) -> PipelineConfig:
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGTACGTACGT\n>s2\nACGTACGAACGT\n", encoding="utf-8")
    return PipelineConfig(target_name="target", input_fasta=fasta)


def test_plan_without_outdir_runs_every_stage(tmp_path):
    config = _config(tmp_path)

    plan = plan_pipeline(config, None)

    assert [(item.stage, item.action) for item in plan.decisions] == [
        (stage, "RUN") for stage in STAGE_ORDER
    ]
    assert plan.reused_results == {}


def test_resume_plan_reuses_every_valid_stage(tmp_path):
    config = _config(tmp_path)
    outdir = tmp_path / "run"
    run_pipeline(config, outdir)

    plan = plan_pipeline(config, outdir, execution=ExecutionPolicy(resume=True))

    assert [item.action for item in plan.decisions] == ["REUSE"] * len(STAGE_ORDER)
    assert set(plan.reused_results) == set(STAGE_ORDER)
