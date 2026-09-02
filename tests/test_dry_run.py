from pathlib import Path
from unittest.mock import patch

from qpcr_pipeline.config import AlignmentConfig, PipelineConfig
from qpcr_pipeline.diagnostics import ComponentReport, EnvironmentReport
from qpcr_pipeline.dry_run import dry_run_pipeline
from qpcr_pipeline.execution import ExecutionPolicy
from qpcr_pipeline.pipeline import run_pipeline
from pipeline_checkpoint_fixtures import checkpoint_alignment


class FakeInspector:
    def __init__(self, report):
        self.report = report
        self.calls = []

    def inspect(self, config):
        self.calls.append(config)
        return self.report


class StaticToolIdentityProvider:
    def identity(self, tool_name):
        assert tool_name == "mafft"
        return {"name": "mafft", "version": "7.526"}


def _environment() -> EnvironmentReport:
    python = ComponentReport("Python", "USED", True, True, "3.12")
    geison = ComponentReport("Geison", "USED", True, True, "0.1")
    git = ComponentReport("Git", "UNAVAILABLE", False, False, None)
    tools = {
        "cd-hit-est": ComponentReport("cd-hit-est", "UNAVAILABLE", False, False, None),
        "mafft": ComponentReport("mafft", "UNAVAILABLE", False, False, None),
        "primer3_core": ComponentReport("primer3_core", "UNAVAILABLE", False, False, None),
        "blast+": ComponentReport("blast+", "NOT_USED", False, False, None),
    }
    return EnvironmentReport(python, geison, git, tools)


def _environment_with_mafft() -> EnvironmentReport:
    report = _environment()
    tools = dict(report.tools)
    tools["mafft"] = ComponentReport("mafft", "AVAILABLE", True, True, "7.526")
    return EnvironmentReport(report.python, report.geison, report.git, tools)


def _config(tmp_path: Path) -> PipelineConfig:
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGTACGTACGT\n>s2\nACGTACGAACGT\n", encoding="utf-8")
    return PipelineConfig(target_name="target", input_fasta=fasta)


def test_dry_run_does_not_create_output_directory(tmp_path):
    config = _config(tmp_path)
    outdir = tmp_path / "absent"
    inspector = FakeInspector(_environment())

    report = dry_run_pipeline(config, outdir, inspector=inspector)

    assert report.target_name == "target"
    assert all(item.action == "RUN" for item in report.decisions)
    assert not outdir.exists()
    assert inspector.calls == [config]


def test_dry_run_does_not_modify_existing_output(tmp_path):
    config = _config(tmp_path)
    outdir = tmp_path / "existing"
    outdir.mkdir()
    sentinel = outdir / "sentinel.bin"
    sentinel.write_bytes(b"unchanged")
    before = {path.relative_to(outdir): path.read_bytes() for path in outdir.rglob("*") if path.is_file()}

    dry_run_pipeline(config, outdir, inspector=FakeInspector(_environment()))

    after = {path.relative_to(outdir): path.read_bytes() for path in outdir.rglob("*") if path.is_file()}
    assert after == before


def test_resume_dry_run_reuses_environment_tool_identity_without_probing_binaries(tmp_path):
    config = PipelineConfig(
        target_name="target",
        input_fasta=_config(tmp_path).input_fasta,
        alignment=AlignmentConfig(enabled=True),
    )
    outdir = tmp_path / "run"

    with patch(
        "qpcr_pipeline.pipeline.align_discovery",
        side_effect=lambda *args, **kwargs: checkpoint_alignment(outdir),
    ):
        run_pipeline(
            config,
            outdir,
            tool_identity_provider=StaticToolIdentityProvider(),
        )

    inspector = FakeInspector(_environment_with_mafft())
    with patch(
        "qpcr_pipeline.checkpoint_stages.SubprocessToolIdentityProvider.identity",
        side_effect=AssertionError("dry-run must not probe scientific binaries during planning"),
    ) as identity:
        report = dry_run_pipeline(
            config,
            outdir,
            execution=ExecutionPolicy(resume=True),
            inspector=inspector,
        )

    identity.assert_not_called()
    actions = {item.stage: item.action for item in report.decisions}
    assert actions["alignment"] == "REUSE"
