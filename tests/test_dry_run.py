from pathlib import Path

from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.diagnostics import ComponentReport, EnvironmentReport
from qpcr_pipeline.dry_run import dry_run_pipeline


class FakeInspector:
    def __init__(self, report):
        self.report = report
        self.calls = []

    def inspect(self, config):
        self.calls.append(config)
        return self.report


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
