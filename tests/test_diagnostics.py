from qpcr_pipeline.config import AlignmentConfig, ClusteringConfig, PipelineConfig, PrimerDesignConfig
from qpcr_pipeline.diagnostics import (
    CommandResult,
    ComponentReport,
    EnvironmentInspector,
    EnvironmentReport,
    doctor_exit_code,
    render_environment_report,
)


class FakeRunner:
    def __init__(self, responses):
        self.responses = responses
        self.calls = []

    def run(self, argv):
        self.calls.append(argv)
        return self.responses.get(argv, CommandResult(127, "", "not found"))


def test_inspector_marks_enabled_tools_required_and_blast_not_used(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    config = PipelineConfig(
        target_name="target",
        input_fasta=fasta,
        clustering=ClusteringConfig(enabled=True),
        alignment=AlignmentConfig(enabled=True),
        primer_design=PrimerDesignConfig(enabled=True),
    )
    runner = FakeRunner({
        ("cd-hit-est", "-h"): CommandResult(0, "CD-HIT version 4.8.1", ""),
        ("mafft", "--version"): CommandResult(0, "v7.526", ""),
        ("primer3_core", "--version"): CommandResult(0, "primer3 release 2.6.1", ""),
        ("blastn", "-version"): CommandResult(127, "", "not found"),
    })

    report = EnvironmentInspector(runner=runner).inspect(config)

    assert report.tools["cd-hit-est"].required is True
    assert report.tools["cd-hit-est"].version == "4.8.1"
    assert report.tools["mafft"].required is True
    assert report.tools["mafft"].version == "v7.526"
    assert report.tools["primer3_core"].required is True
    assert report.tools["primer3_core"].version == "primer3 release 2.6.1"
    assert report.tools["blast+"].status == "NOT_USED"
    assert report.tools["blast+"].required is False
    assert report.tools["blast+"].installed is False


def test_inspector_reports_missing_required_tool(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    config = PipelineConfig(
        target_name="target",
        input_fasta=fasta,
        alignment=AlignmentConfig(enabled=True),
    )
    report = EnvironmentInspector(runner=FakeRunner({})).inspect(config)
    assert report.tools["mafft"].required is True
    assert report.tools["mafft"].installed is False
    assert report.missing_required_tools == ("mafft",)


def test_doctor_context_treats_external_tools_as_optional():
    report = EnvironmentInspector(runner=FakeRunner({})).inspect(None)
    assert report.missing_required_tools == ()
    assert all(not item.required for item in report.tools.values())


def test_doctor_rendering_names_missing_optional_tools_without_failure():
    missing = ComponentReport("mafft", "UNAVAILABLE", False, False, None)
    used = ComponentReport("Python", "USED", True, True, "3.12")
    report = EnvironmentReport(
        used,
        used,
        ComponentReport("Git", "UNAVAILABLE", False, False, None),
        {"mafft": missing},
    )

    rendered = render_environment_report(report)

    assert "mafft" in rendered
    assert "UNAVAILABLE" in rendered
    assert doctor_exit_code(report) == 0
