from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from qpcr_pipeline.config import (
    AlignmentConfig,
    ConservationConfig,
    InclusivityConfig,
    OffTargetConfig,
    PipelineConfig,
    PrimerDesignConfig,
    RankingConfig,
    SpecificityConfig,
)
from qpcr_pipeline.pipeline import run_pipeline
from qpcr_pipeline.run_recording import RunRecorder
from panel_fixtures import approved_panel_config
from pipeline_checkpoint_fixtures import (
    checkpoint_alignment,
    checkpoint_clustering,
    checkpoint_conservation,
    checkpoint_inclusivity,
    checkpoint_primer,
    checkpoint_specificity,
)
from ranking_fixtures import (
    make_inclusivity_result,
    make_primer_result,
    make_specificity_result,
)


FIXTURE_FASTA = Path("tests/fixtures/target_small.fasta")


class RecordingRunRecorder(RunRecorder):
    def __init__(self, output: Path, events: list[str]) -> None:
        super().__init__(output)
        self.events = events

    def complete(self, *args, **kwargs):
        self.events.append("complete")
        return super().complete(*args, **kwargs)


class PipelineResearcherReportTests(unittest.TestCase):
    def _config(self, root: Path) -> PipelineConfig:
        off_target = root / "off.fa"
        off_target.write_text(">off\nACGT\n", encoding="utf-8")
        return PipelineConfig(
            target_name="target",
            input_fasta=FIXTURE_FASTA,
            panel=approved_panel_config(root, "target"),
            alignment=AlignmentConfig(enabled=True),
            conservation=ConservationConfig(enabled=True),
            primer_design=PrimerDesignConfig(enabled=True),
            inclusivity=InclusivityConfig(enabled=True),
            off_targets=(OffTargetConfig(name="off", fasta=off_target),),
            specificity=SpecificityConfig(enabled=True),
            ranking=RankingConfig(enabled=True),
        )

    @staticmethod
    def _patches(output: Path):
        primer = make_primer_result()
        inclusivity = make_inclusivity_result(primer)
        specificity = make_specificity_result(primer)
        return (
            patch(
                "qpcr_pipeline.pipeline.cluster_sequences",
                side_effect=lambda *args, **kwargs: checkpoint_clustering(output),
            ),
            patch(
                "qpcr_pipeline.pipeline.align_discovery",
                side_effect=lambda *args, **kwargs: checkpoint_alignment(output),
            ),
            patch(
                "qpcr_pipeline.pipeline.analyze_conservation",
                side_effect=lambda *args, **kwargs: checkpoint_conservation(output),
            ),
            patch(
                "qpcr_pipeline.pipeline.design_primers",
                side_effect=lambda *args, **kwargs: checkpoint_primer(output, primer),
            ),
            patch(
                "qpcr_pipeline.pipeline.evaluate_inclusivity",
                side_effect=lambda *args, **kwargs: checkpoint_inclusivity(output, inclusivity),
            ),
            patch(
                "qpcr_pipeline.pipeline.evaluate_specificity",
                side_effect=lambda *args, **kwargs: checkpoint_specificity(output, specificity),
            ),
        )

    def test_researcher_report_runs_after_recorder_has_persisted_final_status(self) -> None:
        events: list[str] = []
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "out"
            patches = self._patches(output)

            def generate_report(path: Path) -> Path:
                manifest = json.loads((path / "run_manifest.json").read_text(encoding="utf-8"))
                self.assertEqual(manifest["status"], "COMPLETED")
                events.append("report")
                report = path / "report.html"
                report.write_text("final researcher report", encoding="utf-8")
                return report

            with (
                patches[0], patches[1], patches[2], patches[3], patches[4], patches[5],
                patch("qpcr_pipeline.pipeline.generate_researcher_report", side_effect=generate_report),
            ):
                summary = run_pipeline(
                    self._config(root),
                    output,
                    recorder_factory=lambda path: RecordingRunRecorder(path, events),
                )

            self.assertEqual(summary.status, "COMPLETED")
            self.assertEqual(events[-2:], ["complete", "report"])
            self.assertEqual((output / "report.html").read_text(encoding="utf-8"), "final researcher report")

    def test_report_rendering_failure_does_not_change_scientific_run_status(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "out"
            patches = self._patches(output)
            with (
                patches[0], patches[1], patches[2], patches[3], patches[4], patches[5],
                patch(
                    "qpcr_pipeline.pipeline.generate_researcher_report",
                    side_effect=RuntimeError("render boom"),
                ),
            ):
                summary = run_pipeline(self._config(root), output)

            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            report_error = json.loads((output / "report_error.json").read_text(encoding="utf-8"))

        self.assertEqual(summary.status, "COMPLETED")
        self.assertEqual(manifest["status"], "COMPLETED")
        self.assertEqual(report_error["error_type"], "RuntimeError")
        self.assertIn("render boom", report_error["message"])


if __name__ == "__main__":
    unittest.main()
