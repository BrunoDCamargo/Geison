import json
import tempfile
import unittest
from pathlib import Path
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
from qpcr_pipeline.models import EvaluationSet, TargetSequenceSet
from qpcr_pipeline.pipeline import run_pipeline
from qpcr_pipeline.qc import QCResult
from panel_fixtures import proposal_panel_config
from qpcr_pipeline.ranking_guard import (
    evaluate_ranking_with_execution_guard as real_evaluate_ranking,
)
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


class PipelineRankingTests(unittest.TestCase):
    def _enabled_config(self, off_target_path: Path) -> PipelineConfig:
        return PipelineConfig(
            target_name="target",
            input_fasta=FIXTURE_FASTA,
            panel=proposal_panel_config("target"),
            alignment=AlignmentConfig(enabled=True),
            conservation=ConservationConfig(enabled=True),
            primer_design=PrimerDesignConfig(enabled=True),
            inclusivity=InclusivityConfig(enabled=True),
            off_targets=(OffTargetConfig(name="off", fasta=off_target_path),),
            specificity=SpecificityConfig(enabled=True),
            ranking=RankingConfig(enabled=True),
        )

    @staticmethod
    def _base_upstream_patches(output, primer, inclusivity, specificity):
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

    def test_default_pipeline_publishes_skipped_ranking_summary_and_report(self):
        config = PipelineConfig(target_name="target", input_fasta=FIXTURE_FASTA)
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "out"
            run_pipeline(config, output)
            qc = json.loads((output / "qc_report.json").read_text(encoding="utf-8"))
            ranking_report = json.loads(
                (output / "ranking" / "ranking_report.json").read_text(encoding="utf-8")
            )
        self.assertEqual(
            qc["ranking"],
            {
                "status": "SKIPPED",
                "assay_count": 0,
                "in_silico_pass_count": 0,
                "review_count": 0,
                "high_risk_count": 0,
                "complete_score_count": 0,
                "incomplete_score_count": 0,
                "top_recommended_assay_id": None,
            },
        )
        self.assertEqual(ranking_report["status"], "SKIPPED")

    def test_enabled_pipeline_publishes_pass_ranking_and_top_recommendation(self):
        primer = make_primer_result()
        inclusivity = make_inclusivity_result(primer)
        specificity = make_specificity_result(primer)
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "out"
            off_target = root / "off.fa"
            off_target.write_text(">off\nACGT\n", encoding="utf-8")
            patches = self._base_upstream_patches(output, primer, inclusivity, specificity)
            with patches[0], patches[1], patches[2], patches[3], patches[4], patches[5]:
                summary = run_pipeline(self._enabled_config(off_target), output)
            qc = json.loads((output / "qc_report.json").read_text(encoding="utf-8"))
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(summary.status, "COMPLETED")
            self.assertEqual(manifest["status"], "COMPLETED")
            self.assertEqual(manifest["scientific_completeness"]["missing_evidence"], [])
            self.assertEqual(qc["ranking"]["status"], "COMPLETE")
            self.assertEqual(qc["ranking"]["assay_count"], 1)
            self.assertEqual(qc["ranking"]["in_silico_pass_count"], 1)
            self.assertEqual(qc["ranking"]["review_count"], 0)
            self.assertEqual(qc["ranking"]["high_risk_count"], 0)
            self.assertEqual(qc["ranking"]["complete_score_count"], 1)
            self.assertEqual(qc["ranking"]["incomplete_score_count"], 0)
            self.assertEqual(qc["ranking"]["top_recommended_assay_id"], "a1")
            self.assertTrue((output / "ranking" / "assay_ranking.tsv").exists())
            self.assertTrue((output / "ranking" / "ranking_report.json").exists())
            self.assertTrue((output / "report.html").exists())

    def test_incomplete_run_evidence_blocks_pass_before_ranking_artifacts(self):
        primer = make_primer_result()
        inclusivity = make_inclusivity_result(primer)
        specificity = make_specificity_result(primer)
        empty_qc = QCResult(
            records=(),
            target_sequence_set=TargetSequenceSet(sequence_ids=()),
            evaluation_set=EvaluationSet(sequence_ids=()),
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "out"
            off_target = root / "off.fa"
            off_target.write_text(">off\nACGT\n", encoding="utf-8")
            patches = self._base_upstream_patches(output, primer, inclusivity, specificity)
            with (
                patch("qpcr_pipeline.pipeline.evaluate_sequences", return_value=empty_qc),
                patches[0],
                patches[1],
                patches[2],
                patches[3],
                patches[4],
                patches[5],
            ):
                summary = run_pipeline(self._enabled_config(off_target), output)
            ranking_report = json.loads(
                (output / "ranking" / "ranking_report.json").read_text(encoding="utf-8")
            )
            self.assertEqual(summary.status, "PARTIAL")
            self.assertEqual(ranking_report["counts"]["in_silico_pass"], 0)
            self.assertEqual(ranking_report["counts"]["review"], 1)
            self.assertEqual(ranking_report["counts"]["incomplete_score"], 1)

    def test_review_only_pipeline_never_promotes_review_to_recommended(self):
        primer = make_primer_result()
        sequence_ids = tuple(f"s{i}" for i in range(10))
        compatibility = {
            ("a1", sequence_id): index < 9
            for index, sequence_id in enumerate(sequence_ids)
        }
        inclusivity = make_inclusivity_result(
            primer,
            compatibility=compatibility,
            sequence_ids=sequence_ids,
        )
        specificity = make_specificity_result(primer)
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "out"
            off_target = root / "off.fa"
            off_target.write_text(">off\nACGT\n", encoding="utf-8")
            patches = self._base_upstream_patches(output, primer, inclusivity, specificity)
            with patches[0], patches[1], patches[2], patches[3], patches[4], patches[5]:
                run_pipeline(self._enabled_config(off_target), output)
            qc = json.loads((output / "qc_report.json").read_text(encoding="utf-8"))
        self.assertEqual(qc["ranking"]["review_count"], 1)
        self.assertEqual(qc["ranking"]["in_silico_pass_count"], 0)
        self.assertIsNone(qc["ranking"]["top_recommended_assay_id"])

    def test_specificity_completes_before_ranking_is_invoked(self):
        primer = make_primer_result()
        inclusivity = make_inclusivity_result(primer)
        specificity = make_specificity_result(primer)
        events = []

        def ranking_side_effect(*args, **kwargs):
            events.append("ranking")
            return real_evaluate_ranking(*args, **kwargs)

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "out"
            off_target = root / "off.fa"
            off_target.write_text(">off\nACGT\n", encoding="utf-8")
            patches = self._base_upstream_patches(output, primer, inclusivity, specificity)

            def specificity_side_effect(*args, **kwargs):
                events.append("specificity")
                return checkpoint_specificity(output, specificity)

            with (
                patches[0],
                patches[1],
                patches[2],
                patches[3],
                patches[4],
                patch(
                    "qpcr_pipeline.pipeline.evaluate_specificity",
                    side_effect=specificity_side_effect,
                ),
                patch(
                    "qpcr_pipeline.pipeline.evaluate_ranking_with_execution_guard",
                    side_effect=ranking_side_effect,
                ),
            ):
                run_pipeline(self._enabled_config(off_target), output)
        self.assertEqual(events, ["specificity", "ranking"])


if __name__ == "__main__":
    unittest.main()
