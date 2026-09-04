from pathlib import Path

from qpcr_pipeline.contrastive_conservation import (
    ContrastiveConservationResult,
)


def test_contrastive_codec_round_trips_skipped_result(tmp_path):
    from qpcr_pipeline.checkpoint_codecs import CONTRASTIVE_CONSERVATION_CODEC

    outdir = tmp_path / "out"
    report = outdir / "contrastive_conservation" / "contrastive_conservation_report.json"
    result = ContrastiveConservationResult(
        status="SKIPPED",
        reference_id=None,
        windows=(),
        dataset_evidence=(),
        candidates=(),
        challenge_datasets=(),
        window_metrics_path=None,
        dataset_metrics_path=None,
        candidate_regions_path=None,
        report_path=report,
        html_report_path=None,
    )
    payload = CONTRASTIVE_CONSERVATION_CODEC.encode(result, outdir)
    decoded = CONTRASTIVE_CONSERVATION_CODEC.decode(payload, outdir)
    assert decoded == result


def test_contrastive_stage_definition_and_pipeline_import():
    from qpcr_pipeline.checkpoint_stages import STAGE_DEFINITIONS
    import qpcr_pipeline.pipeline as pipeline

    definition = STAGE_DEFINITIONS["contrastive_conservation"]
    assert definition.dependencies == ("panel", "conservation")
    assert pipeline._run_stage is not None


def test_contrastive_stage_declares_all_published_artifacts(tmp_path):
    from qpcr_pipeline.checkpoint_stages import stage_outputs

    outdir = tmp_path / "out"
    stage = outdir / "contrastive_conservation"
    result = ContrastiveConservationResult(
        status="COMPLETE",
        reference_id="ref",
        windows=(),
        dataset_evidence=(),
        candidates=(),
        challenge_datasets=(),
        window_metrics_path=stage / "window_metrics.tsv",
        dataset_metrics_path=stage / "dataset_metrics.tsv",
        candidate_regions_path=stage / "candidate_regions.tsv",
        report_path=stage / "contrastive_conservation_report.json",
        html_report_path=stage / "report.html",
    )
    assert stage_outputs("contrastive_conservation", result, outdir) == (
        stage / "window_metrics.tsv",
        stage / "dataset_metrics.tsv",
        stage / "candidate_regions.tsv",
        stage / "contrastive_conservation_report.json",
        stage / "report.html",
    )
