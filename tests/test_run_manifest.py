import json

from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.pipeline import run_pipeline


def test_successful_incomplete_fixture_is_partial_and_has_manifest(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGTACGTACGT\n>s2\nACGTACGAACGT\n", encoding="utf-8")
    config = PipelineConfig(target_name="target", input_fasta=fasta)
    outdir = tmp_path / "run"

    summary = run_pipeline(config, outdir)

    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    assert summary.status == "PARTIAL"
    assert manifest["status"] == "PARTIAL"
    assert "RANKING_NOT_COMPLETE" in manifest["scientific_completeness"]["missing_evidence"]
    assert (outdir / "run.log.jsonl").is_file()
