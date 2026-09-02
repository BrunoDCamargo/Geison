import json

from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.pipeline import run_pipeline


def _minimal_config(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGTACGTACGT\n>s2\nACGTACGAACGT\n", encoding="utf-8")
    return PipelineConfig(target_name="target", input_fasta=fasta)


def test_successful_incomplete_fixture_is_partial_and_has_manifest(tmp_path):
    config = _minimal_config(tmp_path)
    outdir = tmp_path / "run"

    summary = run_pipeline(config, outdir)

    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    assert summary.status == "PARTIAL"
    assert manifest["status"] == "PARTIAL"
    assert "RANKING_NOT_COMPLETE" in manifest["scientific_completeness"]["missing_evidence"]
    assert (outdir / "run.log.jsonl").is_file()


def test_local_run_records_effective_config_hash_counts_and_reference(tmp_path):
    config = _minimal_config(tmp_path)
    outdir = tmp_path / "run"

    run_pipeline(config, outdir)

    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    provenance = manifest["input_provenance"]
    assert manifest["effective_config"]["target_name"] == "target"
    assert provenance["kind"] == "fasta"
    assert provenance["source_sha256"].startswith("sha256:")
    assert provenance["accepted_count"] == 2
    assert provenance["rejected_count"] == 0
    assert "records" not in provenance
    assert set(manifest["reference"]) == {"id", "mode"}
    serialized = json.dumps(manifest)
    assert "ACGTACGTACGT" not in serialized
