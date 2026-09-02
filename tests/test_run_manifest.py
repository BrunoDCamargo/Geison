import json

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import NcbiInputConfig, PipelineConfig
from qpcr_pipeline.ncbi import NcbiFetchedRecord, ResolvedNcbiQuery, acquire_ncbi_dataset
from qpcr_pipeline.pipeline import run_pipeline


def _minimal_config(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGTACGTACGT\n>s2\nACGTACGAACGT\n", encoding="utf-8")
    return PipelineConfig(target_name="target", input_fasta=fasta)


class ProvenanceNcbiClient:
    def resolve_query(self, query, max_records):
        raise AssertionError("accession mode must not resolve a query")

    def fetch_records(self, identifiers, *, identifier_kind):
        assert identifier_kind == "accession"
        rows = []
        for identifier in identifiers:
            record = SeqRecord(Seq("ACGTACGTACGT"), id=identifier, name=identifier)
            record.annotations["molecule_type"] = "DNA"
            rows.append(NcbiFetchedRecord(request_id=identifier, record=record))
        return tuple(rows)


class QueryProvenanceNcbiClient:
    def resolve_query(self, query, max_records):
        assert query == "synthetic target[Title]"
        assert max_records == 1
        return ResolvedNcbiQuery(
            uids=("12345",),
            reported_count=1,
            query_translation="synthetic target[Title]",
        )

    def fetch_records(self, identifiers, *, identifier_kind):
        assert identifiers == ("12345",)
        assert identifier_kind == "uid"
        record = SeqRecord(
            Seq("ACGTACGTACGT"),
            id="NC_000002.2",
            name="NC_000002",
        )
        record.annotations["molecule_type"] = "DNA"
        return (NcbiFetchedRecord(request_id="12345", record=record),)


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


def test_ncbi_run_records_request_and_resolved_accession_versions(tmp_path):
    config = PipelineConfig(
        target_name="target",
        input_ncbi=NcbiInputConfig(accessions=("NC_000001.1",)),
    )
    outdir = tmp_path / "run"

    run_pipeline(config, outdir, ncbi_client=ProvenanceNcbiClient())

    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    provenance = manifest["input_provenance"]
    assert provenance["kind"] == "ncbi"
    assert provenance["mode"] == "accessions"
    assert provenance["requested_accessions"] == ["NC_000001.1"]
    assert provenance["resolved_accession_versions"] == ["NC_000001.1"]
    assert provenance["dataset_sha256"].startswith("sha256:")
    serialized = json.dumps(provenance)
    assert "completed_batches" not in serialized
    assert "record_ids" not in serialized
    assert "NCBI_API_KEY" not in serialized


def test_query_ncbi_run_records_query_without_internal_resolution_details(tmp_path):
    config = PipelineConfig(
        target_name="target",
        input_ncbi=NcbiInputConfig(
            query="synthetic target[Title]",
            max_records=1,
            retries=0,
        ),
    )
    outdir = tmp_path / "run"

    run_pipeline(config, outdir, ncbi_client=QueryProvenanceNcbiClient())

    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    provenance = manifest["input_provenance"]
    assert provenance["kind"] == "ncbi"
    assert provenance["mode"] == "query"
    assert provenance["query"] == "synthetic target[Title]"
    assert provenance["resolved_accession_versions"] == ["NC_000002.2"]
    assert provenance["dataset_sha256"].startswith("sha256:")
    serialized = json.dumps(provenance)
    assert "resolved_uids" not in serialized
    assert "query_translation" not in serialized
    assert "completed_batches" not in serialized
    assert "record_ids" not in serialized


def test_frozen_ncbi_run_records_dataset_identity_without_internal_batches(tmp_path):
    frozen = tmp_path / "frozen"
    acquire_ncbi_dataset(
        NcbiInputConfig(accessions=("NC_000001.1",)),
        frozen,
        client=ProvenanceNcbiClient(),
    )
    config = PipelineConfig(
        target_name="target",
        input_ncbi=NcbiInputConfig(frozen_dataset=frozen),
    )
    outdir = tmp_path / "run"

    run_pipeline(config, outdir)

    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    provenance = manifest["input_provenance"]
    assert provenance["kind"] == "ncbi"
    assert provenance["mode"] == "frozen_dataset"
    assert provenance["source_dataset_mode"] == "accessions"
    assert provenance["dataset_sha256"].startswith("sha256:")
    assert provenance["source_manifest_sha256"].startswith("sha256:")
    assert provenance["resolved_accession_versions"] == ["NC_000001.1"]
    serialized = json.dumps(provenance)
    assert "completed_batches" not in serialized
    assert "record_ids" not in serialized
