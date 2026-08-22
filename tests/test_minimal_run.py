import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import NcbiInputConfig, PipelineConfig
from qpcr_pipeline.ncbi import NcbiFetchedRecord, acquire_ncbi_dataset
from qpcr_pipeline.pipeline import run_pipeline


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_FASTA = REPO_ROOT / "tests" / "fixtures" / "target_small.fasta"
LOCAL_PACKAGE_ENV = {
    **os.environ,
    "PYTHONPATH": str(REPO_ROOT) + os.pathsep + os.environ.get("PYTHONPATH", ""),
}


class MinimalPipelineRunTests(unittest.TestCase):
    @staticmethod
    def _ncbi_record(accession_version: str, sequence: str) -> SeqRecord:
        record = SeqRecord(
            Seq(sequence), id=accession_version, description=f"{accession_version} record"
        )
        record.annotations["molecule_type"] = "DNA"
        return record

    @staticmethod
    def _directory_bytes(directory: Path) -> dict[Path, bytes]:
        return {
            path.relative_to(directory): path.read_bytes()
            for path in directory.rglob("*")
            if path.is_file()
        }

    def test_run_routes_acquired_ncbi_records_through_existing_qc(self):
        valid_accession = "NC_VALID.1"
        invalid_accession = "NC_INVALID.1"
        records = {
            valid_accession: self._ncbi_record(valid_accession, "ACGTACGT"),
            invalid_accession: self._ncbi_record(invalid_accession, "ACGTXCGT"),
        }

        class FakeNcbiClient:
            def resolve_query(self, query, max_records):
                raise AssertionError("query resolution was not expected")

            def fetch_records(self, identifiers, *, identifier_kind):
                if identifier_kind != "accession":
                    raise AssertionError("accession records were expected")
                return tuple(
                    NcbiFetchedRecord(request_id=identifier, record=records[identifier])
                    for identifier in reversed(identifiers)
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            summary = run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_ncbi=NcbiInputConfig(
                        accessions=(valid_accession, invalid_accession), batch_size=1
                    ),
                ),
                outdir,
                ncbi_client=FakeNcbiClient(),
            )
            qc_report = json.loads((outdir / "qc_report.json").read_text(encoding="utf-8"))
            effective_manifest = json.loads(
                (outdir / "ncbi_dataset_manifest.json").read_text(encoding="utf-8")
            )

        self.assertEqual(summary.sequence_ids, [valid_accession])
        self.assertEqual(
            qc_report["records"],
            [
                {"sequence_id": valid_accession, "status": "ACCEPTED", "reason_codes": []},
                {
                    "sequence_id": invalid_accession,
                    "status": "REJECTED",
                    "reason_codes": ["INVALID_NUCLEOTIDE"],
                },
            ],
        )
        self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], [valid_accession])
        self.assertEqual(effective_manifest["status"], "COMPLETE")

    def test_run_uses_frozen_ncbi_dataset_without_mutating_its_source(self):
        accession = "NC_FROZEN.1"
        record = self._ncbi_record(accession, "ACGTACGT")

        class FakeNcbiClient:
            def resolve_query(self, query, max_records):
                raise AssertionError("frozen datasets must not resolve queries")

            def fetch_records(self, identifiers, *, identifier_kind):
                raise AssertionError("frozen datasets must not fetch records")

        class DatasetWriter:
            def resolve_query(self, query, max_records):
                raise AssertionError("query resolution was not expected")

            def fetch_records(self, identifiers, *, identifier_kind):
                return tuple(
                    NcbiFetchedRecord(request_id=identifier, record=record)
                    for identifier in identifiers
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            frozen_dir = tmp_path / "frozen"
            acquire_ncbi_dataset(
                NcbiInputConfig(accessions=(accession,)),
                frozen_dir,
                client=DatasetWriter(),
                clock=lambda: "2026-08-21T00:00:00+00:00",
            )
            before = self._directory_bytes(frozen_dir)
            outdir = tmp_path / "run"

            summary = run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_ncbi=NcbiInputConfig(frozen_dataset=frozen_dir),
                ),
                outdir,
                ncbi_client=FakeNcbiClient(),
            )

            after = self._directory_bytes(frozen_dir)
            source_manifest = (frozen_dir / "dataset_manifest.json").read_bytes()
            effective_manifest = (outdir / "ncbi_dataset_manifest.json").read_bytes()

        self.assertEqual(summary.sequence_ids, [accession])
        self.assertEqual(after, before)
        self.assertEqual(effective_manifest, source_manifest)

    def test_run_creates_completed_summary_for_fixture(self):
        executable = shutil.which("qpcr-pipeline")
        self.assertIsNotNone(executable, "qpcr-pipeline console command is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            config_path = tmp_path / "config.yaml"
            outdir = tmp_path / "run"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                f"  fasta: {FIXTURE_FASTA.as_posix()}\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [executable, "run", str(config_path), "--outdir", str(outdir)],
                capture_output=True,
                text=True,
                check=False,
                env=LOCAL_PACKAGE_ENV,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("COMPLETED", result.stdout)

            summary_path = outdir / "run_summary.json"
            self.assertTrue(summary_path.exists(), "run_summary.json was not created")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))

        self.assertEqual(summary["status"], "COMPLETED")
        self.assertEqual(summary["target_name"], "synthetic-target")
        self.assertEqual(summary["sequence_count"], 3)
        self.assertEqual(summary["sequence_ids"], ["seq-1", "seq-2", "seq-3"])

    def test_run_writes_traceable_qc_report_for_approved_sequences(self):
        executable = shutil.which("qpcr-pipeline")
        self.assertIsNotNone(executable, "qpcr-pipeline console command is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "target.fasta"
            config_path = tmp_path / "config.yaml"
            outdir = tmp_path / "run"
            fasta_path.write_text(
                ">accepted-1\n"
                "ACGTACGT\n"
                ">invalid\n"
                "ACGTXCGT\n"
                ">accepted-2\n"
                "ACGTACGA\n",
                encoding="utf-8",
            )
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                f"  fasta: {fasta_path.as_posix()}\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [executable, "run", str(config_path), "--outdir", str(outdir)],
                capture_output=True,
                text=True,
                check=False,
                env=LOCAL_PACKAGE_ENV,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            qc_report_path = outdir / "qc_report.json"
            self.assertTrue(qc_report_path.exists(), "qc_report.json was not created")
            qc_report = json.loads(qc_report_path.read_text(encoding="utf-8"))

        self.assertEqual(summary["sequence_ids"], ["accepted-1", "accepted-2"])
        self.assertEqual(
            qc_report["records"],
            [
                {"sequence_id": "accepted-1", "status": "ACCEPTED", "reason_codes": []},
                {"sequence_id": "invalid", "status": "REJECTED", "reason_codes": ["INVALID_NUCLEOTIDE"]},
                {"sequence_id": "accepted-2", "status": "ACCEPTED", "reason_codes": []},
            ],
        )
        self.assertEqual(qc_report["target_sequence_set"]["sequence_ids"], ["accepted-1", "accepted-2"])
        self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["accepted-1", "accepted-2"])

    def test_run_applies_configured_minimum_length_qc_threshold(self):
        executable = shutil.which("qpcr-pipeline")
        self.assertIsNotNone(executable, "qpcr-pipeline console command is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "target.fasta"
            config_path = tmp_path / "config.yaml"
            outdir = tmp_path / "run"
            fasta_path.write_text(
                ">accepted\n"
                "ACGTACGT\n"
                ">too-short\n"
                "ACGT\n",
                encoding="utf-8",
            )
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                f"  fasta: {fasta_path.as_posix()}\n"
                "qc:\n"
                "  min_length: 8\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [executable, "run", str(config_path), "--outdir", str(outdir)],
                capture_output=True,
                text=True,
                check=False,
                env=LOCAL_PACKAGE_ENV,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            qc_report = json.loads((outdir / "qc_report.json").read_text(encoding="utf-8"))

        self.assertEqual(summary["sequence_ids"], ["accepted"])
        self.assertEqual(
            qc_report["records"],
            [
                {"sequence_id": "accepted", "status": "ACCEPTED", "reason_codes": []},
                {"sequence_id": "too-short", "status": "REJECTED", "reason_codes": ["TOO_SHORT"]},
            ],
        )
        self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["accepted"])


if __name__ == "__main__":
    unittest.main()
