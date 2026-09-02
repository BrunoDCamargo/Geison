import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.alignment import MafftError
from qpcr_pipeline.config import (
    AlignmentConfig,
    ClusteringConfig,
    ConservationConfig,
    InclusivityConfig,
    NcbiInputConfig,
    PipelineConfig,
    PrimerDesignConfig,
)
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

    def test_enabled_clustering_publishes_discovery_set_from_approved_records(self):
        class FakeCdHitRunner:
            def __init__(self):
                self.input_records = []

            def run(self, input_path, output_path, config):
                with Path(input_path).open(encoding="utf-8") as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
                self.input_records = [
                    (record.id, str(record.seq)) for record in records
                ]
                records_by_id = {record.id: record for record in records}
                SeqIO.write(
                    [
                        records_by_id["geison-00000000"],
                        records_by_id["geison-00000002"],
                    ],
                    output_path,
                    "fasta",
                )
                Path(str(output_path) + ".clstr").write_text(
                    ">Cluster 0\n"
                    "0 12nt, >geison-00000000... *\n"
                    "1 12nt, >geison-00000001... at +/99.00%\n"
                    ">Cluster 1\n"
                    "0 12nt, >geison-00000002... *\n",
                    encoding="utf-8",
                )

        runner = FakeCdHitRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "target.fasta"
            outdir = tmp_path / "run"
            fasta_path.write_text(
                ">rejected\n"
                "ACGTXCGTACGT\n"
                ">s1\n"
                "ACGTACGTACGT\n"
                ">s2\n"
                "ACGTACGAACGT\n"
                ">s3\n"
                "ACGTACCCACGT\n",
                encoding="utf-8",
            )

            summary = run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_fasta=fasta_path,
                    clustering=ClusteringConfig(enabled=True),
                ),
                outdir,
                cdhit_runner=runner,
            )
            qc_report = json.loads(
                (outdir / "qc_report.json").read_text(encoding="utf-8")
            )
            clustering_report = json.loads(
                (outdir / "clustering_report.json").read_text(encoding="utf-8")
            )

        self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["s1", "s2", "s3"])
        self.assertEqual(qc_report["discovery_set"]["sequence_ids"], ["s1", "s3"])
        self.assertEqual(summary.sequence_ids, ["s1", "s2", "s3"])
        self.assertEqual(clustering_report["counts"], {"evaluation": 3, "discovery": 2})
        self.assertEqual(
            runner.input_records,
            [
                ("geison-00000000", "ACGTACGTACGT"),
                ("geison-00000001", "ACGTACGAACGT"),
                ("geison-00000002", "ACGTACCCACGT"),
            ],
        )

    def test_enabled_empty_clustering_publishes_empty_raw_cluster_artifact(self):
        class FailingCdHitRunner:
            def __init__(self):
                self.calls = []

            def run(self, input_path, output_path, config):
                self.calls.append((input_path, output_path, config))
                raise AssertionError("empty clustering must not call the runner")

        runner = FailingCdHitRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "target.fasta"
            outdir = tmp_path / "run"
            fasta_path.write_text(
                ">rejected-1\n"
                "ACGTXCGT\n"
                ">rejected-2\n"
                "ACGTXCGC\n",
                encoding="utf-8",
            )

            summary = run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_fasta=fasta_path,
                    clustering=ClusteringConfig(enabled=True),
                ),
                outdir,
                cdhit_runner=runner,
            )
            qc_report = json.loads(
                (outdir / "qc_report.json").read_text(encoding="utf-8")
            )
            clustering_report = json.loads(
                (outdir / "clustering_report.json").read_text(encoding="utf-8")
            )
            discovery_fasta = (outdir / "discovery_set.fasta").read_text(
                encoding="utf-8"
            )
            raw_cluster = (outdir / "clustering" / "cd-hit-est.clstr").read_text(
                encoding="utf-8"
            )

        self.assertEqual(runner.calls, [])
        self.assertEqual(summary.sequence_ids, [])
        self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], [])
        self.assertEqual(qc_report["discovery_set"]["sequence_ids"], [])
        self.assertEqual(discovery_fasta, "")
        self.assertEqual(raw_cluster, "")
        self.assertEqual(clustering_report["status"], "COMPLETE")
        self.assertEqual(
            clustering_report["artifacts"]["raw_cluster"],
            "clustering/cd-hit-est.clstr",
        )

    def test_enabled_alignment_uses_discovery_records_and_publishes_traceability(self):
        class FakeCdHitRunner:
            def run(self, input_path, output_path, config):
                del config
                with Path(input_path).open(encoding="utf-8") as handle:
                    records = {record.id: record for record in SeqIO.parse(handle, "fasta")}
                SeqIO.write(
                    [records["geison-00000000"], records["geison-00000002"]],
                    output_path,
                    "fasta",
                )
                Path(str(output_path) + ".clstr").write_text(
                    ">Cluster 0\n"
                    "0 12nt, >geison-00000000... *\n"
                    "1 12nt, >geison-00000001... at +/99.00%\n"
                    ">Cluster 1\n"
                    "0 12nt, >geison-00000002... *\n",
                    encoding="utf-8",
                )

        class FakeMafftRunner:
            def __init__(self):
                self.input_records = []

            def run(self, input_path, output_path, config):
                del config
                with Path(input_path).open(encoding="utf-8") as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
                self.input_records = [(record.id, str(record.seq)) for record in records]
                SeqIO.write(records, output_path, "fasta")

        mafft_runner = FakeMafftRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "target.fasta"
            outdir = tmp_path / "run"
            fasta_path.write_text(
                ">rejected\nACGTXCGTACGT\n"
                ">s1\nACGTACGTACGT\n"
                ">s2\nACGTACGAACGT\n"
                ">s3\nACGTACCCACGT\n",
                encoding="utf-8",
            )

            summary = run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_fasta=fasta_path,
                    clustering=ClusteringConfig(enabled=True),
                    alignment=AlignmentConfig(enabled=True, reference_id="s3"),
                ),
                outdir,
                cdhit_runner=FakeCdHitRunner(),
                mafft_runner=mafft_runner,
            )
            qc_report = json.loads((outdir / "qc_report.json").read_text(encoding="utf-8"))
            with (outdir / "alignment" / "discovery_alignment.fasta").open(
                encoding="utf-8"
            ) as handle:
                aligned_ids = [record.id for record in SeqIO.parse(handle, "fasta")]
            coordinate_columns = (
                outdir / "alignment" / "coordinate_map.tsv"
            ).read_text(encoding="utf-8").splitlines()[0].split("\t")

        self.assertEqual(
            mafft_runner.input_records,
            [
                ("geison-00000001", "ACGTACCCACGT"),
                ("geison-00000000", "ACGTACGTACGT"),
            ],
        )
        self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["s1", "s2", "s3"])
        self.assertEqual(qc_report["discovery_set"]["sequence_ids"], ["s1", "s3"])
        self.assertEqual(qc_report["alignment"], {
            "status": "COMPLETE",
            "reference_id": "s3",
            "reference_mode": "explicit",
        })
        self.assertEqual(summary.sequence_ids, ["s1", "s2", "s3"])
        self.assertEqual(aligned_ids, ["s1", "s3"])
        self.assertEqual(
            coordinate_columns,
            ["alignment_position", "reference_position", "reference_base"],
        )

    def test_enabled_conservation_uses_only_aligned_discovery_records(self):
        class FakeCdHitRunner:
            def run(self, input_path, output_path, config):
                del config
                with Path(input_path).open(encoding="utf-8") as handle:
                    records = {record.id: record for record in SeqIO.parse(handle, "fasta")}
                SeqIO.write(
                    [records["geison-00000000"], records["geison-00000002"]],
                    output_path,
                    "fasta",
                )
                Path(str(output_path) + ".clstr").write_text(
                    ">Cluster 0\n"
                    "0 4nt, >geison-00000000... *\n"
                    "1 4nt, >geison-00000001... at +/99.00%\n"
                    ">Cluster 1\n"
                    "0 4nt, >geison-00000002... *\n",
                    encoding="utf-8",
                )

        class FakeMafftRunner:
            def run(self, input_path, output_path, config):
                del config
                with Path(input_path).open(encoding="utf-8") as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
                SeqIO.write(records, output_path, "fasta")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "target.fasta"
            outdir = tmp_path / "run"
            fasta_path.write_text(
                ">rejected\nACGX\n"
                ">s1\nACGT\n"
                ">s2\nACAT\n"
                ">s3\nACCT\n",
                encoding="utf-8",
            )

            with patch("qpcr_pipeline.clustering.derive_word_length", return_value=4):
                summary = run_pipeline(
                    PipelineConfig(
                        target_name="synthetic-target",
                        input_fasta=fasta_path,
                        clustering=ClusteringConfig(enabled=True, identity=0.80),
                        alignment=AlignmentConfig(enabled=True, reference_id="s3"),
                        conservation=ConservationConfig(
                            enabled=True, window_size=3, step_size=2
                        ),
                    ),
                    outdir,
                    cdhit_runner=FakeCdHitRunner(),
                    mafft_runner=FakeMafftRunner(),
                )

            qc_report = json.loads(
                (outdir / "qc_report.json").read_text(encoding="utf-8")
            )
            conservation_report = json.loads(
                (outdir / "conservation" / "conservation_report.json").read_text(
                    encoding="utf-8"
                )
            )
            position_lines = (
                outdir / "conservation" / "position_metrics.tsv"
            ).read_text(encoding="utf-8").splitlines()
            position_header = position_lines[0].split("\t")
            position_rows = [
                dict(zip(position_header, line.split("\t"), strict=True))
                for line in position_lines[1:]
            ]
            conservation_files = {
                path.name for path in (outdir / "conservation").iterdir()
            }
            major_consensus = (
                outdir / "conservation" / "consensus_major.fasta"
            ).read_text(encoding="utf-8")
            iupac_consensus = (
                outdir / "conservation" / "consensus_iupac.fasta"
            ).read_text(encoding="utf-8")
            html = (outdir / "report.html").read_text(encoding="utf-8")

        self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["s1", "s2", "s3"])
        self.assertEqual(qc_report["discovery_set"]["sequence_ids"], ["s1", "s3"])
        self.assertEqual(
            qc_report["conservation"],
            {
                "status": "COMPLETE",
                "reference_id": "s3",
                "position_count": 4,
                "window_count": 2,
            },
        )
        self.assertEqual(summary.sequence_ids, ["s1", "s2", "s3"])
        self.assertEqual(conservation_report["discovery_set_ids"], ["s1", "s3"])
        self.assertEqual(
            conservation_files,
            {
                "conservation_report.json",
                "position_metrics.tsv",
                "window_metrics.tsv",
                "consensus_major.fasta",
                "consensus_iupac.fasta",
            },
        )
        self.assertEqual(len(position_rows), 4)
        self.assertTrue(all(row["depth"] == "2" for row in position_rows))
        self.assertEqual(position_rows[2]["frequency_a"], "0.0")
        self.assertEqual(position_rows[2]["frequency_c"], "0.5")
        self.assertEqual(position_rows[2]["frequency_g"], "0.5")
        self.assertEqual(major_consensus, ">geison-major-consensus\nACCT\n")
        self.assertEqual(iupac_consensus, ">geison-iupac-consensus\nACST\n")
        self.assertIn("<canvas", html)

    def test_run_atomically_replaces_summary_and_qc_hardlinks_to_frozen_artifacts(self):
        accession = "NC_FROZEN.1"
        record = self._ncbi_record(accession, "ACGTACGT")

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
            outdir.mkdir()
            summary_path = outdir / "run_summary.json"
            qc_path = outdir / "qc_report.json"
            source_records = frozen_dir / "records.gb"
            source_batch = frozen_dir / "batches" / "batch-00000.gb"
            os.link(source_records, summary_path)
            os.link(source_batch, qc_path)

            run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_ncbi=NcbiInputConfig(frozen_dataset=frozen_dir),
                ),
                outdir,
            )

            self.assertEqual(self._directory_bytes(frozen_dir), before)
            self.assertFalse(os.path.samefile(source_records, summary_path))
            self.assertFalse(os.path.samefile(source_batch, qc_path))
            self.assertEqual(json.loads(summary_path.read_text())["status"], "PARTIAL")
            self.assertEqual(
                json.loads(qc_path.read_text())["evaluation_set"]["sequence_ids"],
                [accession],
            )

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
            self.assertIn("PARTIAL", result.stdout)

            summary_path = outdir / "run_summary.json"
            self.assertTrue(summary_path.exists(), "run_summary.json was not created")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))

        self.assertEqual(summary["status"], "PARTIAL")
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
