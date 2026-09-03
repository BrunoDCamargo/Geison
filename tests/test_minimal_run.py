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
from panel_fixtures import approved_panel_config


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

            # This fixture isolates pipeline routing from CD-HIT's five-base
            # executable limit while keeping the real clustering parser active.
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

    def test_enabled_primer_design_uses_major_consensus_and_publishes_assays(self):
        class FakeMafftRunner:
            def run(self, input_path, output_path, config):
                del config
                shutil.copyfile(input_path, output_path)

        class FakePrimer3Runner:
            def __init__(self):
                self.inputs = []

            def run(self, input_text):
                self.inputs.append(input_text)
                return (
                    "SEQUENCE_ID=region-001\n"
                    "PRIMER_LEFT_NUM_RETURNED=2\n"
                    "PRIMER_INTERNAL_NUM_RETURNED=2\n"
                    "PRIMER_RIGHT_NUM_RETURNED=2\n"
                    "PRIMER_PAIR_NUM_RETURNED=2\n"
                    "PRIMER_LEFT_0=10,20\n"
                    "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
                    "PRIMER_LEFT_0_TM=60.0\n"
                    "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
                    "PRIMER_INTERNAL_0=35,25\n"
                    "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
                    "PRIMER_INTERNAL_0_TM=70.0\n"
                    "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
                    "PRIMER_RIGHT_0=89,20\n"
                    "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
                    "PRIMER_RIGHT_0_TM=60.0\n"
                    "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
                    "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
                    "PRIMER_LEFT_1=100,20\n"
                    "PRIMER_LEFT_1_SEQUENCE=AAAACCCCGGGGTTTTAAAA\n"
                    "PRIMER_LEFT_1_TM=59.0\n"
                    "PRIMER_LEFT_1_GC_PERCENT=40.0\n"
                    "PRIMER_INTERNAL_1=130,20\n"
                    "PRIMER_INTERNAL_1_SEQUENCE=CCCCAAAATTTTGGGGCCCC\n"
                    "PRIMER_INTERNAL_1_TM=69.0\n"
                    "PRIMER_INTERNAL_1_GC_PERCENT=60.0\n"
                    "PRIMER_RIGHT_1=169,20\n"
                    "PRIMER_RIGHT_1_SEQUENCE=TTTTGGGGCCCCAAAATTTT\n"
                    "PRIMER_RIGHT_1_TM=61.0\n"
                    "PRIMER_RIGHT_1_GC_PERCENT=40.0\n"
                    "PRIMER_PAIR_1_PRODUCT_SIZE=70\n"
                    "=\n"
                )

        consensus = "ACGT" * 75
        primer3_runner = FakePrimer3Runner()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "target.fasta"
            outdir = tmp_path / "run"
            fasta_path.write_text(
                f">seq-1\n{consensus}\n>seq-2\n{consensus}\n",
                encoding="utf-8",
            )

            run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_fasta=fasta_path,
                    panel=approved_panel_config(tmp_path, "synthetic-target"),
                    alignment=AlignmentConfig(enabled=True, reference_id="seq-1"),
                    conservation=ConservationConfig(
                        enabled=True, window_size=100, step_size=50
                    ),
                    primer_design=PrimerDesignConfig(
                        enabled=True,
                        max_candidate_regions=1,
                        assays_per_region=2,
                    ),
                ),
                outdir,
                mafft_runner=FakeMafftRunner(),
                primer3_runner=primer3_runner,
            )

            qc_report = json.loads(
                (outdir / "qc_report.json").read_text(encoding="utf-8")
            )
            assay_ids = [
                line.split("\t", 1)[0]
                for line in (outdir / "primer_design" / "assays.tsv")
                .read_text(encoding="utf-8")
                .splitlines()[1:]
            ]

        templates = [
            line.removeprefix("SEQUENCE_TEMPLATE=")
            for input_text in primer3_runner.inputs
            for line in input_text.splitlines()
            if line.startswith("SEQUENCE_TEMPLATE=")
        ]
        self.assertEqual(templates, [consensus])
        self.assertEqual(
            assay_ids,
            ["region-001-assay-001", "region-001-assay-002"],
        )
        self.assertEqual(qc_report["primer_design"], {
            "status": "COMPLETE",
            "reference_id": "seq-1",
            "candidate_region_count": 1,
            "assay_count": 2,
        })
        self.assertEqual(qc_report["inclusivity"], {
            "status": "SKIPPED",
            "evaluation_sequence_count": 0,
            "assay_count": 0,
            "assay_evaluation_count": 0,
            "original_compatible_count": 0,
            "proposed_compatible_count": 0,
        })

    def test_enabled_inclusivity_evaluates_full_evaluation_set_beyond_discovery(self):
        class FakeCdHitRunner:
            def run(self, input_path, output_path, config):
                del config
                with Path(input_path).open(encoding="utf-8") as handle:
                    records = {record.id: record for record in SeqIO.parse(handle, "fasta")}
                SeqIO.write(
                    [records["geison-00000000"], records["geison-00000001"]],
                    output_path,
                    "fasta",
                )
                Path(str(output_path) + ".clstr").write_text(
                    ">Cluster 0\n"
                    "0 100nt, >geison-00000000... *\n"
                    "1 100nt, >geison-00000002... at +/99.00%\n"
                    ">Cluster 1\n"
                    "0 100nt, >geison-00000001... *\n",
                    encoding="utf-8",
                )

        class FakeMafftRunner:
            def run(self, input_path, output_path, config):
                del config
                shutil.copyfile(input_path, output_path)

        class FakePrimer3Runner:
            def run(self, input_text):
                del input_text
                return (
                    "SEQUENCE_ID=region-001\n"
                    "PRIMER_LEFT_NUM_RETURNED=1\n"
                    "PRIMER_INTERNAL_NUM_RETURNED=1\n"
                    "PRIMER_RIGHT_NUM_RETURNED=1\n"
                    "PRIMER_PAIR_NUM_RETURNED=1\n"
                    "PRIMER_LEFT_0=10,20\n"
                    "PRIMER_LEFT_0_SEQUENCE=AAAAAAAAAAAAAAAAAAAA\n"
                    "PRIMER_LEFT_0_TM=60.0\n"
                    "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
                    "PRIMER_INTERNAL_0=35,25\n"
                    "PRIMER_INTERNAL_0_SEQUENCE=CCCCCCCCCCCCCCCCCCCCCCCCC\n"
                    "PRIMER_INTERNAL_0_TM=70.0\n"
                    "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
                    "PRIMER_RIGHT_0=89,20\n"
                    "PRIMER_RIGHT_0_SEQUENCE=GGGGGGGGGGGGGGGGGGGG\n"
                    "PRIMER_RIGHT_0_TM=60.0\n"
                    "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
                    "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
                    "=\n"
                )

        sequence = "T" * 10 + "A" * 20 + "T" * 5 + "C" * 25 + "T" * 10 + "C" * 20 + "T" * 10
        second_sequence = "G" + sequence[1:]
        varied_sequence = "C" + sequence[1:14] + "G" + sequence[15:]
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "target.fasta"
            outdir = tmp_path / "run"
            fasta_path.write_text(
                f">s1\n{sequence}\n>s2\n{second_sequence}\n>s3\n{varied_sequence}\n",
                encoding="utf-8",
            )

            run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_fasta=fasta_path,
                    panel=approved_panel_config(tmp_path, "synthetic-target"),
                    clustering=ClusteringConfig(enabled=True),
                    alignment=AlignmentConfig(enabled=True, reference_id="s1"),
                    conservation=ConservationConfig(
                        enabled=True, window_size=100, step_size=100
                    ),
                    primer_design=PrimerDesignConfig(
                        enabled=True,
                        max_candidate_regions=1,
                        assays_per_region=1,
                        min_minimum_conservation=0.5,
                    ),
                    inclusivity=InclusivityConfig(enabled=True, search_flank=0),
                ),
                outdir,
                cdhit_runner=FakeCdHitRunner(),
                mafft_runner=FakeMafftRunner(),
                primer3_runner=FakePrimer3Runner(),
            )

            qc = json.loads((outdir / "qc_report.json").read_text(encoding="utf-8"))
            report = json.loads(
                (outdir / "inclusivity" / "inclusivity_report.json").read_text(
                    encoding="utf-8"
                )
            )

        self.assertTrue(report["variations"], report)
        variation = next(item for item in report["variations"] if item["role"] == "FORWARD")
        proposal = next(item for item in report["proposals"] if item["role"] == "FORWARD")
        self.assertEqual(report["evaluation_sequence_ids"], ["s1", "s2", "s3"])
        self.assertEqual(qc["inclusivity"]["evaluation_sequence_count"], 3)
        self.assertEqual(qc["inclusivity"]["assay_evaluation_count"], 3)
        self.assertIn("s3", variation["affected_sequence_ids"])
        self.assertEqual(proposal["status"], "ACCEPTED")
        self.assertNotEqual(proposal["original_sequence"], proposal["proposed_sequence"])

    def test_disabled_local_clustering_does_not_call_runner_and_keeps_evaluation_set(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            with patch("qpcr_pipeline.clustering.SubprocessCdHitRunner") as runner_factory, patch(
                "qpcr_pipeline.alignment.SubprocessMafftRunner"
            ) as mafft_runner_factory:
                summary = run_pipeline(
                    PipelineConfig(
                        target_name="synthetic-target", input_fasta=FIXTURE_FASTA
                    ),
                    outdir,
                )
            qc_report = json.loads(
                (outdir / "qc_report.json").read_text(encoding="utf-8")
            )
            alignment_report = json.loads(
                (outdir / "alignment" / "alignment_report.json").read_text(
                    encoding="utf-8"
                )
            )
            conservation_report = json.loads(
                (outdir / "conservation" / "conservation_report.json").read_text(
                    encoding="utf-8"
                )
            )

        runner_factory.assert_not_called()
        mafft_runner_factory.assert_not_called()
        self.assertEqual(summary.sequence_ids, ["seq-1", "seq-2", "seq-3"])
        self.assertEqual(
            qc_report["evaluation_set"]["sequence_ids"], ["seq-1", "seq-2", "seq-3"]
        )
        self.assertEqual(
            qc_report["discovery_set"]["sequence_ids"], ["seq-1", "seq-2", "seq-3"]
        )
        self.assertEqual(
            qc_report["alignment"],
            {"status": "SKIPPED", "reference_id": None, "reference_mode": None},
        )
        self.assertEqual(alignment_report["status"], "SKIPPED")
        self.assertEqual(
            qc_report["conservation"],
            {
                "status": "SKIPPED",
                "reference_id": None,
                "position_count": 0,
                "window_count": 0,
            },
        )
        self.assertEqual(conservation_report["status"], "SKIPPED")
        self.assertFalse((outdir / "conservation" / "position_metrics.tsv").exists())
        self.assertFalse((outdir / "conservation" / "window_metrics.tsv").exists())
        self.assertFalse((outdir / "conservation" / "consensus_major.fasta").exists())
        self.assertFalse((outdir / "conservation" / "consensus_iupac.fasta").exists())
        self.assertFalse((outdir / "report.html").exists())

    def test_run_rejects_invalid_clustering_config_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"

            with self.assertRaisesRegex(ValueError, "identity.*0.80 and 1.0"):
                run_pipeline(
                    PipelineConfig(
                        target_name="synthetic-target",
                        input_fasta=FIXTURE_FASTA,
                        clustering=ClusteringConfig(identity=0.799),
                    ),
                    outdir,
                )

            self.assertFalse(outdir.exists())

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
            with patch("qpcr_pipeline.clustering.SubprocessCdHitRunner") as runner_factory, patch(
                "qpcr_pipeline.alignment.SubprocessMafftRunner"
            ) as mafft_runner_factory:
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
            alignment_report = json.loads(
                (outdir / "alignment" / "alignment_report.json").read_text(
                    encoding="utf-8"
                )
            )
            conservation_report = json.loads(
                (outdir / "conservation" / "conservation_report.json").read_text(
                    encoding="utf-8"
                )
            )

        runner_factory.assert_not_called()
        mafft_runner_factory.assert_not_called()
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
        self.assertEqual(qc_report["discovery_set"]["sequence_ids"], [valid_accession])
        self.assertEqual(
            qc_report["alignment"],
            {"status": "SKIPPED", "reference_id": None, "reference_mode": None},
        )
        self.assertEqual(alignment_report["status"], "SKIPPED")
        self.assertEqual(
            qc_report["conservation"],
            {
                "status": "SKIPPED",
                "reference_id": None,
                "position_count": 0,
                "window_count": 0,
            },
        )
        self.assertEqual(conservation_report["status"], "SKIPPED")
        self.assertFalse((outdir / "conservation" / "position_metrics.tsv").exists())
        self.assertFalse((outdir / "conservation" / "window_metrics.tsv").exists())
        self.assertFalse((outdir / "conservation" / "consensus_major.fasta").exists())
        self.assertFalse((outdir / "conservation" / "consensus_iupac.fasta").exists())
        self.assertFalse((outdir / "report.html").exists())
        self.assertEqual(effective_manifest["status"], "COMPLETE")

    def test_run_rejects_invalid_primer_design_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            config = PipelineConfig(
                target_name="invalid-primer-design",
                input_fasta=FIXTURE_FASTA,
                primer_design=PrimerDesignConfig(candidate_region_length=0),
            )

            with self.assertRaisesRegex(ValueError, "candidate_region_length"):
                run_pipeline(config, outdir)

            self.assertFalse(outdir.exists())

    def test_run_rejects_enabled_primer_design_without_conservation_before_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            config = PipelineConfig(
                target_name="primer-design-without-conservation",
                input_fasta=FIXTURE_FASTA,
                primer_design=PrimerDesignConfig(enabled=True),
            )

            with self.assertRaisesRegex(ValueError, "requires enabled conservation"):
                run_pipeline(config, outdir)

            self.assertFalse(outdir.exists())

    def test_run_rejects_enabled_inclusivity_without_primer_design_before_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            config = PipelineConfig(
                target_name="inclusivity-without-primer-design",
                input_fasta=FIXTURE_FASTA,
                inclusivity=InclusivityConfig(enabled=True),
            )

            with self.assertRaisesRegex(ValueError, "requires enabled primer design"):
                run_pipeline(config, outdir)

            self.assertFalse(outdir.exists())

    def test_disabled_primer_design_publishes_only_skipped_report_without_runner(self):
        class FailingPrimer3Runner:
            def __init__(self):
                self.calls = []

            def run(self, input_text):
                self.calls.append(input_text)
                raise AssertionError("disabled primer design must not call Primer3")

        runner = FailingPrimer3Runner()
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_fasta=FIXTURE_FASTA,
                ),
                outdir,
                primer3_runner=runner,
            )

            primer_design_files = {
                path.name for path in (outdir / "primer_design").iterdir()
            }
            inclusivity_files = {
                path.name for path in (outdir / "inclusivity").iterdir()
            }
            primer_design_report = json.loads(
                (outdir / "primer_design" / "primer_design_report.json").read_text(
                    encoding="utf-8"
                )
            )
            qc_report = json.loads(
                (outdir / "qc_report.json").read_text(encoding="utf-8")
            )

        self.assertEqual(runner.calls, [])
        self.assertEqual(primer_design_files, {"primer_design_report.json"})
        self.assertEqual(inclusivity_files, {"inclusivity_report.json"})
        self.assertEqual(primer_design_report["status"], "SKIPPED")
        self.assertEqual(qc_report["primer_design"], {
            "status": "SKIPPED",
            "reference_id": None,
            "candidate_region_count": 0,
            "assay_count": 0,
        })
        self.assertEqual(qc_report["inclusivity"], {
            "status": "SKIPPED",
            "evaluation_sequence_count": 0,
            "assay_count": 0,
            "assay_evaluation_count": 0,
            "original_compatible_count": 0,
            "proposed_compatible_count": 0,
        })

    def test_run_rejects_invalid_conservation_before_output_or_analysis(self):
        cases = (
            (
                PipelineConfig(
                    target_name="invalid-conservation",
                    input_fasta=FIXTURE_FASTA,
                    conservation=ConservationConfig(window_size=0),
                ),
                "window_size",
            ),
            (
                PipelineConfig(
                    target_name="conservation-without-alignment",
                    input_fasta=FIXTURE_FASTA,
                    conservation=ConservationConfig(enabled=True),
                ),
                "requires enabled alignment",
            ),
        )

        for config, message in cases:
            with self.subTest(config=config), tempfile.TemporaryDirectory() as tmpdir:
                outdir = Path(tmpdir) / "run"
                with patch("qpcr_pipeline.pipeline.analyze_conservation") as analyze:
                    with self.assertRaisesRegex(ValueError, message):
                        run_pipeline(config, outdir)

                analyze.assert_not_called()
                self.assertFalse(outdir.exists())

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

    def _assert_frozen_manifest_destination_does_not_mutate_source(self, link):
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
            source_manifest = (frozen_dir / "dataset_manifest.json").read_bytes()
            outdir = tmp_path / "run"
            outdir.mkdir()
            destination = outdir / "ncbi_dataset_manifest.json"
            link(frozen_dir / "records.gb", destination)

            run_pipeline(
                PipelineConfig(
                    target_name="synthetic-target",
                    input_ncbi=NcbiInputConfig(frozen_dataset=frozen_dir),
                ),
                outdir,
            )

            after = self._directory_bytes(frozen_dir)
            effective_manifest = destination.read_bytes()

        self.assertEqual(after, before)
        self.assertEqual(effective_manifest, source_manifest)
        self.assertEqual(json.loads(effective_manifest)["status"], "COMPLETE")

    def test_run_replaces_hardlinked_frozen_manifest_destination(self):
        self._assert_frozen_manifest_destination_does_not_mutate_source(
            lambda source, destination: os.link(source, destination)
        )

    def test_run_replaces_symlinked_frozen_manifest_destination_when_supported(self):
        def symlink(source, destination):
            try:
                destination.symlink_to(source)
            except OSError as error:
                self.skipTest(f"symlink creation is unavailable: {error}")

        self._assert_frozen_manifest_destination_does_not_mutate_source(symlink)

    def test_run_rejects_invalid_direct_config_before_creating_output(self):
        cases = (
            PipelineConfig(
                target_name="ambiguous",
                input_fasta=FIXTURE_FASTA,
                input_genbank=Path("unused.gb"),
            ),
            PipelineConfig(
                target_name="ambiguous-ncbi",
                input_ncbi=NcbiInputConfig(
                    query="example[Organism]", accessions=("NC_1",)
                ),
            ),
            PipelineConfig(
                target_name="invalid-alignment",
                input_fasta=FIXTURE_FASTA,
                alignment=AlignmentConfig(threads=0),
            ),
        )

        for config in cases:
            with self.subTest(config=config), tempfile.TemporaryDirectory() as tmpdir:
                outdir = Path(tmpdir) / "run"

                with self.assertRaisesRegex(ValueError, "Exactly one|exactly one|Alignment threads"):
                    run_pipeline(config, outdir)

                self.assertFalse(outdir.exists())

    def test_enabled_alignment_rejects_reference_outside_discovery_without_calling_mafft(self):
        class FailingMafftRunner:
            def __init__(self):
                self.calls = []

            def run(self, input_path, output_path, config):
                self.calls.append((input_path, output_path, config))
                raise AssertionError("MAFFT must not run for a reference outside Discovery")

        runner = FailingMafftRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            with self.assertRaisesRegex(MafftError, "not in the Discovery Set"):
                run_pipeline(
                    PipelineConfig(
                        target_name="synthetic-target",
                        input_fasta=FIXTURE_FASTA,
                        alignment=AlignmentConfig(enabled=True, reference_id="missing"),
                    ),
                    outdir,
                    mafft_runner=runner,
                )

            report_path = outdir / "alignment" / "alignment_report.json"
            report_status = (
                json.loads(report_path.read_text(encoding="utf-8"))["status"]
                if report_path.exists()
                else None
            )

        self.assertEqual(runner.calls, [])
        self.assertNotEqual(report_status, "COMPLETE")

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

    def test_run_rejects_output_inside_frozen_dataset_without_mutation(self):
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

        for name in ("same directory", "nested directory"):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                frozen_dir = Path(tmpdir) / "frozen"
                acquire_ncbi_dataset(
                    NcbiInputConfig(accessions=(accession,)),
                    frozen_dir,
                    client=DatasetWriter(),
                    clock=lambda: "2026-08-21T00:00:00+00:00",
                )
                before = self._directory_bytes(frozen_dir)
                outdir = (
                    frozen_dir if name == "same directory" else frozen_dir / "nested-run"
                )

                with self.assertRaisesRegex(ValueError, "output directory.*frozen dataset"):
                    run_pipeline(
                        PipelineConfig(
                            target_name="synthetic-target",
                            input_ncbi=NcbiInputConfig(frozen_dataset=frozen_dir),
                        ),
                        outdir,
                    )

                if outdir != frozen_dir:
                    self.assertFalse(outdir.exists())
                self.assertEqual(self._directory_bytes(frozen_dir), before)

    def test_atomic_summary_write_cleans_temporary_file_after_replace_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            outdir.mkdir()

            with patch.object(Path, "replace", side_effect=OSError("replace failed")):
                with self.assertRaisesRegex(OSError, "replace failed"):
                    run_pipeline(
                        PipelineConfig(
                            target_name="synthetic-target", input_fasta=FIXTURE_FASTA
                        ),
                        outdir,
                    )

            self.assertFalse((outdir / "run_summary.json").exists())
            self.assertFalse((outdir / "qc_report.json").exists())
            self.assertEqual(list(outdir.iterdir()), [])

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
