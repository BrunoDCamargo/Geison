import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from Bio import SeqIO

from qpcr_pipeline import clustering
from qpcr_pipeline.clustering import CdHitError, cluster_sequences, derive_word_length
from qpcr_pipeline.config import ClusteringConfig
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet, EvaluationSet


VALID_CLUSTER_TEXT = (
    ">Cluster 9\n"
    "0 8nt, >geison-00000002... *\n"
    "1 8nt, >geison-00000000... at +/99.00%\n"
    ">Cluster 2\n"
    "0 8nt, >geison-00000001... *\n"
)


class FakeCdHitRunner:
    def __init__(self, representative_internal_ids, cluster_text):
        self.representative_internal_ids = representative_internal_ids
        self.cluster_text = cluster_text
        self.calls = []

    def run(self, input_path, output_path, config):
        self.calls.append((input_path, output_path, config))
        with Path(input_path).open(encoding="utf-8") as handle:
            records = {record.id: record for record in SeqIO.parse(handle, "fasta")}
        SeqIO.write(
            [records[record_id] for record_id in self.representative_internal_ids],
            output_path,
            "fasta",
        )
        Path(str(output_path) + ".clstr").write_text(
            self.cluster_text, encoding="utf-8"
        )


class FailingCdHitRunner:
    def __init__(self):
        self.calls = []

    def run(self, input_path, output_path, config):
        self.calls.append((input_path, output_path, config))
        raise CdHitError("fake CD-HIT failure")


class SubprocessCdHitRunnerTests(unittest.TestCase):
    def setUp(self):
        self.config = ClusteringConfig(
            enabled=True, identity=0.95, threads=4, memory_mb=2048
        )

    def test_invokes_cdhit_with_structured_arguments(self):
        # Omitting or changing a CD-HIT option would make clustering unsafe or incorrect.
        with tempfile.TemporaryDirectory(prefix="input path ") as tmpdir:
            root = Path(tmpdir)
            input_path = root / "input sequences.fasta"
            output_path = root / "representatives output.fasta"
            input_path.write_text(">sequence\nACGT\n", encoding="utf-8")
            output_path.write_text(">sequence\nACGT\n", encoding="utf-8")
            Path(str(output_path) + ".clstr").write_text(">Cluster 0\n", encoding="utf-8")

            with mock.patch("qpcr_pipeline.clustering.shutil.which", return_value="C:/tools/cd-hit-est.exe") as which, mock.patch(
                "qpcr_pipeline.clustering.subprocess.run",
                return_value=__import__("subprocess").CompletedProcess([], 0, "", ""),
            ) as run:
                clustering.SubprocessCdHitRunner().run(input_path, output_path, self.config)

        which.assert_called_once_with("cd-hit-est")
        run.assert_called_once_with(
            [
                "C:/tools/cd-hit-est.exe",
                "-i", str(input_path),
                "-o", str(output_path),
                "-c", "0.95",
                "-n", "10",
                "-d", "0",
                "-g", "1",
                "-r", "0",
                "-T", "4",
                "-M", "2048",
            ],
            capture_output=True,
            text=True,
            check=False,
        )

    def test_missing_executable_fails_before_starting_a_process(self):
        # Running with an unresolved executable would produce an opaque OS failure.
        with mock.patch("qpcr_pipeline.clustering.shutil.which", return_value=None), mock.patch(
            "qpcr_pipeline.clustering.subprocess.run"
        ) as run:
            with self.assertRaisesRegex(CdHitError, "not found on PATH"):
                clustering.SubprocessCdHitRunner().run(
                    Path("input.fasta"), Path("output.fasta"), self.config
                )

        run.assert_not_called()

    def test_nonzero_exit_reports_the_code_and_bounded_normalized_stderr(self):
        # Returning raw, unlimited stderr can hide the failure reason or exhaust reports.
        stderr = "  diagnostic   details\n" + ("x" * 2_100)
        with mock.patch("qpcr_pipeline.clustering.shutil.which", return_value="cd-hit-est"), mock.patch(
            "qpcr_pipeline.clustering.subprocess.run",
            return_value=__import__("subprocess").CompletedProcess([], 23, "", stderr),
        ):
            with self.assertRaises(CdHitError) as raised:
                clustering.SubprocessCdHitRunner().run(
                    Path("input.fasta"), Path("output.fasta"), self.config
                )

        message = str(raised.exception)
        self.assertIn("exit code 23", message)
        self.assertIn("diagnostic details", message)
        self.assertNotIn("\n", message)
        self.assertNotIn("x" * 2_001, message)

    def test_success_without_both_expected_output_files_fails(self):
        # Accepting a partial CD-HIT result would publish a corrupt clustering result.
        for missing_path in ("representatives.fasta", "representatives.fasta.clstr"):
            with self.subTest(missing_path=missing_path), tempfile.TemporaryDirectory() as tmpdir:
                output_path = Path(tmpdir) / "representatives.fasta"
                available_path = (
                    Path(str(output_path) + ".clstr")
                    if missing_path == output_path.name
                    else output_path
                )
                available_path.write_text(">sequence\nACGT\n", encoding="utf-8")
                with mock.patch("qpcr_pipeline.clustering.shutil.which", return_value="cd-hit-est"), mock.patch(
                    "qpcr_pipeline.clustering.subprocess.run",
                    return_value=__import__("subprocess").CompletedProcess([], 0, "", ""),
                ):
                    with self.assertRaisesRegex(CdHitError, "representative FASTA and .clstr"):
                        clustering.SubprocessCdHitRunner().run(
                            Path(tmpdir) / "input.fasta", output_path, self.config
                        )


class DefaultRunnerTests(unittest.TestCase):
    def test_default_runner_is_constructed_only_for_enabled_nonempty_clustering(self):
        # Constructing the external-tool runner for bypass paths would make those paths depend on CD-HIT.
        records = (LocalSequenceRecord("seq-1", "ACGTACGT"),)
        evaluation_set = EvaluationSet(("seq-1",))
        fake_runner = FakeCdHitRunner(
            ("geison-00000000",), ">Cluster 0\n0 8nt, >geison-00000000... *\n"
        )
        with mock.patch(
            "qpcr_pipeline.clustering.SubprocessCdHitRunner", return_value=fake_runner
        ) as runner_factory, tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            cluster_sequences(records, evaluation_set, ClusteringConfig(), output_dir)
            cluster_sequences((), EvaluationSet(()), ClusteringConfig(enabled=True), output_dir)
            cluster_sequences(
                records, evaluation_set, ClusteringConfig(enabled=True), output_dir
            )

        runner_factory.assert_called_once_with()
        self.assertEqual(len(fake_runner.calls), 1)


class ClusteringTests(unittest.TestCase):
    def setUp(self):
        self.records = (
            LocalSequenceRecord("seq-0", "ACGTACGT"),
            LocalSequenceRecord("seq-1", "GGGGCCCC"),
            LocalSequenceRecord("seq-2", "TTAATTAA"),
        )
        self.evaluation_set = EvaluationSet(("seq-0", "seq-1", "seq-2"))

    def test_derives_compatible_cdhit_word_lengths_at_every_boundary(self):
        # A wrong identity range must select an incompatible CD-HIT word length.
        for identity, expected in (
            (1.0, 10),
            (0.95, 10),
            (0.949, 8),
            (0.90, 8),
            (0.88, 7),
            (0.85, 6),
            (0.80, 5),
            (0.75, 4),
        ):
            with self.subTest(identity=identity):
                self.assertEqual(derive_word_length(identity), expected)

    def test_enabled_clustering_restores_ids_and_stabilizes_cluster_order(self):
        runner = FakeCdHitRunner(
            ("geison-00000002", "geison-00000001"), VALID_CLUSTER_TEXT
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            result = cluster_sequences(
                self.records,
                self.evaluation_set,
                ClusteringConfig(enabled=True),
                Path(tmpdir),
                runner=runner,
            )

            report = json.loads(result.report_path.read_text(encoding="utf-8"))
            with result.discovery_fasta_path.open(encoding="utf-8") as handle:
                fasta_records = list(SeqIO.parse(handle, "fasta"))
            raw_text = result.raw_cluster_path.read_text(encoding="utf-8")

        self.assertEqual(len(runner.calls), 1)
        self.assertEqual(result.discovery_set, DiscoverySet(("seq-1", "seq-2")))
        self.assertEqual(
            tuple(
                (
                    cluster.cluster_id,
                    cluster.representative_id,
                    tuple(
                        (member.sequence_id, member.representative, member.identity, member.strand)
                        for member in cluster.members
                    ),
                )
                for cluster in result.clusters
            ),
            (
                (
                    "cluster-00000",
                    "seq-2",
                    (("seq-2", True, None, None), ("seq-0", False, 99.0, "+")),
                ),
                (
                    "cluster-00001",
                    "seq-1",
                    (("seq-1", True, None, None),),
                ),
            ),
        )
        self.assertEqual(
            [(record.id, str(record.seq)) for record in fasta_records],
            [("seq-1", "GGGGCCCC"), ("seq-2", "TTAATTAA")],
        )
        self.assertEqual(raw_text, VALID_CLUSTER_TEXT)
        self.assertEqual(report["clusters"][0]["cluster_id"], "cluster-00000")
        self.assertEqual(report["clusters"][0]["members"][1], {
            "sequence_id": "seq-0", "representative": False, "identity": 99.0, "strand": "+"
        })

    def test_disabled_clustering_never_requires_or_calls_a_runner(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            result = cluster_sequences(
                self.records, self.evaluation_set, ClusteringConfig(), output_dir
            )
            report = json.loads(result.report_path.read_text(encoding="utf-8"))

        self.assertEqual(result.discovery_set, DiscoverySet(self.evaluation_set.sequence_ids))
        self.assertIsNone(result.raw_cluster_path)
        self.assertEqual(report["clusters"], [])
        self.assertFalse(report["clustering_enabled"])

    def test_empty_enabled_clustering_never_calls_runner(self):
        runner = FailingCdHitRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            result = cluster_sequences(
                (),
                EvaluationSet(()),
                ClusteringConfig(enabled=True),
                Path(tmpdir),
                runner=runner,
            )

        self.assertEqual(result.discovery_set, DiscoverySet(()))
        self.assertEqual(runner.calls, [])
        self.assertIsNone(result.raw_cluster_path)

    def test_publishes_complete_traceable_artifacts_without_internal_ids(self):
        runner = FakeCdHitRunner(
            ("geison-00000002", "geison-00000001"), VALID_CLUSTER_TEXT
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            result = cluster_sequences(
                self.records,
                self.evaluation_set,
                ClusteringConfig(enabled=True, identity=0.95, threads=4, memory_mb=2048),
                output_dir,
                runner=runner,
            )
            report_text = result.report_path.read_text(encoding="utf-8")
            report = json.loads(report_text)

            self.assertEqual(result.discovery_fasta_path, output_dir / "discovery_set.fasta")
            self.assertEqual(result.report_path, output_dir / "clustering_report.json")
            self.assertEqual(result.raw_cluster_path, output_dir / "clustering" / "cd-hit-est.clstr")
            self.assertEqual(report["schema_version"], 1)
            self.assertEqual(report["status"], "COMPLETE")
            self.assertTrue(report["clustering_enabled"])
            self.assertEqual(report["tool"], "cd-hit-est")
            self.assertEqual(report["parameters"], {
                "identity": 0.95, "word_length": 10, "threads": 4, "memory_mb": 2048
            })
            self.assertEqual(report["evaluation_set"], {"sequence_ids": ["seq-0", "seq-1", "seq-2"]})
            self.assertEqual(report["discovery_set"], {"sequence_ids": ["seq-1", "seq-2"]})
            self.assertEqual(report["counts"], {"evaluation": 3, "discovery": 2})
            self.assertEqual(report["artifacts"], {
                "discovery_fasta": "discovery_set.fasta",
                "raw_cluster": "clustering/cd-hit-est.clstr",
            })
            self.assertNotIn("geison-", report_text)
            self.assertNotIn("geison-", result.discovery_fasta_path.read_text(encoding="utf-8"))

    def test_replaces_a_preexisting_hardlinked_artifact_without_mutating_its_source(self):
        runner = FakeCdHitRunner(
            ("geison-00000002", "geison-00000001"), VALID_CLUSTER_TEXT
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "input-artifact.fasta"
            target = root / "discovery_set.fasta"
            source.write_text("original input artifact", encoding="utf-8")
            os.link(source, target)

            cluster_sequences(
                self.records,
                self.evaluation_set,
                ClusteringConfig(enabled=True),
                root,
                runner=runner,
            )

            self.assertEqual(source.read_text(encoding="utf-8"), "original input artifact")
            self.assertIn(">seq-1", target.read_text(encoding="utf-8"))

    def test_failing_runner_does_not_publish_a_complete_report(self):
        runner = FailingCdHitRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            with self.assertRaisesRegex(CdHitError, "fake CD-HIT failure"):
                cluster_sequences(
                    self.records,
                    self.evaluation_set,
                    ClusteringConfig(enabled=True),
                    output_dir,
                    runner=runner,
                )

            self.assertFalse((output_dir / "clustering_report.json").exists())

    def test_invalid_cluster_compositions_do_not_publish_complete_reports(self):
        cases = (
            ("duplicate member", ">Cluster 0\n0 8nt, >geison-00000000... *\n1 8nt, >geison-00000000... at +/99.00%\n", ("geison-00000000",)),
            ("missing member", ">Cluster 0\n0 8nt, >geison-00000000... *\n", ("geison-00000000",)),
            ("unknown member", ">Cluster 0\n0 8nt, >geison-00000000... *\n1 8nt, >unknown... at +/99.00%\n", ("geison-00000000",)),
            ("zero representatives", ">Cluster 0\n0 8nt, >geison-00000000... at +/99.00%\n1 8nt, >geison-00000001... at +/99.00%\n2 8nt, >geison-00000002... at +/99.00%\n", ("geison-00000000",)),
            ("multiple representatives", ">Cluster 0\n0 8nt, >geison-00000000... *\n1 8nt, >geison-00000001... *\n2 8nt, >geison-00000002... at +/99.00%\n", ("geison-00000000",)),
            ("representative output disagrees", VALID_CLUSTER_TEXT, ("geison-00000000", "geison-00000001")),
            ("malformed header", "Cluster 0\n0 8nt, >geison-00000000... *\n", ("geison-00000000",)),
            ("malformed member", ">Cluster 0\n0 8nt >geison-00000000... *\n", ("geison-00000000",)),
        )
        for name, cluster_text, representative_ids in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                output_dir = Path(tmpdir)
                runner = FakeCdHitRunner(representative_ids, cluster_text)
                with self.assertRaises(CdHitError):
                    cluster_sequences(
                        self.records,
                        self.evaluation_set,
                        ClusteringConfig(enabled=True),
                        output_dir,
                        runner=runner,
                    )
                self.assertFalse((output_dir / "clustering_report.json").exists())

    def test_duplicate_approved_original_ids_fail_before_runner_or_artifacts(self):
        records = (
            LocalSequenceRecord("duplicate", "ACGTACGT"),
            LocalSequenceRecord("duplicate", "GGGGCCCC"),
        )
        evaluation_set = EvaluationSet(("duplicate", "duplicate"))
        runner = FailingCdHitRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            with self.assertRaisesRegex(CdHitError, "Duplicate"):
                cluster_sequences(
                    records,
                    evaluation_set,
                    ClusteringConfig(enabled=True),
                    output_dir,
                    runner=runner,
                )
            self.assertEqual(runner.calls, [])
            self.assertFalse((output_dir / "clustering_report.json").exists())


if __name__ == "__main__":
    unittest.main()
