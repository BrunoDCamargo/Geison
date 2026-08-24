import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.alignment import (
    AlignedSequence,
    AlignmentCoordinate,
    MafftError,
    align_discovery,
)
from qpcr_pipeline.config import AlignmentConfig
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet


class FakeMafftRunner:
    def __init__(self, output_records):
        self.output_records = output_records
        self.calls = []
        self.input_records = []

    def run(self, input_path, output_path, config):
        self.calls.append((input_path, output_path, config))
        with Path(input_path).open(encoding="utf-8") as handle:
            self.input_records = [
                (record.id, str(record.seq))
                for record in SeqIO.parse(handle, "fasta")
            ]
        SeqIO.write(
            [
                SeqRecord(Seq(sequence), id=sequence_id, description="")
                for sequence_id, sequence in self.output_records
            ],
            output_path,
            "fasta",
        )


class FailingRunner:
    def run(self, input_path, output_path, config):
        raise AssertionError("runner should not be called")


class RaisingRunner:
    def run(self, input_path, output_path, config):
        raise MafftError("MAFFT failed")


class EmptyOutputRunner:
    def run(self, input_path, output_path, config):
        Path(output_path).write_text("", encoding="utf-8")


def record(sequence_id, sequence):
    return LocalSequenceRecord(sequence_id=sequence_id, sequence=sequence)


class AlignmentTests(unittest.TestCase):
    def setUp(self):
        self.records = (
            record("ref", "ACGT"),
            record("reverse", "ACGT"),
            record("other", "ACGA"),
        )
        self.discovery = DiscoverySet(("ref", "reverse", "other"))

    def _run(self, records=None, discovery=None, config=None, runner=None, directory=None):
        return align_discovery(
            self.records if records is None else records,
            self.discovery if discovery is None else discovery,
            AlignmentConfig(enabled=True, reference_id="ref")
            if config is None
            else config,
            Path(directory),
            runner=runner,
        )

    def _assert_no_complete_report(self, directory):
        report = Path(directory) / "alignment" / "alignment_report.json"
        if report.exists():
            self.assertNotEqual(json.loads(report.read_text(encoding="utf-8"))["status"], "COMPLETE")

    def test_runner_alignment_restores_ids_order_and_orientations(self):
        runner = FakeMafftRunner(
            (
                ("geison-00000002", "ACGA"),
                ("_R_geison-00000001", "ACGT"),
                ("geison-00000000", "ACGT"),
            )
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self._run(runner=runner, directory=tmpdir)
            fasta = (Path(tmpdir) / "alignment" / "discovery_alignment.fasta").read_text(
                encoding="utf-8"
            )
            report = (Path(tmpdir) / "alignment" / "alignment_report.json").read_text(
                encoding="utf-8"
            )

        self.assertEqual(
            runner.input_records,
            [
                ("geison-00000000", "ACGT"),
                ("geison-00000001", "ACGT"),
                ("geison-00000002", "ACGA"),
            ],
        )
        self.assertEqual(result.discovery_set, self.discovery)
        self.assertEqual([item.sequence_id for item in result.sequences], ["ref", "reverse", "other"])
        self.assertEqual(
            [item.orientation for item in result.sequences],
            ["forward", "reverse_complemented", "forward"],
        )
        self.assertNotIn("geison-", fasta + report)
        self.assertNotIn("_R_", fasta + report)

    def test_automatic_reference_uses_ambiguity_length_then_discovery_order(self):
        cases = (
            (
                (record("ambiguous", "ACGN"), record("clean", "ACGT")),
                "clean",
            ),
            (
                (record("short", "ACG"), record("long", "ACGT")), "long"),
            (
                (record("first", "ACGT"), record("second", "TGCA")), "first"),
        )
        for records, expected in cases:
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as tmpdir:
                discovery = DiscoverySet(tuple(item.sequence_id for item in records))
                runner = FakeMafftRunner(
                    tuple(
                        (f"geison-{position:08d}", item.sequence + ("-" if item.sequence_id == "short" else ""))
                        for position, item in enumerate(records)
                    )
                )
                result = self._run(
                    records=records,
                    discovery=discovery,
                    config=AlignmentConfig(enabled=True),
                    runner=runner,
                    directory=tmpdir,
                )
                report = json.loads(result.report_path.read_text(encoding="utf-8"))
                self.assertEqual(result.reference_id, expected)
                self.assertEqual(report["reference"]["mode"], "automatic")
                self.assertEqual(
                    report["reference"]["automatic_selection_rule"],
                    "lowest_ambiguity_fraction_then_longest_then_discovery_order",
                )

    def test_coordinates_and_tsv_represent_reference_gaps(self):
        records = (record("ref", "ACGT"), record("other", "ATCGT"))
        runner = FakeMafftRunner(
            (("geison-00000000", "A-CGT"), ("geison-00000001", "ATCGT"))
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self._run(
                records=records,
                discovery=DiscoverySet(("ref", "other")),
                runner=runner,
                directory=tmpdir,
            )
            rows = result.coordinate_map_path.read_text(encoding="utf-8").splitlines()

        self.assertEqual(
            result.coordinates,
            (
                AlignmentCoordinate(1, 1, "A"),
                AlignmentCoordinate(2, None, None),
                AlignmentCoordinate(3, 2, "C"),
                AlignmentCoordinate(4, 3, "G"),
                AlignmentCoordinate(5, 4, "T"),
            ),
        )
        self.assertEqual(rows[0], "alignment_position\treference_position\treference_base")
        self.assertEqual(rows[2], "2\t\t")

    def test_validation_errors_do_not_publish_complete_report(self):
        valid_output = (("geison-00000000", "ACGT"), ("geison-00000001", "ACGA"))
        cases = (
            ("duplicate discovery", self.records, DiscoverySet(("ref", "ref")), None, None),
            ("duplicate records", (record("ref", "ACGT"), record("ref", "ACGA")), DiscoverySet(("ref", "ref")), None, None),
            ("membership mismatch", self.records[:2], self.discovery, None, None),
            ("reference absent", self.records, self.discovery, None, AlignmentConfig(enabled=True, reference_id="missing")),
            ("empty input", (record("ref", ""),), DiscoverySet(("ref",)), None, AlignmentConfig(enabled=True)),
            ("missing output", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner(valid_output[:1]), None),
            ("unknown output", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner((("unknown", "ACGT"), ("geison-00000001", "ACGT"))), None),
            ("duplicate output", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner((("geison-00000000", "ACGT"), ("geison-00000000", "ACGT"))), None),
            ("double reverse", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner((("_R__R_geison-00000000", "ACGT"), ("geison-00000001", "ACGT"))), None),
            ("reversed reference", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner((("_R_geison-00000000", "ACGT"), ("geison-00000001", "ACGT"))), None),
            ("unequal lengths", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner((("geison-00000000", "ACGT"), ("geison-00000001", "ACG-T"))), None),
            ("all gap", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner((("geison-00000000", "ACGT"), ("geison-00000001", "----"))), None),
            ("invalid output symbol", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner((("geison-00000000", "ACGT"), ("geison-00000001", "ACGZ"))), None),
            ("forward mutation", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner((("geison-00000000", "ACGT"), ("geison-00000001", "ACGG"))), None),
            ("reverse mutation", self.records[:2], DiscoverySet(("ref", "reverse")), FakeMafftRunner((("geison-00000000", "ACGT"), ("_R_geison-00000001", "ACGG"))), None),
            ("runner error", self.records[:2], DiscoverySet(("ref", "reverse")), RaisingRunner(), None),
            ("empty output", self.records[:2], DiscoverySet(("ref", "reverse")), EmptyOutputRunner(), None),
        )
        for name, records, discovery, runner, config in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                with self.assertRaises((MafftError, ValueError)):
                    self._run(records, discovery, config, runner, tmpdir)
                self._assert_no_complete_report(tmpdir)

    def test_disabled_removes_only_stale_data_and_publishes_skipped_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment_dir = Path(tmpdir) / "alignment"
            alignment_dir.mkdir()
            stale_fasta = alignment_dir / "discovery_alignment.fasta"
            stale_tsv = alignment_dir / "coordinate_map.tsv"
            sibling = alignment_dir / "keep.txt"
            stale_fasta.write_text("stale", encoding="utf-8")
            stale_tsv.write_text("stale", encoding="utf-8")
            sibling.write_text("keep", encoding="utf-8")
            result = self._run(
                config=AlignmentConfig(enabled=False), runner=FailingRunner(), directory=tmpdir
            )
            report = json.loads(result.report_path.read_text(encoding="utf-8"))

            self.assertEqual(result.status, "SKIPPED")
            self.assertEqual(set(alignment_dir.iterdir()), {sibling, result.report_path})
            self.assertIsNone(result.alignment_fasta_path)
            self.assertIsNone(result.coordinate_map_path)
            self.assertEqual(report["artifacts"], {"alignment_fasta": None, "coordinate_map": None})

    def test_enabled_empty_publishes_empty_artifacts_without_runner(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self._run(
                records=(), discovery=DiscoverySet(()), config=AlignmentConfig(enabled=True),
                runner=FailingRunner(), directory=tmpdir,
            )
            self.assertEqual(result.status, "COMPLETE")
            self.assertEqual(result.alignment_fasta_path.read_text(encoding="utf-8"), "")
            self.assertEqual(
                result.coordinate_map_path.read_text(encoding="utf-8"),
                "alignment_position\treference_position\treference_base\n",
            )
            self.assertIsNone(result.reference_id)
            self.assertEqual(result.sequences, ())
            self.assertEqual(result.coordinates, ())

    def test_enabled_singleton_is_identity_without_runner(self):
        singleton = (record("only", "ACGT"),)
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self._run(
                records=singleton, discovery=DiscoverySet(("only",)),
                config=AlignmentConfig(enabled=True), runner=FailingRunner(), directory=tmpdir,
            )
            self.assertEqual(result.reference_id, "only")
            self.assertEqual(result.reference_mode, "automatic")
            self.assertEqual(result.sequences, (AlignedSequence("only", "ACGT", "forward"),))
            self.assertEqual([item.reference_position for item in result.coordinates], [1, 2, 3, 4])

    def test_report_contains_traceable_schema_and_relative_artifacts(self):
        runner = FakeMafftRunner(
            (("geison-00000000", "ACGT"), ("geison-00000001", "ACGT"), ("geison-00000002", "ACGA"))
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self._run(runner=runner, directory=tmpdir)
            report = json.loads(result.report_path.read_text(encoding="utf-8"))

        self.assertEqual(report["schema_version"], 1)
        self.assertEqual(report["enabled"], True)
        self.assertEqual(report["tool"], {"name": "mafft", "parameters": {"threads": 1}})
        self.assertEqual(report["discovery_set_ids"], ["ref", "reverse", "other"])
        self.assertEqual(report["reference"], {"id": "ref", "mode": "explicit", "automatic_selection_rule": None})
        self.assertEqual([item["orientation"] for item in report["orientations"]], ["forward", "forward", "forward"])
        self.assertEqual(report["counts"], {"discovery": 3, "alignment_length": 4, "reverse_complemented": 0})
        self.assertEqual(report["artifacts"], {"alignment_fasta": "alignment/discovery_alignment.fasta", "coordinate_map": "alignment/coordinate_map.tsv"})

    def test_publication_replaces_hardlink_without_modifying_source(self):
        runner = FakeMafftRunner(
            (("geison-00000000", "ACGT"), ("geison-00000001", "ACGT"), ("geison-00000002", "ACGA"))
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            source = base / "source.fasta"
            source.write_text("do not change", encoding="utf-8")
            alignment_dir = base / "alignment"
            alignment_dir.mkdir()
            os.link(source, alignment_dir / "discovery_alignment.fasta")
            self._run(runner=runner, directory=tmpdir)
            self.assertEqual(source.read_text(encoding="utf-8"), "do not change")

    def test_data_publication_failure_does_not_publish_complete_report(self):
        runner = FakeMafftRunner(
            (("geison-00000000", "ACGT"), ("geison-00000001", "ACGT"), ("geison-00000002", "ACGA"))
        )
        with tempfile.TemporaryDirectory() as tmpdir, patch.object(Path, "replace", side_effect=OSError("replace failed")):
            with self.assertRaises(OSError):
                self._run(runner=runner, directory=tmpdir)
            self._assert_no_complete_report(tmpdir)

    def test_enabled_multiple_records_requires_an_injected_runner_for_now(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(MafftError, "runner.*Task 3|inject"):
                self._run(directory=tmpdir)
            self._assert_no_complete_report(tmpdir)


if __name__ == "__main__":
    unittest.main()
