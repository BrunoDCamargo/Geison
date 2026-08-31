import gc
import tempfile
import unittest
import warnings
from pathlib import Path

from qpcr_pipeline.local_input import load_fasta
from qpcr_pipeline.qc import QCStatus, evaluate_sequences


class FastaQcTests(unittest.TestCase):
    def test_load_fasta_closes_its_input_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "target.fasta"
            fasta_path.write_text(">seq-1\nACGT\n", encoding="utf-8")

            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter("always", ResourceWarning)
                records = load_fasta(fasta_path)
                gc.collect()

        self.assertEqual(tuple(record.sequence_id for record in records), ("seq-1",))
        self.assertEqual(
            [str(item.message) for item in captured if issubclass(item.category, ResourceWarning)],
            [],
        )

    def test_invalid_nucleotide_is_rejected_from_validated_sets(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "target.fasta"
            fasta_path.write_text(
                ">accepted\n"
                "ACGTACGT\n"
                ">invalid\n"
                "ACGTXCGT\n",
                encoding="utf-8",
            )

            records = load_fasta(fasta_path)
            result = evaluate_sequences(records)

        self.assertEqual(result.records[0].sequence_id, "accepted")
        self.assertEqual(result.records[0].status, QCStatus.ACCEPTED)
        self.assertEqual(result.records[0].reason_codes, ())

        self.assertEqual(result.records[1].sequence_id, "invalid")
        self.assertEqual(result.records[1].status, QCStatus.REJECTED)
        self.assertEqual(result.records[1].reason_codes, ("INVALID_NUCLEOTIDE",))

        self.assertEqual(result.target_sequence_set.sequence_ids, ("accepted",))
        self.assertEqual(result.evaluation_set.sequence_ids, ("accepted",))

    def test_sequence_shorter_than_configured_minimum_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "target.fasta"
            fasta_path.write_text(
                ">accepted\n"
                "ACGTACGT\n"
                ">too-short\n"
                "ACGT\n",
                encoding="utf-8",
            )

            records = load_fasta(fasta_path)
            result = evaluate_sequences(records, min_length=5)

        self.assertEqual(result.records[0].status, QCStatus.ACCEPTED)
        self.assertEqual(result.records[1].status, QCStatus.REJECTED)
        self.assertEqual(result.records[1].reason_codes, ("TOO_SHORT",))
        self.assertEqual(result.target_sequence_set.sequence_ids, ("accepted",))
        self.assertEqual(result.evaluation_set.sequence_ids, ("accepted",))

    def test_sequence_above_configured_ambiguity_fraction_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "target.fasta"
            fasta_path.write_text(
                ">accepted\n"
                "ACGTACGT\n"
                ">too-ambiguous\n"
                "ACGTNNNN\n",
                encoding="utf-8",
            )

            records = load_fasta(fasta_path)
            result = evaluate_sequences(records, max_ambiguous_fraction=0.25)

        self.assertEqual(result.records[0].status, QCStatus.ACCEPTED)
        self.assertEqual(result.records[1].status, QCStatus.REJECTED)
        self.assertEqual(result.records[1].reason_codes, ("EXCESSIVE_AMBIGUITY",))
        self.assertEqual(result.target_sequence_set.sequence_ids, ("accepted",))
        self.assertEqual(result.evaluation_set.sequence_ids, ("accepted",))

    def test_exact_duplicate_after_first_occurrence_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "target.fasta"
            fasta_path.write_text(
                ">first\n"
                "ACGTACGT\n"
                ">unique\n"
                "ACGTACGA\n"
                ">duplicate\n"
                "ACGTACGT\n",
                encoding="utf-8",
            )

            records = load_fasta(fasta_path)
            result = evaluate_sequences(records)

        self.assertEqual(result.records[0].status, QCStatus.ACCEPTED)
        self.assertEqual(result.records[1].status, QCStatus.ACCEPTED)
        self.assertEqual(result.records[2].status, QCStatus.REJECTED)
        self.assertEqual(result.records[2].reason_codes, ("DUPLICATE_SEQUENCE",))
        self.assertEqual(result.target_sequence_set.sequence_ids, ("first", "unique"))
        self.assertEqual(result.evaluation_set.sequence_ids, ("first", "unique"))

    def test_sequence_outside_expected_length_tolerance_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "target.fasta"
            fasta_path.write_text(
                ">expected\n"
                "ACGTACGT\n"
                ">length-outlier\n"
                "ACGTACGTACGT\n",
                encoding="utf-8",
            )

            records = load_fasta(fasta_path)
            result = evaluate_sequences(
                records,
                expected_length=8,
                length_tolerance_fraction=0.25,
            )

        self.assertEqual(result.records[0].status, QCStatus.ACCEPTED)
        self.assertEqual(result.records[1].status, QCStatus.REJECTED)
        self.assertEqual(result.records[1].reason_codes, ("INCONSISTENT_LENGTH",))
        self.assertEqual(result.target_sequence_set.sequence_ids, ("expected",))
        self.assertEqual(result.evaluation_set.sequence_ids, ("expected",))


if __name__ == "__main__":
    unittest.main()
