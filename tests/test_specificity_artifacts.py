import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import OffTargetConfig, SpecificityConfig
from qpcr_pipeline.off_targets import load_off_target_dataset, load_off_target_datasets
from qpcr_pipeline.primer_design import AssayCandidate, DesignedOligo, PrimerDesignResult
from qpcr_pipeline.specificity import evaluate_specificity


class SpecificityHelpers:
    @staticmethod
    def _oligo(sequence: str, start: int) -> DesignedOligo:
        return DesignedOligo(
            sequence=sequence,
            reference_start=start,
            reference_end=start + len(sequence) - 1,
            length=len(sequence),
            tm=60.0,
            gc_percent=50.0,
            penalty=None,
            metrics=(),
        )

    @classmethod
    def _assay(cls, assay_id: str = "a1") -> AssayCandidate:
        f = cls._oligo("ACGT", 1)
        p = cls._oligo("TTAA", 9)
        r = cls._oligo("AGTC", 17)
        return AssayCandidate(
            assay_id=assay_id,
            region_id="r1",
            primer3_index=0,
            forward_primer=f,
            probe=p,
            reverse_primer=r,
            product_size=20,
            pair_penalty=None,
            metrics=(),
        )

    @classmethod
    def _primer_result(cls, assays=None, status="COMPLETE") -> PrimerDesignResult:
        return PrimerDesignResult(
            status=status,
            reference_id="ref" if status == "COMPLETE" else None,
            candidates=(),
            assays=tuple(assays if assays is not None else (cls._assay(),)),
            candidate_regions_path=None,
            assays_path=None,
            primer3_input_path=None,
            primer3_output_path=None,
            report_path=Path("primer-design-report.json"),
        )


class OffTargetDatasetTests(unittest.TestCase):
    def test_loads_fasta_with_stable_sha_and_sequence_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "off.fasta"
            raw = b">z-first\nACGT\n>a-second\nTTTT\n"
            path.write_bytes(raw)
            dataset = load_off_target_dataset(OffTargetConfig(name="human", fasta=path))
        self.assertEqual(dataset.name, "human")
        self.assertEqual(dataset.source_type, "FASTA")
        self.assertEqual(dataset.sequence_ids, ("z-first", "a-second"))
        self.assertEqual(dataset.sha256, hashlib.sha256(raw).hexdigest())
        self.assertIsNone(dataset.frozen_manifest_path)

    def test_empty_fasta_is_valid_but_nonempty_malformed_fasta_is_not(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            empty = Path(tmpdir) / "empty.fasta"
            malformed = Path(tmpdir) / "bad.fasta"
            empty.write_text("", encoding="utf-8")
            malformed.write_text("ACGT without header\n", encoding="utf-8")
            self.assertEqual(
                load_off_target_dataset(OffTargetConfig(name="empty", fasta=empty)).records,
                (),
            )
            with self.assertRaisesRegex(ValueError, "bad.*FASTA"):
                load_off_target_dataset(OffTargetConfig(name="bad", fasta=malformed))

    def test_duplicate_fasta_ids_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "dup.fasta"
            path.write_text(">x\nACGT\n>x\nTTTT\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unique"):
                load_off_target_dataset(OffTargetConfig(name="dup", fasta=path))

    def test_frozen_dataset_uses_existing_validator_and_preserves_manifest(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            records_path = directory / "records.gb"
            manifest_path = directory / "dataset_manifest.json"
            record = SeqRecord(Seq("ACGT"), id="NC_1.2")
            record.annotations["molecule_type"] = "DNA"
            SeqIO.write((record,), records_path, "genbank")
            manifest = {
                "source": {"mode": "query", "query": "example[Organism]"},
                "resolved_entries": [
                    {"accession": "NC_1", "accession_version": "NC_1.2"}
                ],
            }
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with patch(
                "qpcr_pipeline.off_targets.validate_frozen_dataset",
                return_value=SimpleNamespace(
                    records_path=records_path, manifest_path=manifest_path
                ),
            ) as validate:
                dataset = load_off_target_dataset(
                    OffTargetConfig(name="neighbors", frozen_dataset=directory)
                )
        validate.assert_called_once_with(directory)
        self.assertEqual(dataset.source_type, "NCBI_FROZEN")
        self.assertEqual(dataset.sequence_ids, ("NC_1.2",))
        self.assertEqual(dataset.frozen_manifest["source"]["query"], "example[Organism]")
        self.assertEqual(
            dataset.frozen_manifest["resolved_entries"][0]["accession_version"],
            "NC_1.2",
        )

    def test_dataset_order_is_configuration_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            z = Path(tmpdir) / "z.fa"
            a = Path(tmpdir) / "a.fa"
            z.write_text(">z\nACGT\n", encoding="utf-8")
            a.write_text(">a\nACGT\n", encoding="utf-8")
            loaded = load_off_target_datasets(
                (
                    OffTargetConfig(name="z-first", fasta=z),
                    OffTargetConfig(name="a-second", fasta=a),
                )
            )
        self.assertEqual(tuple(item.name for item in loaded), ("z-first", "a-second"))


class SpecificityArtifactTests(SpecificityHelpers, unittest.TestCase):
    def test_disabled_does_not_read_datasets_and_publishes_only_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir)
            with patch(
                "qpcr_pipeline.specificity.load_off_target_datasets",
                side_effect=AssertionError("datasets must not be read"),
            ):
                result = evaluate_specificity(
                    self._primer_result(status="SKIPPED", assays=()),
                    (OffTargetConfig(name="bad", fasta=Path("missing.fa")),),
                    SpecificityConfig(enabled=False),
                    output,
                )
            self.assertEqual(result.status, "SKIPPED")
            self.assertEqual(
                {path.name for path in (output / "specificity").iterdir()},
                {"specificity_report.json"},
            )

    def test_complete_run_publishes_hits_amplicons_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "off.fa"
            fasta.write_text(">off-1\nACGTCCCCTTAACCCCGACT\n", encoding="utf-8")
            result = evaluate_specificity(
                self._primer_result(),
                (OffTargetConfig(name="human", fasta=fasta),),
                SpecificityConfig(
                    enabled=True,
                    max_primer_mismatches=0,
                    max_probe_mismatches=0,
                    max_amplicon_size=20,
                ),
                root / "out",
            )
            report = json.loads(result.report_path.read_text(encoding="utf-8"))
            hit_header = result.off_target_hits_path.read_text(encoding="utf-8").splitlines()[0]
            amp_header = result.plausible_amplicons_path.read_text(encoding="utf-8").splitlines()[0]
        self.assertEqual(result.status, "COMPLETE")
        self.assertTrue(any(item.detectable_off_target for item in result.amplicons))
        self.assertEqual(report["datasets"][0]["name"], "human")
        self.assertEqual(report["datasets"][0]["sequence_ids"], ["off-1"])
        self.assertIn("mismatch_positions", hit_header)
        self.assertIn("detectable_off_target", amp_header)

    def test_hit_truncation_never_removes_real_amplicon_or_its_coordinates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "off.fa"
            # Forward at 1 is retained but too far from reverse; forward at 15 is
            # truncated from individual hits and forms the only valid amplicon.
            sequence = "ACGT" + "C" * 10 + "ACGT" + "CC" + "TTAA" + "CC" + "GACT"
            fasta.write_text(f">off-1\n{sequence}\n", encoding="utf-8")
            result = evaluate_specificity(
                self._primer_result(),
                (OffTargetConfig(name="human", fasta=fasta),),
                SpecificityConfig(
                    enabled=True,
                    max_hits_per_oligo_per_dataset=1,
                    max_primer_mismatches=0,
                    max_probe_mismatches=0,
                    max_amplicon_size=16,
                ),
                root / "out",
            )
        risky = [item for item in result.amplicons if item.detectable_off_target]
        self.assertTrue(risky)
        later = next(item for item in risky if item.forward_source_start == 15)
        self.assertIsNone(later.forward_hit_rank)
        self.assertEqual((later.forward_source_start, later.forward_source_end), (15, 18))
        self.assertEqual((later.reverse_source_start, later.reverse_source_end), (27, 30))
        summary = next(
            item for item in result.retention
            if item.dataset_name == "human" and item.assay_id == "a1" and item.role == "FORWARD"
        )
        self.assertGreater(summary.total_hit_count, summary.retained_hit_count)
        self.assertTrue(summary.truncated)

    def test_enabled_empty_assays_produce_header_only_tsvs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "empty.fa"
            fasta.write_text("", encoding="utf-8")
            result = evaluate_specificity(
                self._primer_result(assays=()),
                (OffTargetConfig(name="empty", fasta=fasta),),
                SpecificityConfig(enabled=True),
                root / "out",
            )
            self.assertEqual(
                len(result.off_target_hits_path.read_text(encoding="utf-8").splitlines()), 1
            )
            self.assertEqual(
                len(result.plausible_amplicons_path.read_text(encoding="utf-8").splitlines()), 1
            )


if __name__ == "__main__":
    unittest.main()
