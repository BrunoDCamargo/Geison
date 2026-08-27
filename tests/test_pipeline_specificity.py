import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import (
    AlignmentConfig,
    ConservationConfig,
    InclusivityConfig,
    OffTargetConfig,
    PipelineConfig,
    PrimerDesignConfig,
    SpecificityConfig,
)
from qpcr_pipeline.inclusivity import evaluate_inclusivity
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet, EvaluationSet
from qpcr_pipeline.pipeline import run_pipeline
from qpcr_pipeline.primer_design import AssayCandidate, DesignedOligo, PrimerDesignResult
from qpcr_pipeline.specificity_matching import enumerate_compatible_hits


FIXTURE_FASTA = Path("tests/fixtures/target_small.fasta")


class PipelineSpecificityTests(unittest.TestCase):
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
    def _assay(cls, forward: str = "ACGT") -> AssayCandidate:
        f = cls._oligo(forward, 1)
        p = cls._oligo("TTAA", 9)
        r = cls._oligo("AGTC", 17)
        return AssayCandidate(
            assay_id="a1",
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
    def _primer_result(cls, assay=None) -> PrimerDesignResult:
        return PrimerDesignResult(
            status="COMPLETE",
            reference_id="seq-1",
            candidates=(),
            assays=(assay or cls._assay(),),
            candidate_regions_path=None,
            assays_path=None,
            primer3_input_path=None,
            primer3_output_path=None,
            report_path=Path("unused-primer-report.json"),
        )

    @staticmethod
    def _inclusivity_skipped():
        return SimpleNamespace(
            status="SKIPPED",
            evaluation_sequence_ids=(),
            assay_results=(),
        )

    def _upstream_patches(self, primer_result):
        return (
            patch(
                "qpcr_pipeline.pipeline.cluster_sequences",
                return_value=SimpleNamespace(discovery_set=DiscoverySet(("seq-1",))),
            ),
            patch(
                "qpcr_pipeline.pipeline.align_discovery",
                return_value=SimpleNamespace(
                    status="COMPLETE",
                    reference_id="seq-1",
                    reference_mode="automatic",
                ),
            ),
            patch(
                "qpcr_pipeline.pipeline.analyze_conservation",
                return_value=SimpleNamespace(
                    status="COMPLETE",
                    reference_id="seq-1",
                    positions=(),
                    windows=(),
                ),
            ),
            patch("qpcr_pipeline.pipeline.design_primers", return_value=primer_result),
            patch(
                "qpcr_pipeline.pipeline.evaluate_inclusivity",
                return_value=self._inclusivity_skipped(),
            ),
        )

    def test_default_pipeline_publishes_skipped_specificity_summary(self):
        config = PipelineConfig(target_name="target", input_fasta=FIXTURE_FASTA)
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "out"
            run_pipeline(config, output)
            qc = json.loads((output / "qc_report.json").read_text(encoding="utf-8"))
        self.assertEqual(
            qc["specificity"],
            {
                "status": "SKIPPED",
                "dataset_count": 0,
                "sequence_count": 0,
                "assay_count": 0,
                "retained_hit_count": 0,
                "plausible_amplicon_count": 0,
                "detectable_off_target_count": 0,
            },
        )

    def test_enabled_specificity_runs_after_primer_design_even_when_inclusivity_disabled(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            off = root / "off.fa"
            off.write_text(">off-1\nACGTCCCCTTAACCCCGACT\n", encoding="utf-8")
            config = PipelineConfig(
                target_name="target",
                input_fasta=FIXTURE_FASTA,
                alignment=AlignmentConfig(enabled=True),
                conservation=ConservationConfig(enabled=True),
                primer_design=PrimerDesignConfig(enabled=True),
                inclusivity=InclusivityConfig(enabled=False),
                off_targets=(OffTargetConfig(name="human", fasta=off),),
                specificity=SpecificityConfig(
                    enabled=True,
                    max_primer_mismatches=0,
                    max_probe_mismatches=0,
                    max_amplicon_size=20,
                ),
            )
            primer_result = self._primer_result()
            patches = self._upstream_patches(primer_result)
            with patches[0], patches[1], patches[2], patches[3], patches[4]:
                output = root / "out"
                run_pipeline(config, output)
            qc = json.loads((output / "qc_report.json").read_text(encoding="utf-8"))
        self.assertEqual(qc["specificity"]["status"], "COMPLETE")
        self.assertEqual(qc["specificity"]["dataset_count"], 1)
        self.assertGreaterEqual(qc["specificity"]["plausible_amplicon_count"], 1)
        self.assertGreaterEqual(qc["specificity"]["detectable_off_target_count"], 1)

    def test_frozen_specificity_path_never_calls_live_ncbi_acquisition(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            frozen = root / "frozen"
            frozen.mkdir()
            records_path = frozen / "records.gb"
            manifest_path = frozen / "dataset_manifest.json"
            record = SeqRecord(Seq("ACGTCCCCTTAACCCCGACT"), id="NC_1.2")
            record.annotations["molecule_type"] = "DNA"
            SeqIO.write((record,), records_path, "genbank")
            manifest_path.write_text(
                json.dumps(
                    {
                        "source": {"mode": "query", "query": "example[Organism]"},
                        "resolved_entries": [
                            {"accession": "NC_1", "accession_version": "NC_1.2"}
                        ],
                    }
                ),
                encoding="utf-8",
            )
            from qpcr_pipeline.specificity import evaluate_specificity

            with patch(
                "qpcr_pipeline.off_targets.validate_frozen_dataset",
                return_value=SimpleNamespace(
                    records_path=records_path, manifest_path=manifest_path
                ),
            ), patch(
                "qpcr_pipeline.ncbi.acquire_ncbi_dataset",
                side_effect=AssertionError("live acquisition must not run"),
            ), patch(
                "qpcr_pipeline.ncbi.BioEntrezClient.from_environment",
                side_effect=AssertionError("Entrez client must not be created"),
            ):
                result = evaluate_specificity(
                    self._primer_result(),
                    (OffTargetConfig(name="frozen", frozen_dataset=frozen),),
                    SpecificityConfig(
                        enabled=True,
                        max_primer_mismatches=0,
                        max_probe_mismatches=0,
                        max_amplicon_size=20,
                    ),
                    root / "out",
                )
        self.assertEqual(result.status, "COMPLETE")

    def test_specificity_matches_public_inclusivity_iupac_semantics(self):
        cases = (
            ("ARGT", "AAGT", 0, True),
            ("ARGT", "AGGT", 0, True),
            ("ARGT", "ANGT", 0, False),
            ("ACGT", "ATGT", 1, True),
            ("ACGT", "ACGA", 1, False),
        )
        for oligo, target, max_mismatches, expected in cases:
            with self.subTest(
                oligo=oligo,
                target=target,
                max_mismatches=max_mismatches,
            ):
                assay = self._assay(forward=oligo)
                specificity_config = SpecificityConfig(
                    max_primer_mismatches=max_mismatches,
                    reject_primer_3_prime_mismatch=True,
                    primer_3_prime_bases=2,
                )
                specificity_ok = bool(
                    enumerate_compatible_hits(
                        "d",
                        LocalSequenceRecord("s1", target),
                        assay,
                        "FORWARD",
                        specificity_config,
                    )
                )
                with tempfile.TemporaryDirectory() as tmpdir:
                    inclusivity = evaluate_inclusivity(
                        (LocalSequenceRecord("s1", target),),
                        EvaluationSet(("s1",)),
                        self._primer_result(assay),
                        InclusivityConfig(
                            enabled=True,
                            search_flank=0,
                            max_hits_per_oligo=20,
                            max_primer_mismatches=max_mismatches,
                            reject_primer_3_prime_mismatch=True,
                            primer_3_prime_bases=2,
                        ),
                        Path(tmpdir),
                    )
                inclusivity_ok = any(
                    hit.role == "FORWARD" and hit.compatible
                    for hit in inclusivity.oligo_matches
                )
                self.assertEqual(specificity_ok, expected)
                self.assertEqual(specificity_ok, inclusivity_ok)


if __name__ == "__main__":
    unittest.main()
