import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

from qpcr_pipeline.alignment import (
    AlignedSequence,
    AlignmentCoordinate,
    AlignmentResult,
)
from qpcr_pipeline.config import ConservationConfig
from qpcr_pipeline.conservation import (
    ConservationError,
    ReferenceAnnotation,
    WindowConservation,
    analyze_conservation,
)
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet


POSITION_HEADER = (
    "alignment_position\treference_position\treference_base\tdepth\tcoverage\t"
    "frequency_a\tfrequency_c\tfrequency_g\tfrequency_t\tgap_frequency\t"
    "major_allele_frequency\tentropy_bits\tmajor_consensus\tiupac_consensus"
)
WINDOW_HEADER = (
    "reference_start\treference_end\tposition_count\tmean_conservation\t"
    "minimum_conservation\tmean_coverage\tmean_gap_frequency\tmean_entropy_bits"
)


def record(sequence_id, sequence, *, features=()):
    metadata = {} if features is None else {"features": tuple(features)}
    return LocalSequenceRecord(sequence_id, sequence, metadata)


def alignment(*sequence_items, reference_id="ref", status="COMPLETE", coordinates=None):
    if len(sequence_items) == 1 and isinstance(sequence_items[0], tuple) and (
        not sequence_items[0] or isinstance(sequence_items[0][0], tuple)
    ):
        sequence_items = sequence_items[0]
    sequence_items = tuple(sequence_items)
    ids = tuple(sequence_id for sequence_id, _ in sequence_items)
    sequences = tuple(
        AlignedSequence(sequence_id, sequence, "forward")
        for sequence_id, sequence in sequence_items
    )
    if coordinates is None:
        if not sequence_items:
            coordinates = ()
        else:
            reference_sequence = dict(sequence_items)[reference_id]
            reference_position = 0
            built = []
            for alignment_position, base in enumerate(reference_sequence, 1):
                if base == "-":
                    built.append(AlignmentCoordinate(alignment_position, None, None))
                else:
                    reference_position += 1
                    built.append(
                        AlignmentCoordinate(alignment_position, reference_position, base)
                    )
            coordinates = tuple(built)
    return AlignmentResult(
        status=status,
        discovery_set=DiscoverySet(ids),
        reference_id=reference_id if sequence_items else None,
        reference_mode="explicit" if sequence_items else None,
        sequences=sequences,
        coordinates=tuple(coordinates),
        alignment_fasta_path=None,
        coordinate_map_path=None,
        report_path=Path("unused-alignment-report.json"),
    )


def matching_records(alignment_result, *, features=()):
    if not alignment_result.sequences:
        return ()
    records = []
    for item in alignment_result.sequences:
        sequence = item.aligned_sequence.replace("-", "")
        records.append(
            record(
                item.sequence_id,
                sequence,
                features=features if item.sequence_id == alignment_result.reference_id else (),
            )
        )
    return tuple(records)


class _ConservationTestCase(unittest.TestCase):
    def run_analysis(self, alignment_result, *, records=None, config=None, directory=None):
        if records is None:
            records = matching_records(alignment_result)
        if config is None:
            config = ConservationConfig(enabled=True, window_size=5, step_size=4)
        if directory is None:
            self._temporary_directory = tempfile.TemporaryDirectory()
            self.addCleanup(self._temporary_directory.cleanup)
            directory = self._temporary_directory.name
        return analyze_conservation(
            records, alignment_result, config, Path(directory), target_name="target <one>"
        )


class ConservationCalculationTests(_ConservationTestCase):

    def test_known_frequencies_entropy_and_reference_aware_consensus(self):
        result = self.run_analysis(
            alignment(("ref", "AAA"), ("s2", "AAC"), ("s3", "ACG"), ("s4", "ACT"))
        )

        self.assertEqual(
            [(p.frequency_a, p.frequency_c, p.frequency_g, p.frequency_t) for p in result.positions],
            [(1.0, 0.0, 0.0, 0.0), (0.5, 0.5, 0.0, 0.0), (0.25, 0.25, 0.25, 0.25)],
        )
        self.assertEqual([p.entropy_bits for p in result.positions], [0.0, 1.0, 2.0])
        self.assertEqual([p.major_allele_frequency for p in result.positions], [1.0, 0.5, 0.25])
        self.assertEqual(result.major_consensus, "AAA")
        self.assertEqual(result.iupac_consensus, "AMN")

    def test_reference_gap_column_is_measured_but_omitted_from_consensus_artifacts(self):
        result = self.run_analysis(alignment(("ref", "A-C"), ("s2", "ATC"), ("s3", "A-C")))

        insertion = result.positions[1]
        self.assertIsNone(insertion.reference_position)
        self.assertEqual(insertion.depth, 1)
        self.assertEqual(insertion.coverage, 1 / 3)
        self.assertEqual(insertion.gap_frequency, 2 / 3)
        self.assertEqual(insertion.frequency_t, 1.0)
        self.assertEqual(result.major_consensus, "AC")
        self.assertEqual(result.iupac_consensus, "AC")
        self.assertEqual(result.major_consensus_path.read_text(encoding="utf-8"), ">geison-major-consensus\nAC\n")
        self.assertEqual(result.iupac_consensus_path.read_text(encoding="utf-8"), ">geison-iupac-consensus\nAC\n")

    def test_fractional_iupac_observations_have_exact_frequencies(self):
        result = self.run_analysis(alignment(("ref", "ACG"), ("s2", "RBN")))
        expected = (
            (0.75, 0.0, 0.25, 0.0, "R"),
            (0.0, 2 / 3, 1 / 6, 1 / 6, "B"),
            (1 / 8, 1 / 8, 5 / 8, 1 / 8, "N"),
        )
        self.assertEqual(
            tuple(
                (p.frequency_a, p.frequency_c, p.frequency_g, p.frequency_t, p.iupac_consensus)
                for p in result.positions
            ),
            expected,
        )

    def test_exact_fractional_tie_prefers_the_reference_base(self):
        result = self.run_analysis(
            alignment(
                ("ref", "C"),
                ("s2", "A"),
                ("s3", "A"),
                ("s4", "B"),
                ("s5", "B"),
                ("s6", "B"),
                ("s7", "H"),
            )
        )

        position = result.positions[0]
        self.assertEqual(position.frequency_a, 1 / 3)
        self.assertEqual(position.frequency_c, 1 / 3)
        self.assertEqual(position.major_allele_frequency, 1 / 3)
        self.assertEqual(position.major_consensus, "C")
        self.assertEqual(result.major_consensus, "C")

    def test_every_iupac_symbol_expands_to_its_canonical_support_set(self):
        expected_support = {
            "A": "A", "C": "C", "G": "G", "T": "T",
            "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
            "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
        }
        for symbol, support in expected_support.items():
            with self.subTest(symbol=symbol):
                result = self.run_analysis(alignment(("ref", "A"), ("s2", symbol)))
                observed = {
                    base
                    for base, frequency in zip(
                        "ACGT",
                        (
                            result.positions[0].frequency_a,
                            result.positions[0].frequency_c,
                            result.positions[0].frequency_g,
                            result.positions[0].frequency_t,
                        ),
                    )
                    if frequency > 0
                }
                self.assertEqual(observed, set(support) | {"A"})

    def test_gaps_are_excluded_from_entropy(self):
        result = self.run_analysis(alignment(("ref", "A"), ("s2", "C"), ("s3", "-"), ("s4", "-")))
        self.assertEqual(result.positions[0].entropy_bits, 1.0)

    def test_all_gap_column_fails_before_publication(self):
        coordinates = (AlignmentCoordinate(1, None, None),)
        bad_alignment = alignment(("ref", "-"), ("s2", "-"), coordinates=coordinates)
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ConservationError, "all-gap"):
                self.run_analysis(
                    bad_alignment,
                    records=(record("ref", ""), record("s2", "")),
                    directory=tmpdir,
                )
            self.assertFalse((Path(tmpdir) / "conservation" / "conservation_report.json").exists())

    def test_windows_use_regular_starts_and_one_anchored_tail(self):
        result = self.run_analysis(
            alignment(
                ("ref", "AAAAAAAAAAAA"),
                ("s2", "AACAAAAAAAAA"),
                ("s3", "ACG-AAAAAAAA"),
                ("s4", "ACT-AAAAAAAA"),
            )
        )
        self.assertEqual([(w.reference_start, w.reference_end) for w in result.windows], [(1, 5), (5, 9), (8, 12)])
        self.assertEqual(
            result.windows[0],
            WindowConservation(1, 5, 5, 0.75, 0.25, 0.9, 0.1, 0.6),
        )

    def test_empty_short_and_step_one_window_boundaries(self):
        empty = self.run_analysis(alignment(()))
        short = self.run_analysis(alignment(("ref", "AAAA"), ("s2", "AAAA")))
        step_one = self.run_analysis(
            alignment(("ref", "AAAAAAA"), ("s2", "AAAAAAA")),
            config=ConservationConfig(enabled=True, window_size=5, step_size=1),
        )
        self.assertEqual(empty.windows, ())
        self.assertEqual([(w.reference_start, w.reference_end) for w in short.windows], [(1, 4)])
        self.assertEqual(
            [(w.reference_start, w.reference_end) for w in step_one.windows],
            [(1, 5), (2, 6), (3, 7)],
        )


class ConservationAnnotationTests(_ConservationTestCase):
    def test_simple_and_compound_features_are_clipped_labeled_and_sorted(self):
        features = (
            SeqFeature(SimpleLocation(1, 5, strand=1), type="gene", qualifiers={"gene": ["abc"], "locus_tag": ["ignored"]}),
            SeqFeature(
                CompoundLocation((SimpleLocation(10, 14, strand=-1), SimpleLocation(6, 8, strand=-1))),
                type="CDS",
                qualifiers={"product": ["protein X"]},
            ),
        )
        aligned = alignment(("ref", "A" * 12), ("s2", "A" * 12))
        result = self.run_analysis(aligned, records=matching_records(aligned, features=features))
        self.assertEqual(
            result.annotations,
            (
                ReferenceAnnotation("gene", 2, 5, 1, "abc"),
                ReferenceAnnotation("CDS", 7, 8, -1, "protein X"),
                ReferenceAnnotation("CDS", 11, 12, -1, "protein X"),
            ),
        )

    def test_annotation_label_precedence_and_skipped_counts_are_reported(self):
        external = SimpleLocation(0, 2, ref="external")
        features = (
            SeqFeature(SimpleLocation(0, 1), type="misc", qualifiers={"gene": ["g"], "locus_tag": ["l"], "product": ["p"]}),
            SeqFeature(SimpleLocation(1, 2), type="misc", qualifiers={"locus_tag": ["l"], "product": ["p"]}),
            SeqFeature(SimpleLocation(2, 3), type="misc", qualifiers={"product": ["p"]}),
            SeqFeature(SimpleLocation(3, 4), type="repeat_region"),
            SeqFeature(SimpleLocation(0, 4), type="source"),
            SeqFeature(external, type="gene", qualifiers={"gene": ["external"]}),
            SeqFeature(SimpleLocation(2, 2), type="gene"),
            SeqFeature(SimpleLocation(10, 12), type="gene"),
            "malformed",
        )
        aligned = alignment(("ref", "AAAA"), ("s2", "AAAA"))
        result = self.run_analysis(aligned, records=matching_records(aligned, features=features))
        report = json.loads(result.report_path.read_text(encoding="utf-8"))

        self.assertEqual([item.label for item in result.annotations], ["g", "l", "p", "repeat_region"])
        self.assertEqual(report["annotation_counts"], {"published": 4, "skipped": 5})

    def test_reference_without_feature_metadata_has_no_annotations(self):
        aligned = alignment(("ref", "AAAA"), ("s2", "AAAA"))
        records = (record("ref", "AAAA", features=None), record("s2", "AAAA", features=None))
        self.assertEqual(self.run_analysis(aligned, records=records).annotations, ())


class ConservationValidationTests(_ConservationTestCase):
    def assert_invalid(self, alignment_result, *, records=None, message=None):
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.assertRaisesRegex(ConservationError, message) if message else self.assertRaises(ConservationError)
            with context:
                self.run_analysis(alignment_result, records=records, directory=tmpdir)
            self.assertFalse((Path(tmpdir) / "conservation" / "conservation_report.json").exists())

    def test_enabled_requires_complete_alignment(self):
        self.assert_invalid(alignment((), status="SKIPPED"), message="COMPLETE")

    def test_alignment_ids_must_equal_discovery_once_and_in_order(self):
        base = alignment(("ref", "A"), ("s2", "A"))
        cases = (
            ("duplicate", (AlignedSequence("ref", "A", "forward"), AlignedSequence("ref", "A", "forward"))),
            ("missing", (AlignedSequence("ref", "A", "forward"),)),
            ("reordered", tuple(reversed(base.sequences))),
        )
        for name, sequences in cases:
            with self.subTest(name=name):
                bad = AlignmentResult(base.status, base.discovery_set, base.reference_id, base.reference_mode, sequences, base.coordinates, None, None, base.report_path)
                self.assert_invalid(bad, records=(record("ref", "A"), record("s2", "A")), message="IDs")

    def test_record_membership_must_exactly_match_discovery(self):
        base = alignment(("ref", "A"), ("s2", "A"))
        for name, records in (
            ("duplicate", (record("ref", "A"), record("ref", "A"))),
            ("missing", (record("ref", "A"),)),
            ("unknown", (record("ref", "A"), record("other", "A"))),
        ):
            with self.subTest(name=name):
                self.assert_invalid(base, records=records, message="records")

    def test_rejects_unequal_lengths_unknown_symbols_and_invalid_coordinates(self):
        unequal = alignment(("ref", "AA"), ("s2", "A"))
        unknown = alignment(("ref", "A"), ("s2", "Z"))
        coordinate_base = alignment(("ref", "AA"), ("s2", "AA"))
        invalid_coordinates = (
            (AlignmentCoordinate(2, 1, "A"), AlignmentCoordinate(2, 2, "A")),
            (AlignmentCoordinate(1, 1, "C"), AlignmentCoordinate(2, 2, "A")),
            (AlignmentCoordinate(1, 2, "A"), AlignmentCoordinate(2, 1, "A")),
        )
        self.assert_invalid(unequal, message="length")
        self.assert_invalid(unknown, message="symbol")
        for coordinates in invalid_coordinates:
            bad = AlignmentResult(
                coordinate_base.status, coordinate_base.discovery_set, coordinate_base.reference_id,
                coordinate_base.reference_mode, coordinate_base.sequences, coordinates,
                None, None, coordinate_base.report_path,
            )
            self.assert_invalid(bad, message="coordinate")

    def test_rejects_missing_reversed_or_mismatched_reference_records(self):
        base = alignment(("ref", "AC"), ("s2", "AC"))
        missing_reference = AlignmentResult(base.status, base.discovery_set, None, base.reference_mode, base.sequences, base.coordinates, None, None, base.report_path)
        reversed_reference = AlignmentResult(
            base.status, base.discovery_set, base.reference_id, base.reference_mode,
            (AlignedSequence("ref", "GT", "reverse_complemented"), base.sequences[1]),
            base.coordinates, None, None, base.report_path,
        )
        self.assert_invalid(missing_reference, message="reference")
        self.assert_invalid(reversed_reference, message="reference")
        self.assert_invalid(base, records=(record("ref", "AG"), record("s2", "AC")), message="reference")


class ConservationArtifactTests(_ConservationTestCase):
    def test_enabled_empty_publishes_empty_artifacts_and_complete_report(self):
        result = self.run_analysis(alignment(()))
        self.assertEqual(result.position_metrics_path.read_text(encoding="utf-8"), POSITION_HEADER + "\n")
        self.assertEqual(result.window_metrics_path.read_text(encoding="utf-8"), WINDOW_HEADER + "\n")
        self.assertEqual(result.major_consensus_path.read_text(encoding="utf-8"), "")
        self.assertEqual(result.iupac_consensus_path.read_text(encoding="utf-8"), "")
        self.assertIn("No conservation windows are available.", result.html_report_path.read_text(encoding="utf-8"))
        self.assertEqual(json.loads(result.report_path.read_text(encoding="utf-8"))["status"], "COMPLETE")

    def test_disabled_does_not_read_sequences_removes_only_exact_stale_paths(self):
        class UnreadableAlignment:
            discovery_set = DiscoverySet(("ref",))

            @property
            def sequences(self):
                raise AssertionError("disabled conservation read alignment sequences")

        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            conservation_dir = base / "conservation"
            conservation_dir.mkdir()
            stale_names = ("position_metrics.tsv", "window_metrics.tsv", "consensus_major.fasta", "consensus_iupac.fasta")
            for name in stale_names:
                (conservation_dir / name).write_text("stale", encoding="utf-8")
            (base / "report.html").write_text("stale", encoding="utf-8")
            sibling = conservation_dir / "keep.txt"
            sibling.write_text("keep", encoding="utf-8")

            result = analyze_conservation(
                (), UnreadableAlignment(), ConservationConfig(enabled=False), base, target_name="target"
            )
            report = json.loads(result.report_path.read_text(encoding="utf-8"))

            self.assertEqual(result.status, "SKIPPED")
            self.assertEqual(set(conservation_dir.iterdir()), {sibling, result.report_path})
            self.assertFalse((base / "report.html").exists())
            self.assertEqual(report["status"], "SKIPPED")
            self.assertTrue(all(value is None for value in report["artifacts"].values()))

    def test_tsv_headers_null_fields_and_report_schema_are_traceable(self):
        aligned = alignment(("ref", "A-C"), ("s2", "ATC"), ("s3", "A-C"))
        result = self.run_analysis(aligned)
        position_lines = result.position_metrics_path.read_text(encoding="utf-8").splitlines()
        window_lines = result.window_metrics_path.read_text(encoding="utf-8").splitlines()
        report = json.loads(result.report_path.read_text(encoding="utf-8"))

        self.assertEqual(position_lines[0], POSITION_HEADER)
        self.assertTrue(position_lines[2].startswith("2\t\t\t1\t"))
        self.assertEqual(window_lines[0], WINDOW_HEADER)
        self.assertEqual(report["schema_version"], 1)
        self.assertEqual(report["enabled"], True)
        self.assertEqual(report["window_parameters"], {"window_size": 5, "step_size": 4})
        self.assertEqual(
            report["metric_definitions"],
            {
                "base_frequencies": "fractional IUPAC support normalized by non-gap depth",
                "conservation": "major allele frequency",
                "coverage": "non-gap depth divided by sequence count",
                "entropy_bits": "Shannon entropy of A/C/G/T frequencies in bits",
                "gap_frequency": "gap count divided by sequence count",
            },
        )
        self.assertEqual(
            report["counts"],
            {"sequences": 3, "alignment_columns": 3, "reference_positions": 2, "windows": 1, "annotations": 0},
        )
        self.assertEqual(report["annotation_counts"], {"published": 0, "skipped": 0})
        self.assertEqual(report["consensus_lengths"], {"major": 2, "iupac": 2})
        self.assertEqual(
            report["artifacts"],
            {
                "position_metrics": "conservation/position_metrics.tsv",
                "window_metrics": "conservation/window_metrics.tsv",
                "major_consensus": "conservation/consensus_major.fasta",
                "iupac_consensus": "conservation/consensus_iupac.fasta",
                "html_report": "report.html",
            },
        )

    def test_consensus_headers_do_not_leak_original_ids_and_html_escapes_supplied_text(self):
        feature = SeqFeature(SimpleLocation(0, 1), type="gene", qualifiers={"gene": ["</script><b>bad</b>"]})
        aligned = alignment((("reference-secret", "A"), ("member-secret", "A")), reference_id="reference-secret")
        result = analyze_conservation(
            matching_records(aligned, features=(feature,)), aligned,
            ConservationConfig(enabled=True, window_size=1, step_size=1),
            Path(self._new_tempdir()), target_name="<img src=x onerror=bad>",
        )
        fasta = result.major_consensus_path.read_text(encoding="utf-8")
        html = result.html_report_path.read_text(encoding="utf-8")
        self.assertEqual(fasta, ">geison-major-consensus\nA\n")
        self.assertNotIn("reference-secret", fasta)
        self.assertNotIn("member-secret", fasta)
        self.assertNotIn("</script><b>bad</b>", html)
        self.assertNotIn("<img src=x onerror=bad>", html)
        self.assertIn("textContent", html)

    def _new_tempdir(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        return temporary.name

    def test_each_enabled_artifact_replacement_is_hardlink_safe(self):
        aligned = alignment(("ref", "AAAAAA"), ("s2", "AACAAA"))
        destinations = (
            Path("conservation/position_metrics.tsv"),
            Path("conservation/window_metrics.tsv"),
            Path("conservation/consensus_major.fasta"),
            Path("conservation/consensus_iupac.fasta"),
            Path("report.html"),
            Path("conservation/conservation_report.json"),
        )
        for relative in destinations:
            with self.subTest(path=str(relative)), tempfile.TemporaryDirectory() as tmpdir:
                base = Path(tmpdir)
                source = base / "source.txt"
                source.write_text("do not mutate", encoding="utf-8")
                destination = base / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                os.link(source, destination)
                self.run_analysis(aligned, directory=base)
                self.assertEqual(source.read_text(encoding="utf-8"), "do not mutate")

    def test_prior_report_is_invalidated_before_every_enabled_data_replacement(self):
        aligned = alignment(("ref", "AAAAAA"), ("s2", "AACAAA"))
        failing_names = (
            "position_metrics.tsv", "window_metrics.tsv", "consensus_major.fasta",
            "consensus_iupac.fasta", "report.html",
        )
        for failing_name in failing_names:
            with self.subTest(failing_name=failing_name), tempfile.TemporaryDirectory() as tmpdir:
                base = Path(tmpdir)
                prior = self.run_analysis(aligned, directory=base)
                failing_path = base / (failing_name if failing_name == "report.html" else f"conservation/{failing_name}")
                original_replace = Path.replace

                def fail_selected(source, destination):
                    if Path(destination) == failing_path:
                        raise OSError(f"failed {failing_name}")
                    return original_replace(source, destination)

                with patch.object(Path, "replace", new=fail_selected):
                    with self.assertRaisesRegex(OSError, "failed"):
                        self.run_analysis(aligned, directory=base)
                self.assertFalse(prior.report_path.exists())

    def test_prior_report_is_invalidated_during_every_disabled_publication_failure(self):
        aligned = alignment(("ref", "AAAAAA"), ("s2", "AACAAA"))
        stale_names = (
            "position_metrics.tsv", "window_metrics.tsv", "consensus_major.fasta",
            "consensus_iupac.fasta", "report.html",
        )
        cases = tuple(("unlink", name) for name in stale_names) + (("replace", "conservation_report.json"),)
        for operation, failing_name in cases:
            with self.subTest(operation=operation, failing_name=failing_name), tempfile.TemporaryDirectory() as tmpdir:
                base = Path(tmpdir)
                prior = self.run_analysis(aligned, directory=base)
                failing_path = base / (
                    failing_name if failing_name == "report.html" else f"conservation/{failing_name}"
                )
                if operation == "unlink":
                    original = Path.unlink

                    def fail_selected(path, *args, **kwargs):
                        if Path(path) == failing_path:
                            raise OSError(f"failed {failing_name}")
                        return original(path, *args, **kwargs)

                    failure_patch = patch.object(Path, "unlink", new=fail_selected)
                else:
                    original = Path.replace

                    def fail_selected(source, destination):
                        if Path(destination) == failing_path:
                            raise OSError(f"failed {failing_name}")
                        return original(source, destination)

                    failure_patch = patch.object(Path, "replace", new=fail_selected)
                with failure_patch:
                    with self.assertRaisesRegex(OSError, "failed"):
                        analyze_conservation(
                            (), aligned, ConservationConfig(enabled=False), base, target_name="target"
                        )
                self.assertFalse(prior.report_path.exists())

    def test_validation_and_renderer_failures_publish_no_new_complete_report(self):
        aligned = alignment(("ref", "AA"), ("s2", "AA"))
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ConservationError):
                self.run_analysis(aligned, records=(record("ref", "AA"),), directory=tmpdir)
            self.assertFalse((Path(tmpdir) / "conservation" / "conservation_report.json").exists())
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch("qpcr_pipeline.conservation.render_conservation_html", side_effect=RuntimeError("renderer failed")):
                with self.assertRaisesRegex(RuntimeError, "renderer failed"):
                    self.run_analysis(aligned, directory=tmpdir)
            self.assertFalse((Path(tmpdir) / "conservation" / "conservation_report.json").exists())


if __name__ == "__main__":
    unittest.main()
