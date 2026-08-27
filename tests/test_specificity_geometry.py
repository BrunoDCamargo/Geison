import unittest

from qpcr_pipeline.config import SpecificityConfig
from qpcr_pipeline.specificity_matching import (
    GeometryAmplicon,
    MatchHit,
    find_plausible_amplicons,
)


class SpecificityGeometryTests(unittest.TestCase):
    @staticmethod
    def _hit(
        role: str,
        start: int,
        end: int,
        *,
        sequence_id: str = "s1",
        orientation: str = "FORWARD",
        dataset: str = "d1",
        assay: str = "a1",
        source_start: int | None = None,
        source_end: int | None = None,
    ) -> MatchHit:
        return MatchHit(
            dataset_name=dataset,
            assay_id=assay,
            sequence_id=sequence_id,
            role=role,
            orientation=orientation,
            source_start=start if source_start is None else source_start,
            source_end=end if source_end is None else source_end,
            oriented_start=start,
            oriented_end=end,
            mismatch_positions=(),
            mismatch_count=0,
            exact_match=True,
            three_prime_mismatch=False,
            compatible=True,
        )

    def test_isolated_or_incompatible_contexts_do_not_create_amplicons(self):
        config = SpecificityConfig(max_amplicon_size=100)
        cases = (
            (self._hit("FORWARD", 1, 4),),
            (
                self._hit("FORWARD", 1, 4, sequence_id="s1"),
                self._hit("REVERSE", 20, 23, sequence_id="s2"),
            ),
            (
                self._hit("FORWARD", 1, 4, orientation="FORWARD"),
                self._hit("REVERSE", 20, 23, orientation="REVERSE_COMPLEMENT"),
            ),
            (
                self._hit("REVERSE", 1, 4),
                self._hit("FORWARD", 20, 23),
            ),
            (
                self._hit("FORWARD", 1, 4),
                self._hit("REVERSE", 150, 153),
            ),
        )
        for hits in cases:
            with self.subTest(hits=hits):
                self.assertEqual(find_plausible_amplicons(hits, config), ())

    def test_plausible_primers_without_probe_are_not_detectable(self):
        result = find_plausible_amplicons(
            (self._hit("FORWARD", 3, 6), self._hit("REVERSE", 15, 18)),
            SpecificityConfig(max_amplicon_size=20),
        )
        self.assertEqual(len(result), 1)
        self.assertTrue(result[0].primer_amplicon_plausible)
        self.assertFalse(result[0].detectable_off_target)
        self.assertEqual(result[0].amplicon_size, 16)

    def test_probe_must_be_strictly_inside_primer_geometry(self):
        base = (self._hit("FORWARD", 3, 6), self._hit("REVERSE", 15, 18))
        outside = find_plausible_amplicons(
            base + (self._hit("PROBE", 19, 22),),
            SpecificityConfig(max_amplicon_size=30),
        )[0]
        inside = find_plausible_amplicons(
            base + (self._hit("PROBE", 9, 12),),
            SpecificityConfig(max_amplicon_size=30),
        )[0]
        self.assertFalse(outside.detectable_off_target)
        self.assertTrue(inside.detectable_off_target)
        self.assertEqual(len(inside.probes), 1)

    def test_preserves_multiple_valid_amplicons(self):
        hits = (
            self._hit("FORWARD", 1, 4),
            self._hit("FORWARD", 5, 8),
            self._hit("REVERSE", 20, 23),
            self._hit("REVERSE", 24, 27),
            self._hit("PROBE", 12, 15),
        )
        result = find_plausible_amplicons(hits, SpecificityConfig(max_amplicon_size=40))
        self.assertEqual(len(result), 4)
        self.assertTrue(all(item.detectable_off_target for item in result))

    def test_reverse_complement_geometry_keeps_source_interval_normalized(self):
        forward = self._hit(
            "FORWARD", 3, 6, orientation="REVERSE_COMPLEMENT",
            source_start=15, source_end=18,
        )
        probe = self._hit(
            "PROBE", 9, 12, orientation="REVERSE_COMPLEMENT",
            source_start=9, source_end=12,
        )
        reverse = self._hit(
            "REVERSE", 15, 18, orientation="REVERSE_COMPLEMENT",
            source_start=3, source_end=6,
        )
        result = find_plausible_amplicons(
            (forward, probe, reverse), SpecificityConfig(max_amplicon_size=30)
        )
        self.assertEqual(len(result), 1)
        self.assertEqual((result[0].source_start, result[0].source_end), (3, 18))
        self.assertTrue(result[0].detectable_off_target)

    def test_exact_max_amplicon_boundary_is_accepted(self):
        result = find_plausible_amplicons(
            (self._hit("FORWARD", 1, 4), self._hit("REVERSE", 17, 20)),
            SpecificityConfig(max_amplicon_size=20),
        )
        self.assertEqual(result[0].amplicon_size, 20)


if __name__ == "__main__":
    unittest.main()
