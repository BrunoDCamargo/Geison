import unittest

from qpcr_pipeline.iupac import (
    IupacError,
    IUPAC_COMPLEMENT,
    IUPAC_SUPPORT,
    iupac_support,
    minimal_iupac_symbol,
    mismatch_positions,
    normalize_iupac,
    reverse_complement_iupac,
    sequence_degeneracy,
)


class IupacTests(unittest.TestCase):
    def test_support_and_minimal_symbols_cover_all_nonempty_base_sets(self):
        self.assertEqual(iupac_support("R"), frozenset(("A", "G")))
        self.assertEqual(minimal_iupac_symbol(("A", "G")), "R")
        self.assertEqual(minimal_iupac_symbol(("A", "C", "G", "T")), "N")

    def test_normalizes_lowercase_and_rejects_invalid_contextually(self):
        self.assertEqual(normalize_iupac("acgtryn", context="record 's1'"), "ACGTRYN")
        with self.assertRaisesRegex(IupacError, "record 's1'.*position 3.*X"):
            normalize_iupac("ACX", context="record 's1'")

    def test_all_fifteen_symbols_have_support_and_complements(self):
        self.assertEqual(set(IUPAC_SUPPORT), set("ACGTRYSWKMBDHVN"))
        self.assertEqual(set(IUPAC_COMPLEMENT), set("ACGTRYSWKMBDHVN"))
        for symbol, support in IUPAC_SUPPORT.items():
            self.assertEqual(iupac_support(symbol.lower()), support)
            complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
            self.assertEqual(iupac_support(IUPAC_COMPLEMENT[symbol]),
                             frozenset(complement[base] for base in support))

    def test_minimal_symbol_rejects_empty_and_invalid_support(self):
        with self.assertRaises(IupacError):
            minimal_iupac_symbol(())
        with self.assertRaises(IupacError):
            minimal_iupac_symbol(("A", "X"))

    def test_mismatch_positions_use_target_subset_semantics(self):
        self.assertEqual(mismatch_positions("ARNT", "AGCT"), ())
        self.assertEqual(mismatch_positions("AANT", "ARCT"), (2,))
        self.assertEqual(mismatch_positions("ACGT", "ACGA"), (4,))

    def test_reverse_complement_and_total_degeneracy(self):
        self.assertEqual(reverse_complement_iupac("ARYKMBVDHN"), "NDHBVKMRYT")
        self.assertEqual(sequence_degeneracy("ARYN"), 1 * 2 * 2 * 4)

    def test_mismatch_positions_rejects_different_lengths(self):
        with self.assertRaisesRegex(IupacError, "lengths must match"):
            mismatch_positions("AC", "A")

    def test_comparisons_report_invalid_symbol_context(self):
        with self.assertRaisesRegex(IupacError, "target.*position 2.*X"):
            mismatch_positions("AA", "AX")
        with self.assertRaisesRegex(IupacError, "sequence.*position 2.*X"):
            sequence_degeneracy("AX")
