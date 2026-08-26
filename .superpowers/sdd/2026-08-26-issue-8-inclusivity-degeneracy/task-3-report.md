# Task 3 report: local oligo search and mismatch classification

## Changes

- Added `qpcr_pipeline/inclusivity.py` with the exact immutable public result models, `InclusivityError`, role/orientation aliases, private `_Hit`, role lookup, and bounded `_enumerate_hits()`.
- Local matching normalizes IUPAC input, applies conservative target-support subset semantics, retains a per-orientation cap independently of compatibility, and retains raw target sites only in private `_Hit`.
- Search results use synthesis-orientation mismatch positions, reverse-primer target reverse complementation, configured 3-prime/probe compatibility flags, deterministic rank ordering, and reverse-complement coordinate conversion.
- Added `tests/test_inclusivity.py` with local builders and coverage for exact, mismatch, IUPAC, orientation, clipping, short-record, ordering, cap, and invalid-input behavior.
- Preserved the unrelated modified Task 2 report and untracked `.tmp` / `__pycache__` paths.

## TDD evidence

Initial RED command:

```text
python -m unittest tests.test_inclusivity.InclusivitySearchTests.test_enumerates_exact_forward_hit_with_source_coordinates -v
ERROR: ModuleNotFoundError: No module named 'qpcr_pipeline.inclusivity'
```

This is the expected failure before the new module and enumerator existed.

Focused GREEN command:

```text
python -m unittest tests.test_inclusivity.InclusivitySearchTests -q
Ran 15 tests in 0.001s
OK
```

Full-suite GREEN command:

```text
python -m unittest discover -s tests
Ran 268 tests in 26.989s
OK (skipped=1)
```

## Self-review

- Forward coordinates remain 1-based inclusive; reverse-complement intervals use `[L-end+1, L-start+1]`.
- Reverse primers compare the reverse-complemented segment, preserving 1-based synthesis-orientation positions.
- Ranking is mismatch count, 3-prime mismatch count, absolute displacement, oriented start, then raw normalized target segment; ranks are assigned after sorting and before the cap is returned.
- The final-five predicate uses the required strict `position > length - bases` boundary, and compatibility exactly follows the prescribed primer/probe predicates.
- Invalid IUPAC values are surfaced as `InclusivityError` with assay, record, role, and helper position context, never the full sequence.
- Public dataclasses contain no raw target site; it is available solely in `_Hit` for later tasks. No geometry, Evaluation Set orchestration, proposal, or publication behavior was added.

## Concerns

None.

Implementation commit: `25e7d19 feat: find and classify oligo matches`
