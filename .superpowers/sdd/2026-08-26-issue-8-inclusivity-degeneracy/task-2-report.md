# Task 2 report: reusable IUPAC primitives

## Changes

- Added `qpcr_pipeline/iupac.py` with all 15 canonical IUPAC symbols, support and complement tables, `IupacError`, normalization, minimal-symbol inversion, conservative mismatch matching, reverse complement, and sequence degeneracy.
- Added `tests/test_iupac.py` covering the brief examples, complete symbol coverage, invalid supports, length validation, conservative subset semantics, and contextual invalid-symbol errors.
- Existing untracked `.tmp` and `__pycache__` directories were preserved.

## TDD evidence

Initial RED:

```text
python -m unittest tests.test_iupac -v
ERROR: ModuleNotFoundError: No module named 'qpcr_pipeline.iupac'
```

After the initial minimal implementation, the literal normalization test exposed the brief's ordering conflict. Its regex requires context, position, then symbol, while the brief's sample implementation places symbol first. The error was adjusted to satisfy the literal test while retaining all context.

Second focused RED command:

```text
python -m unittest tests.test_iupac -v
FAIL: complete-complement test assertion had an incorrect test-side complement expression
FAIL: sequence_degeneracy('ARYN') was 16, while the brief's expression 2*2*2*4 evaluates to 32
```

The complement assertion was corrected. The degeneracy expectation was corrected to canonical support semantics (`1*2*2*4 == 16`), consistent with the prescribed production algorithm.

Focused GREEN:

```text
python -m unittest tests.test_iupac -v
Ran 8 tests in 0.001s
OK
```

Full suite GREEN:

```text
python -m unittest discover -s tests -q
Ran 253 tests in 29.126s
OK (skipped=1)
```

## Self-review

- All 15 symbols A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N are explicitly represented in both support and complement tables.
- Matching uses conservative target-subset semantics and validates both sequences with contextual 1-based positions.
- Normalization is uppercase-only and does not include the full input sequence in errors.
- Degeneracy multiplies canonical support cardinalities; reverse complements preserve ambiguity semantics.
- Implementation is dependency-free, deterministic, typed at public collection boundaries, and contains no unrelated refactoring.

## Concerns

The brief contains two literal inconsistencies: the normalization sample error ordering conflicts with its required regex, and the requested `ARYN` degeneracy expression evaluates to 32 although canonical IUPAC support cardinalities yield 16. This implementation follows the explicit test regex and the prescribed canonical degeneracy algorithm, respectively.

Commit: `27e3d82 feat: add conservative IUPAC operations`
