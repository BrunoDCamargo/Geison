# Issue #8 Inclusivity and IUPAC Degeneracy Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evaluate every Primer3 assay against the complete Evaluation Set and publish auditable mismatch, compatibility, variation, and bounded IUPAC-degeneracy results without replacing original oligos.

**Architecture:** Add validated pipeline configuration, reusable pure IUPAC operations, and one pure-Python inclusivity service. The service performs directed local hit search in both record orientations, chooses deterministic assay geometry, aggregates selected binding sites, proposes bounded IUPAC expansions, compares original and proposed compatibility, and atomically publishes normalized artifacts.

**Tech Stack:** Python 3.11+, standard library (`dataclasses`, `itertools`, `json`, `pathlib`, `typing`, `uuid`), existing Biopython-backed `LocalSequenceRecord`, PyYAML, `unittest`, and pytest as the aggregate runner. No new dependency or external executable.

**Spec:** `docs/superpowers/specs/2026-08-26-issue-8-inclusivity-degeneracy-design.md`

## Global Constraints

- Inclusivity evaluates every and only QC-approved member of the complete `EvaluationSet`, in its existing order.
- Accepted IUPAC symbols are `ACGTRYSWKMBDHVN`; comparison is conservative target-support subset of oligo-support.
- All reported mismatch positions are 1-based in oligo synthesis orientation.
- Search both supplied and reverse-complement record orientations inside expected coordinates plus `search_flank`.
- A valid assay geometry is `forward_end < probe_start <= probe_end < reverse_start` with the configured amplicon-size tolerance.
- Primer compatibility defaults to at most 2 mismatches and no mismatch in the final five 3-prime bases; probe compatibility defaults to at most 1 mismatch.
- Total sequence degeneracy defaults to at most 16 for primers and 4 for probes; primer 3-prime degeneracy is disabled by default.
- Never mutate `AssayCandidate` or silently substitute a proposed sequence.
- Proposal evaluation uses the same selected geometric binding site as the original assay.
- Raw Evaluation sequences must not appear in JSON or TSV artifacts.
- Data artifacts are written atomically before `inclusivity_report.json`; disabled runs remove only exact issue-#8 TSVs and publish `SKIPPED` last.
- Preserve all existing behavior when `inclusivity.enabled` is false.
- Follow red-green-refactor TDD and commit after every independently reviewed task.

## File map

- Modify `qpcr_pipeline/config.py`: define, parse, and validate `InclusivityConfig`; attach it to `PipelineConfig` and enforce the primer-design dependency.
- Create `qpcr_pipeline/iupac.py`: reusable IUPAC normalization, support-set comparison, reverse complement, mismatch positions, minimal symbol, and total degeneracy.
- Create `qpcr_pipeline/inclusivity.py`: public immutable models, local hit search, geometry selection, variation/proposal logic, result construction, TSV/JSON rendering, and atomic publication.
- Modify `qpcr_pipeline/pipeline.py`: call inclusivity after Primer3 and add the top-level QC summary.
- Create `tests/test_iupac.py`: pure IUPAC behavior.
- Create `tests/test_inclusivity.py`: search, compatibility, proposals, validation, artifacts, and publication failure behavior.
- Modify `tests/test_config.py`: configuration parsing and direct validation.
- Modify `tests/test_minimal_run.py`: pipeline dependency, Evaluation Set routing, disabled behavior, and QC integration.
- Modify `README.md`: configuration, matching semantics, artifacts, and correct issue responsibility boundaries.

---

### Task 1: Inclusivity configuration and dependency validation

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Modify: `tests/test_config.py`

**Interfaces:**
- Consumes: existing `PipelineConfig`, `PrimerDesignConfig`, `load_config()`, and `validate_pipeline_config()`.
- Produces: `InclusivityConfig`, `validate_inclusivity_config(config: InclusivityConfig) -> None`, `PipelineConfig.inclusivity`, and YAML section `inclusivity`.

- [ ] **Step 1: Write failing default and YAML parsing tests**

Add imports for `InclusivityConfig` and `validate_inclusivity_config`, then add:

```python
def test_loads_inclusivity_configuration_and_defaults_when_omitted(self):
    base = (
        "target:\n  name: target\n"
        f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
        "alignment:\n  enabled: true\n"
        "conservation:\n  enabled: true\n"
        "primer_design:\n  enabled: true\n"
    )
    self.assertEqual(self._load_yaml(base).inclusivity, InclusivityConfig())

    loaded = self._load_yaml(
        base
        + "inclusivity:\n"
        "  enabled: true\n"
        "  search_flank: 40\n"
        "  max_hits_per_oligo: 7\n"
        "  max_primer_mismatches: 1\n"
        "  max_probe_mismatches: 0\n"
        "  reject_primer_3_prime_mismatch: false\n"
        "  primer_3_prime_bases: 4\n"
        "  max_primer_degeneracy: 8\n"
        "  max_probe_degeneracy: 2\n"
        "  allow_primer_3_prime_degeneracy: true\n"
        "  max_amplicon_size_delta: 5\n"
    )
    self.assertEqual(
        loaded.inclusivity,
        InclusivityConfig(
            enabled=True,
            search_flank=40,
            max_hits_per_oligo=7,
            max_primer_mismatches=1,
            max_probe_mismatches=0,
            reject_primer_3_prime_mismatch=False,
            primer_3_prime_bases=4,
            max_primer_degeneracy=8,
            max_probe_degeneracy=2,
            allow_primer_3_prime_degeneracy=True,
            max_amplicon_size_delta=5,
        ),
    )
```

- [ ] **Step 2: Run the focused test and verify red**

Run: `python -m unittest tests.test_config.PipelineConfigTests.test_loads_inclusivity_configuration_and_defaults_when_omitted -v`

Expected: import or attribute failure for `InclusivityConfig`.

- [ ] **Step 3: Add the immutable model, parser, and pipeline field**

Implement the exact model:

```python
@dataclass(frozen=True, slots=True)
class InclusivityConfig:
    enabled: bool = False
    search_flank: int = 250
    max_hits_per_oligo: int = 20
    max_primer_mismatches: int = 2
    max_probe_mismatches: int = 1
    reject_primer_3_prime_mismatch: bool = True
    primer_3_prime_bases: int = 5
    max_primer_degeneracy: int = 16
    max_probe_degeneracy: int = 4
    allow_primer_3_prime_degeneracy: bool = False
    max_amplicon_size_delta: int = 20
```

Add `inclusivity: InclusivityConfig = field(default_factory=InclusivityConfig)` to
`PipelineConfig`, parse `raw.get("inclusivity", {})`, reject every unknown field,
and pass the parsed object into the `PipelineConfig` constructor.

- [ ] **Step 4: Add failing invalid-value and dependency tests**

Use table-driven cases covering a non-mapping section, unknown field, integer
boolean, negative zero-allowed fields, zero positive-only fields, fractional
integers, direct construction, and enabled inclusivity without enabled primer
design. Include this dependency assertion:

```python
with self.assertRaisesRegex(ValueError, "inclusivity.*requires enabled primer design"):
    _ = PipelineConfig(
        target_name="target",
        input_fasta=FIXTURE_FASTA,
        inclusivity=InclusivityConfig(enabled=True),
    ).selected_input
```

- [ ] **Step 5: Implement shared exact-type validation**

Implement:

```python
def validate_inclusivity_config(config: InclusivityConfig) -> None:
    if not isinstance(config, InclusivityConfig):
        raise ValueError("Inclusivity configuration must be an InclusivityConfig.")
    for name in (
        "enabled",
        "reject_primer_3_prime_mismatch",
        "allow_primer_3_prime_degeneracy",
    ):
        if not isinstance(getattr(config, name), bool):
            raise ValueError(f"Inclusivity {name} must be a boolean.")
    for name in (
        "search_flank",
        "max_primer_mismatches",
        "max_probe_mismatches",
        "max_amplicon_size_delta",
    ):
        value = getattr(config, name)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"Inclusivity {name} must be a non-negative integer.")
    for name in (
        "max_hits_per_oligo",
        "primer_3_prime_bases",
        "max_primer_degeneracy",
        "max_probe_degeneracy",
    ):
        value = getattr(config, name)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"Inclusivity {name} must be a positive integer.")
```

Call it from `validate_pipeline_config()` and enforce the enabled dependency after
both configs have been validated.

- [ ] **Step 6: Run configuration tests and commit**

Run: `python -m unittest tests.test_config -q`

Expected: all configuration tests pass.

Commit:

```bash
git add qpcr_pipeline/config.py tests/test_config.py
git commit -m "feat: configure inclusivity evaluation"
```

---

### Task 2: Reusable IUPAC primitives

**Files:**
- Create: `qpcr_pipeline/iupac.py`
- Create: `tests/test_iupac.py`

**Interfaces:**
- Consumes: plain strings and iterables of canonical bases.
- Produces: `IupacError`, `normalize_iupac()`, `iupac_support()`, `minimal_iupac_symbol()`, `reverse_complement_iupac()`, `mismatch_positions()`, and `sequence_degeneracy()`.

- [ ] **Step 1: Write failing mapping and normalization tests**

```python
class IupacTests(unittest.TestCase):
    def test_support_and_minimal_symbols_cover_all_nonempty_base_sets(self):
        self.assertEqual(iupac_support("R"), frozenset(("A", "G")))
        self.assertEqual(minimal_iupac_symbol(("A", "G")), "R")
        self.assertEqual(minimal_iupac_symbol(("A", "C", "G", "T")), "N")

    def test_normalizes_lowercase_and_rejects_invalid_contextually(self):
        self.assertEqual(normalize_iupac("acgtryn", context="record 's1'"), "ACGTRYN")
        with self.assertRaisesRegex(IupacError, "record 's1'.*position 3.*X"):
            normalize_iupac("ACX", context="record 's1'")
```

- [ ] **Step 2: Verify red**

Run: `python -m unittest tests.test_iupac -v`

Expected: import failure for `qpcr_pipeline.iupac`.

- [ ] **Step 3: Implement the canonical support and complement tables**

Define all 15 symbols explicitly and construct the inverse table from sorted
support sets. Implement `IupacError(ValueError)` with these behaviors:

```python
def normalize_iupac(sequence: str, *, context: str) -> str:
    normalized = sequence.upper()
    for position, symbol in enumerate(normalized, 1):
        if symbol not in IUPAC_SUPPORT:
            raise IupacError(
                f"Invalid IUPAC symbol {symbol!r} in {context} at position {position}."
            )
    return normalized

def iupac_support(symbol: str) -> frozenset[str]:
    normalized = symbol.upper()
    if len(normalized) != 1 or normalized not in IUPAC_SUPPORT:
        raise IupacError(f"Invalid IUPAC symbol {symbol!r}.")
    return IUPAC_SUPPORT[normalized]

def minimal_iupac_symbol(bases: Iterable[str]) -> str:
    support = frozenset(base.upper() for base in bases)
    try:
        return SUPPORT_TO_IUPAC[support]
    except KeyError as error:
        raise IupacError("IUPAC support must be a non-empty subset of A, C, G, T.") from error

def reverse_complement_iupac(sequence: str) -> str:
    normalized = normalize_iupac(sequence, context="sequence")
    return "".join(IUPAC_COMPLEMENT[symbol] for symbol in reversed(normalized))
```

Every error reports context, symbol, and 1-based position, never the full input
sequence.

- [ ] **Step 4: Write failing conservative-match and degeneracy tests**

```python
def test_mismatch_positions_use_target_subset_semantics(self):
    self.assertEqual(mismatch_positions("ARNT", "AGCT"), ())
    self.assertEqual(mismatch_positions("AANT", "ARCT"), (2,))
    self.assertEqual(mismatch_positions("ACGT", "ACGA"), (4,))

def test_reverse_complement_and_total_degeneracy(self):
    self.assertEqual(reverse_complement_iupac("ARYKMBVDHN"), "NDHBVKMRYT")
    self.assertEqual(sequence_degeneracy("ARYN"), 2 * 2 * 2 * 4)
```

- [ ] **Step 5: Implement comparison and degeneracy**

```python
def mismatch_positions(oligo: str, target: str) -> tuple[int, ...]:
    if len(oligo) != len(target):
        raise IupacError("Oligo and target lengths must match.")
    return tuple(
        index
        for index, (oligo_symbol, target_symbol) in enumerate(zip(oligo, target), 1)
        if not iupac_support(target_symbol) <= iupac_support(oligo_symbol)
    )

def sequence_degeneracy(sequence: str) -> int:
    result = 1
    for symbol in sequence:
        result *= len(iupac_support(symbol))
    return result
```

- [ ] **Step 6: Run tests and commit**

Run: `python -m unittest tests.test_iupac -q`

Expected: all IUPAC tests pass.

Commit:

```bash
git add qpcr_pipeline/iupac.py tests/test_iupac.py
git commit -m "feat: add conservative IUPAC operations"
```

---

### Task 3: Local oligo search and mismatch classification

**Files:**
- Create: `qpcr_pipeline/inclusivity.py`
- Create: `tests/test_inclusivity.py`

**Interfaces:**
- Consumes: `AssayCandidate`, `DesignedOligo`, `InclusivityConfig`, and IUPAC helpers.
- Produces: `InclusivityError`, role/orientation aliases, all public result dataclasses from the spec, private `_Hit`, `_oligo_for_role()`, and `_enumerate_hits()`.

- [ ] **Step 1: Write public-model and exact-hit tests with assay builders**

Create local test helpers that build `DesignedOligo` and `AssayCandidate` with
fixed coordinates. Then test a forward oligo at expected coordinates:

```python
def test_enumerates_exact_forward_hit_with_source_coordinates(self):
    hits = _enumerate_hits(
        assay_id="a1",
        sequence_id="s1",
        oriented_sequence="TTACGTACTT",
        orientation="FORWARD",
        role="FORWARD",
        oligo=self._oligo("ACGT", 3, 6),
        config=InclusivityConfig(search_flank=0, max_hits_per_oligo=3),
    )
    self.assertEqual(hits[0].public.source_start, 3)
    self.assertEqual(hits[0].public.source_end, 6)
    self.assertEqual(hits[0].public.mismatch_positions, ())
    self.assertTrue(hits[0].public.exact_match)
    self.assertTrue(hits[0].public.compatible)
```

- [ ] **Step 2: Verify red**

Run: `python -m unittest tests.test_inclusivity.InclusivitySearchTests.test_enumerates_exact_forward_hit_with_source_coordinates -v`

Expected: import failure for `qpcr_pipeline.inclusivity`.

- [ ] **Step 3: Add the exact public and internal models**

Implement `OligoMatch`, `ProposedOligoCompatibility`, `AssayInclusivity`,
`OligoVariation`, `DegeneracyProposal`, and `InclusivityResult` exactly as defined
in the spec. Add this private carrier because target sites must never be published:

```python
@dataclass(frozen=True, slots=True)
class _Hit:
    public: OligoMatch
    oriented_start: int
    oriented_end: int
    target_in_synthesis_orientation: str
```

Add role lookup:

```python
def _oligo_for_role(assay: AssayCandidate, role: OligoRole) -> DesignedOligo:
    return {
        "FORWARD": assay.forward_primer,
        "PROBE": assay.probe,
        "REVERSE": assay.reverse_primer,
    }[role]
```

- [ ] **Step 4: Implement bounded window enumeration minimally**

`_enumerate_hits()` must normalize the oriented record and oligo, clip the
expanded 1-based interval, reverse-complement reverse-primer target segments,
calculate 1-based mismatch positions, mark the configured 3-prime suffix, convert
reverse-oriented intervals back to source coordinates, and sort by the exact spec
key before assigning `hit_rank` and truncating.

Use these exact compatibility predicates:

```python
three_prime_mismatch = role != "PROBE" and any(
    position > len(oligo_sequence) - config.primer_3_prime_bases
    for position in positions
)
compatible = (
    len(positions) <= (
        config.max_probe_mismatches if role == "PROBE"
        else config.max_primer_mismatches
    )
    and not (
        role != "PROBE"
        and config.reject_primer_3_prime_mismatch
        and three_prime_mismatch
    )
)
```

- [ ] **Step 5: Add synthetic mismatch, orientation, and deterministic-cap tests**

Cover exact, internal, first-base 5-prime, every position in the final-five
3-prime suffix, probe mismatch, IUPAC ambiguity, reverse primer synthesis
orientation, reverse-complement source coordinates, clipped windows, short
records, repeated equally scoring hits, and `max_hits_per_oligo`. Assert the
default 3-prime hit is incompatible while an internal two-mismatch primer remains
compatible.

Use one explicit reverse-coordinate assertion:

```python
self.assertEqual((hit.public.source_start, hit.public.source_end), (7, 10))
self.assertEqual(hit.public.orientation, "REVERSE_COMPLEMENT")
```

for an oriented interval `[3, 6]` in a length-12 source record.

- [ ] **Step 6: Run search tests and commit**

Run: `python -m unittest tests.test_inclusivity.InclusivitySearchTests -q`

Expected: all local search tests pass.

Commit:

```bash
git add qpcr_pipeline/inclusivity.py tests/test_inclusivity.py
git commit -m "feat: find and classify oligo matches"
```

---

### Task 4: Complete-assay geometry across the Evaluation Set

**Files:**
- Modify: `qpcr_pipeline/inclusivity.py`
- Modify: `tests/test_inclusivity.py`

**Interfaces:**
- Consumes: `_enumerate_hits()`, `LocalSequenceRecord`, `EvaluationSet`, `AssayCandidate`, and `InclusivityConfig`.
- Produces: private `_SelectedBinding`, `_validate_enabled_inputs()`, `_select_binding()`, and `_evaluate_original()`.

- [ ] **Step 1: Write failing geometry and tie-break tests**

Use an assay with forward `[3, 6]`, probe `[9, 12]`, reverse `[15, 18]`, product
size 16, and a record containing matching sites. Assert strict ordering,
amplicon-size tolerance, chosen orientation, and complete compatibility. Add
negative cases where the probe overlaps a primer and product size is outside the
delta.

```python
selected = _select_binding(record, assay, InclusivityConfig(search_flank=2))
self.assertTrue(selected.geometry_found)
self.assertEqual(selected.orientation, "FORWARD")
self.assertEqual(selected.amplicon_size, 16)
self.assertTrue(selected.original_compatible)
```

- [ ] **Step 2: Verify red**

Run: `python -m unittest tests.test_inclusivity.InclusivityGeometryTests -v`

Expected: import failure for `_select_binding`.

- [ ] **Step 3: Implement selected-binding state and deterministic triplets**

```python
@dataclass(frozen=True, slots=True)
class _SelectedBinding:
    assay: AssayCandidate
    sequence_id: str
    orientation: Orientation | None
    forward: _Hit | None
    probe: _Hit | None
    reverse: _Hit | None
    geometry_found: bool
    source_amplicon_start: int | None
    source_amplicon_end: int | None
    amplicon_size: int | None
    original_compatible: bool
    retained_hits: tuple[_Hit, ...]
```

Enumerate triplets only within one orientation. Accept strict geometry and an
amplicon size in
`[assay.product_size - max_amplicon_size_delta, assay.product_size + max_amplicon_size_delta]`.
Sort with the spec key. Mark only the chosen triplet's public matches selected.
If no triplet exists, select one independent best hit per role across both
orientations, keep `orientation=None`, and force complete compatibility false.

- [ ] **Step 4: Write failing Evaluation Set validation and ordering tests**

Test blank/duplicate Evaluation IDs, missing/extra/reordered records, invalid
IUPAC with contextual errors, a reverse-complement record, empty Evaluation Set,
and two assays. Assert output work items are assay-major then Evaluation-order:

```python
self.assertEqual(
    [(item.assay.assay_id, item.sequence_id) for item in selected],
    [("a1", "s1"), ("a1", "s2"), ("a2", "s1"), ("a2", "s2")],
)
```

- [ ] **Step 5: Implement enabled-input validation and whole-set evaluation**

```python
def _evaluate_original(
    records: tuple[LocalSequenceRecord, ...],
    evaluation_set: EvaluationSet,
    assays: tuple[AssayCandidate, ...],
    config: InclusivityConfig,
) -> tuple[_SelectedBinding, ...]:
    validated = _validate_enabled_inputs(records, evaluation_set, assays)
    return tuple(
        _select_binding(record, assay, config)
        for assay in assays
        for record in validated
    )
```

Validation requires exact IDs and order, unique nonblank assay IDs, valid oligo
length/coordinates, strict issue-#7 forward/probe/reverse assay geometry, and
positive product sizes. Error text identifies assay, record, role, or position but
never contains complete sequences.

- [ ] **Step 6: Run geometry tests and commit**

Run: `python -m unittest tests.test_inclusivity.InclusivityGeometryTests -q`

Expected: all geometry and Evaluation Set tests pass.

Commit:

```bash
git add qpcr_pipeline/inclusivity.py tests/test_inclusivity.py
git commit -m "feat: evaluate complete assay geometry"
```

---

### Task 5: Variation aggregation and bounded proposals

**Files:**
- Modify: `qpcr_pipeline/inclusivity.py`
- Modify: `tests/test_inclusivity.py`

**Interfaces:**
- Consumes: `_SelectedBinding`, IUPAC helpers, `EvaluationSet`, and `InclusivityConfig`.
- Produces: `_variation_rows()`, `_proposal_for_role()`, `_proposals()`, and `_assay_results_with_proposals()`.

- [ ] **Step 1: Write failing positional-variation tests**

Build selected geometric sites for three Evaluation records where one internal
base differs, one 5-prime base differs, and one reverse-primer base differs.
Assert positions remain in synthesis orientation, affected IDs preserve
Evaluation order, and affected fraction uses the full Evaluation Set denominator.

```python
self.assertEqual(variation.oligo_position, 2)
self.assertEqual(variation.original_symbol, "A")
self.assertEqual(variation.observed_symbol, "R")
self.assertEqual(variation.affected_sequence_ids, ("s2",))
self.assertEqual(variation.affected_fraction, 1 / 3)
```

- [ ] **Step 2: Verify red**

Run: `python -m unittest tests.test_inclusivity.InclusivityDegeneracyTests.test_aggregates_variations_in_evaluation_order -v`

Expected: import failure for `_variation_rows`.

- [ ] **Step 3: Implement variation aggregation**

Aggregate only selected geometric sites. Emit a row when at least one target
support differs from original support. Calculate observed union and minimal IUPAC
symbol. Mark final `primer_3_prime_bases` only for forward/reverse roles. Return
rows ordered assay order, role order `FORWARD`, `PROBE`, `REVERSE`, then position.

- [ ] **Step 4: Write failing proposal-status and optimization tests**

Add separate tests for:

- internal `A/G` variation proposing `R` and increasing exact count;
- proposal preserving `AssayCandidate.forward_primer.sequence`;
- no geometric sites using status `UNCHANGED`, reason `NO_GEOMETRIC_SITES`;
- invariant sites using `UNCHANGED`, `NO_VARIATION`;
- useful variation with no strict gain using `UNCHANGED`, `NO_IMPROVEMENT`;
- final-five primer variation using `REJECTED`, `REJECTED_3_PRIME` by default;
- the same variation accepted when `allow_primer_3_prime_degeneracy=True`;
- primer and probe proposals rejected independently by degeneracy limits;
- an already-degenerate original whose total degeneracy is included in the cap;
- deterministic tie-break by exact count, degeneracy, changed count, sequence;
- zero binding-site fraction serialized as `None`.

Use this preservation assertion:

```python
original = assay.forward_primer.sequence
proposal = _proposal_for_role(assay, "FORWARD", selected, evaluation_set, config)
self.assertEqual(assay.forward_primer.sequence, original)
self.assertNotEqual(proposal.proposed_sequence, "")
```

- [ ] **Step 5: Implement deterministic bounded search**

For each varying position, offer only the original symbol and the minimal symbol
covering original plus observed support. Enumerate with a depth-first search in
ascending position order, prune when `sequence_degeneracy(candidate)` exceeds the
role limit, score exact coverage over all selected binding sites, and choose:

```python
key = (
    -exact_count,
    sequence_degeneracy(candidate),
    len(changed_positions),
    candidate,
)
```

Map outcomes exactly:

```python
("ACCEPTED", "ACCEPTED_IMPROVEMENT")
("UNCHANGED", "NO_VARIATION")
("UNCHANGED", "NO_GEOMETRIC_SITES")
("UNCHANGED", "NO_IMPROVEMENT")
("REJECTED", "REJECTED_3_PRIME")
("REJECTED", "REJECTED_LIMIT")
```

- [ ] **Step 6: Recompute proposed mismatch details at fixed selected sites**

Implement `ProposedOligoCompatibility` for each selected role using the accepted
proposal sequence or original sequence. Preserve geometry and coordinates.
`AssayInclusivity.proposed_compatible` is true only when geometry exists and all
three proposed role results are compatible. Missing geometry yields three `None`
proposed role results and false complete compatibility.

- [ ] **Step 7: Run proposal tests and commit**

Run: `python -m unittest tests.test_inclusivity.InclusivityDegeneracyTests -q`

Expected: all variation and proposal tests pass.

Commit:

```bash
git add qpcr_pipeline/inclusivity.py tests/test_inclusivity.py
git commit -m "feat: propose bounded IUPAC degeneracy"
```

---

### Task 6: Public service and atomic artifacts

**Files:**
- Modify: `qpcr_pipeline/inclusivity.py`
- Modify: `tests/test_inclusivity.py`

**Interfaces:**
- Consumes: all Task 3-5 evaluation functions and models.
- Produces: public `evaluate_inclusivity(records, evaluation_set, primer_design, config, output_dir) -> InclusivityResult` and five fixed artifacts under `inclusivity/`.

- [ ] **Step 1: Write failing disabled and complete service tests**

For disabled mode, pass deliberately invalid records and a `SKIPPED`
`PrimerDesignResult`; assert the service does not inspect them, publishes only
`inclusivity_report.json`, and returns empty tuples and `None` data paths. For
enabled mode, assert four data TSVs, schema-version-1 JSON, relative paths, exact
Evaluation IDs, assay-major ordering, nullable fractions, and no full Evaluation
sequence in any artifact.

```python
result = evaluate_inclusivity(
    records=(),
    evaluation_set=EvaluationSet(()),
    primer_design=self._primer_result(status="SKIPPED"),
    config=InclusivityConfig(enabled=False),
    output_dir=output_dir,
)
self.assertEqual(result.status, "SKIPPED")
self.assertEqual(
    {path.name for path in (output_dir / "inclusivity").iterdir()},
    {"inclusivity_report.json"},
)
```

- [ ] **Step 2: Verify red**

Run: `python -m unittest tests.test_inclusivity.InclusivityArtifactTests.test_disabled_publishes_only_skipped_report -v`

Expected: import failure for `evaluate_inclusivity`.

- [ ] **Step 3: Implement service orchestration and stable renderers**

Use this exact public signature:

```python
def evaluate_inclusivity(
    records: tuple[LocalSequenceRecord, ...],
    evaluation_set: EvaluationSet,
    primer_design: PrimerDesignResult,
    config: InclusivityConfig,
    output_dir: Path,
) -> InclusivityResult:
```

Disabled mode validates only config, removes exact stale TSV names, and publishes
the final report. Enabled mode requires `primer_design.status == "COMPLETE"`, then
evaluates, proposes, builds results, and renders:

```text
oligo_matches.tsv
assay_inclusivity.tsv
oligo_variations.tsv
degeneracy_proposals.tsv
inclusivity_report.json
```

Use these exact TSV columns and explicit row field ordering; do not rely on
`asdict().values()`:

```python
OLIGO_MATCH_COLUMNS = (
    "assay_id", "sequence_id", "role", "orientation", "hit_rank",
    "source_start", "source_end", "expected_start", "expected_end",
    "displacement", "mismatch_positions", "mismatch_count", "exact_match",
    "three_prime_mismatch", "probe_mismatch", "compatible", "selected",
)
ASSAY_COLUMNS = (
    "assay_id", "sequence_id", "orientation", "geometry_found",
    "source_amplicon_start", "source_amplicon_end", "amplicon_size",
    "forward_hit_rank", "probe_hit_rank", "reverse_hit_rank",
    "original_forward_compatible", "original_probe_compatible",
    "original_reverse_compatible", "original_compatible",
    "proposed_forward_compatible", "proposed_probe_compatible",
    "proposed_reverse_compatible", "proposed_compatible",
)
VARIATION_COLUMNS = (
    "assay_id", "role", "oligo_position", "original_symbol",
    "original_support", "observed_symbol", "observed_support",
    "affected_sequence_ids", "affected_sequence_count", "affected_fraction",
    "primer_3_prime_position",
)
PROPOSAL_COLUMNS = (
    "assay_id", "role", "original_sequence", "proposed_sequence", "status",
    "reason", "original_degeneracy", "proposed_degeneracy",
    "changed_positions", "binding_site_count", "original_exact_count",
    "original_exact_fraction", "proposed_exact_count", "proposed_exact_fraction",
)
```

Serialize tuples as comma-separated values in TSV and JSON arrays in the report.
Render absent scalar values as empty TSV fields. The JSON top level has exactly
`schema_version`, `status`, `enabled`, `configuration`,
`evaluation_sequence_ids`, `counts`, `oligo_matches`, `assay_results`,
`variations`, `proposals`, and `artifacts`. Use
`json.dumps(report, indent=2, sort_keys=True, allow_nan=False)`.

The `counts` mapping has exactly `evaluation_sequences`, `assays`,
`assay_evaluations`, `retained_oligo_hits`, `variations`, `proposals`,
`accepted_proposals`, `original_compatible`, and `proposed_compatible`.
The `artifacts` mapping has exactly `oligo_matches`, `assay_inclusivity`,
`oligo_variations`, and `degeneracy_proposals`, using relative POSIX paths or
`null` when skipped.

- [ ] **Step 4: Write failing empty, stale cleanup, and atomic-failure tests**

Cover enabled empty assays, enabled empty Evaluation Set, header-only TSVs,
disabled exact stale cleanup, preservation of an unrelated sibling, normalized
LF newlines, old-report invalidation before the first data replacement, report
written last, and simulated replacement failure for each of the five destinations.

Patch `Path.replace` and record destination names:

```python
self.assertEqual(attempted_destinations[-1], "inclusivity_report.json")
self.assertFalse((artifact_dir / "inclusivity_report.json").exists())
```

when a data write fails after invalidation.

- [ ] **Step 5: Implement exact cleanup and atomic JSON-last publication**

Use unique siblings and cleanup in `finally`:

```python
def _atomic_write_text(destination: Path, content: str) -> None:
    temporary = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.tmp"
    try:
        with temporary.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
```

Create the directory, unlink the prior report, publish four TSVs in the documented
order, then publish JSON. Disabled mode unlinks only the four exact TSV paths.

- [ ] **Step 6: Run the complete module tests and commit**

Run: `python -m unittest tests.test_inclusivity -q`

Expected: all inclusivity tests pass.

Commit:

```bash
git add qpcr_pipeline/inclusivity.py tests/test_inclusivity.py
git commit -m "feat: publish inclusivity artifacts"
```

---

### Task 7: Pipeline integration, documentation, and regression verification

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `tests/test_minimal_run.py`
- Modify: `README.md`

**Interfaces:**
- Consumes: `PipelineConfig.inclusivity`, approved records, exact `QCResult.evaluation_set`, `PrimerDesignResult`, and `evaluate_inclusivity()`.
- Produces: pipeline issue-#8 execution, top-level `qc_report.json.inclusivity`, README configuration/artifact contract, and corrected issue boundaries.

- [ ] **Step 1: Write failing pre-output dependency and disabled pipeline tests**

Add an invalid direct configuration test proving enabled inclusivity without
primer design fails before output creation. Extend the minimal disabled run to
assert only the skipped report exists and QC contains:

```python
{
    "status": "SKIPPED",
    "evaluation_sequence_count": 0,
    "assay_count": 0,
    "assay_evaluation_count": 0,
    "original_compatible_count": 0,
    "proposed_compatible_count": 0,
}
```

- [ ] **Step 2: Verify red**

Run: `python -m unittest tests.test_minimal_run.MinimalRunTests.test_disabled_primer_design_publishes_only_skipped_report_without_runner -v`

Expected: missing `inclusivity` QC key or artifact directory.

- [ ] **Step 3: Call the service and publish QC counts**

Immediately after `design_primers()`, call:

```python
inclusivity = evaluate_inclusivity(
    approved_records,
    result.evaluation_set,
    primer_design,
    config.inclusivity,
    output_dir,
)
```

Add the exact QC mapping, where `assay_count` is the number of upstream assays and
the last two counts sum `AssayInclusivity.original_compatible` and
`.proposed_compatible`.

- [ ] **Step 4: Write and pass an enabled end-to-end Evaluation Set test**

Extend the existing fake-MAFFT/fake-Primer3 pattern with three QC-approved target
records. Make the clustering runner retain only two Discovery representatives,
while the third Evaluation member carries an internal primer variation. Enable
inclusivity and assert:

```python
self.assertEqual(report["evaluation_sequence_ids"], ["s1", "s2", "s3"])
self.assertEqual(qc["inclusivity"]["evaluation_sequence_count"], 3)
self.assertEqual(qc["inclusivity"]["assay_evaluation_count"], 3)
self.assertIn("s3", variation["affected_sequence_ids"])
self.assertEqual(proposal["status"], "ACCEPTED")
self.assertNotEqual(proposal["original_sequence"], proposal["proposed_sequence"])
```

Run: `python -m unittest tests.test_minimal_run -q`

Expected: all pipeline tests pass.

- [ ] **Step 5: Document the effective configuration and artifacts**

Add the complete default YAML block from the spec, describe conservative IUPAC
matching, both orientations, 1-based synthesis mismatch positions, complete-assay
geometry, original/proposed side-by-side behavior, four TSVs plus JSON, and the
scientific limitation that proposals are computational candidates.

Replace README lines claiming issue #8 uses the Discovery Set and issue #9 uses
the Evaluation Set with this responsibility statement:

```text
A issue #8 avalia inclusividade contra todo o Evaluation Set. A issue #9 avalia
especificidade contra conjuntos off-target; propostas IUPAC nunca substituem
silenciosamente os oligos originais.
```

- [ ] **Step 6: Run focused and full regression suites**

Run in order:

```bash
python -m unittest tests.test_config tests.test_iupac tests.test_inclusivity tests.test_minimal_run -q
python -m unittest discover -s tests -q
python -m unittest discover -s integration_tests -q
python -m pytest -q
git diff --check
```

Expected: all offline tests pass; integration tests may skip only documented
missing external executables; `git diff --check` prints nothing.

- [ ] **Step 7: Commit integration and documentation**

```bash
git add qpcr_pipeline/pipeline.py tests/test_minimal_run.py README.md
git commit -m "feat: integrate Evaluation Set inclusivity"
```

---

## Final review and delivery

- Run `superpowers:requesting-code-review` against the issue #8 spec, this plan,
  and the complete issue #8 commit range.
- Resolve every confirmed finding with `superpowers:receiving-code-review` and a
  focused red-green regression test where behavior changes.
- Run `superpowers:verification-before-completion` with fresh root commands:
  `python -m unittest discover -s tests -q`,
  `python -m unittest discover -s integration_tests -q`,
  `python -m pytest -q`, and `git diff --check`.
- Use `superpowers:finishing-a-development-branch` to present or carry out the
  integration action authorized by the user.
