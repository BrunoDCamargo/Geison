# Issue #7 Primer3 Assays Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rank conserved reference regions and generate multiple auditable forward/probe/reverse qPCR assays through the official Primer3 executable.

**Architecture:** Extend configuration first, then build pure candidate selection in `primer_design.py`, isolate Boulder-IO and subprocess behavior in `primer3.py`, and compose both behind an atomic `design_primers` service. Integrate that service immediately after conservation without changing earlier set semantics or the issue-#6 HTML.

**Tech Stack:** Python 3.10+, frozen/slotted dataclasses, standard-library subprocess/JSON/TSV, `unittest`, Biopython/PyYAML already present, external `primer3_core` only when enabled.

**Spec:** `docs/superpowers/specs/2026-08-25-issue-7-primer3-assay-design.md`

## Global Constraints

- Follow strict RED, GREEN, REFACTOR TDD and record the expected failure before production edits.
- Add no Python dependency and never invoke Primer3 through a shell.
- Keep primer design disabled by default; enabled primer design requires enabled conservation.
- Use the complete issue-#6 majority consensus as every Primer3 template.
- Preserve 1-based inclusive reference coordinates in all public models and artifacts.
- Publish data atomically and publish `primer_design_report.json` last after invalidating any prior report.
- Do not modify `report.html`; assay visualization belongs to issues #10 and #17.
- Preserve unrelated files and all existing sequence-set, alignment, and conservation semantics.

---

### Task 1: Primer-design configuration contract

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Modify: `tests/test_config.py`

**Interfaces:**
- Consumes: existing `PipelineConfig`, `load_config`, and `validate_pipeline_config`.
- Produces: `OligoConstraints`, `PrimerDesignConfig`, `validate_primer_design_config`, and `PipelineConfig.primer_design` with the exact fields/defaults in the spec.

- [ ] **Step 1: Write failing default and YAML tests**

Add imports and tests that assert these hand-written values, including nested overrides:

```python
self.assertEqual(config.primer_design.max_candidate_regions, 10)
self.assertEqual(config.primer_design.primer.opt_tm, 60.0)
self.assertEqual(config.primer_design.probe.opt_tm, 70.0)

loaded = self._load_yaml(
    minimal_yaml
    + "primer_design:\n"
      "  enabled: true\n"
      "  max_candidate_regions: 3\n"
      "  assays_per_region: 4\n"
      "  candidate_region_length: 220\n"
      "  product_size_max: 180\n"
      "  primer:\n"
      "    opt_tm: 59.5\n"
      "  probe:\n"
      "    min_size: 20\n"
)
self.assertTrue(loaded.primer_design.enabled)
self.assertEqual(loaded.primer_design.max_candidate_regions, 3)
self.assertEqual(loaded.primer_design.primer.opt_tm, 59.5)
self.assertEqual(loaded.primer_design.probe.min_size, 20)
```

The fixture must enable alignment and conservation so only the intended configuration is under test.

- [ ] **Step 2: Run the focused tests and confirm RED**

Run: `python -m unittest tests.test_config.PipelineConfigTests.test_loads_primer_design_configuration_and_defaults_when_omitted -v`

Expected: import/attribute failure because the new configuration types do not exist.

- [ ] **Step 3: Implement immutable models, parsing, and defaults**

Use immutable default factories and merge nested YAML values over each preset:

```python
@dataclass(frozen=True, slots=True)
class OligoConstraints:
    min_size: int
    opt_size: int
    max_size: int
    min_tm: float
    opt_tm: float
    max_tm: float
    min_gc_percent: float
    max_gc_percent: float

def _primer_defaults() -> OligoConstraints:
    return OligoConstraints(18, 20, 25, 58.0, 60.0, 62.0, 40.0, 60.0)

def _probe_defaults() -> OligoConstraints:
    return OligoConstraints(18, 25, 30, 68.0, 70.0, 72.0, 30.0, 80.0)
```

Add `primer_design` to `PipelineConfig`, load it from YAML, and validate it before output creation.

- [ ] **Step 4: Add failing validation tables, then implement each rule**

Cover booleans masquerading as numbers, non-finite floats, bounds, unknown top-level/nested fields, unordered size/Tm triples, invalid product range, candidate region shorter than product maximum, and enabled design without conservation:

```python
invalid = PrimerDesignConfig(
    candidate_region_length=150,
    product_size_max=200,
)
with self.assertRaisesRegex(ValueError, "candidate_region_length"):
    validate_primer_design_config(invalid)
```

Run each new table first to observe its intended failure, implement the minimal shared validators, then rerun.

- [ ] **Step 5: Verify configuration regressions**

Run: `python -m unittest tests.test_config -v`

Expected: all configuration tests pass.

- [ ] **Step 6: Commit**

```powershell
git add -- qpcr_pipeline/config.py tests/test_config.py
git commit -m "feat: configure Primer3 assay design"
```

---

### Task 2: Pure candidate-region selection

**Files:**
- Create: `qpcr_pipeline/primer_design.py`
- Create: `tests/test_primer_design.py`

**Interfaces:**
- Consumes: `ConservationResult`, `PositionConservation`, `WindowConservation`, and validated `PrimerDesignConfig`.
- Produces: frozen/slotted `CandidateRegion`, `DesignedOligo`, `AssayCandidate`, `PrimerDesignResult`, `PrimerDesignError`, and `_select_candidate_regions(conservation, config) -> tuple[CandidateRegion, ...]`.

- [ ] **Step 1: Write a failing boundary-expansion test**

Build literal conservation positions/windows and assert intervals shift at both ends instead of shrinking:

```python
regions = _select_candidate_regions(conservation, permissive_config)
self.assertEqual(
    [(item.reference_start, item.reference_end) for item in regions],
    [(1, 300), (201, 500)],
)
```

Expected RED: module/function missing.

- [ ] **Step 2: Implement models, input validation, and interval expansion**

Define the exact dataclasses from the spec. Validate `COMPLETE` status, reference ID, contiguous reference positions, consensus length, and in-bounds windows before constructing candidate intervals.

- [ ] **Step 3: Add failing metric and eligibility tests**

Use hand-calculated positions to assert exact mean/minimum conservation, coverage, gap, entropy, usable length, and usable fraction. Parameterize one mutation per threshold so each threshold independently excludes the region.

```python
self.assertEqual(region.usable_length, 3)
self.assertAlmostEqual(region.mean_conservation, 0.925)
self.assertAlmostEqual(region.mean_entropy_bits, 0.25)
```

Observe RED for missing aggregation/eligibility, then implement with `math.fsum` and the configured thresholds.

- [ ] **Step 4: Add failing ranking, deduplication, overlap, and cap tests**

Construct literal tied windows that isolate each ordering key. Assert identical intervals retain the best peak, overlap is `overlap_length / min(lengths)`, equality at the configured maximum is allowed, excess overlap is rejected, and rank/IDs are reassigned only after selection.

```python
self.assertEqual([item.region_id for item in regions], ["region-001", "region-002"])
self.assertEqual([item.rank for item in regions], [1, 2])
```

Observe RED, implement deterministic sorting/greedy selection, and refactor helpers only while green.

- [ ] **Step 5: Add empty and malformed-input tests**

Assert an empty complete conservation result returns `()`. Assert non-complete status, coordinate gaps/duplicates, consensus mismatch, and out-of-range windows raise `PrimerDesignError` before any filesystem interaction.

- [ ] **Step 6: Verify candidate selection**

Run: `python -m unittest tests.test_primer_design -v`

Expected: all candidate/model tests pass.

- [ ] **Step 7: Commit**

```powershell
git add -- qpcr_pipeline/primer_design.py tests/test_primer_design.py
git commit -m "feat: rank conserved candidate regions"
```

---

### Task 3: Safe Primer3 Boulder-IO boundary

**Files:**
- Create: `qpcr_pipeline/primer3.py`
- Create: `tests/test_primer3.py`

**Interfaces:**
- Consumes: `CandidateRegion`, `DesignedOligo`, `AssayCandidate`, `PrimerDesignConfig`, and the full majority consensus.
- Produces: `Primer3Runner` protocol, `SubprocessPrimer3Runner`, `build_primer3_input(consensus, candidates, config) -> str`, and `parse_primer3_output(text, candidates, consensus) -> tuple[tuple[AssayCandidate, ...], dict[str, dict[str, str]]]` where the mapping preserves record errors/warnings/explanations.

- [ ] **Step 1: Write an exact failing Boulder writer test**

Assert a complete literal record, including:

```text
SEQUENCE_ID=region-001
SEQUENCE_TEMPLATE=<full consensus>
SEQUENCE_INCLUDED_REGION=0,300
PRIMER_TASK=generic
PRIMER_FIRST_BASE_INDEX=0
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_NUM_RETURN=5
PRIMER_PRODUCT_SIZE_RANGE=70-200
...
=
```

Include every configured `PRIMER_*` and `PRIMER_INTERNAL_*` size/Tm/GC tag in a fixed order. Observe module/function failure, then implement deterministic serialization. Reject IDs or sequences containing Boulder delimiters/newlines and noncanonical consensus symbols.

- [ ] **Step 2: Write failing parser tests for one complete assay**

Use a literal output record with left `(10,20)`, internal `(35,25)`, right `(89,20)`, and product size `80`. Independently assert 1-based intervals `11..30`, `36..60`, and `71..90`, synthesis-orientation sequences, Tm/GC/penalties, pair penalty, stable assay ID, and sorted extra metrics.

Observe RED, then implement record splitting, scalar parsing, and exact coordinate conversion.

- [ ] **Step 3: Extend parser tests across records and failures**

Write separate failing cases for out-of-order records, multiple assays, a valid zero-pair record with explanations, global error tags, duplicate/missing/unknown IDs, partial assays, malformed numeric values, invalid lengths, out-of-region coordinates, sequence-length mismatch, and product-coordinate disagreement. Implement one validation branch at a time.

- [ ] **Step 4: Write failing subprocess runner tests using real temporary executables**

Create executable `.cmd` fixtures on Windows and executable scripts on POSIX that read stdin and emit controlled stdout/stderr. Assert observable stdout and errors rather than mocked subprocess calls. Cover fixed arguments, stdin exchange, missing executable, nonzero exit with a 2,000-character bounded stderr excerpt, and empty output.

- [ ] **Step 5: Implement the non-shell runner**

Resolve the executable with `shutil.which` and call:

```python
subprocess.run(
    [executable, "--strict_tags", "--io_version=4"],
    input=input_text,
    capture_output=True,
    text=True,
    encoding="utf-8",
    errors="strict",
    check=False,
    shell=False,
)
```

Normalize unexpected runner exceptions into `PrimerDesignError` without including the full consensus.

- [ ] **Step 6: Verify the Primer3 boundary**

Run: `python -m unittest tests.test_primer3 -v`

Expected: all Boulder writer/parser/runner tests pass.

- [ ] **Step 7: Commit**

```powershell
git add -- qpcr_pipeline/primer3.py tests/test_primer3.py
git commit -m "feat: add safe Primer3 Boulder IO boundary"
```

---

### Task 4: Atomic primer-design service and artifacts

**Files:**
- Modify: `qpcr_pipeline/primer_design.py`
- Modify: `tests/test_primer_design.py`

**Interfaces:**
- Consumes: Task 2 selection/models and Task 3 writer/runner/parser.
- Produces: `design_primers(conservation, config, output_dir, *, runner=None) -> PrimerDesignResult`, stable issue-#7 TSV/raw/JSON artifacts, and skip/empty/failure semantics.

- [ ] **Step 1: Write a failing enabled-service happy-path test**

Use a deterministic fake runner that returns a complete literal Primer3 record. Assert the real service result and filesystem artifacts, including candidate/assay counts, paths, exact TSV headers, raw input/output equality, typed JSON values, preserved extra metrics, and relative POSIX artifact paths.

Observe RED because `design_primers` does not exist, then implement the minimal orchestration and serializers.

- [ ] **Step 2: Write failing disabled and empty-candidate tests**

Assert disabled design never reads malformed conservation data or calls the runner, removes only exact stale issue-#7 artifacts, preserves a sibling file, and publishes `SKIPPED`. Assert enabled/no-candidate design writes both header-only TSVs, omits raw files, avoids the runner, and publishes `COMPLETE` with zero counts.

Implement exact cleanup and the two early-return branches.

- [ ] **Step 3: Write failing zero-assay and runner-error tests**

Assert candidates plus a zero-pair response publishes raw exchanges and a header-only assay TSV. Assert runner/parsing errors invalidate a prior `COMPLETE` report and leave no new report.

Implement runner construction only when candidates exist and preserve explanation mappings in JSON.

- [ ] **Step 4: Write failing atomic-publication tests**

For each destination replacement, patch `Path.replace` to fail and assert the prior report is absent, temporary siblings are removed, unrelated siblings remain, and no hardlinked source is modified. Include disabled cleanup and enabled data/report ordering.

Implement unique temporary siblings with `uuid.uuid4`, `open("x")`, and `Path.replace`; publish report last.

- [ ] **Step 5: Verify the service**

Run: `python -m unittest tests.test_primer_design -v`

Expected: all candidate, service, artifact, skip, empty, and failure tests pass.

- [ ] **Step 6: Commit**

```powershell
git add -- qpcr_pipeline/primer_design.py tests/test_primer_design.py
git commit -m "feat: publish auditable Primer3 assays"
```

---

### Task 5: Pipeline, documentation, and real-binary integration

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `tests/test_minimal_run.py`
- Create: `integration_tests/test_primer3_core.py`
- Modify: `README.md`

**Interfaces:**
- Consumes: `PipelineConfig.primer_design`, `Primer3Runner`, and `design_primers`.
- Produces: `run_pipeline(..., primer3_runner=None)`, the exact `qc_report.json.primer_design` summary, user-facing configuration/artifact documentation, and conditional real `primer3_core` verification.

- [ ] **Step 1: Write a failing pipeline integration test**

Configure alignment, conservation, and primer design against a deterministic alignment/Primer3 runner pair. Assert the fake Primer3 runner receives the entire majority consensus in each `SEQUENCE_TEMPLATE`, multiple returned assays are persisted, and QC contains exactly:

```python
self.assertEqual(qc_report["primer_design"], {
    "status": "COMPLETE",
    "reference_id": "seq-1",
    "candidate_region_count": 1,
    "assay_count": 2,
})
```

Observe RED because `run_pipeline` has no `primer3_runner` or design stage, then add the injected parameter, service call after conservation, and summary fields.

- [ ] **Step 2: Add failing pipeline guard and disabled tests**

Assert invalid direct primer configuration fails before output creation, enabled design without conservation fails before output creation, and the default disabled stage publishes only a `SKIPPED` report without invoking a failing runner. Implement only the necessary integration/validation corrections.

- [ ] **Step 3: Write the conditional real-binary integration test**

Create `integration_tests/test_primer3_core.py` with `@unittest.skipUnless(shutil.which("primer3_core"), ...)`. Feed a fixed canonical synthetic consensus and one literal candidate through `SubprocessPrimer3Runner`, parse it, and assert record identity, all returned oligos lie inside the candidate, product-size coordinate consistency, and at least one complete forward/probe/reverse assay. Do not assert version-sensitive penalty values.

- [ ] **Step 4: Document configuration, artifacts, coordinates, and runtime requirement**

Extend `README.md` with the complete default YAML shape, candidate ordering and overlap rule, Primer3 `PATH` requirement, output files, 1-based inclusive coordinates, disabled/empty behavior, and the explicit boundary that inclusivity/specificity/final risk/UI are later stages.

- [ ] **Step 5: Run focused and full verification**

Run:

```powershell
python -m unittest tests.test_config tests.test_primer_design tests.test_primer3 tests.test_minimal_run -v
python -m unittest discover -s tests -v
python -m unittest discover -s integration_tests -p 'test_*.py' -v
python -m pytest -q
```

Expected: all unit tests pass; existing Windows symlink skips may remain; external-tool integration tests skip only when their executables are absent.

- [ ] **Step 6: Commit**

```powershell
git add -- qpcr_pipeline/pipeline.py tests/test_minimal_run.py integration_tests/test_primer3_core.py README.md
git commit -m "feat: integrate Primer3 assay design pipeline"
```

---

### Task 6: Whole-issue verification and acceptance audit

**Files:**
- Modify only files required by verified defects.

**Interfaces:**
- Consumes: complete issue-#7 branch and its progress ledger.
- Produces: a clean whole-branch review, reproducible verification evidence, and no unresolved Critical/Important findings.

- [ ] **Step 1: Run formatting and repository checks**

Run `git diff --check HEAD~5..HEAD` and inspect `git status --short`. Confirm only intentional issue-#7 files are tracked and pre-existing `.tmp`/`__pycache__` paths remain untouched.

- [ ] **Step 2: Audit every acceptance criterion against tests and artifacts**

Map region ranking, configurable cap, complete majority template, configurable hydrolysis-probe preset, multiple assays, persisted metrics, and real-binary fixture to named passing tests. Any missing criterion receives a new failing test before a fix.

- [ ] **Step 3: Run the final full suite from a clean process**

Run:

```powershell
python -m unittest discover -s tests -v
python -m unittest discover -s integration_tests -p 'test_*.py' -v
python -m pytest -q
```

Record exact pass/skip/fail counts and verify no new warnings or errors.

- [ ] **Step 4: Complete final code review and fix wave**

Package the full issue diff from the pre-issue base, obtain a Superpowers whole-branch review, dispatch one TDD fix wave for all real findings, and run one scoped re-review plus the affected/full suites.

- [ ] **Step 5: Commit only verified fixes**

```powershell
git add -- qpcr_pipeline/config.py qpcr_pipeline/primer_design.py qpcr_pipeline/primer3.py qpcr_pipeline/pipeline.py tests/test_config.py tests/test_primer_design.py tests/test_primer3.py tests/test_minimal_run.py integration_tests/test_primer3_core.py README.md
git commit -m "fix: resolve Primer3 assay review findings"
```

Skip the commit when review is clean and the worktree has no tracked changes.
