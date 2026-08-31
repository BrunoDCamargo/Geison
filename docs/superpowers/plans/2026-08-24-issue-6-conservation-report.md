# Interactive Conservation Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Calculate traceable position and window conservation metrics from the Discovery Set alignment and publish an offline interactive genome report with annotations, zoom, and hover.

**Architecture:** Add validated conservation configuration and a focused `qpcr_pipeline.conservation` scientific service that produces immutable metrics and atomic TSV/FASTA/JSON artifacts. Keep browser presentation isolated in `qpcr_pipeline.report_html`, which renders a deterministic Canvas report from compact window and annotation data; route the completed analysis through the pipeline after alignment.

**Tech Stack:** Python 3.10+, standard-library `math`/`json`/`html`/`pathlib`/`uuid`, Biopython `SeqFeature` location models, PyYAML, `unittest`, pytest, inline HTML/CSS/JavaScript Canvas, Codex in-app browser verification.

**Spec:** `docs/superpowers/specs/2026-08-24-issue-6-conservation-report-design.md`

## Global Constraints

- Conservation defaults to disabled; existing runs must not require a completed alignment or a plotting dependency.
- Enabled conservation requires enabled alignment and a `COMPLETE` `AlignmentResult`.
- `window_size` and `step_size` are non-boolean integers from 1 through 1,000,000, and `step_size <= window_size`.
- IUPAC ambiguity contributes fractional A/C/G/T observations; gaps affect depth, coverage, and gap frequency but never base frequencies or entropy.
- Position conservation is `major_allele_frequency`; coverage, gaps, and entropy remain separate metrics.
- Consensus FASTA omits reference-gap columns and uses IDs `geison-major-consensus` and `geison-iupac-consensus`.
- Window coordinates are 1-based reference coordinates; short references get one partial window and longer references get full stepped windows plus one anchored terminal window when required.
- GenBank feature coordinates convert from Biopython zero-based/end-exclusive locations to 1-based inclusive spans; source and external-reference features are omitted.
- `report.html` is deterministic and self-contained, with no network request, CDN, remote font/image, or external JavaScript dependency.
- Artifact publication invalidates the exact previous conservation report before fixed-name data mutation, publishes data/HTML first, and publishes JSON last.
- Standard tests remain offline; browser verification uses a deterministic generated fixture.
- No candidate-region selection, Primer3, assay design, inclusivity, specificity, ranking, checkpoints, or desktop UI is added.

---

### Task 1: Parse and validate conservation configuration

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Modify: `tests/test_config.py`

**Interfaces:**
- Produces: `ConservationConfig(enabled: bool = False, window_size: int = 100, step_size: int = 10)`
- Produces: `validate_conservation_config(config: ConservationConfig) -> None`
- Extends: `PipelineConfig.conservation: ConservationConfig`
- Consumes later: Task 2 conservation service and Task 4 pipeline integration

- [ ] **Step 1: Add RED defaults and explicit YAML tests**

Import `ConservationConfig`. Prove omitted configuration uses exact defaults and
explicit YAML parses exactly:

```python
config = self._load_yaml(
    "target:\n  name: target\n"
    f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
    "alignment:\n  enabled: true\n"
    "conservation:\n"
    "  enabled: true\n"
    "  window_size: 120\n"
    "  step_size: 12\n"
)
self.assertEqual(
    config.conservation,
    ConservationConfig(enabled=True, window_size=120, step_size=12),
)
self.assertEqual(
    self._load_yaml(minimal_yaml).conservation,
    ConservationConfig(),
)
```

- [ ] **Step 2: Add RED YAML and direct-construction validation tables**

Cover a non-mapping section, unknown key, integer `enabled`, boolean/fractional/
zero/1,000,001 sizes, and `step_size > window_size`. Direct configurations must
fail through `selected_input`, including the cross-stage dependency:

```python
invalid = (
    ConservationConfig(enabled=1),
    ConservationConfig(window_size=True),
    ConservationConfig(window_size=0),
    ConservationConfig(step_size=1_000_001),
    ConservationConfig(window_size=10, step_size=11),
)
for conservation in invalid:
    with self.subTest(conservation=conservation):
        config = PipelineConfig(
            target_name="target",
            input_fasta=FIXTURE_FASTA,
            conservation=conservation,
        )
        with self.assertRaises(ValueError):
            _ = config.selected_input

config = PipelineConfig(
    target_name="target",
    input_fasta=FIXTURE_FASTA,
    alignment=AlignmentConfig(enabled=False),
    conservation=ConservationConfig(enabled=True),
)
with self.assertRaisesRegex(ValueError, "requires enabled alignment"):
    _ = config.selected_input
```

- [ ] **Step 3: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_config.PipelineConfigTests -v
```

Expected: import failure because `ConservationConfig` does not exist.

Commit:

```powershell
git add -- tests/test_config.py
git commit -m "test: define conservation configuration contract"
```

- [ ] **Step 4: Implement parsing and shared validation**

Add the frozen/slotted dataclass and defaulted `PipelineConfig` field. Parse only
the three allowed keys from an optional mapping. Use this validator verbatim:

```python
def validate_conservation_config(config: ConservationConfig) -> None:
    if not isinstance(config, ConservationConfig):
        raise ValueError("Conservation configuration must be a ConservationConfig.")
    if not isinstance(config.enabled, bool):
        raise ValueError("Conservation enabled must be a boolean.")
    for field_name, value in (
        ("window_size", config.window_size),
        ("step_size", config.step_size),
    ):
        if (
            isinstance(value, bool)
            or not isinstance(value, int)
            or not 1 <= value <= 1_000_000
        ):
            raise ValueError(
                f"Conservation {field_name} must be an integer between 1 and 1000000."
            )
    if config.step_size > config.window_size:
        raise ValueError("Conservation step_size cannot exceed window_size.")
```

Call it from `validate_pipeline_config`, followed by:

```python
if config.conservation.enabled and not config.alignment.enabled:
    raise ValueError("Enabled conservation requires enabled alignment.")
```

- [ ] **Step 5: Verify GREEN and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_config -v
python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v
```

Expected: all pass with only the existing Windows symlink skip and no warnings.

Commit:

```powershell
git add -- qpcr_pipeline/config.py
git commit -m "feat: configure conservation analysis"
```

---

### Task 2: Build the scientific conservation service and atomic artifacts

**Files:**
- Create: `qpcr_pipeline/conservation.py`
- Create: `qpcr_pipeline/report_html.py`
- Create: `tests/test_conservation.py`

**Interfaces:**
- Consumes: `ConservationConfig`, `LocalSequenceRecord`, `AlignmentResult`, `AlignedSequence`, `AlignmentCoordinate`
- Produces: `ConservationError`, `PositionConservation`, `WindowConservation`, `ReferenceAnnotation`, `ConservationResult`
- Produces: `analyze_conservation(records, alignment, config, output_dir, *, target_name) -> ConservationResult`
- Produces for Task 3: `render_conservation_html(*, target_name, reference_id, sequence_count, config, windows, annotations) -> str`
- Consumes later: Task 4 pipeline integration

- [ ] **Step 1: Add RED known-frequency and entropy tests**

Create alignment helpers that construct real immutable alignment models. Use four
aligned sequences with columns `AAAA`, `AACC`, and `ACGT` and assert:

```python
self.assertEqual(
    [(p.frequency_a, p.frequency_c, p.frequency_g, p.frequency_t)
     for p in result.positions],
    [(1.0, 0.0, 0.0, 0.0),
     (0.5, 0.5, 0.0, 0.0),
     (0.25, 0.25, 0.25, 0.25)],
)
self.assertEqual([p.entropy_bits for p in result.positions], [0.0, 1.0, 2.0])
self.assertEqual(
    [p.major_allele_frequency for p in result.positions],
    [1.0, 0.5, 0.25],
)
self.assertEqual(result.major_consensus, "AAA")
self.assertEqual(result.iupac_consensus, "AMN")
```

The reference is `AAA`, so reference-aware ties choose A.

- [ ] **Step 2: Add RED fractional-IUPAC, gap, insertion, and consensus tests**

Use reference `A-C` plus aligned members `ATC` and `A-C`. Assert the insertion
column has null reference coordinate, depth 1, coverage `1/3`, gap frequency
`2/3`, T frequency 1, and does not appear in either consensus FASTA sequence.

Use pairs `A/R`, `C/B`, and `G/N` to assert exact fractional frequencies:

```python
expected = (
    (0.75, 0.0, 0.25, 0.0, "R"),
    (0.0, 2 / 3, 1 / 6, 1 / 6, "B"),
    (1 / 8, 1 / 8, 5 / 8, 1 / 8, "N"),
)
```

Add a table containing every IUPAC symbol and its expected canonical support set.
Prove gaps are excluded from entropy and an all-gap column raises
`ConservationError` before publication.

- [ ] **Step 3: Add RED sliding-window tests**

Use four aligned sequences:

```python
("ref", "AAAAAAAAAAAA")
("s2",  "AACAAAAAAAAA")
("s3",  "ACG-AAAAAAAA")
("s4",  "ACT-AAAAAAAA")
```

With `window_size=5`, `step_size=4`, assert exact windows `(1, 5)`,
`(5, 9)`, `(8, 12)`. The first window must equal:

```python
WindowConservation(
    reference_start=1,
    reference_end=5,
    position_count=5,
    mean_conservation=0.75,
    minimum_conservation=0.25,
    mean_coverage=0.9,
    mean_gap_frequency=0.1,
    mean_entropy_bits=0.6,
)
```

Add:

- zero positions produces zero windows;
- length 4 with window size 5 produces exactly `(1, 4)`;
- step 1 produces every fitting full window without duplicating the anchored tail.

- [ ] **Step 4: Add RED reference annotation tests**

Create `SeqFeature` values using `SimpleLocation` and `CompoundLocation`. Assert:

```python
self.assertEqual(
    result.annotations,
    (
        ReferenceAnnotation("gene", 2, 5, 1, "abc"),
        ReferenceAnnotation("CDS", 7, 8, -1, "protein X"),
        ReferenceAnnotation("CDS", 11, 12, -1, "protein X"),
    ),
)
```

Cover label precedence `gene`, `locus_tag`, `product`, then type; omit `source`,
external-reference, empty, and out-of-range parts and assert published/skipped
counts in JSON. FASTA metadata with no features produces an empty tuple.

- [ ] **Step 5: Add RED input, empty, disabled, and artifact tests**

Cover:

- enabled requires `AlignmentResult(status="COMPLETE")`;
- duplicate/missing/reordered alignment IDs and record membership mismatches;
- unequal alignment lengths, unknown symbols, invalid coordinates, missing or
  reversed reference records;
- enabled empty publishes header-only TSVs, empty consensus FASTA files, an
  empty-state HTML, and COMPLETE JSON;
- disabled never reads alignment sequences, removes only exact stale TSV/FASTA/
  HTML paths, preserves unrelated siblings, and publishes SKIPPED;
- exact TSV headers and blank null coordinate fields;
- schema version 1, metric definitions, effective window parameters, counts,
  annotation counts, consensus lengths, and relative artifact paths;
- original IDs do not leak into consensus FASTA headers;
- hardlinked final artifacts are replaced without mutating their sources;
- a failure at every enabled data/HTML replacement point and every disabled
  removal/report point starts with an old COMPLETE report and leaves no report;
- no new COMPLETE report exists after any validation or renderer failure.

- [ ] **Step 6: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_conservation -v
```

Expected: import failure because `qpcr_pipeline.conservation` does not exist.

Commit:

```powershell
git add -- tests/test_conservation.py
git commit -m "test: define conservation analysis contract"
```

- [ ] **Step 7: Implement models, validation, metrics, windows, and annotations**

Define the public models exactly as the spec. Use explicit tables:

```python
_IUPAC_BASES = {
    "A": frozenset("A"), "C": frozenset("C"),
    "G": frozenset("G"), "T": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"),
    "S": frozenset("CG"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"),
    "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
_BASES_TO_IUPAC = {bases: code for code, bases in _IUPAC_BASES.items()}
```

For each non-gap symbol, add `1 / len(support)` to each supported base, then
divide totals by depth. Compute entropy with `math.log2`. Determine major base by
maximum frequency, reference-base tie preference, then `ACGT` order.

Build windows only from non-null reference coordinates using the exact start rules
in the spec. Use `math.fsum` for deterministic means.

Extract local `SeqFeature` parts from the selected reference metadata. Convert and
clip coordinates, use qualifier precedence, skip unsupported parts while counting
them, and sort annotations deterministically.

- [ ] **Step 8: Implement atomic scientific artifacts and the basic renderer boundary**

Implement `analyze_conservation` with all preflight validation before mutation.
Create a focused `qpcr_pipeline.report_html` containing the final function
signature and a basic deterministic, self-contained HTML renderer sufficient for
artifact and empty-state tests. Task 3 replaces its plot body with the full Canvas
interaction.

Use exact artifact filenames from the spec. Invalidate the prior JSON report before
the first fixed-name data mutation. Write through unique sibling files and
`Path.replace`; publish the JSON last. Pass only compact window and annotation
models to the renderer.

Use these exact TSV headers:

```text
alignment_position\treference_position\treference_base\tdepth\tcoverage\tfrequency_a\tfrequency_c\tfrequency_g\tfrequency_t\tgap_frequency\tmajor_allele_frequency\tentropy_bits\tmajor_consensus\tiupac_consensus
```

```text
reference_start\treference_end\tposition_count\tmean_conservation\tminimum_conservation\tmean_coverage\tmean_gap_frequency\tmean_entropy_bits
```

- [ ] **Step 9: Verify GREEN and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_conservation -v
python -W default::ResourceWarning -m unittest tests.test_config tests.test_alignment tests.test_conservation -v
```

Expected: all pass offline with no warnings.

Commit:

```powershell
git add -- qpcr_pipeline/conservation.py qpcr_pipeline/report_html.py
git commit -m "feat: calculate traceable conservation metrics"
```

---

### Task 3: Build and browser-verify the interactive genome report

**Files:**
- Modify: `qpcr_pipeline/report_html.py`
- Create: `tests/test_report_html.py`
- Create: `integration_tests/generate_conservation_report.py`

**Interfaces:**
- Consumes: `ConservationConfig`, `WindowConservation`, `ReferenceAnnotation`
- Produces: `render_conservation_html(*, target_name, reference_id, sequence_count, config, windows, annotations) -> str`
- Preserves: Task 2 artifact/service interface and self-contained HTML contract

- [ ] **Step 1: Add RED self-contained and safe-embedding tests**

Render a target name containing `</script><img src=x onerror=alert(1)>`, an
annotation label containing markup, and two windows. Assert:

```python
self.assertNotIn("</script><img", html)
self.assertNotIn("https://", html)
self.assertNotIn("http://", html)
self.assertNotIn("src=\"//", html)
self.assertIn("<canvas", html)
self.assertIn('id="reset-zoom"', html)
self.assertIn('id="hover-details"', html)
self.assertIn('id="top-windows"', html)
self.assertIn("wheel", html)
self.assertIn("pointermove", html)
self.assertIn("pointerdown", html)
self.assertIn("textContent", html)
```

Parse the exact embedded marker
`<script id="geison-report-data" type="application/json">` and prove the
original strings round-trip as data while every literal `</script` sequence is
escaped in source.

- [ ] **Step 2: Add RED semantic, ranking, and empty-state tests**

Use 12 windows with deliberate ties. Assert the embedded top-ten order is mean
conservation descending, minimum descending, coverage descending, entropy
ascending, start ascending. Assert the report contains:

- target/reference/sequence/window summary;
- visible legend labels for conservation and coverage;
- zoom/pan/reset instructions;
- annotation lane label;
- hover field labels for interval, mean/minimum conservation, coverage, gaps,
  entropy, and overlapping annotations;
- click handlers on top-window rows;
- an accessible empty-state element when windows are empty;
- Canvas resizing for device pixel ratio and no one-node-per-window SVG.

- [ ] **Step 3: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_report_html -v
```

Expected: failures for missing Canvas controls and interaction code.

Commit:

```powershell
git add -- tests/test_report_html.py
git commit -m "test: define interactive conservation report"
```

- [ ] **Step 4: Implement the Canvas report**

Replace the basic renderer body with one deterministic template containing:

```javascript
canvas.addEventListener("wheel", onWheel, {passive: false});
canvas.addEventListener("pointerdown", onPointerDown);
canvas.addEventListener("pointermove", onPointerMove);
canvas.addEventListener("pointerup", onPointerUp);
resetButton.addEventListener("click", resetZoom);
```

Maintain `viewStart`/`viewEnd` in 1-based reference coordinates. Wheel zoom is
pointer-centered and clamped to the genome bounds; drag changes the interval while
preserving width. Binary-search the nearest visible window for hover. Draw
conservation, coverage, axes, top-window emphasis, and annotation spans on Canvas.
Populate all textual values with `textContent`.

Embed data as JSON arrays and escape `<`, `>`, `&`, U+2028, and U+2029 before
placing it in a script block. Sort top windows in Python so UI ranking is
deterministic and directly testable.

- [ ] **Step 5: Run unit GREEN**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_report_html tests.test_conservation -v
```

Expected: all pass with no warnings.

- [ ] **Step 6: Generate and inspect a deterministic browser fixture**

Create `integration_tests/generate_conservation_report.py` as a reusable, manually
invoked fixture generator. Its core setup is:

```python
import sys
from pathlib import Path

from Bio.SeqFeature import SeqFeature, SimpleLocation

from qpcr_pipeline.alignment import (
    AlignedSequence,
    AlignmentCoordinate,
    AlignmentResult,
)
from qpcr_pipeline.config import ConservationConfig
from qpcr_pipeline.conservation import analyze_conservation
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet

if len(sys.argv) != 2:
    raise SystemExit("usage: generate_conservation_report.py OUTPUT_DIRECTORY")
output_dir = Path(sys.argv[1]).resolve()

reference = "ACGT" * 100
variant_one = list(reference)
variant_two = list(reference)
for position in range(140, 161):
    variant_one[position] = "A" if reference[position] != "A" else "C"
for position in range(270, 291):
    variant_two[position] = "T" if reference[position] != "T" else "G"

records = (
    LocalSequenceRecord(
        "ref",
        reference,
        metadata={
            "features": (
                SeqFeature(SimpleLocation(20, 90, strand=1), type="gene",
                           qualifiers={"gene": ["alpha"]}),
                SeqFeature(SimpleLocation(120, 190, strand=1), type="CDS",
                           qualifiers={"product": ["beta protein"]}),
                SeqFeature(SimpleLocation(240, 320, strand=-1), type="gene",
                           qualifiers={"gene": ["gamma"]}),
                SeqFeature(SimpleLocation(340, 390, strand=1), type="regulatory",
                           qualifiers={"note": ["terminal region"]}),
            )
        },
    ),
    LocalSequenceRecord("variant-1", "".join(variant_one)),
    LocalSequenceRecord("variant-2", "".join(variant_two)),
)
discovery = DiscoverySet(tuple(record.sequence_id for record in records))
alignment = AlignmentResult(
    status="COMPLETE",
    discovery_set=discovery,
    reference_id="ref",
    reference_mode="explicit",
    sequences=tuple(
        AlignedSequence(record.sequence_id, record.sequence, "forward")
        for record in records
    ),
    coordinates=tuple(
        AlignmentCoordinate(position, position, base)
        for position, base in enumerate(reference, start=1)
    ),
    alignment_fasta_path=None,
    coordinate_map_path=None,
    report_path=output_dir / "alignment" / "alignment_report.json",
)
result = analyze_conservation(
    records,
    alignment,
    ConservationConfig(enabled=True, window_size=50, step_size=10),
    output_dir,
    target_name="Geison browser fixture",
)
print(result.html_report_path.resolve())
```

The script accepts exactly one output-directory argument and prints the absolute
`report.html` path. Generate it with:

```powershell
python integration_tests/generate_conservation_report.py .tmp/conservation-browser-report
```

Serve `.tmp/conservation-browser-report` on localhost and use the Codex in-app
browser skill to verify:

- the page loads without console errors or network requests;
- Canvas and top-ten table are visible;
- wheel zoom changes the displayed reference interval;
- reset returns the full interval;
- pointer hover populates the detail panel;
- clicking a top-window row narrows the interval;
- annotation labels appear when their spans overlap the hover window.

Capture one screenshot for visual inspection. The fixture is temporary and is not
committed.

- [ ] **Step 7: Commit production**

Run `git diff --check`, then:

```powershell
git add -- qpcr_pipeline/report_html.py integration_tests/generate_conservation_report.py
git commit -m "feat: render interactive genome conservation report"
```

---

### Task 4: Integrate conservation into pipeline outputs and documentation

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `tests/test_minimal_run.py`
- Modify: `README.md`

**Interfaces:**
- Consumes: `analyze_conservation` and `ConservationResult`
- Produces: conservation artifacts and root `report.html` through `run_pipeline`
- Extends: `qc_report.json.conservation` with status, reference ID, position count, and window count
- Preserves: `RunSummary`, Evaluation Set, Discovery Set, clustering, and alignment contracts

- [ ] **Step 1: Add RED enabled pipeline integration test**

Reuse fake CD-HIT and MAFFT runners so Discovery is `s1`, `s3` with reference
`s3`. Enable conservation with `window_size=3`, `step_size=2`. Assert:

```python
self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["s1", "s2", "s3"])
self.assertEqual(qc_report["discovery_set"]["sequence_ids"], ["s1", "s3"])
self.assertEqual(qc_report["conservation"], {
    "status": "COMPLETE",
    "reference_id": "s3",
    "position_count": 4,
    "window_count": 2,
})
self.assertEqual(summary.sequence_ids, ["s1", "s2", "s3"])
```

Assert exact conservation artifacts, `report.html`, consensus sequences, and that
position TSV uses only aligned Discovery members rather than rejected or
non-representative Evaluation members.

- [ ] **Step 2: Add RED disabled and early-validation pipeline tests**

Extend default local and NCBI tests to assert SKIPPED conservation report, no
scientific data/HTML artifacts, and top-level QC traceability. Add direct invalid
`ConservationConfig` and enabled-conservation/disabled-alignment cases proving
failure before output-directory creation.

Patch `analyze_conservation` only in the invalid-configuration cases and prove it
is never called. Do not mock it in valid pipeline tests.

- [ ] **Step 3: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_minimal_run -v
```

Expected: missing conservation routing/artifact/QC failures.

Commit:

```powershell
git add -- tests/test_minimal_run.py
git commit -m "test: define conservation pipeline integration"
```

- [ ] **Step 4: Implement pipeline routing**

Immediately after alignment:

```python
conservation = analyze_conservation(
    discovery_records,
    alignment,
    config.conservation,
    output_dir,
    target_name=config.target_name,
)
```

Add this exact top-level QC object:

```python
"conservation": {
    "status": conservation.status,
    "reference_id": conservation.reference_id,
    "position_count": len(conservation.positions),
    "window_count": len(conservation.windows),
},
```

Do not change summary or any existing set/artifact semantics.

- [ ] **Step 5: Document configuration, metrics, artifacts, and interaction**

Extend `README.md` with:

- default-disabled YAML configuration and alignment dependency;
- precise coverage, gap, fractional IUPAC, major-frequency, entropy, consensus,
  and window definitions;
- reference-coordinate and insertion-column behavior;
- annotation extraction from GenBank references;
- all conservation artifact paths;
- self-contained `report.html`, Canvas zoom/pan/reset/hover/top-window controls;
- the boundary that issue #6 visualizes peaks but issue #7 selects candidates.

- [ ] **Step 6: Verify GREEN and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_minimal_run tests.test_config tests.test_alignment tests.test_conservation tests.test_report_html -v
python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v
python -m pytest -q
```

Expected: all standard tests pass with only the existing Windows symlink skip and
no real CD-HIT, MAFFT, or network access.

Commit:

```powershell
git add -- qpcr_pipeline/pipeline.py README.md
git commit -m "feat: publish conservation analysis from pipeline"
```

---

### Task 5: Review, verify, and publish issue #6

**Files:**
- Review: all changes from issue #6 base `ab89ba1` through HEAD
- No production changes unless review demonstrates a defect

**Interfaces:**
- Produces: reviewed issue #6 implementation and fast-forwarded `origin/develop`

- [ ] **Step 1: Run fresh final verification**

Run:

```powershell
python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v
python -m pytest -q
python -W default::ResourceWarning -m unittest discover -s integration_tests -p "test_*.py" -v
git diff --check ab89ba1..HEAD
```

Expected: standard suites pass; existing conditional CD-HIT/MAFFT integrations
pass when installed and otherwise skip; no whitespace errors.

- [ ] **Step 2: Browser-verify the final generated report**

Regenerate the deterministic report fixture from Task 3 against final HEAD and
repeat page load, console/network, wheel zoom, reset, hover, top-window click, and
annotation checks. Record the screenshot path and observed interaction results in
the execution ledger.

- [ ] **Step 3: Review the complete issue diff**

Review against GitHub issue #6 and the approved design. Critical and Important
findings return through focused RED/GREEN fixes and scoped re-review. Minor
observations are recorded and triaged by the whole-issue review.

- [ ] **Step 4: Confirm remote fast-forward and publish**

Run:

```powershell
git fetch origin develop
git merge-base --is-ancestor origin/develop HEAD
git rev-list --left-right --count origin/develop...HEAD
git push origin develop
git ls-remote origin refs/heads/develop
```

Expected: zero commits behind, fast-forward push, and remote `develop` equal to
local HEAD. Standing authorization is direct `origin/develop` publication only.
Do not open a PR, modify `main`, or close/mutate GitHub issues.
