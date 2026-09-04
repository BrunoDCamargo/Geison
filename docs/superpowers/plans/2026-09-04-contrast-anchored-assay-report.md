# Contrast-Anchored Assay Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Constrain Primer3 designs to the contrastive evidence that selected each candidate region, then publish a self-contained researcher-facing HTML report and evidence bundle that clearly separate technical completion from scientific outcome.

**Architecture:** Preserve the existing stage order and scientific boundaries. `CandidateRegion.reference_start/reference_end` remains the broad Primer3 design interval, while `peak_start/peak_end` becomes an active contrast anchor only when `PrimerDesignResult.candidate_source == CONTRASTIVE_CONSERVATION`. The final researcher report is generated after run-state recording from already-published Geison artifacts, not by recalculating science in the notebook or HTML layer.

**Tech Stack:** Python 3, unittest/pytest-compatible tests, Primer3 Boulder-IO, existing Geison JSON/TSV artifacts, static offline HTML, `zipfile`, Google Colab notebook.

**Spec:** `docs/superpowers/specs/2026-09-04-contrast-anchored-assay-report-design.md`

## Global Constraints

- Keep `CONSERVATION_ONLY` design behavior unchanged.
- A `COMPLETE` contrastive result is authoritative; never silently fall back to conservation-only regions.
- `SEQUENCE_INCLUDED_REGION` remains the broad candidate region.
- `SEQUENCE_TARGET` represents the complete contrast anchor for contrastive candidates only.
- Final assay-level specificity remains independent and authoritative.
- `COMPLETED` is a technical run state, not a synonym for `IN SILICO PASS`.
- Primary report is `output/report.html`, English, single-file, static, offline, escaped safely, and printable.
- Report rendering must consume published artifacts and must not recalculate scientific metrics.
- Report-generation failure must preserve underlying scientific run state and artifacts.
- Evidence bundle packaging must not recalculate scientific results.
- Do not weaken specificity thresholds or synthetic challenge rules merely to force a passing demo.

---

### Task 1: Enforce the contrast anchor at the Primer3 boundary

**Files:**
- Modify: `qpcr_pipeline/primer3.py`
- Modify: `qpcr_pipeline/primer_design/__init__.py`
- Test: `tests/test_primer3.py`
- Test: `tests/test_primer_design.py`

**Interfaces:**
- Consumes: existing `CandidateRegion.peak_start`, `CandidateRegion.peak_end`, and `PrimerDesignResult.candidate_source`.
- Produces: `build_primer3_input(..., require_contrast_anchor: bool = False) -> str` and `parse_primer3_output(..., require_contrast_anchor: bool = False) -> tuple[tuple[AssayCandidate, ...], dict[str, dict[str, str]]]`.
- Contract: when `require_contrast_anchor=True`, every candidate emits `SEQUENCE_TARGET` and every accepted amplicon contains `peak_start..peak_end` completely.

- [ ] **Step 1: Add failing Boulder-IO serialization tests**

Add tests proving a contrastive candidate with region `501..800` and peak `601..700` emits the broad included region and zero-based target coordinates:

```python
def test_contrastive_input_serializes_peak_as_sequence_target(self) -> None:
    region = candidate(reference_start=501, reference_end=800)
    region = replace(region, peak_start=601, peak_end=700)

    text = build_primer3_input(
        "A" * 1200,
        (region,),
        PrimerDesignConfig(),
        require_contrast_anchor=True,
    )

    self.assertIn("SEQUENCE_INCLUDED_REGION=500,300\n", text)
    self.assertIn("SEQUENCE_TARGET=600,100\n", text)
```

Also add a conservation-only test asserting `SEQUENCE_TARGET` is absent when `require_contrast_anchor=False`.

- [ ] **Step 2: Run the focused serialization tests and confirm RED**

Run:

```bash
python -m pytest tests/test_primer3.py -k "sequence_target or serializes_complete_record" -q
```

Expected: new contrast-anchor test fails because `build_primer3_input` has no anchor mode / target tag yet; legacy serialization still passes.

- [ ] **Step 3: Implement minimal Primer3 target serialization**

Update `build_primer3_input` to accept a keyword-only flag:

```python
def build_primer3_input(
    consensus: str,
    candidates: tuple[CandidateRegion, ...],
    config: PrimerDesignConfig,
    *,
    require_contrast_anchor: bool = False,
) -> str:
```

For each candidate, keep:

```python
included_start = candidate.reference_start - 1
included_length = candidate.reference_end - candidate.reference_start + 1
```

and, only when required, insert after `SEQUENCE_INCLUDED_REGION`:

```python
anchor_start = candidate.peak_start - 1
anchor_length = candidate.peak_end - candidate.peak_start + 1
lines.append(f"SEQUENCE_TARGET={anchor_start},{anchor_length}")
```

Validate that the anchor lies wholly inside the candidate region before serialization; raise `PrimerDesignError` if not.

- [ ] **Step 4: Add the regression test for an amplicon that escapes the anchor**

Use the observed failure geometry: candidate `501..800`, peak `601..700`, assay product approximately `714..800`. The parser must reject it when the anchor contract is enabled:

```python
def test_rejects_contrastive_amplicon_outside_anchor(self) -> None:
    region = replace(
        candidate(reference_start=501, reference_end=800),
        peak_start=601,
        peak_end=700,
    )
    output = _single_assay_output(
        forward=(713, 20),
        probe=(743, 25),
        reverse=(799, 20),
        product_size=87,
    )

    with self.assertRaisesRegex(PrimerDesignError, "contrast anchor"):
        parse_primer3_output(
            output,
            (region,),
            "A" * 1200,
            require_contrast_anchor=True,
        )
```

Add a positive case where forward begins before `601`, reverse ends after `700`, and parsing succeeds.

- [ ] **Step 5: Run the parser regression and confirm RED**

Run:

```bash
python -m pytest tests/test_primer3.py -k "contrastive_amplicon or contrast_anchor" -q
```

Expected: new escaping-amplicon test fails because parser currently validates only broad region bounds and oligo geometry.

- [ ] **Step 6: Implement post-parse anchor validation**

Extend `parse_primer3_output` with keyword-only `require_contrast_anchor: bool = False`. After product coordinates are established, require:

```python
if require_contrast_anchor and not (
    forward.reference_start <= candidate.peak_start
    and candidate.peak_end <= reverse.reference_end
):
    raise PrimerDesignError(
        f"Primer3 output pair {index} does not contain the required contrast anchor "
        f"{candidate.peak_start}..{candidate.peak_end}."
    )
```

Do not require any one oligo itself to overlap the anchor; require the complete amplicon to contain it.

- [ ] **Step 7: Wire candidate source into the Primer3 boundary**

In `design_primers`, derive:

```python
require_contrast_anchor = candidate_source == "CONTRASTIVE_CONSERVATION"
```

and pass that flag to both `build_primer3_input` and `parse_primer3_output`. Preserve the existing no-fallback behavior in `_candidates_and_source`.

- [ ] **Step 8: Add primer-design source regression tests**

In `tests/test_primer_design.py`, verify a `COMPLETE` contrastive result passes anchor mode through a recording/fake runner and that a `COMPLETE` contrastive result with zero candidates remains `CONTRASTIVE_CONSERVATION` with zero assays and does not call legacy candidate selection.

- [ ] **Step 9: Run focused tests and commit**

Run:

```bash
python -m pytest tests/test_primer3.py tests/test_primer_design.py -q
```

Expected: PASS.

Commit:

```bash
git add qpcr_pipeline/primer3.py qpcr_pipeline/primer_design/__init__.py tests/test_primer3.py tests/test_primer_design.py
git commit -m "fix: anchor contrastive assays to discriminant windows"
```

---

### Task 2: Build an artifact-driven researcher report model

**Files:**
- Create: `qpcr_pipeline/researcher_report.py`
- Test: `tests/test_researcher_report.py`

**Interfaces:**
- Produces: `scientific_outcome(run_status: str, ranking_report: dict[str, object] | None) -> tuple[str, str]`.
- Produces: `load_researcher_report_data(output_dir: Path) -> ResearcherReportData`.
- `ResearcherReportData` contains only parsed published artifacts / provenance required for rendering; it does not recompute scientific metrics.

- [ ] **Step 1: Write failing outcome-state tests**

Cover the required states with minimal ranking-report fixtures:

```python
def test_completed_with_pass_has_positive_outcome(self):
    title, _ = scientific_outcome(
        "COMPLETED",
        {"assays": [{"classification": "IN SILICO PASS"}]},
    )
    self.assertEqual(title, "In-silico candidate(s) identified")


def test_completed_all_high_risk_has_negative_outcome(self):
    title, body = scientific_outcome(
        "COMPLETED",
        {"assays": [{"classification": "HIGH_RISK", "reason_codes": ["DETECTABLE_OFF_TARGET"]}]},
    )
    self.assertEqual(title, "No in-silico acceptable assay candidates identified")
    self.assertIn("HIGH_RISK", body)


def test_partial_is_inconclusive(self):
    title, _ = scientific_outcome("PARTIAL", None)
    self.assertEqual(title, "Inconclusive - insufficient evidence")
```

Also cover `REVIEW`, zero assays, and `FAILED`.

- [ ] **Step 2: Run outcome tests and confirm RED**

Run:

```bash
python -m pytest tests/test_researcher_report.py -k outcome -q
```

Expected: FAIL because module/functions do not exist.

- [ ] **Step 3: Implement the report-data loader and explicit outcome mapping**

Create a focused module with a frozen dataclass:

```python
@dataclass(frozen=True, slots=True)
class ResearcherReportData:
    output_dir: Path
    run_manifest: dict[str, object]
    panel: dict[str, object] | None
    conservation: dict[str, object] | None
    contrastive: dict[str, object] | None
    primer_design: dict[str, object] | None
    inclusivity: dict[str, object] | None
    specificity: dict[str, object] | None
    ranking: dict[str, object] | None
```

Load only known files when present:

```text
run_manifest.json
panel/approved_panel.json
conservation/conservation_report.json
contrastive_conservation/contrastive_conservation_report.json
primer_design/primer_design_report.json
inclusivity/inclusivity_report.json
specificity/specificity_report.json
ranking/ranking_report.json
```

Missing stage files must map to `None`, not fabricated empty evidence.

- [ ] **Step 4: Add loader tests for missing and present artifacts**

Use a temporary output directory with only a run manifest, then add one artifact at a time. Assert missing sections remain `None` and that malformed JSON raises a clear report-data error rather than being silently ignored.

- [ ] **Step 5: Run focused tests and commit**

Run:

```bash
python -m pytest tests/test_researcher_report.py -q
```

Expected: PASS.

Commit:

```bash
git add qpcr_pipeline/researcher_report.py tests/test_researcher_report.py
git commit -m "feat: model researcher report from published artifacts"
```

---

### Task 3: Render the complete self-contained study report

**Files:**
- Create: `qpcr_pipeline/researcher_report_html.py`
- Modify: `qpcr_pipeline/assay_report_html.py` only to reuse small safe formatting helpers if needed; do not duplicate ranking logic.
- Test: `tests/test_researcher_report_html.py`
- Keep existing: `tests/test_assay_report_html.py`

**Interfaces:**
- Consumes: `ResearcherReportData` and the explicit outcome from Task 2.
- Produces: `render_researcher_report_html(data: ResearcherReportData) -> str`.
- Renderer may transform existing numeric arrays into inline SVG coordinates for display, but must not recompute scientific classifications, similarity, conservation, specificity, or ranking scores.

- [ ] **Step 1: Write failing summary and safety tests**

Construct artifact dictionaries containing one PASS assay and assert the HTML includes:

```text
Geison Researcher Report
Run summary
Scientific outcome
Approved panel and study context
Target conservation
Target vs non-target contrast
Assay design
Target coverage / inclusivity
Specificity
Final candidates
Interpretation and limitations
Reproducibility
```

Also assert:

```python
self.assertNotIn("https://", html.lower())
self.assertNotIn("http://", html.lower())
self.assertNotIn("<script", html.lower())
```

and verify hostile target/dataset names are escaped.

- [ ] **Step 2: Run the HTML tests and confirm RED**

Run:

```bash
python -m pytest tests/test_researcher_report_html.py -q
```

Expected: FAIL because the renderer does not exist.

- [ ] **Step 3: Implement the report shell and summary hierarchy**

Render a single `<!doctype html>` document with embedded CSS only. Put the outcome ahead of aggregate scores, for example:

```html
<section class="outcome">
  <p class="eyebrow">Scientific outcome</p>
  <h2>...</h2>
  <p>...</p>
</section>
```

Summary cards must read counts already stored in report artifacts; do not derive new biological metrics from sequence payloads.

- [ ] **Step 4: Render panel, conservation, and contrast evidence**

Panel section reads the approved manifest definition. Conservation and contrast sections read their report JSON arrays/candidates. Render compact inline SVG plots using recorded coordinates and metric values only. Clearly draw candidate-region span separately from `peak_start..peak_end` / contrast anchor.

If a section is unavailable, render a visible `Evidence unavailable` block that names the missing artifact rather than inventing values.

- [ ] **Step 5: Render assay, inclusivity, specificity, and ranking evidence**

Reuse the existing final assay detail semantics: primer/probe sequences, coordinates, Tm, GC, penalty, product size, inclusivity counts, specificity evidence, classification, score components, and reason codes. For contrastive assays, show candidate source and anchor coordinates and state whether recorded assay coordinates contain the anchor.

Do not mark `REVIEW` or `HIGH_RISK` as recommended. Mark the top `IN SILICO PASS` only when one exists.

- [ ] **Step 6: Render limitations and reproducibility**

Always include the in-silico boundary. Render sanitized run IDs, timestamps, effective config summary, panel provenance/hash, reference provenance, and environment/tool versions from `run_manifest.json`. Never dump diagnostic logs or long raw sequence fields into the report.

- [ ] **Step 7: Add state-coverage tests**

Test at least PASS, REVIEW-only, HIGH_RISK-only, zero assays/PARTIAL, and FAILED. Verify the visible outcome differs correctly and missing evidence remains explicit.

- [ ] **Step 8: Run report tests plus legacy report tests and commit**

Run:

```bash
python -m pytest tests/test_researcher_report.py tests/test_researcher_report_html.py tests/test_assay_report_html.py -q
```

Expected: PASS.

Commit:

```bash
git add qpcr_pipeline/researcher_report.py qpcr_pipeline/researcher_report_html.py qpcr_pipeline/assay_report_html.py tests/test_researcher_report.py tests/test_researcher_report_html.py tests/test_assay_report_html.py
git commit -m "feat: publish complete researcher study report"
```

---

### Task 4: Publish the final report after run state is authoritative

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `qpcr_pipeline/ranking.py`
- Modify: `qpcr_pipeline/ranking_guard.py`
- Test: `tests/test_pipeline.py`
- Test: `tests/test_ranking.py`

**Interfaces:**
- Produces: `generate_researcher_report(output_dir: Path) -> Path` in `qpcr_pipeline/researcher_report.py`.
- Ranking continues to publish `ranking/assay_ranking.tsv` and `ranking/ranking_report.json`.
- `output/report.html` becomes the run-level researcher report generated after `RunRecorder.complete(...)` has written authoritative status/provenance.

- [ ] **Step 1: Write a failing pipeline ordering test**

Use a recorder/test double to verify researcher report generation happens only after `recorder.complete(...)`, so the renderer can read final `COMPLETED` / `PARTIAL` state and provenance from `run_manifest.json`.

- [ ] **Step 2: Run the focused test and confirm RED**

Run:

```bash
python -m pytest tests/test_pipeline.py -k researcher_report -q
```

Expected: FAIL because pipeline does not yet generate the run-level report after recorder completion.

- [ ] **Step 3: Stop ranking from owning the root report path**

Keep ranking's structured artifacts authoritative, but move any ranking-only HTML to an internal ranking path if retained, for example:

```python
"html": output_dir / "ranking" / "report.html"
```

Update its artifact JSON accordingly. The root `output/report.html` is reserved for the researcher report.

- [ ] **Step 4: Generate the researcher report after recorder completion**

After:

```python
recorder.complete(...)
```

call:

```python
try:
    generate_researcher_report(output_dir)
except Exception as report_error:
    _write_json_atomic(
        output_dir / "report_error.json",
        {"error_type": type(report_error).__name__, "message": str(report_error)[:1000]},
    )
```

Do not call `recorder.fail` for this presentation-only error, and do not alter assay classifications or the already-recorded scientific run state.

- [ ] **Step 5: Add tests for rendering failure isolation**

Patch `generate_researcher_report` to raise. Assert the scientific run manifest remains `COMPLETED` or `PARTIAL` as previously determined, scientific artifacts remain present, and `report_error.json` is written.

- [ ] **Step 6: Run pipeline/ranking tests and commit**

Run:

```bash
python -m pytest tests/test_pipeline.py tests/test_ranking.py tests/test_ranking_guard.py -q
```

Expected: PASS.

Commit:

```bash
git add qpcr_pipeline/pipeline.py qpcr_pipeline/ranking.py qpcr_pipeline/ranking_guard.py tests/test_pipeline.py tests/test_ranking.py tests/test_ranking_guard.py
git commit -m "feat: publish report after authoritative run completion"
```

---

### Task 5: Package the evidence bundle and expose report actions in Colab

**Files:**
- Create: `qpcr_pipeline/evidence_bundle.py`
- Modify: `notebooks/geison_guided_colab.ipynb`
- Modify: `docs/guided-colab.md`
- Test: `tests/test_evidence_bundle.py`
- Test: notebook static tests in the existing notebook-test file.

**Interfaces:**
- Produces: `create_evidence_bundle(output_dir: Path, destination: Path, *, extra_files: tuple[Path, ...] = ()) -> Path`.
- Bundle contains only an allowlist of published artifacts plus explicitly supplied approved config/panel files.

- [ ] **Step 1: Write failing allowlist packaging tests**

Create a temporary output tree with representative files and an unrelated secret-like file. Assert the ZIP contains the approved artifact set and excludes unrelated files:

```python
bundle = create_evidence_bundle(output_dir, tmp / "evidence.zip", extra_files=(approved_config,))
with zipfile.ZipFile(bundle) as archive:
    names = set(archive.namelist())
self.assertIn("report.html", names)
self.assertIn("run_manifest.json", names)
self.assertIn("ranking/ranking_report.json", names)
self.assertNotIn("unrelated-secret.txt", names)
```

- [ ] **Step 2: Run bundle tests and confirm RED**

Run:

```bash
python -m pytest tests/test_evidence_bundle.py -q
```

Expected: FAIL because helper does not exist.

- [ ] **Step 3: Implement deterministic evidence-bundle packaging**

Use `zipfile.ZipFile(..., ZIP_DEFLATED)` and an explicit allowlist covering:

```text
report.html
run_manifest.json
run_summary.json
qc_report.json
panel/approved_panel.json
contrastive_conservation/**
primer_design/**
inclusivity/**
specificity/**
ranking/**
.checkpoints/**/manifest.json
```

Include explicit `extra_files` by basename under `inputs/` (for example `config-approved.yaml`). Do not include arbitrary recursive files outside this list.

- [ ] **Step 4: Add Colab section `11. Researcher report`**

The final notebook cell must:

1. check `output/report.html` exists;
2. display a clear link/preview action;
3. build `evidence_bundle.zip` using `create_evidence_bundle` with the approved config/panel paths;
4. expose Colab downloads using `google.colab.files.download` for `report.html` and the ZIP.

Keep notebook code as presentation/packaging only; do not calculate scientific metrics there.

- [ ] **Step 5: Extend notebook static tests**

Assert the notebook contains:

```text
## 11. Researcher report
report.html
evidence_bundle.zip
create_evidence_bundle
files.download
```

and still contains all existing approval/run-manifest sections.

- [ ] **Step 6: Update guided-Colab documentation and commit**

Document the difference between the notebook and final report, the evidence ZIP, and browser-print-to-PDF behavior.

Run:

```bash
python -m pytest tests/test_evidence_bundle.py tests/test_guided_colab_notebook.py -q
```

Use the actual existing notebook static-test filename if it differs.

Commit:

```bash
git add qpcr_pipeline/evidence_bundle.py notebooks/geison_guided_colab.ipynb docs/guided-colab.md tests/test_evidence_bundle.py tests
git commit -m "feat: add researcher report downloads and evidence bundle"
```

---

### Task 6: Make the synthetic demo prove scientific coherence end to end

**Files:**
- Modify if needed: `examples/guided_demo/generate_demo_data.py`
- Modify: `integration_tests/test_guided_contrastive_demo.py`
- Modify if fixture expectations require it: `tests/test_guided_demo.py`

**Interfaces:**
- The synthetic fixture must retain a target-stable, challenge-shared region and a target-stable, challenge-altered discriminant region.
- At least one valid anchored assay should be possible without weakening specificity matching rules.

- [ ] **Step 1: Strengthen the integration test before touching fixture data**

Extend the E2E test to load `primer_design_report.json`, `ranking/ranking_report.json`, and `report.html`. Assert:

```python
self.assertEqual(primer_report["candidate_source"], "CONTRASTIVE_CONSERVATION")
self.assertTrue(primer_report["assays"])
```

For every accepted assay, map its region and assert the product contains `peak_start..peak_end`. Assert final specificity/ranking completes and the final classification is coherent with the synthetic challenge design. Require `output/report.html` and evidence-bundle creation.

- [ ] **Step 2: Run the guided integration test and observe the real failure**

Run:

```bash
python -m pytest integration_tests/test_guided_contrastive_demo.py -q
```

Expected before fixture adjustment: either Primer3 returns zero anchored pairs or the resulting anchored assays remain scientifically high-risk. Record which condition actually occurs; do not guess.

- [ ] **Step 3: Adjust only the synthetic fixture if necessary**

If the anchor constraint leaves no valid specific assay, alter only the synthetic challenge mutations / local sequence structure so the discriminant interval contains realistic Primer3-compatible target oligos while challenges disrupt one or more assay members sufficiently for the existing specificity engine to reject off-target amplification.

Preserve:

- deterministic seed `20260904`;
- target variants remaining inclusive;
- shared interval that remains highly challenge-similar;
- discriminant interval that remains target-conserved and challenge-different;
- no real pathogen assay sequence or network dependency.

Do not change specificity mismatch thresholds just to make the demo pass.

- [ ] **Step 4: Run demo unit tests plus real-tool E2E**

Run:

```bash
python -m pytest tests/test_guided_demo.py integration_tests/test_guided_contrastive_demo.py -q
```

Expected: PASS when MAFFT and Primer3 are installed; unit portion must pass everywhere.

- [ ] **Step 5: Commit the verified demo behavior**

```bash
git add examples/guided_demo/generate_demo_data.py integration_tests/test_guided_contrastive_demo.py tests/test_guided_demo.py
git commit -m "test: require scientifically coherent guided demo outcome"
```

---

### Task 7: Full regression, documentation, and acceptance verification

**Files:**
- Modify as needed: `README.md`
- Modify as needed: `docs/guided-colab.md`
- No unrelated refactors.

**Interfaces:**
- Final root artifacts: `output/report.html`, `output/run_manifest.json`, scientific stage artifacts, optional user-created `evidence_bundle.zip`.

- [ ] **Step 1: Run the complete unit/regression suite**

Run:

```bash
python -m pytest -q
```

Expected: all existing and new tests pass; only pre-existing intentional skips remain.

- [ ] **Step 2: Run the real-tool guided integration test explicitly**

Run:

```bash
python -m pytest integration_tests/test_guided_contrastive_demo.py -q
```

Expected: PASS when MAFFT and Primer3 are present.

- [ ] **Step 3: Inspect generated report and evidence bundle artifacts**

Generate the synthetic demo once and verify:

```text
output/report.html
output/run_manifest.json
output/ranking/ranking_report.json
evidence_bundle.zip
```

Open/search the HTML and confirm the visible scientific outcome matches ranking classifications, the contrast anchor is visible, and no external network resource is referenced.

- [ ] **Step 4: Update user-facing docs**

Document:

- contrast anchor semantics;
- technical completion vs scientific outcome;
- researcher report location;
- evidence bundle download;
- in-silico/wet-lab boundary.

- [ ] **Step 5: Final commit**

```bash
git add README.md docs/guided-colab.md
git commit -m "docs: explain anchored design and researcher report"
```

- [ ] **Step 6: Record verification evidence before integration**

Capture the exact final test counts and Git commit SHA. Do not claim completion until the full suite and guided integration evidence are available.
