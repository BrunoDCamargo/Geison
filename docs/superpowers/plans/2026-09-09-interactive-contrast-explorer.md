# Interactive Contrast Explorer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the existing contrastive-conservation HTML report interactive for zooming, panning, filtering challenge datasets, and inspecting candidate regions without changing scientific results.

**Architecture:** Extend the existing self-contained `render_contrastive_html()` output rather than creating a second report. Python continues to serialize the same scientific evidence plus compact candidate metadata; browser-native Canvas and JavaScript handle view state, redraws, hit-testing, dataset visibility, and synchronized candidate selection. No external dependencies or write-back paths are added.

**Tech Stack:** Python stdlib (`html`, `json`), existing Geison dataclasses, HTML5 Canvas, browser-native JavaScript/CSS, pytest.

**Spec:** `docs/superpowers/specs/2026-09-09-interactive-contrast-explorer-design.md`

## Global Constraints

- Keep the report one deterministic offline HTML file.
- Add no network requests, CDN assets, external JavaScript libraries, fonts, or dependencies.
- Keep the explorer read-only: no edit handles, write-back, Primer3 rerun, or mutation of evidence.
- Default rendering must show target windows against their worst challenge similarity rather than plotting every per-dataset evidence row.
- Individual challenge datasets are opt-in overlays controlled in the report.
- Existing complete evidence tables remain available.
- Preserve current escaping guarantees for HTML and embedded JSON.

---

### Task 1: Lock the interaction contract with focused failing tests

**Files:**
- Modify: `tests/test_contrastive_report_html.py`

**Interfaces:**
- Consumes: `render_contrastive_html(...) -> str` from `qpcr_pipeline.contrastive_report_html`.
- Produces: regression expectations for controls, candidate metadata, offline behavior, read-only behavior, and determinism.

- [ ] **Step 1: Extend the existing report fixture so the candidate has metrics useful to the detail card**

Keep `_candidate()` as the single candidate fixture and reuse the existing synthetic window/dataset evidence. Do not add another fixture module.

- [ ] **Step 2: Add a failing test for interactive controls and candidate serialization**

Add this focused test:

```python
def test_report_exposes_read_only_interactive_explorer_controls():
    windows = (
        ContrastWindowEvidence(
            150, 200, 0.99, 0.95, 1.0, 0.0, 0.02, True,
            "challenge-a", "CRITICAL", 0.12, 0.12, None, 0.87,
        ),
    )
    dataset = (
        DatasetWindowEvidence(
            150, 200, "challenge-a", "CRITICAL", 2,
            "seq-1", "forward", 0.12,
        ),
    )
    rendered = render_contrastive_html(
        target_name="Synthetic target",
        reference_id="ref-1",
        windows=windows,
        dataset_evidence=dataset,
        candidates=(_candidate(),),
    )

    for marker in (
        'id="zoom-in"',
        'id="zoom-out"',
        'id="reset-view"',
        'id="show-all"',
        'id="hide-all"',
        'id="dataset-filters"',
        'id="candidate-detail"',
        "const candidateRegions=",
        "const datasetNames=",
    ):
        assert marker in rendered

    assert 'type="range"' not in rendered
    assert 'contenteditable="true"' not in rendered
    assert "http://" not in rendered
    assert "https://" not in rendered
```

- [ ] **Step 3: Add a failing test that the default chart does not automatically render every challenge row**

The HTML may still embed `challengePoints` for opt-in overlays, but its JavaScript should initialize the visible dataset set empty and draw the overview from `targetPoints` first:

```python
def test_report_defaults_to_worst_challenge_overview_and_opt_in_dataset_overlays():
    rendered = render_contrastive_html(
        target_name="target",
        reference_id=None,
        windows=(),
        dataset_evidence=(),
        candidates=(),
    )
    assert "const visibleDatasets=new Set();" in rendered
    assert "drawOverview" in rendered
    assert "visibleDatasets.has" in rendered
```

- [ ] **Step 4: Run the focused tests and verify RED**

Run:

```bash
pytest tests/test_contrastive_report_html.py -q
```

Expected: FAIL because the controls, `candidateRegions`, `datasetNames`, and view-state JavaScript do not exist yet.

- [ ] **Step 5: Commit the RED test**

```bash
git add tests/test_contrastive_report_html.py
git commit -m "test: define interactive contrast explorer contract"
```

---

### Task 2: Serialize compact explorer data and render native controls

**Files:**
- Modify: `qpcr_pipeline/contrastive_report_html.py`
- Test: `tests/test_contrastive_report_html.py`

**Interfaces:**
- Consumes: existing `windows`, `dataset_evidence`, and `candidates` arguments.
- Produces: embedded JavaScript constants `candidateRegions` and `datasetNames`, plus stable HTML control IDs used by Task 3.

- [ ] **Step 1: Build candidate-region JSON from existing dataclasses**

Inside `render_contrastive_html`, create a compact list before rendering:

```python
candidate_regions = [
    {
        "id": item.region.region_id,
        "rank": item.region.rank,
        "start": item.region.reference_start,
        "end": item.region.reference_end,
        "peakStart": item.region.peak_start,
        "peakEnd": item.region.peak_end,
        "contributing": item.contributing_windows,
        "worstDataset": item.worst_dataset_name,
        "criticality": item.worst_dataset_criticality,
        "similarity": item.worst_similarity,
        "margin": item.contrast_margin,
        "meanConservation": item.region.mean_conservation,
        "minimumConservation": item.region.minimum_conservation,
    }
    for item in candidates
]
dataset_names = sorted({row.dataset_name for row in dataset_evidence}, key=str.casefold)
```

Use `_safe_json()` for both arrays. Do not add a new serializer helper.

- [ ] **Step 2: Make candidate rows selectable without introducing edit semantics**

Give each candidate table row a stable read-only selector:

```python
f'<tr class="candidate-row" data-region-id="{_cell(item.region.region_id)}" tabindex="0">'
```

Keep the existing cells and evidence columns.

- [ ] **Step 3: Add native controls and a detail card to the HTML**

Immediately before the canvases, add:

```html
<div class="toolbar" aria-label="Explorer controls">
  <button id="zoom-in" type="button">Zoom in</button>
  <button id="zoom-out" type="button">Zoom out</button>
  <button id="reset-view" type="button">Reset view</button>
  <button id="show-all" type="button">Show all datasets</button>
  <button id="hide-all" type="button">Hide all datasets</button>
</div>
<div id="dataset-filters" class="filters" aria-label="Challenge dataset overlays"></div>
```

After the canvases, add:

```html
<div id="candidate-detail" class="card" aria-live="polite">
  Select a candidate region to inspect its recorded evidence.
</div>
```

Use CSS already present in the file. Add only the small rules needed for `.toolbar`, `.filters`, `.candidate-row`, `.selected`, and canvas cursor states.

- [ ] **Step 4: Embed the new JSON constants**

At the top of the existing `<script>` block add:

```javascript
const candidateRegions=...;
const datasetNames=...;
const visibleDatasets=new Set();
```

Keep `targetPoints` and `challengePoints` unchanged as the scientific source arrays.

- [ ] **Step 5: Run focused tests**

Run:

```bash
pytest tests/test_contrastive_report_html.py -q
```

Expected: the controls/serialization test passes; the behavior test may still fail until redraw/filter logic is implemented.

- [ ] **Step 6: Commit the data/control layer**

```bash
git add qpcr_pipeline/contrastive_report_html.py tests/test_contrastive_report_html.py
git commit -m "feat: add contrast explorer controls"
```

---

### Task 3: Add zoom, pan, filtering, hover, and synchronized candidate selection

**Files:**
- Modify: `qpcr_pipeline/contrastive_report_html.py`
- Test: `tests/test_contrastive_report_html.py`

**Interfaces:**
- Consumes: `targetPoints`, `challengePoints`, `candidateRegions`, `datasetNames`, and stable control IDs from Task 2.
- Produces: browser-local explorer behavior only; no Python or pipeline API changes.

- [ ] **Step 1: Replace one-shot drawing with two redraw functions**

Keep Canvas and browser APIs. Introduce browser-local view state:

```javascript
const fullStart=Math.min(0,...targetPoints.map(p=>p.start));
const fullEnd=Math.max(1,...targetPoints.map(p=>p.end));
let viewStart=fullStart;
let viewEnd=fullEnd;
let selectedRegionId=null;
```

Add coordinate helpers that map reference positions using `viewStart`/`viewEnd` and skip points outside the visible interval.

`drawOverview()` must:
- clear and redraw axes;
- draw one point per `targetPoints` item using `x = worst similarity` and `y = target conservation`;
- then draw only `challengePoints` whose `dataset` exists in `visibleDatasets`.

`drawTrack()` must:
- clear and redraw axes;
- draw target conservation as an ordered line or compact points by reference position;
- draw worst challenge similarity from `targetPoints` as the comparison series;
- draw selected challenge overlays only when enabled;
- draw each candidate region as a reference-coordinate band, with the selected band visually distinct.

- [ ] **Step 2: Build dataset checkboxes from `datasetNames`**

For each dataset, append a label with a checkbox. On change, add/remove the dataset name from `visibleDatasets` and call both redraw functions.

`Show all` fills the set with every `datasetNames` entry and checks all boxes. `Hide all` clears the set and unchecks all boxes.

- [ ] **Step 3: Implement zoom and reset with one shared function**

Use a single helper:

```javascript
function zoomAt(referencePosition,factor) {
  const width=viewEnd-viewStart;
  const next=Math.max(100,Math.min(fullEnd-fullStart,width*factor));
  const ratio=(referencePosition-viewStart)/width;
  viewStart=Math.max(fullStart,referencePosition-next*ratio);
  viewEnd=Math.min(fullEnd,viewStart+next);
  viewStart=Math.max(fullStart,viewEnd-next);
  drawTrack();
}
```

- `zoom-in`: factor `0.7` centered on the current view midpoint.
- `zoom-out`: factor `1 / 0.7` centered on the current view midpoint.
- mouse wheel on the reference track: same helper centered on the mouse reference coordinate; call `preventDefault()`.
- `reset-view`: restore `fullStart`/`fullEnd` and redraw.

Do not add a slider or third-party zoom package.

- [ ] **Step 4: Implement drag-to-pan only on the reference track**

On pointer down, store the pointer x and the starting interval. While dragging, translate pixel delta into reference-coordinate delta, clamp to `fullStart`/`fullEnd`, and redraw. On pointer up/leave, stop dragging.

Set `track.style.cursor` between `grab` and `grabbing`; do not make the overview draggable.

- [ ] **Step 5: Synchronize candidate selection across band, table, and detail card**

Implement:

```javascript
function selectRegion(regionId) {
  selectedRegionId=regionId;
  document.querySelectorAll('.candidate-row').forEach(row => {
    row.classList.toggle('selected', row.dataset.regionId===regionId);
  });
  const region=candidateRegions.find(item=>item.id===regionId);
  candidateDetail.textContent = region
    ? `${region.id} · rank ${region.rank} · ${region.start}-${region.end} · anchor ${region.peakStart}-${region.peakEnd} · worst ${region.worstDataset} · similarity ${region.similarity} · margin ${region.margin}`
    : 'Select a candidate region to inspect its recorded evidence.';
  drawTrack();
}
```

Candidate row click and Enter/Space keyboard activation call `selectRegion`.

Track click hit-tests candidate bands first; clicking a band calls `selectRegion`.

Keep the detail card read-only text. Do not introduce form inputs for region coordinates.

- [ ] **Step 6: Extend hover evidence**

Retain the existing hover panel. For the overview, show dataset, window, similarity, target conservation, sequence id, and orientation when the closest visible point is a challenge overlay. For target/worst-challenge points, show the worst dataset and recorded window metrics.

For the reference track, show the nearest visible target/worst-challenge point and its reference interval.

- [ ] **Step 7: Run focused tests until GREEN**

Run:

```bash
pytest tests/test_contrastive_report_html.py -q
```

Expected: PASS.

- [ ] **Step 8: Commit the interactive behavior**

```bash
git add qpcr_pipeline/contrastive_report_html.py tests/test_contrastive_report_html.py
git commit -m "feat: make contrast report interactive"
```

---

### Task 4: Regression verification

**Files:**
- No new files.
- Verify: `qpcr_pipeline/contrastive_report_html.py`
- Verify: `tests/test_contrastive_report_html.py`

**Interfaces:**
- Consumes: completed explorer implementation.
- Produces: evidence that report-only changes did not regress the scientific pipeline.

- [ ] **Step 1: Run the focused report tests**

```bash
pytest tests/test_contrastive_report_html.py -q
```

Expected: PASS.

- [ ] **Step 2: Run the normal repository test suite**

```bash
pytest -q
```

Expected: PASS.

- [ ] **Step 3: Run the repository integration gate used for `main`**

Use the existing CI configuration unchanged. Confirm the branch CI completes successfully; do not modify `.circleci/config.yml` just to make the feature branch run unless the repository already permits it.

- [ ] **Step 4: Review the final diff against `main`**

The changed implementation surface should remain limited to:

```text
qpcr_pipeline/contrastive_report_html.py
tests/test_contrastive_report_html.py
docs/superpowers/specs/2026-09-09-interactive-contrast-explorer-design.md
docs/superpowers/plans/2026-09-09-interactive-contrast-explorer.md
```

Reject unrelated refactors or scientific-algorithm changes.

- [ ] **Step 5: Commit any test-only correction if required**

Only if verification revealed a legitimate issue:

```bash
git add qpcr_pipeline/contrastive_report_html.py tests/test_contrastive_report_html.py
git commit -m "test: tighten contrast explorer regression coverage"
```
