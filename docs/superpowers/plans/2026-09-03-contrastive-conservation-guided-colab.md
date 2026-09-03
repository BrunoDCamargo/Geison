# Contrastive Conservation + Guided Colab Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a first-class `contrastive_conservation` stage between target conservation and primer design, then expose the workflow through a guided Colab notebook that remains CLI-first, reproducible, and synthetic by default.

**Architecture:** Keep Geison as the scientific source of truth. The core resolves approved CHALLENGE datasets from the Panel Model, measures region-level target/non-target contrast, consolidates overlapping evidence into deterministic candidate regions, checkpoints the new stage, and lets `primer_design` consume those regions while retaining the existing conservation-only fallback. A new guided notebook only builds configuration, calls the installed CLI, and renders published artifacts; it never reimplements conservation, contrast, primer design, inclusivity, specificity, or ranking.

**Tech Stack:** Python 3.10+, Biopython 1.85+, PyYAML 6+, stdlib dataclasses/pathlib/json/hashlib, pytest, unittest integration tests, existing Geison checkpoint/resume infrastructure, Google Colab forms, pandas/matplotlib only for notebook rendering.

**Spec:** `docs/superpowers/specs/2026-09-03-contrastive-conservation-guided-colab-design.md`

## Global Constraints

- `contrastive_conservation.enabled` defaults to `false`.
- Existing configurations with the stage omitted/disabled preserve the current conservation-only primer candidate selection behavior.
- Enabled contrastive analysis requires a frozen approved panel and at least one resolvable `CHALLENGE` non-target dataset.
- Panel criticality (`CRITICAL`, `IMPORTANT`, `BACKGROUND`) affects interpretation and deterministic ordering, not raw similarity measurement.
- The core must not embed organism-specific biological cutoffs.
- The same `off_targets` datasets are reused by contrastive analysis and final assay specificity.
- Per-dataset evidence is preserved in artifacts; no aggregate score may be the only explanation for ordering.
- The contrastive checkpoint includes only candidate-selection fields from `PrimerDesignConfig`, not primer/probe thermodynamic constraints.
- The existing `specificity` stage remains an independent assay-level check after oligo design.
- The current `notebooks/geison_colab.ipynb` remains available until the guided notebook is accepted.
- The bundled West Nile-like demonstration uses synthetic sequences only and contains no organism-specific assay sequence or hidden biological cutoff.
- The guided notebook may render artifacts with pandas/matplotlib, but must not import Geison scientific internals or BioPython scientific functions.
- `ACTION_REQUIRED`, `PARTIAL`, and `FAILED` must remain visually distinct from successful completion.

---

## File Map

### New production files

- `qpcr_pipeline/region_selection.py`: shared target-region geometry, target eligibility, overlap logic, and legacy conservation-only candidate selection.
- `qpcr_pipeline/challenge_panel.py`: deterministic mapping from approved panel CHALLENGE entries to existing `off_targets` datasets.
- `qpcr_pipeline/contrastive_similarity.py`: focused region-similarity interface and deterministic Biopython implementation.
- `qpcr_pipeline/contrastive_conservation.py`: typed result model, per-window/per-dataset evidence, deterministic ordering, region consolidation, TSV/JSON artifacts.
- `qpcr_pipeline/contrastive_report_html.py`: self-contained stage-specific contrast report.

### Modified production files

- `qpcr_pipeline/config.py`: add `ContrastiveConservationConfig`, YAML parsing, and validation.
- `qpcr_pipeline/execution.py`: insert `contrastive_conservation` after `conservation`.
- `qpcr_pipeline/checkpoint_codecs.py`: add `CONTRASTIVE_CONSERVATION_CODEC`.
- `qpcr_pipeline/checkpoint_stages.py`: parameters, panel/challenge identities, outputs, and Primer3 tool identity boundary.
- `qpcr_pipeline/pipeline.py`: execute the new stage, feed it into primer design, and expose summary counts.
- `qpcr_pipeline/primer_design.py`: consume contrastive regions when available and report candidate source.
- `tests/pipeline_checkpoint_fixtures.py`: add typed checkpoint fixture for the new stage.
- `README.md`: describe the new stage and legacy fallback.
- `docs/colab.md`: link the guided notebook while keeping the operational notebook documented.

### New tests

- `tests/test_region_selection.py`
- `tests/test_challenge_panel.py`
- `tests/test_contrastive_similarity.py`
- `tests/test_contrastive_conservation.py`
- `tests/test_contrastive_report_html.py`
- `tests/test_guided_colab_notebook.py`
- `tests/test_guided_demo.py`
- `integration_tests/test_guided_contrastive_demo.py`

### New guided demo files

- `examples/guided_demo/generate_demo_data.py`
- `notebooks/geison_guided_colab.ipynb`
- `docs/guided-colab.md`

### Existing tests expected to change

- `tests/test_config.py`
- `tests/test_execution.py`
- `tests/test_execution_plan.py`
- `tests/test_checkpoint_codecs.py`
- `tests/test_checkpoint_stages.py`
- `tests/test_pipeline_resume.py`
- `tests/test_minimal_run.py`
- `tests/test_primer_design.py`
- `tests/test_colab_notebook.py` only if shared documentation assertions need a second notebook link; do not weaken its existing low-level notebook checks.

---

### Task 1: Add contrastive configuration and execution-graph node

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Modify: `qpcr_pipeline/execution.py`
- Modify: `tests/test_config.py`
- Modify: `tests/test_execution.py`
- Modify: `tests/test_execution_plan.py`

**Interfaces:**
- Produces `ContrastiveConservationConfig(enabled: bool = False)`.
- Adds `PipelineConfig.contrastive_conservation: ContrastiveConservationConfig`.
- Adds `contrastive_conservation` to `StageName`, `STAGE_ORDER`, and `STAGE_DEPENDENCIES`.
- `primer_design` depends on `contrastive_conservation`; disabled mode still executes a checkpointed SKIPPED contrastive stage, so the graph is structurally stable.

- [ ] **Step 1: Write failing configuration tests**

Add tests that require the new dataclass, default-disabled parsing, enabled validation, and frozen-panel requirement:

```python
from qpcr_pipeline.config import (
    ContrastiveConservationConfig,
    PipelineConfig,
    ConservationConfig,
)


def test_contrastive_conservation_defaults_disabled():
    config = PipelineConfig(target_name="synthetic", input_fasta=Path("target.fa"))
    assert config.contrastive_conservation == ContrastiveConservationConfig(enabled=False)


def test_enabled_contrast_requires_enabled_conservation(tmp_path):
    fasta = tmp_path / "target.fa"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    config = PipelineConfig(
        target_name="synthetic",
        input_fasta=fasta,
        contrastive_conservation=ContrastiveConservationConfig(enabled=True),
    )
    with pytest.raises(ValueError, match="requires enabled conservation"):
        validate_pipeline_config(config)
```

Also add YAML coverage for:

```yaml
contrastive_conservation:
  enabled: true
```

Expected enabled validation rules:

```python
if config.contrastive_conservation.enabled and not config.conservation.enabled:
    raise ValueError("Enabled contrastive conservation requires enabled conservation.")
if config.contrastive_conservation.enabled:
    if config.panel is None or config.panel.frozen_manifest is None:
        raise ValueError("Enabled contrastive conservation requires an approved frozen panel.")
    if not config.off_targets:
        raise ValueError("Enabled contrastive conservation requires off-target datasets.")
```

- [ ] **Step 2: Run the focused config tests and confirm RED**

Run:

```bash
python -m pytest -q tests/test_config.py
```

Expected: failures because `ContrastiveConservationConfig` and parsing do not exist.

- [ ] **Step 3: Implement the minimal config parser/validator**

Add:

```python
@dataclass(frozen=True, slots=True)
class ContrastiveConservationConfig:
    enabled: bool = False
```

Parse `raw.get("contrastive_conservation", {})` with the same strict-boolean pattern used by other enabled stages. Insert the new field between `conservation` and `primer_design` in `PipelineConfig`.

- [ ] **Step 4: Write the failing graph test**

Require:

```python
assert STAGE_ORDER == (
    "panel", "input", "qc", "clustering", "alignment", "conservation",
    "contrastive_conservation", "primer_design", "inclusivity", "specificity", "ranking",
)
assert STAGE_DEPENDENCIES["contrastive_conservation"] == ("conservation",)
assert STAGE_DEPENDENCIES["primer_design"] == ("contrastive_conservation",)
```

- [ ] **Step 5: Run graph tests and confirm RED**

Run:

```bash
python -m pytest -q tests/test_execution.py tests/test_execution_plan.py
```

Expected: stage-order assertions fail.

- [ ] **Step 6: Update `execution.py` and re-run focused tests**

Run:

```bash
python -m pytest -q tests/test_config.py tests/test_execution.py tests/test_execution_plan.py
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add qpcr_pipeline/config.py qpcr_pipeline/execution.py tests/test_config.py tests/test_execution.py tests/test_execution_plan.py
git commit -m "feat: add contrastive conservation stage config"
```

---

### Task 2: Extract shared target-region selection primitives without behavior change

**Files:**
- Create: `qpcr_pipeline/region_selection.py`
- Modify: `qpcr_pipeline/primer_design.py`
- Create: `tests/test_region_selection.py`
- Modify: `tests/test_primer_design.py`

**Interfaces:**
- Produces/re-exports `CandidateRegion`.
- Produces `candidate_region_from_window(conservation, window, config) -> CandidateRegion`.
- Produces `is_target_eligible(region, config) -> bool`.
- Produces `overlap_fraction(left, right) -> float`.
- Produces `select_conservation_candidate_regions(conservation, config) -> tuple[CandidateRegion, ...]`.
- `qpcr_pipeline.primer_design.CandidateRegion` remains import-compatible by importing the class into that module namespace.

- [ ] **Step 1: Write characterization tests around the current selector**

Use a small synthetic `ConservationResult` already shaped like existing primer-design tests and assert exact candidate IDs, coordinates, ordering, deduplication, and overlap suppression.

```python
regions = select_conservation_candidate_regions(conservation, config)
assert [(r.region_id, r.reference_start, r.reference_end) for r in regions] == [
    ("region-001", 1, 300),
    ("region-002", 201, 500),
]
```

The fixture values must be chosen so these intervals are deterministic under the existing rules.

- [ ] **Step 2: Run the new test and confirm RED because the module does not exist**

```bash
python -m pytest -q tests/test_region_selection.py
```

- [ ] **Step 3: Move, do not rewrite, the existing candidate geometry/eligibility logic**

Move the current `CandidateRegion`, interval expansion, metric aggregation, target eligibility, overlap calculation, and conservation-only selection functions into `region_selection.py`. Keep the current ranking tuple unchanged for legacy selection.

`primer_design.py` should import:

```python
from qpcr_pipeline.region_selection import (
    CandidateRegion,
    select_conservation_candidate_regions,
)
```

and replace the private selector call with the new public helper.

- [ ] **Step 4: Run the full primer-design test file**

```bash
python -m pytest -q tests/test_region_selection.py tests/test_primer_design.py
```

Expected: PASS with no changed public behavior.

- [ ] **Step 5: Commit**

```bash
git add qpcr_pipeline/region_selection.py qpcr_pipeline/primer_design.py tests/test_region_selection.py tests/test_primer_design.py
git commit -m "refactor: share candidate region selection"
```

---

### Task 3: Resolve approved CHALLENGE panel entries to existing off-target datasets

**Files:**
- Create: `qpcr_pipeline/challenge_panel.py`
- Create: `tests/test_challenge_panel.py`
- Modify: `tests/panel_fixtures.py`

**Interfaces:**

```python
@dataclass(frozen=True, slots=True)
class ChallengeDatasetBinding:
    name: str
    criticality: Criticality
    dataset: OffTargetDataset


def resolve_challenge_datasets(
    manifest: ApprovedPanelManifest,
    configs: tuple[OffTargetConfig, ...],
) -> tuple[ChallengeDatasetBinding, ...]: ...
```

Normalization is `value.strip().casefold()`. Returned order follows approved panel non-target order, not off-target configuration order.

- [ ] **Step 1: Write failing resolver tests**

Cover all required rules:

```python
bindings = resolve_challenge_datasets(manifest, configs)
assert [(item.name, item.criticality) for item in bindings] == [
    ("critical-a", "CRITICAL"),
    ("important-b", "IMPORTANT"),
]
```

Add tests that:

- ignore panel non-targets without `CHALLENGE` role;
- reject missing mapping with `ValueError("Challenge dataset missing ...")`;
- reject duplicate normalized off-target names;
- preserve panel criticality even if file/config order differs;
- load datasets through `load_off_target_dataset()` so specificity and contrast share one physical data path.

- [ ] **Step 2: Run and verify RED**

```bash
python -m pytest -q tests/test_challenge_panel.py
```

- [ ] **Step 3: Implement the resolver**

Core algorithm:

```python
by_name = {}
for config in configs:
    key = config.name.strip().casefold()
    if key in by_name:
        raise ValueError("Off-target dataset names must be unique after normalization.")
    by_name[key] = config

resolved = []
for item in manifest.definition.non_targets:
    if "CHALLENGE" not in item.dataset_roles:
        continue
    key = item.name.strip().casefold()
    config = by_name.get(key)
    if config is None:
        raise ValueError(f"Challenge dataset missing for approved panel entry {item.name!r}.")
    resolved.append(
        ChallengeDatasetBinding(
            name=item.name,
            criticality=item.criticality,
            dataset=load_off_target_dataset(config),
        )
    )
return tuple(resolved)
```

- [ ] **Step 4: Re-run focused tests**

```bash
python -m pytest -q tests/test_challenge_panel.py tests/test_panel.py tests/test_panel_manifest.py
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add qpcr_pipeline/challenge_panel.py tests/test_challenge_panel.py tests/panel_fixtures.py
git commit -m "feat: resolve challenge panel datasets"
```

---

### Task 4: Add a focused, swappable region-similarity engine

**Files:**
- Create: `qpcr_pipeline/contrastive_similarity.py`
- Create: `tests/test_contrastive_similarity.py`

**Interfaces:**

```python
@dataclass(frozen=True, slots=True)
class RegionSimilarity:
    sequence_id: str
    similarity: float
    orientation: Literal["forward", "reverse"]


class RegionSimilarityEngine(Protocol):
    def best_match(
        self,
        query: str,
        records: tuple[LocalSequenceRecord, ...],
    ) -> RegionSimilarity | None: ...


class BiopythonLocalSimilarityEngine:
    def best_match(...) -> RegionSimilarity | None: ...
```

The v1 measurement is a continuous normalized local similarity, not a biological PASS threshold.

- [ ] **Step 1: Write failing pure-engine tests**

Require these relational properties:

```python
engine = BiopythonLocalSimilarityEngine()
identical = engine.best_match("ACGTACGT", (record("same", "TTACGTACGTGG"),))
divergent = engine.best_match("ACGTACGT", (record("different", "TTTTTTTTTTTT"),))
assert identical is not None and identical.similarity == 1.0
assert divergent is not None and divergent.similarity < identical.similarity
```

Also assert that reverse-complement evidence can win and that an empty dataset returns `None`.

- [ ] **Step 2: Run and confirm RED**

```bash
python -m pytest -q tests/test_contrastive_similarity.py
```

- [ ] **Step 3: Implement with `Bio.Align.PairwiseAligner`**

Use one configured aligner:

```python
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 1.0
aligner.mismatch_score = 0.0
aligner.open_gap_score = -1.0
aligner.extend_gap_score = -0.5
```

For each record, score forward and reverse-complement sequence and normalize:

```python
similarity = max(0.0, min(1.0, aligner.score(query, subject) / len(query)))
```

Tie-break deterministically by:

1. larger similarity;
2. lexicographically smaller `sequence_id`;
3. `forward` before `reverse`.

Reject blank query strings with `ValueError`.

- [ ] **Step 4: Re-run focused tests**

```bash
python -m pytest -q tests/test_contrastive_similarity.py
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add qpcr_pipeline/contrastive_similarity.py tests/test_contrastive_similarity.py
git commit -m "feat: add contrastive region similarity engine"
```

---

### Task 5: Implement typed contrastive evidence, transparent ordering, and artifacts

**Files:**
- Create: `qpcr_pipeline/contrastive_conservation.py`
- Create: `tests/test_contrastive_conservation.py`

**Interfaces:**

```python
@dataclass(frozen=True, slots=True)
class DatasetWindowEvidence:
    reference_start: int
    reference_end: int
    dataset_name: str
    criticality: Criticality
    sequence_count: int
    best_sequence_id: str | None
    best_orientation: str | None
    similarity: float


@dataclass(frozen=True, slots=True)
class ContrastWindowEvidence:
    reference_start: int
    reference_end: int
    target_mean_conservation: float
    target_minimum_conservation: float
    target_mean_coverage: float
    target_mean_gap_frequency: float
    target_mean_entropy_bits: float
    target_eligible: bool
    worst_dataset_name: str | None
    worst_dataset_criticality: Criticality | None
    worst_similarity: float | None
    worst_critical_similarity: float | None
    worst_important_similarity: float | None
    contrast_margin: float | None


@dataclass(frozen=True, slots=True)
class ContrastCandidateRegion:
    region: CandidateRegion
    contributing_windows: tuple[tuple[int, int], ...]
    worst_dataset_name: str | None
    worst_dataset_criticality: Criticality | None
    worst_similarity: float | None
    worst_critical_similarity: float | None
    worst_important_similarity: float | None
    contrast_margin: float | None


@dataclass(frozen=True, slots=True)
class ContrastiveConservationResult:
    status: Literal["SKIPPED", "COMPLETE"]
    reference_id: str | None
    windows: tuple[ContrastWindowEvidence, ...]
    dataset_evidence: tuple[DatasetWindowEvidence, ...]
    candidates: tuple[ContrastCandidateRegion, ...]
    challenge_dataset_names: tuple[str, ...]
    window_metrics_path: Path | None
    dataset_metrics_path: Path | None
    candidate_regions_path: Path | None
    report_path: Path
    html_report_path: Path | None
```

Primary entry point:

```python
def analyze_contrastive_conservation(
    conservation: ConservationResult,
    approved_panel: ApprovedPanelManifest | None,
    off_target_configs: tuple[OffTargetConfig, ...],
    config: ContrastiveConservationConfig,
    primer_config: PrimerDesignConfig,
    output_dir: Path,
    *,
    similarity_engine: RegionSimilarityEngine | None = None,
) -> ContrastiveConservationResult: ...
```

- [ ] **Step 1: Write the disabled-path test**

```python
result = analyze_contrastive_conservation(
    conservation,
    None,
    (),
    ContrastiveConservationConfig(enabled=False),
    primer_config,
    outdir,
)
assert result.status == "SKIPPED"
assert result.windows == ()
assert result.candidates == ()
assert result.report_path.is_file()
```

The disabled path must not load panel or off-target datasets.

- [ ] **Step 2: Write the deterministic synthetic contrast test**

Create a 1,200-base synthetic conservation fixture with two target-stable areas:

- one area whose challenge records remain similar;
- one area whose challenge records are deliberately less similar.

Use an injected fake `RegionSimilarityEngine` keyed by query coordinates so the unit test tests stage policy rather than Biopython behavior.

Assert:

```python
assert result.status == "COMPLETE"
assert result.candidates
assert result.candidates[0].region.peak_start >= 650
assert result.candidates[0].region.peak_end <= 950
```

Also include overlapping peak windows and assert they consolidate into one selected region with multiple `contributing_windows`.

- [ ] **Step 3: Write the criticality-ordering test**

Use equal target metrics and challenge similarities that differ by criticality. Required deterministic ordering key for target-eligible windows is:

```python
(
    worst_critical_similarity_or_minus_one,
    worst_important_similarity_or_minus_one,
    worst_similarity_or_minus_one,
    -target_mean_conservation,
    target_mean_entropy_bits,
    reference_start,
    reference_end,
)
```

Lower similarity sorts earlier. Missing criticality-class evidence uses `-1.0`, so a region is not penalized for a class absent from the approved panel. This is an ordering policy only; raw values remain unchanged.

- [ ] **Step 4: Run the new tests and confirm RED**

```bash
python -m pytest -q tests/test_contrastive_conservation.py
```

- [ ] **Step 5: Implement target-window measurement**

For each `conservation.windows` entry:

1. slice `conservation.major_consensus[reference_start - 1:reference_end]`;
2. build the expanded `CandidateRegion` with `candidate_region_from_window()`;
3. compute `target_eligible = is_target_eligible(region, primer_config)`;
4. call the similarity engine once per resolved challenge dataset;
5. preserve one `DatasetWindowEvidence` per window × dataset;
6. compute worst-overall, worst-CRITICAL, and worst-IMPORTANT summaries;
7. compute `contrast_margin = target_mean_conservation - worst_similarity` only as an explanatory value, not the sole ranking policy.

- [ ] **Step 6: Implement deterministic consolidation**

Group target-eligible windows by expanded `(reference_start, reference_end)` interval. For each exact interval, retain all contributing peak-window coordinates and choose the best aggregate evidence by the ordering key. Sort unique intervals by the same key, then apply `max_region_overlap_fraction` using `overlap_fraction()`. Assign stable IDs:

```python
contrast-region-001
contrast-region-002
```

Do not apply a hidden contrast threshold.

- [ ] **Step 7: Publish TSV and JSON artifacts atomically**

Paths:

```text
contrastive_conservation/window_metrics.tsv
contrastive_conservation/dataset_metrics.tsv
contrastive_conservation/candidate_regions.tsv
contrastive_conservation/contrastive_conservation_report.json
```

The JSON report must include configuration, approved panel identity, challenge dataset provenance, counts, criticality summaries, artifact paths, and the explicit ordering policy fields.

- [ ] **Step 8: Re-run focused tests**

```bash
python -m pytest -q tests/test_contrastive_conservation.py tests/test_region_selection.py tests/test_challenge_panel.py tests/test_contrastive_similarity.py
```

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add qpcr_pipeline/contrastive_conservation.py tests/test_contrastive_conservation.py
git commit -m "feat: add contrastive conservation analysis"
```

---

### Task 6: Make primer design consume contrastive candidate regions without breaking legacy behavior

**Files:**
- Modify: `qpcr_pipeline/primer_design.py`
- Modify: `tests/test_primer_design.py`

**Interfaces:**

Change:

```python
def design_primers(
    conservation: ConservationResult,
    config: PrimerDesignConfig,
    output_dir: Path,
    *,
    contrastive: ContrastiveConservationResult | None = None,
    runner: Primer3Runner | None = None,
) -> PrimerDesignResult: ...
```

Add to `PrimerDesignResult`:

```python
candidate_source: Literal["CONSERVATION_ONLY", "CONTRASTIVE_CONSERVATION"]
```

Change `primer3_required()` to accept the same optional contrastive result.

- [ ] **Step 1: Write the failing legacy compatibility test**

Call `design_primers()` with `contrastive=None` and assert the same candidates/assays as the existing fixture plus:

```python
assert result.candidate_source == "CONSERVATION_ONLY"
```

- [ ] **Step 2: Write the failing contrastive-source test**

Create a `ContrastiveConservationResult(status="COMPLETE", candidates=(...))` with one `CandidateRegion` whose ID is `contrast-region-001`. Assert:

```python
assert result.candidate_source == "CONTRASTIVE_CONSERVATION"
assert [candidate.region_id for candidate in result.candidates] == ["contrast-region-001"]
```

A COMPLETE contrastive result with zero candidates must yield zero primer candidates; it must not silently fall back to conservation-only selection.

- [ ] **Step 3: Run and confirm RED**

```bash
python -m pytest -q tests/test_primer_design.py
```

- [ ] **Step 4: Implement the two-path selector**

```python
if contrastive is not None and contrastive.status == "COMPLETE":
    candidates = tuple(item.region for item in contrastive.candidates)
    candidate_source = "CONTRASTIVE_CONSERVATION"
else:
    candidates = select_conservation_candidate_regions(conservation, config)
    candidate_source = "CONSERVATION_ONLY"
```

Keep target consensus input and Primer3 parsing unchanged.

- [ ] **Step 5: Include `candidate_source` in `primer_design_report.json` and typed result**

Do not modify assay oligo calculations.

- [ ] **Step 6: Re-run tests**

```bash
python -m pytest -q tests/test_primer_design.py tests/test_primer3.py
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add qpcr_pipeline/primer_design.py tests/test_primer_design.py
git commit -m "feat: design assays from contrastive regions"
```

---

### Task 7: Add a self-contained contrast report

**Files:**
- Create: `qpcr_pipeline/contrastive_report_html.py`
- Modify: `qpcr_pipeline/contrastive_conservation.py`
- Create: `tests/test_contrastive_report_html.py`

**Interfaces:**

```python
def render_contrastive_html(
    *,
    target_name: str,
    reference_id: str | None,
    windows: tuple[ContrastWindowEvidence, ...],
    dataset_evidence: tuple[DatasetWindowEvidence, ...],
    candidates: tuple[ContrastCandidateRegion, ...],
) -> str: ...
```

- [ ] **Step 1: Write a failing HTML safety/structure test**

Require:

```python
html = render_contrastive_html(...)
assert "Target vs non-target contrast" in html
assert "<canvas" in html
assert "candidate region" in html.lower()
assert "worst challenge" in html.lower()
assert "https://" not in html
assert "http://" not in html
```

Also inject a dataset name containing `<script>` and assert it is JSON/HTML escaped and cannot terminate a script block.

- [ ] **Step 2: Run and confirm RED**

```bash
python -m pytest -q tests/test_contrastive_report_html.py
```

- [ ] **Step 3: Implement a deterministic Canvas report using the existing offline report style**

Render:

1. summary cards: target/reference/windows/challenge datasets/candidate regions;
2. quadrant plot: X = worst challenge similarity, Y = target conservation;
3. reference track: target conservation and worst challenge similarity across reference;
4. read-only hover details;
5. candidate region table/cards including worst dataset and criticality;
6. expandable/per-dataset evidence table.

No click-to-edit behavior.

- [ ] **Step 4: Publish the HTML at**

```text
contrastive_conservation/report.html
```

Set `ContrastiveConservationResult.html_report_path` when enabled; keep it `None` in the disabled path.

- [ ] **Step 5: Run tests**

```bash
python -m pytest -q tests/test_contrastive_report_html.py tests/test_contrastive_conservation.py
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add qpcr_pipeline/contrastive_report_html.py qpcr_pipeline/contrastive_conservation.py tests/test_contrastive_report_html.py tests/test_contrastive_conservation.py
git commit -m "feat: render contrastive conservation report"
```

---

### Task 8: Integrate checkpointing, invalidation, resume, and pipeline execution

**Files:**
- Modify: `qpcr_pipeline/checkpoint_codecs.py`
- Modify: `qpcr_pipeline/checkpoint_stages.py`
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `tests/pipeline_checkpoint_fixtures.py`
- Modify: `tests/test_checkpoint_codecs.py`
- Modify: `tests/test_checkpoint_stages.py`
- Modify: `tests/test_pipeline_resume.py`
- Modify: `tests/test_minimal_run.py`

**Interfaces:**
- Adds `CONTRASTIVE_CONSERVATION_CODEC`.
- Adds contrastive stage parameters and source identities.
- `_run_stage("contrastive_conservation", ...)` calls `analyze_contrastive_conservation()`.
- `_run_stage("primer_design", ...)` passes `results["contrastive_conservation"]`.

- [ ] **Step 1: Write codec RED test**

Construct a `ContrastiveConservationResult` containing nested candidates and paths inside `outdir`. Require strict round-trip equality and path escape rejection.

Run:

```bash
python -m pytest -q tests/test_checkpoint_codecs.py
```

Expected: FAIL because codec is missing.

- [ ] **Step 2: Add the structural codec**

In `checkpoint_codecs.py` import `ContrastiveConservationResult` and add:

```python
CONTRASTIVE_CONSERVATION_CODEC = _StructuralCodec[ContrastiveConservationResult](
    ContrastiveConservationResult
)
```

- [ ] **Step 3: Write checkpoint-stage parameter and identity tests**

Require `stage_parameters("contrastive_conservation", config)` to include:

```python
{
    "config": {"enabled": True},
    "region_selection": {
        "max_candidate_regions": ...,
        "candidate_region_length": ...,
        "max_region_overlap_fraction": ...,
        "min_mean_conservation": ...,
        "min_minimum_conservation": ...,
        "min_mean_coverage": ...,
        "max_mean_gap_frequency": ...,
        "max_mean_entropy_bits": ...,
        "min_usable_fraction": ...,
    },
    "challenge_panel": [
        {"name": ..., "criticality": ..., "dataset_roles": ["CHALLENGE"]},
    ],
}
```

Explicitly assert that `primer.min_tm`, `probe.max_tm`, product-size limits, and `assays_per_region` are absent from the contrastive stage parameter fingerprint.

- [ ] **Step 4: Add contrastive input identities**

When enabled, `stage_input_identities("contrastive_conservation", config)` records:

- approved panel SHA-256;
- each resolved CHALLENGE dataset name/source/record hash;
- frozen manifest hash when applicable.

Reuse the same file-hash logic already used by `specificity`.

- [ ] **Step 5: Add stage outputs**

Declare all enabled artifacts:

```python
window_metrics_path
dataset_metrics_path
candidate_regions_path
report_path
html_report_path
```

- [ ] **Step 6: Add pipeline execution and summary projection**

In `_run_stage`:

```python
if stage == "contrastive_conservation":
    approved_panel = None
    if config.contrastive_conservation.enabled:
        assert config.panel is not None and config.panel.frozen_manifest is not None
        approved_panel = load_approved_panel_manifest(config.panel.frozen_manifest)
    return analyze_contrastive_conservation(
        results["conservation"],
        approved_panel,
        config.off_targets,
        config.contrastive_conservation,
        config.primer_design,
        output_dir,
    )
```

For primer design:

```python
return design_primers(
    conservation,
    config.primer_design,
    output_dir,
    contrastive=results["contrastive_conservation"],
    runner=primer3_runner,
)
```

Add a `contrastive_conservation` block to `qc_report.json` with status, window count, challenge dataset count, candidate region count, and worst-critical evidence availability. Do not duplicate full TSV data into this summary.

- [ ] **Step 7: Update Primer3 tool identity logic**

`primer3_required()` must receive the loaded contrastive result from `stage_context` so checkpoint planning knows whether Primer3 is actually required in contrastive mode.

- [ ] **Step 8: Write invalidation tests before implementation is considered complete**

Use `dataclasses.replace()` and local synthetic off-target files to assert:

1. changing only primer/probe Tm leaves `contrastive_conservation` REUSE and reruns `primer_design` descendants;
2. changing `max_candidate_regions` reruns contrastive and descendants;
3. changing a CHALLENGE FASTA reruns contrastive and descendants but reuses `panel` through `conservation`;
4. changing panel criticality reruns contrastive and descendants but reuses target acquisition/alignment/conservation;
5. changing only specificity mismatch tolerances leaves contrastive reusable;
6. deleting `contrastive_conservation/candidate_regions.tsv` makes the contrastive checkpoint invalid;
7. disabled legacy config still resumes every stage successfully.

- [ ] **Step 9: Run focused checkpoint/pipeline tests**

```bash
python -m pytest -q \
  tests/test_checkpoint_codecs.py \
  tests/test_checkpoint_stages.py \
  tests/test_pipeline_resume.py \
  tests/test_minimal_run.py
```

Expected: PASS.

- [ ] **Step 10: Commit**

```bash
git add qpcr_pipeline/checkpoint_codecs.py qpcr_pipeline/checkpoint_stages.py qpcr_pipeline/pipeline.py tests/pipeline_checkpoint_fixtures.py tests/test_checkpoint_codecs.py tests/test_checkpoint_stages.py tests/test_pipeline_resume.py tests/test_minimal_run.py
git commit -m "feat: checkpoint contrastive conservation stage"
```

---

### Task 9: Build a deterministic West Nile-like synthetic guided-demo dataset

**Files:**
- Create: `examples/guided_demo/generate_demo_data.py`
- Create: `tests/test_guided_demo.py`

**Interfaces:**

CLI:

```bash
python examples/guided_demo/generate_demo_data.py /tmp/geison-guided-demo
```

Published files:

```text
target.fasta
challenge-related-a.fasta
challenge-related-b.fasta
challenge-context-a.fasta
panel-proposal.yaml
config-proposal.yaml
```

The generator creates synthetic sequence data only; it performs no Geison analysis and contains no real assay sequence.

- [ ] **Step 1: Write deterministic generator tests**

Run the generator twice into separate temp directories and assert byte-for-byte equality of every output file. Assert:

- target FASTA has multiple variants;
- at least one panel entry is `CRITICAL` and one is `IMPORTANT`;
- all challenge names exactly match `off_targets[].name` in generated config;
- generated config enables alignment, conservation, contrastive conservation, primer design, inclusivity, specificity, and ranking;
- no network input is configured;
- the text includes `West Nile-like` as narrative only and contains no organism-specific assay sequence field.

- [ ] **Step 2: Run and confirm RED**

```bash
python -m pytest -q tests/test_guided_demo.py
```

- [ ] **Step 3: Implement the generator with a fixed seed**

Use a fixed synthetic reference and deterministic mutations that create:

- target diversity outside a stable candidate area;
- a shared-conserved area retained in challenge datasets;
- a target-stable/challenge-divergent area;
- multiple challenge groups and criticality classes.

The generator may use Python stdlib and BioPython FASTA writing, but it must not call Geison analysis functions.

- [ ] **Step 4: Re-run generator tests**

```bash
python -m pytest -q tests/test_guided_demo.py
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add examples/guided_demo/generate_demo_data.py tests/test_guided_demo.py
git commit -m "test: add deterministic guided demo data"
```

---

### Task 10: Create the guided Colab workbench without duplicating scientific logic

**Files:**
- Create: `notebooks/geison_guided_colab.ipynb`
- Create: `tests/test_guided_colab_notebook.py`
- Create: `docs/guided-colab.md`

**Interfaces:**
- Notebook calls installed `qpcr-pipeline` commands.
- Notebook reads generated YAML, TSV, JSON, and HTML artifacts.
- Primary mode uses Colab form cells (`#@title`, `#@param`).

- [ ] **Step 1: Write static notebook RED tests**

Load the notebook as JSON and assert all required user-facing sections exist:

```python
required_markdown = [
    "Geison - Assay discovery workbench",
    "Project and panel",
    "Data readiness",
    "Target conservation",
    "Target vs non-target contrast",
    "Assay design",
    "Target coverage",
    "Specificity",
    "Final candidates",
    "Reproducibility",
]
```

Require form metadata markers and both modes:

```python
assert "#@param" in code
assert "Demo (synthetic)" in code
assert "Project" in code
```

Require CLI delegation:

```python
assert "qpcr-pipeline doctor" in code
assert "qpcr-pipeline run" in code
assert "qpcr-pipeline panel approve" in code
assert "--resume" in code
```

Reject scientific duplication:

```python
for forbidden in [
    "from qpcr_pipeline", "import qpcr_pipeline", "from Bio", "import Bio",
    "PairwiseAligner", "find_plausible_amplicons", "design_primers(",
    "analyze_conservation(", "analyze_contrastive_conservation(",
]:
    assert forbidden not in code
```

Require `run_manifest.json`, generated config visibility, advanced artifact paths, and explicit handling text for `ACTION_REQUIRED`, `PARTIAL`, and `FAILED`.

- [ ] **Step 2: Run and confirm RED**

```bash
python -m pytest -q tests/test_guided_colab_notebook.py
```

- [ ] **Step 3: Build Section 0-3: setup, mode, project/panel, data readiness**

Normal user flow:

1. clone/update `main` and install Geison/tools;
2. run `qpcr-pipeline doctor`;
3. choose `Demo (synthetic)` or `Project`;
4. demo mode calls `examples/guided_demo/generate_demo_data.py`;
5. project mode collects target FASTA path plus up to five challenge dataset slots (`name`, `path`, `criticality`) through forms;
6. notebook writes generated proposal config under `/content/geison_workbench/`;
7. notebook shows the generated panel in a review block;
8. first run reaches `ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED`;
9. explicit `APROVAR` input triggers `qpcr-pipeline panel approve`;
10. notebook writes frozen-manifest config and resumes.

Generated YAML is shown only in an Advanced cell and persisted as an artifact.

- [ ] **Step 4: Build Section 4-5: conservation and contrast rendering**

Read Geison-published artifacts only.

Use pandas/matplotlib to render:

- target conservation summary from `conservation/conservation_report.json` and `window_metrics.tsv`;
- quadrant plot from `contrastive_conservation/window_metrics.tsv`;
- reference track from the same published TSV;
- candidate region table from `candidate_regions.tsv`;
- per-dataset drill-down from `dataset_metrics.tsv`.

The notebook may transform rows for presentation but must not recompute similarity or candidate ranking.

- [ ] **Step 5: Build Section 6-9: assays, specificity, final candidates, reproducibility**

Render:

- `primer_design/primer_design_report.json` and candidate source;
- inclusivity report as `Target coverage`;
- specificity report and per-dataset evidence;
- ranking report with human-readable reason codes;
- `run_manifest.json` status and `missing_evidence`;
- downloadable/listed artifact paths and Geison commit/version.

Map states exactly:

```text
ACTION_REQUIRED -> Review and approve the panel before scientific execution.
PARTIAL -> show exact missing_evidence.
FAILED -> show failed stage and sanitized diagnostic message.
COMPLETED -> show completed only when manifest says COMPLETED.
```

- [ ] **Step 6: Document the notebook**

`docs/guided-colab.md` must explain:

- demo vs project mode;
- Panel approval gate;
- synthetic-demo scope;
- what contrast means versus final specificity;
- where generated configuration and artifacts live;
- how to restart cleanly;
- how to inspect Advanced evidence.

- [ ] **Step 7: Re-run static notebook tests**

```bash
python -m pytest -q tests/test_guided_colab_notebook.py tests/test_colab_notebook.py tests/test_guided_demo.py
```

Expected: PASS, and the old operational notebook tests remain unchanged.

- [ ] **Step 8: Commit**

```bash
git add notebooks/geison_guided_colab.ipynb tests/test_guided_colab_notebook.py docs/guided-colab.md
git commit -m "feat: add guided assay discovery Colab"
```

---

### Task 11: Add end-to-end synthetic integration coverage and documentation rollout

**Files:**
- Create: `integration_tests/test_guided_contrastive_demo.py`
- Modify: `README.md`
- Modify: `docs/colab.md`

**Interfaces:**
- Integration test uses the installed CLI and real installed MAFFT/Primer3 when available.
- It does not use network access.

- [ ] **Step 1: Write the integration test**

Use `unittest` so existing `python -m unittest discover -s integration_tests -q` picks it up. The test should:

1. generate demo data into a temp directory;
2. run the proposal config through the installed `qpcr-pipeline` executable and assert exit code `3` with `PANEL_APPROVAL_REQUIRED`;
3. approve the generated proposal through `qpcr-pipeline panel approve`;
4. produce a frozen-manifest config from the generated template;
5. run with `--resume`;
6. assert contrastive report and candidate-region artifacts exist;
7. assert `primer_design_report.json` says `candidate_source == "CONTRASTIVE_CONSERVATION"`;
8. assert specificity and ranking reports exist;
9. assert `run_manifest.json` status is `COMPLETED` for the deterministic synthetic fixture;
10. assert the top contrastive region traces to multiple contributing windows.

Decorate with `skipUnless` checks for `mafft` and `primer3_core`; no CD-HIT dependency is required unless the demo explicitly enables clustering.

- [ ] **Step 2: Run the integration test locally with required binaries**

```bash
python -m unittest integration_tests.test_guided_contrastive_demo -q
```

Expected: PASS.

If the fixture does not reach `COMPLETED`, adjust only the synthetic fixture generation or explicit synthetic config parameters; do not weaken completion rules or introduce organism-specific defaults.

- [ ] **Step 3: Update README and Colab docs**

README should present:

```text
Target conservation -> Target vs non-target contrast -> Assay design -> Specificity
```

and link both notebooks:

- guided workbench: `notebooks/geison_guided_colab.ipynb`;
- low-level operational notebook: `notebooks/geison_colab.ipynb`.

`docs/colab.md` should clearly say the existing notebook remains useful for low-level validation and point normal evaluators to the guided notebook.

- [ ] **Step 4: Run focused regression**

```bash
python -m pytest -q \
  tests/test_config.py \
  tests/test_region_selection.py \
  tests/test_challenge_panel.py \
  tests/test_contrastive_similarity.py \
  tests/test_contrastive_conservation.py \
  tests/test_contrastive_report_html.py \
  tests/test_primer_design.py \
  tests/test_checkpoint_codecs.py \
  tests/test_checkpoint_stages.py \
  tests/test_pipeline_resume.py \
  tests/test_guided_demo.py \
  tests/test_guided_colab_notebook.py \
  tests/test_colab_notebook.py
```

Expected: PASS.

- [ ] **Step 5: Run full unit regression**

```bash
python -m pytest -q
```

Expected: zero failures.

- [ ] **Step 6: Run full integration regression**

```bash
python -m unittest discover -s integration_tests -q
```

Expected: zero failures, with skips only for genuinely unavailable external binaries outside CI.

- [ ] **Step 7: Installed CLI smoke**

```bash
qpcr-pipeline doctor
qpcr-pipeline panel approve --help
```

Then execute the synthetic guided demo workflow once through the installed command rather than `python -m` imports.

- [ ] **Step 8: Commit**

```bash
git add integration_tests/test_guided_contrastive_demo.py README.md docs/colab.md
git commit -m "test: cover guided contrastive workflow"
```

---

## Final Verification Checklist

Before calling the subproject complete, verify every acceptance criterion from the spec against evidence:

- [ ] `contrastive_conservation` appears in the typed graph between conservation and primer design.
- [ ] Omitted/disabled configuration preserves conservation-only candidate behavior.
- [ ] Enabled contrast requires a frozen approved panel.
- [ ] Panel CHALLENGE mappings and criticality are preserved and deterministic.
- [ ] Per-dataset metrics are present in TSV/JSON artifacts.
- [ ] No hidden contrast threshold determines candidate acceptance.
- [ ] Candidate regions consolidate overlapping windows deterministically.
- [ ] Primer design reports `CONTRASTIVE_CONSERVATION` source when enabled.
- [ ] Final specificity remains independently executed.
- [ ] Challenge-only changes reuse target acquisition through conservation and invalidate contrastive descendants.
- [ ] Primer thermodynamic-only changes do not invalidate contrastive conservation.
- [ ] Synthetic shared-conserved and discriminant regions are distinguished in tests.
- [ ] Guided notebook contains no Geison/BioPython scientific imports.
- [ ] Normal guided flow requires no manual YAML editing.
- [ ] Generated YAML/manifests/artifacts remain inspectable.
- [ ] `ACTION_REQUIRED`, `PARTIAL`, and `FAILED` are never shown as successful completion.
- [ ] Full pytest regression passes.
- [ ] Full integration regression passes in CircleCI with installed tools.
- [ ] Guided demo remains synthetic and carries no real assay sequence or organism-specific cutoff.
- [ ] Contrast visuals are read-only; policy changes require explicit form/config change and rerun.

## Execution Notes

- Implement on an isolated branch/worktree rather than directly on `main`.
- Use TDD task-by-task: RED -> minimal GREEN -> focused regression -> commit.
- Do not fold the existing synthetic spike notebook into production code; keep it as design evidence under `notebooks/spikes/`.
- Do not redesign the final ranking score in this subproject. The guided notebook may explain regional contrast alongside ranking reasons, but ranking semantics remain unchanged unless a separate approved design changes them.
- Do not add live biological panel construction or live challenge-dataset discovery to the guided demo.
