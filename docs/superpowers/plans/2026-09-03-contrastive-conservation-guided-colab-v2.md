# Contrastive Conservation + Guided Colab Implementation Plan v2

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a first-class `contrastive_conservation` stage between target conservation and primer design, then expose the workflow through a guided Colab notebook that remains CLI-first, reproducible, and synthetic by default.

**Architecture:** Keep Geison as the scientific source of truth. Panel approval remains a pre-execution gate, but target input/conservation checkpoints are independent of panel-challenge semantics so panel criticality changes do not force target reacquisition. `contrastive_conservation` depends on both the approved panel checkpoint and completed target conservation, resolves the same CHALLENGE datasets later used by specificity, preserves per-dataset evidence, consolidates overlapping windows into deterministic regions, and feeds those regions into `primer_design`. The guided notebook only configures, invokes, and renders Geison artifacts.

**Tech Stack:** Python 3.10+, Biopython 1.85+, PyYAML 6+, stdlib dataclasses/pathlib/json/hashlib, pytest, unittest integration tests, existing Geison checkpoint/resume infrastructure, Google Colab forms, pandas/matplotlib for notebook rendering only.

**Specs:**
- `docs/superpowers/specs/2026-09-03-contrastive-conservation-guided-colab-design.md`
- `docs/superpowers/specs/2026-09-03-contrastive-conservation-checkpoint-clarification.md`

## Global Constraints

- `contrastive_conservation.enabled` defaults to `false`.
- Existing configs with the stage omitted/disabled preserve the current conservation-only primer candidate behavior.
- Enabled contrastive analysis requires a frozen approved panel and at least one resolvable `CHALLENGE` non-target dataset.
- A proposal config keeps contrastive analysis disabled until approval; the approved config enables it explicitly.
- Panel criticality (`CRITICAL`, `IMPORTANT`, `BACKGROUND`) affects interpretation and deterministic ordering, not raw similarity measurement.
- No organism-specific biological cutoff is embedded in the core.
- The same `off_targets` datasets are reused by contrastive analysis and final assay specificity.
- Per-dataset evidence is preserved; no opaque aggregate score is the only explanation for candidate ordering.
- The contrastive checkpoint includes only candidate-selection fields from `PrimerDesignConfig`, not primer/probe thermodynamic constraints.
- `specificity` remains an independent assay-level check after oligo design.
- The current `notebooks/geison_colab.ipynb` remains available until the guided notebook is accepted.
- The bundled West Nile-like demonstration uses synthetic sequences only and contains no organism-specific assay sequence or hidden biological cutoff.
- The guided notebook may render artifacts but must not import Geison scientific internals or BioPython scientific functions.
- `ACTION_REQUIRED`, `PARTIAL`, and `FAILED` must never be presented as successful completion.

---

## File Map

### New production files
- `qpcr_pipeline/region_selection.py`: shared target-region geometry, eligibility, overlap, and conservation-only selector.
- `qpcr_pipeline/challenge_panel.py`: approved CHALLENGE panel to `off_targets` resolver.
- `qpcr_pipeline/contrastive_similarity.py`: focused region-similarity interface and Biopython implementation.
- `qpcr_pipeline/contrastive_conservation.py`: typed evidence/result model, ordering, region consolidation, TSV/JSON artifacts.
- `qpcr_pipeline/contrastive_report_html.py`: self-contained read-only contrast report.

### Modified production files
- `qpcr_pipeline/config.py`
- `qpcr_pipeline/execution.py`
- `qpcr_pipeline/checkpoint_codecs.py`
- `qpcr_pipeline/checkpoint_stages.py`
- `qpcr_pipeline/pipeline.py`
- `qpcr_pipeline/primer_design.py`
- `README.md`
- `docs/colab.md`

### New tests and guided-demo files
- `tests/test_region_selection.py`
- `tests/test_challenge_panel.py`
- `tests/test_contrastive_similarity.py`
- `tests/test_contrastive_conservation.py`
- `tests/test_contrastive_report_html.py`
- `tests/test_guided_demo.py`
- `tests/test_guided_colab_notebook.py`
- `integration_tests/test_guided_contrastive_demo.py`
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
- `tests/pipeline_checkpoint_fixtures.py`

---

### Task 1: Add contrastive config and correct checkpoint dependency graph

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Modify: `qpcr_pipeline/execution.py`
- Modify: `tests/test_config.py`
- Modify: `tests/test_execution.py`
- Modify: `tests/test_execution_plan.py`

**Interfaces:**

```python
@dataclass(frozen=True, slots=True)
class ContrastiveConservationConfig:
    enabled: bool = False
```

`PipelineConfig` gains:

```python
contrastive_conservation: ContrastiveConservationConfig = field(
    default_factory=ContrastiveConservationConfig
)
```

Exact checkpoint graph:

```python
STAGE_DEPENDENCIES = {
    "panel": (),
    "input": (),
    "qc": ("input",),
    "clustering": ("qc",),
    "alignment": ("clustering",),
    "conservation": ("alignment",),
    "contrastive_conservation": ("panel", "conservation"),
    "primer_design": ("contrastive_conservation",),
    "inclusivity": ("primer_design", "qc"),
    "specificity": ("primer_design",),
    "ranking": ("primer_design", "inclusivity", "specificity"),
}
```

- [ ] **Step 1: Write config RED tests**

Add tests for default disabled parsing and these exact validation rules:

```python
if config.contrastive_conservation.enabled and not config.conservation.enabled:
    raise ValueError("Enabled contrastive conservation requires enabled conservation.")
if config.contrastive_conservation.enabled:
    if config.panel is None or config.panel.frozen_manifest is None:
        raise ValueError("Enabled contrastive conservation requires an approved frozen panel.")
    if not config.off_targets:
        raise ValueError("Enabled contrastive conservation requires off-target datasets.")
```

Also assert a config with `panel.proposal` plus `contrastive_conservation.enabled: true` is rejected before execution.

- [ ] **Step 2: Run config tests and confirm RED**

```bash
python -m pytest -q tests/test_config.py
```

- [ ] **Step 3: Implement strict parser/validator**

Parse:

```yaml
contrastive_conservation:
  enabled: false
```

with the same strict boolean handling as other stage configs.

- [ ] **Step 4: Write graph RED tests**

Require exact order:

```python
assert STAGE_ORDER == (
    "panel", "input", "qc", "clustering", "alignment", "conservation",
    "contrastive_conservation", "primer_design", "inclusivity", "specificity", "ranking",
)
assert STAGE_DEPENDENCIES["input"] == ()
assert STAGE_DEPENDENCIES["contrastive_conservation"] == ("panel", "conservation")
assert STAGE_DEPENDENCIES["primer_design"] == ("contrastive_conservation",)
```

- [ ] **Step 5: Run graph tests and confirm RED**

```bash
python -m pytest -q tests/test_execution.py tests/test_execution_plan.py
```

- [ ] **Step 6: Implement graph and re-run focused tests**

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

```python
@dataclass(frozen=True, slots=True)
class CandidateRegion:
    region_id: str
    rank: int
    reference_start: int
    reference_end: int
    peak_start: int
    peak_end: int
    position_count: int
    usable_length: int
    usable_fraction: float
    mean_conservation: float
    minimum_conservation: float
    mean_coverage: float
    mean_gap_frequency: float
    mean_entropy_bits: float


def candidate_region_from_window(
    conservation: ConservationResult,
    window: WindowConservation,
    config: PrimerDesignConfig,
) -> CandidateRegion: ...


def is_target_eligible(region: CandidateRegion, config: PrimerDesignConfig) -> bool: ...

def overlap_fraction(left: CandidateRegion, right: CandidateRegion) -> float: ...
def select_conservation_candidate_regions(
    conservation: ConservationResult,
    config: PrimerDesignConfig,
) -> tuple[CandidateRegion, ...]: ...
```

`primer_design.py` imports `CandidateRegion` into its module namespace so existing imports remain valid.

- [ ] **Step 1: Write characterization tests around current candidate behavior**

Assert exact coordinates, target metrics, ordering, exact-interval deduplication, overlap suppression, and stable `region-001` IDs using the existing primer-design fixture style.

- [ ] **Step 2: Run and confirm RED because the new module does not exist**

```bash
python -m pytest -q tests/test_region_selection.py
```

- [ ] **Step 3: Move the current private geometry/eligibility code into `region_selection.py`**

Do not change the legacy ranking tuple or eligibility thresholds during extraction.

`primer_design.py` changes to:

```python
from qpcr_pipeline.region_selection import (
    CandidateRegion,
    select_conservation_candidate_regions,
)
```

- [ ] **Step 4: Run primer-design regression**

```bash
python -m pytest -q tests/test_region_selection.py tests/test_primer_design.py
```

Expected: PASS with unchanged legacy candidates.

- [ ] **Step 5: Commit**

```bash
git add qpcr_pipeline/region_selection.py qpcr_pipeline/primer_design.py tests/test_region_selection.py tests/test_primer_design.py
git commit -m "refactor: share candidate region selection"
```

---

### Task 3: Resolve approved CHALLENGE entries to the same off-target datasets used by specificity

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

Normalization: `name.strip().casefold()`. Return order: approved panel order.

- [ ] **Step 1: Write resolver RED tests**

Cover:
- only non-targets carrying `CHALLENGE` role are selected;
- missing mapping raises `Challenge dataset missing for approved panel entry ...`;
- duplicate normalized off-target names are rejected;
- criticality always comes from the approved panel;
- panel order wins over config order;
- loading goes through `load_off_target_dataset()`.

- [ ] **Step 2: Run and confirm RED**

```bash
python -m pytest -q tests/test_challenge_panel.py
```

- [ ] **Step 3: Implement deterministic resolver**

```python
by_name: dict[str, OffTargetConfig] = {}
for config in configs:
    key = config.name.strip().casefold()
    if key in by_name:
        raise ValueError("Off-target dataset names must be unique after normalization.")
    by_name[key] = config

resolved = []
for item in manifest.definition.non_targets:
    if "CHALLENGE" not in item.dataset_roles:
        continue
    config = by_name.get(item.name.strip().casefold())
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
```

Default implementation: `BiopythonLocalSimilarityEngine`.

- [ ] **Step 1: Write pure-engine RED tests**

Require:

```python
identical = engine.best_match("ACGTACGT", (record("same", "TTACGTACGTGG"),))
divergent = engine.best_match("ACGTACGT", (record("different", "TTTTTTTTTTTT"),))
assert identical is not None and identical.similarity == 1.0
assert divergent is not None and divergent.similarity < identical.similarity
```

Also test reverse-complement winning, blank-query rejection, empty-dataset `None`, and deterministic tie-breaking.

- [ ] **Step 2: Run and confirm RED**

```bash
python -m pytest -q tests/test_contrastive_similarity.py
```

- [ ] **Step 3: Implement with one configured `Bio.Align.PairwiseAligner`**

```python
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 1.0
aligner.mismatch_score = 0.0
aligner.open_gap_score = -1.0
aligner.extend_gap_score = -0.5
```

Normalize each forward/reverse-complement score by query length:

```python
similarity = max(0.0, min(1.0, aligner.score(query, subject) / len(query)))
```

Tie-break: higher similarity, smaller sequence ID, forward before reverse.

- [ ] **Step 4: Re-run tests**

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

### Task 5: Implement typed contrastive evidence, provenance, deterministic ordering, and artifacts

**Files:**
- Create: `qpcr_pipeline/contrastive_conservation.py`
- Create: `tests/test_contrastive_conservation.py`

**Interfaces:**

```python
@dataclass(frozen=True, slots=True)
class ChallengeDatasetSummary:
    name: str
    criticality: Criticality
    source_type: str
    records_sha256: str
    sequence_count: int


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
    challenge_datasets: tuple[ChallengeDatasetSummary, ...]
    window_metrics_path: Path | None
    dataset_metrics_path: Path | None
    candidate_regions_path: Path | None
    report_path: Path
    html_report_path: Path | None
```

Entry point:

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

- [ ] **Step 1: Write disabled-path RED test**

Disabled analysis returns `SKIPPED`, no window/candidate data, writes only a SKIPPED report, and does not load panel/off-target datasets.

- [ ] **Step 2: Write deterministic synthetic contrast RED test**

Use an injected fake similarity engine and a 1,200-base conservation fixture with:
- target-stable/shared challenge area;
- target-stable/challenge-divergent area;
- overlapping peak windows that expand to the same/similar candidate intervals;
- at least one `CRITICAL` and one `IMPORTANT` challenge dataset.

Assert the best selected region falls in the deliberately discriminant zone and contains multiple `contributing_windows`.

- [ ] **Step 3: Write criticality-ordering RED test**

For target-eligible windows, exact ordering key is:

```python
(
    worst_critical_similarity if worst_critical_similarity is not None else -1.0,
    worst_important_similarity if worst_important_similarity is not None else -1.0,
    worst_similarity if worst_similarity is not None else -1.0,
    -target_mean_conservation,
    target_mean_entropy_bits,
    reference_start,
    reference_end,
)
```

Lower challenge similarity sorts earlier. This is transparent ordering, not a PASS threshold.

- [ ] **Step 4: Run and confirm RED**

```bash
python -m pytest -q tests/test_contrastive_conservation.py
```

- [ ] **Step 5: Implement per-window measurement**

For each conservation window:
1. slice `conservation.major_consensus[start - 1:end]`;
2. create the expanded target region with `candidate_region_from_window()`;
3. calculate `target_eligible = is_target_eligible(...)`;
4. evaluate each resolved challenge dataset with the injected/default similarity engine;
5. preserve one `DatasetWindowEvidence` row per window × dataset;
6. summarize worst overall/CRITICAL/IMPORTANT similarity;
7. retain `contrast_margin = target_mean_conservation - worst_similarity` for explanation only.

Build `ChallengeDatasetSummary` directly from resolved dataset provenance and panel criticality.

- [ ] **Step 6: Implement deterministic candidate consolidation**

Group target-eligible windows by expanded interval, preserve all contributing peak-window coordinates, sort unique intervals by the ordering key, then apply `max_region_overlap_fraction` using shared `overlap_fraction()`. Assign stable IDs `contrast-region-001`, `contrast-region-002`, etc. Do not apply a hidden contrast threshold.

- [ ] **Step 7: Publish atomic artifacts**

```text
contrastive_conservation/window_metrics.tsv
contrastive_conservation/dataset_metrics.tsv
contrastive_conservation/candidate_regions.tsv
contrastive_conservation/contrastive_conservation_report.json
```

Report JSON includes configuration, approved-panel identity, typed challenge provenance, counts, criticality summaries, explicit ordering fields, and artifact paths.

- [ ] **Step 8: Run focused tests**

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

### Task 6: Make primer design consume contrastive regions while preserving legacy mode

**Files:**
- Modify: `qpcr_pipeline/primer_design.py`
- Modify: `tests/test_primer_design.py`

**Interfaces:**

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

`PrimerDesignResult` gains:

```python
candidate_source: Literal["CONSERVATION_ONLY", "CONTRASTIVE_CONSERVATION"]
```

`primer3_required()` accepts the same optional contrastive result.

- [ ] **Step 1: Write legacy compatibility RED test**

With `contrastive=None`, assert existing candidate coordinates/assays remain unchanged and `candidate_source == "CONSERVATION_ONLY"`.

- [ ] **Step 2: Write contrastive-source RED test**

Pass a COMPLETE contrastive result with a `contrast-region-001` candidate and assert primer design uses exactly that region and reports `CONTRASTIVE_CONSERVATION`.

A COMPLETE contrastive result with zero candidates must yield zero primer candidates and must not fall back to conservation-only selection.

- [ ] **Step 3: Run and confirm RED**

```bash
python -m pytest -q tests/test_primer_design.py
```

- [ ] **Step 4: Implement the two-path selection**

```python
if contrastive is not None and contrastive.status == "COMPLETE":
    candidates = tuple(item.region for item in contrastive.candidates)
    candidate_source = "CONTRASTIVE_CONSERVATION"
else:
    candidates = select_conservation_candidate_regions(conservation, config)
    candidate_source = "CONSERVATION_ONLY"
```

Keep target consensus and Primer3 oligo calculations unchanged.

- [ ] **Step 5: Persist source in typed result and report JSON**

- [ ] **Step 6: Run regression**

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

### Task 7: Add the self-contained contrast report

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

- [ ] **Step 1: Write HTML RED tests**

Require `Target vs non-target contrast`, `<canvas`, candidate-region details, worst-challenge details, no external `http://`/`https://` resources, and safe escaping of a dataset name containing `<script>`.

- [ ] **Step 2: Run and confirm RED**

```bash
python -m pytest -q tests/test_contrastive_report_html.py
```

- [ ] **Step 3: Implement deterministic read-only report**

Render:
- summary cards;
- quadrant plot: X worst challenge similarity, Y target conservation;
- reference track: target conservation and worst challenge similarity;
- read-only hover details;
- candidate region table/cards;
- per-dataset evidence table.

No click-to-edit selection behavior.

- [ ] **Step 4: Publish `contrastive_conservation/report.html`**

Set `html_report_path` only for COMPLETE enabled analysis.

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

### Task 8: Integrate codec, checkpoint fingerprints, resume invalidation, and pipeline orchestration

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
- `stage_parameters("contrastive_conservation", config)` includes contrast config, challenge-panel semantics, and only the candidate-selection subset of `PrimerDesignConfig`.
- `stage_input_identities("contrastive_conservation", config)` includes approved-panel SHA plus only resolved CHALLENGE dataset identities.
- `_run_stage("contrastive_conservation")` calls the new analyzer.
- `_run_stage("primer_design")` passes the loaded contrastive result.

- [ ] **Step 1: Write codec RED test**

Require strict round-trip equality for nested challenge summaries, window evidence, candidate regions, and relative paths; reject path escape.

```bash
python -m pytest -q tests/test_checkpoint_codecs.py
```

- [ ] **Step 2: Add structural codec**

```python
CONTRASTIVE_CONSERVATION_CODEC = _StructuralCodec[ContrastiveConservationResult](
    ContrastiveConservationResult
)
```

- [ ] **Step 3: Write parameter-boundary RED tests**

Require `region_selection` fingerprint fields:

```python
(
    "max_candidate_regions",
    "candidate_region_length",
    "max_region_overlap_fraction",
    "min_mean_conservation",
    "min_minimum_conservation",
    "min_mean_coverage",
    "max_mean_gap_frequency",
    "max_mean_entropy_bits",
    "min_usable_fraction",
)
```

Explicitly assert these do **not** enter the contrastive fingerprint:
- `assays_per_region`;
- product-size limits;
- primer/probe Tm, size, GC, penalty constraints.

- [ ] **Step 4: Add panel/challenge source identities**

When enabled, hash:
- approved panel manifest;
- each mapped CHALLENGE FASTA record file, or frozen records+manifest;
- include challenge name/source in deterministic panel order.

Do not hash unrelated off-target configs that are not approved `CHALLENGE` entries.

- [ ] **Step 5: Add declared stage outputs**

```python
window_metrics_path
dataset_metrics_path
candidate_regions_path
report_path
html_report_path
```

- [ ] **Step 6: Add pipeline execution**

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

Primer design:

```python
return design_primers(
    results["conservation"],
    config.primer_design,
    output_dir,
    contrastive=results["contrastive_conservation"],
    runner=primer3_runner,
)
```

Add a compact `contrastive_conservation` projection to `qc_report.json`: status, window count, challenge dataset count, candidate region count.

- [ ] **Step 7: Make Primer3 tool identity contrast-aware**

Pass `stage_context.get("contrastive_conservation")` into `primer3_required()` so planning uses the actual candidate source.

- [ ] **Step 8: Write invalidation RED tests**

Required outcomes:

1. primer/probe Tm-only change: `panel` through `contrastive_conservation` REUSE; `primer_design` and descendants rerun.
2. `max_candidate_regions` change: target chain through `conservation` REUSE; `contrastive_conservation` and descendants rerun.
3. CHALLENGE FASTA content change: `panel`, `input`, `qc`, `clustering`, `alignment`, `conservation` REUSE; contrastive descendants rerun.
4. panel criticality or CHALLENGE membership change: `input` through `conservation` REUSE; `panel` and contrastive descendants rerun.
5. specificity tolerance-only change: contrastive and primer design REUSE; specificity and ranking rerun according to existing dependencies.
6. deleted contrastive declared artifact: contrastive checkpoint invalid, target chain reusable.
7. disabled legacy config: resume reuses every valid stage including a SKIPPED contrastive checkpoint.

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

### Task 9: Build deterministic synthetic guided-demo inputs and the proposal/approved config transition

**Files:**
- Create: `examples/guided_demo/generate_demo_data.py`
- Create: `tests/test_guided_demo.py`

**Interface:**

```bash
python examples/guided_demo/generate_demo_data.py /tmp/geison-guided-demo
```

Generated files:

```text
target.fasta
challenge-related-a.fasta
challenge-related-b.fasta
challenge-context-a.fasta
panel-proposal.yaml
config-proposal.yaml
config-approved-template.yaml
```

- [ ] **Step 1: Write deterministic generator RED tests**

Generate twice and compare every file byte-for-byte. Assert:
- multiple synthetic target variants;
- at least one `CRITICAL` and one `IMPORTANT` challenge entry;
- exact name match between panel CHALLENGE entries and `off_targets`;
- no network input;
- `config-proposal.yaml` contains `panel.proposal` and `contrastive_conservation.enabled: false`;
- `config-approved-template.yaml` contains a placeholder/path slot for `panel.frozen_manifest` and `contrastive_conservation.enabled: true`;
- approved template enables alignment, conservation, primer design, inclusivity, specificity, and ranking;
- demo labels are West Nile-like narrative only and contain no real assay sequence field.

- [ ] **Step 2: Run and confirm RED**

```bash
python -m pytest -q tests/test_guided_demo.py
```

- [ ] **Step 3: Implement fixed-seed synthetic fixture generation**

Generate a stable target reference and deterministic variants with:
- target diversity outside a stable region;
- one shared-conserved challenge region;
- one target-stable/challenge-divergent region;
- multiple challenge criticality classes.

The generator writes input/config files only and calls no Geison analysis function.

- [ ] **Step 4: Re-run tests**

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

### Task 10: Create the guided Colab workbench without duplicated scientific logic

**Files:**
- Create: `notebooks/geison_guided_colab.ipynb`
- Create: `tests/test_guided_colab_notebook.py`
- Create: `docs/guided-colab.md`

**Interfaces:**
- Notebook invokes installed `qpcr-pipeline` commands.
- Notebook reads generated YAML/TSV/JSON/HTML artifacts.
- Common inputs use Colab forms (`#@title`, `#@param`).

- [ ] **Step 1: Write static notebook RED tests**

Require user-facing sections:

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

Require CLI delegation and forms:

```python
for marker in [
    "#@param", "Demo (synthetic)", "Project",
    "qpcr-pipeline doctor", "qpcr-pipeline run",
    "qpcr-pipeline panel approve", "--resume",
]:
    assert marker in code
```

Reject duplicated scientific implementation:

```python
for forbidden in [
    "from qpcr_pipeline", "import qpcr_pipeline", "from Bio", "import Bio",
    "PairwiseAligner", "find_plausible_amplicons", "design_primers(",
    "analyze_conservation(", "analyze_contrastive_conservation(",
]:
    assert forbidden not in code
```

Require `run_manifest.json`, generated config visibility, artifact paths, and explicit UI text for `ACTION_REQUIRED`, `PARTIAL`, `FAILED`, and `COMPLETED`.

- [ ] **Step 2: Run and confirm RED**

```bash
python -m pytest -q tests/test_guided_colab_notebook.py
```

- [ ] **Step 3: Build setup + proposal review flow**

Normal flow:
1. clone/update `main`, install Geison/tools, run `qpcr-pipeline doctor`;
2. choose `Demo (synthetic)` or `Project`;
3. demo mode calls the generator;
4. project mode collects target FASTA plus up to five challenge dataset name/path/criticality slots via forms;
5. write/preserve `config-proposal.yaml` with contrastive disabled;
6. run proposal config and assert `ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED`;
7. display `panel_proposal.yaml` in a review card;
8. explicit `APROVAR` input calls `qpcr-pipeline panel approve`;
9. materialize/preserve `config-approved.yaml` from the approved template with the frozen manifest path and contrastive enabled;
10. resume with approved config.

Both configs remain inspectable under Advanced.

- [ ] **Step 4: Build conservation + contrast views from artifacts only**

Use pandas/matplotlib to render:
- conservation summary/report;
- quadrant plot from `contrastive_conservation/window_metrics.tsv`;
- reference track from the same published metrics;
- candidate table from `candidate_regions.tsv`;
- per-dataset detail from `dataset_metrics.tsv`.

No similarity or ranking is recomputed in notebook cells.

- [ ] **Step 5: Build assay, target coverage, specificity, final candidate, and reproducibility views**

Render published report JSON/TSV files and `run_manifest.json`. Map states exactly:

```text
ACTION_REQUIRED -> Review and approve the panel before scientific execution.
PARTIAL -> show exact missing_evidence.
FAILED -> show failed stage and sanitized diagnostic message.
COMPLETED -> show completed only when run_manifest.json says COMPLETED.
```

Advanced view exposes generated YAML, panel manifest, tool versions, checkpoints, raw artifacts, and commit/version.

- [ ] **Step 6: Document normal and advanced use**

`docs/guided-colab.md` explains demo/project mode, proposal/approval transition, synthetic-demo boundary, contrast versus final specificity, clean restart, and artifact locations.

- [ ] **Step 7: Run notebook regression**

```bash
python -m pytest -q tests/test_guided_colab_notebook.py tests/test_colab_notebook.py tests/test_guided_demo.py
```

Expected: PASS; existing low-level notebook checks remain intact.

- [ ] **Step 8: Commit**

```bash
git add notebooks/geison_guided_colab.ipynb tests/test_guided_colab_notebook.py docs/guided-colab.md
git commit -m "feat: add guided assay discovery Colab"
```

---

### Task 11: Add end-to-end synthetic integration coverage and rollout docs

**Files:**
- Create: `integration_tests/test_guided_contrastive_demo.py`
- Modify: `README.md`
- Modify: `docs/colab.md`

- [ ] **Step 1: Write installed-CLI integration test**

Use `unittest` so existing CircleCI discovery finds it. The test:
1. generates demo files in a temp directory;
2. runs `qpcr-pipeline run config-proposal.yaml --outdir ...` and asserts exit `3` + `PANEL_APPROVAL_REQUIRED`;
3. approves the proposal with installed `qpcr-pipeline panel approve`;
4. materializes `config-approved.yaml` from the template;
5. reruns with `--resume`;
6. asserts contrastive TSV/JSON/HTML artifacts exist;
7. asserts multiple raw windows contribute to at least one contrast candidate region;
8. asserts primer report says `candidate_source == "CONTRASTIVE_CONSERVATION"`;
9. asserts specificity and ranking reports exist;
10. asserts `run_manifest.json` is `COMPLETED` for the deterministic synthetic fixture.

Decorate with `skipUnless` for `mafft` and `primer3_core`; do not require network access.

- [ ] **Step 2: Run focused integration**

```bash
python -m unittest integration_tests.test_guided_contrastive_demo -q
```

Expected: PASS.

If completion fails, adjust only synthetic fixture generation or explicit synthetic-demo config values. Do not weaken completion rules or add hidden biological defaults.

- [ ] **Step 3: Update docs**

README shows:

```text
Target conservation -> Target vs non-target contrast -> Assay design -> Specificity
```

and links:
- `notebooks/geison_guided_colab.ipynb` for normal evaluators;
- `notebooks/geison_colab.ipynb` for low-level operational validation.

`docs/colab.md` keeps the existing notebook documented and points normal users to the guided workbench.

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

Expected: zero failures; only environment-dependent binary skips are acceptable outside CI.

- [ ] **Step 7: Installed CLI smoke**

```bash
qpcr-pipeline doctor
qpcr-pipeline panel approve --help
```

Then execute the synthetic guided workflow once through installed `qpcr-pipeline` commands.

- [ ] **Step 8: Commit**

```bash
git add integration_tests/test_guided_contrastive_demo.py README.md docs/colab.md
git commit -m "test: cover guided contrastive workflow"
```

---

## Final Verification Checklist

- [ ] `contrastive_conservation` is typed and ordered between conservation and primer design.
- [ ] `panel` and `input` are independent checkpoint roots after preflight; contrastive depends on both panel and conservation.
- [ ] Disabled/omitted contrast preserves conservation-only candidate behavior.
- [ ] Enabled contrast requires a frozen approved panel.
- [ ] Proposal config keeps contrast disabled; approved config enables it explicitly.
- [ ] CHALLENGE mapping and criticality are deterministic and preserved.
- [ ] Typed result retains challenge provenance, not only names.
- [ ] Per-dataset evidence exists in artifacts.
- [ ] No hidden contrast PASS cutoff controls selection.
- [ ] Overlapping windows consolidate deterministically.
- [ ] Primer design records `CONTRASTIVE_CONSERVATION` source when enabled.
- [ ] Final specificity remains independently evaluated.
- [ ] Panel-only challenge-semantic changes reuse target input through conservation.
- [ ] Challenge file changes reuse panel and target chain but invalidate contrastive descendants.
- [ ] Primer thermodynamic-only changes do not invalidate contrastive.
- [ ] Synthetic shared-conserved and discriminant regions are distinguished.
- [ ] Guided notebook contains no duplicated scientific implementation.
- [ ] Normal guided flow requires no manual YAML editing.
- [ ] Generated YAML/manifests/tool versions/raw artifacts remain inspectable.
- [ ] `ACTION_REQUIRED`, `PARTIAL`, and `FAILED` are never shown as completion.
- [ ] Full pytest regression passes.
- [ ] Full integration regression passes in CircleCI with installed tools.
- [ ] Guided demo remains synthetic and contains no real assay sequence or organism-specific cutoff.
- [ ] Contrast exploration is read-only; policy changes require explicit form/config change and rerun.

## Execution Notes

- Implement in an isolated branch/worktree, not directly on `main`.
- Use TDD for every task: RED -> minimal GREEN -> focused regression -> commit.
- Keep `notebooks/spikes/contrastive_conservation_synthetic.ipynb` as design evidence; do not turn it into production logic.
- Do not redesign final ranking semantics in this subproject.
- Do not add live biological panel construction or live challenge-dataset discovery.
- This v2 plan supersedes `docs/superpowers/plans/2026-09-03-contrastive-conservation-guided-colab.md`.
