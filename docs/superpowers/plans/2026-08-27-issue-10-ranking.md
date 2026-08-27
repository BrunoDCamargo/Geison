# Issue #10 Assay Ranking and Classification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an offline, auditable final assay classification and ranking stage that classifies before scoring, preserves every Primer3 assay, publishes structured reasons plus decomposable scores, and replaces the root report only when ranking is enabled.

**Architecture:** Configuration defines conservative thresholds and score weights. `qpcr_pipeline/ranking.py` first validates and classifies upstream evidence into immutable classified records, then a separate scoring pass calculates five absolute components and performs deterministic class-first ordering. `qpcr_pipeline/assay_report_html.py` renders the final static report, and the ranking stage atomically publishes TSV, JSON, and HTML without parsing upstream artifact files.

**Tech Stack:** Python 3.12, stdlib dataclasses/json/html/pathlib/math/uuid, existing PyYAML configuration, unittest-style tests collected by pytest, existing Geison result dataclasses.

**Spec:** `docs/superpowers/specs/2026-08-27-issue-10-ranking-design.md`

## Global Constraints

- Ranking is disabled by default.
- Default inclusivity thresholds are `min_inclusivity_for_pass: 1.0` and `min_inclusivity_before_high_risk: 0.90`.
- Default score weights are inclusivity `0.35`, specificity `0.25`, conservation `0.20`, Primer3 quality `0.10`, robustness `0.10`; weights must sum to `1.0` within absolute tolerance `1e-9`.
- Classification is finalized before score components are calculated; score never overrides class.
- Class order is `IN SILICO PASS`, `REVIEW`, `HIGH_RISK`.
- A detectable off-target always forces `HIGH_RISK` and cannot be weakened through configuration.
- Only original Primer3 assays are ranked. IUPAC proposals are contextual robustness evidence and never become independent assays.
- Missing evidence is never replaced by zero. Any unavailable component yields `score_status=INCOMPLETE` and `final_score=None`, even when its configured weight is zero.
- Ranking performs no network access and adds no third-party dependency.
- When ranking is disabled, it leaves root `report.html` untouched. When enabled, it owns root `report.html`.
- After valid enabled ranking configuration is known, stale ranking TSV/JSON/root HTML are removed before upstream evidence validation.
- Final text content is fully computed before publication and each file is atomically replaced through a temporary sibling.
- Normal validation uses `python -m pytest -q`; `.circleci/config.yml` is not changed for issue #10.

---

## File Structure

**Create**
- `qpcr_pipeline/ranking.py` — ranking models, integrity checks, reasons/classes, scoring, ordering, artifacts.
- `qpcr_pipeline/assay_report_html.py` — final static offline assay report.
- `tests/ranking_fixtures.py` — deterministic valid upstream result builders.
- `tests/test_ranking_config.py` — ranking YAML/direct config validation.
- `tests/test_ranking_classification.py` — evidence integrity, reason codes, class precedence.
- `tests/test_ranking_scoring.py` — exact five-component formulas and incomplete scores.
- `tests/test_ranking_ordering.py` — class-first and deterministic tie ordering.
- `tests/test_ranking_artifacts.py` — stale cleanup, skipped behavior, TSV/JSON/HTML publication.
- `tests/test_assay_report_html.py` — HTML content, escaping, offline guarantee, empty state.
- `tests/test_pipeline_ranking.py` — pipeline order and QC summary.

**Modify**
- `qpcr_pipeline/config.py` — ranking config dataclasses, parser, validator, pipeline dependency.
- `qpcr_pipeline/pipeline.py` — invoke ranking after specificity and summarize it.
- `README.md` — class/score semantics, artifacts, limitations, report ownership.

**Leave unchanged**
- `qpcr_pipeline/report_html.py` — conservation-specific renderer.
- `.circleci/config.yml` — existing normal `develop` job is sufficient.

---

### Task 1: Add the ranking configuration contract

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Create: `tests/test_ranking_config.py`

**Interfaces:**

```python
@dataclass(frozen=True, slots=True)
class RankingWeights:
    inclusivity: float = 0.35
    specificity: float = 0.25
    conservation: float = 0.20
    primer3_quality: float = 0.10
    robustness: float = 0.10


@dataclass(frozen=True, slots=True)
class RankingConfig:
    enabled: bool = False
    min_inclusivity_for_pass: float = 1.0
    min_inclusivity_before_high_risk: float = 0.90
    weights: RankingWeights = field(default_factory=RankingWeights)
```

Produces `validate_ranking_config(config: RankingConfig) -> None`. `PipelineConfig` gains `ranking: RankingConfig = field(default_factory=RankingConfig)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_ranking_config.py` using the same temporary-YAML pattern as `tests/test_specificity_config.py`. The tests must assert:

```python
self.assertEqual(config.ranking, RankingConfig())
self.assertEqual(config.ranking.weights, RankingWeights())
```

A fully specified YAML must parse to:

```python
RankingConfig(
    enabled=True,
    min_inclusivity_for_pass=0.98,
    min_inclusivity_before_high_risk=0.85,
    weights=RankingWeights(
        inclusivity=0.40,
        specificity=0.30,
        conservation=0.15,
        primer3_quality=0.10,
        robustness=0.05,
    ),
)
```

Use these concrete invalid YAML cases:

```text
ranking:
  enabled: nope

ranking:
  min_inclusivity_for_pass: 1.1

ranking:
  min_inclusivity_before_high_risk: -0.1

ranking:
  min_inclusivity_for_pass: 0.8
  min_inclusivity_before_high_risk: 0.9

ranking:
  weights:
    inclusivity: 0.50
    specificity: 0.50
    conservation: 0.50
    primer3_quality: 0.00
    robustness: 0.00

ranking:
  surprise: 1

ranking:
  weights:
    surprise: 1
```

Also construct a direct `PipelineConfig` with `ranking=RankingConfig(enabled=True)` and `primer_design.enabled=False`; `selected_input` must raise with a message containing `ranking` and `requires enabled primer design`.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_ranking_config.py -q
```

Expected: import/collection failure because the ranking config types do not exist.

- [ ] **Step 3: Implement dataclasses, parser, and validator**

In `config.py`, validate booleans strictly, thresholds as finite numbers in `[0, 1]`, lower threshold not above pass threshold, five weights as finite non-negative numbers, and:

```python
math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1e-9)
```

The ranking parser allows exactly `enabled`, `min_inclusivity_for_pass`, `min_inclusivity_before_high_risk`, and `weights`. The nested weights mapping allows exactly `inclusivity`, `specificity`, `conservation`, `primer3_quality`, and `robustness`; unspecified fields use dataclass defaults.

In `validate_pipeline_config()` add:

```python
validate_ranking_config(config.ranking)
if config.ranking.enabled and not config.primer_design.enabled:
    raise ValueError("Enabled ranking requires enabled primer design.")
```

- [ ] **Step 4: Verify GREEN and config regressions**

```bash
python -m pytest tests/test_ranking_config.py tests/test_specificity_config.py -q
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add qpcr_pipeline/config.py tests/test_ranking_config.py
git commit -m "feat: add ranking configuration"
```

---

### Task 2: Validate upstream evidence and classify assays before scoring

**Files:**
- Create: `qpcr_pipeline/ranking.py`
- Create: `tests/ranking_fixtures.py`
- Create: `tests/test_ranking_classification.py`

**Interfaces:**

```python
ReasonSeverity = Literal["HIGH_RISK", "REVIEW", "ADVISORY"]
AssayClassification = Literal["IN SILICO PASS", "REVIEW", "HIGH_RISK"]


class RankingError(RuntimeError):
    """Raised when ranking evidence is structurally untrustworthy."""


@dataclass(frozen=True, slots=True)
class RankingReason:
    code: str
    severity: ReasonSeverity
    source: str
    message: str
    evidence: tuple[tuple[str, object], ...]


@dataclass(frozen=True, slots=True)
class ClassifiedAssay:
    assay: AssayCandidate
    region: CandidateRegion
    classification: AssayClassification
    reasons: tuple[RankingReason, ...]
    original_compatible_count: int | None
    evaluation_sequence_count: int | None
    compatible_off_target_hit_count: int | None
    plausible_off_target_count: int | None
    detectable_off_target_count: int | None
    proposals: tuple[DegeneracyProposal, ...]
    inclusivity_available: bool
    specificity_available: bool
    missing_components: tuple[str, ...]
```

Produces `classify_assays(primer_design: PrimerDesignResult, inclusivity: InclusivityResult, specificity: SpecificityResult, config: RankingConfig) -> tuple[ClassifiedAssay, ...]`. This function validates evidence, creates reasons, and finalizes class without calculating a score.

- [ ] **Step 1: Add deterministic valid fixtures**

Create `tests/ranking_fixtures.py` with builders for:

```text
make_region(region_id="r1") -> CandidateRegion
make_oligo(sequence, start) -> DesignedOligo
make_assay(assay_id="a1", region_id="r1", primer3_index=0, pair_penalty=1.0) -> AssayCandidate
make_primer_result(assays=None, candidates=None, status="COMPLETE") -> PrimerDesignResult
make_inclusivity_result(primer, compatibility=None, sequence_ids=("s1",), proposals=(), status="COMPLETE") -> InclusivityResult
make_specificity_result(primer, dataset_names=("off",), hit_totals=None, amplicons=(), status="COMPLETE") -> SpecificityResult
```

`make_region()` uses perfect fractions, zero gaps, zero entropy. `make_specificity_result()` creates exactly one retention row for every configured dataset × assay × role in `FORWARD`, `PROBE`, `REVERSE`; `total_hit_count` comes from `hit_totals` and `retained_hit_count=min(total_hit_count, 20)`.

- [ ] **Step 2: Write RED classification tests**

In `tests/test_ranking_classification.py`, use ten Evaluation Set IDs to assert:

```text
10/10 compatible -> IN SILICO PASS
9/10 compatible  -> REVIEW + INCLUSIVITY_BELOW_PASS
8/10 compatible  -> HIGH_RISK + INCLUSIVITY_BELOW_MINIMUM
```

Add a configurable-threshold case with pass `0.95` and high-risk floor `0.80`.

Create a `PlausibleAmplicon` helper and assert:

```text
detectable_off_target=True -> HIGH_RISK + DETECTABLE_OFF_TARGET
primer_amplicon_plausible=True, detectable=False -> REVIEW + PLAUSIBLE_OFF_TARGET_AMPLICON
no amplicon + retention total hits > 0 -> advisory ISOLATED_OFF_TARGET_HITS only
```

With both inclusivity and specificity `SKIPPED`, assert exactly these codes after deduplication:

```python
{
    "EVIDENCE_INCOMPLETE",
    "INCLUSIVITY_EVIDENCE_MISSING",
    "SPECIFICITY_EVIDENCE_MISSING",
}
```

and class `REVIEW`. With skipped inclusivity plus a detectable off-target, class must remain `HIGH_RISK`.

Create one accepted and one rejected `DegeneracyProposal` for different roles. Assert the class stays PASS when other evidence passes and reasons contain one aggregated `IUPAC_PROPOSAL_ACCEPTED` and one aggregated `IUPAC_PROPOSAL_REJECTED` reason. Each reason's evidence must include the affected role names so deduplication never discards proposal context.

- [ ] **Step 3: Verify RED**

```bash
python -m pytest tests/test_ranking_classification.py -q
```

Expected: import failure because `qpcr_pipeline.ranking` does not exist.

- [ ] **Step 4: Implement complete evidence integrity checks**

Before classification, enforce:

```text
Primer design COMPLETE
unique assay_id
unique region_id
every assay.region_id resolves
```

For `InclusivityResult.status == COMPLETE`:

```text
unique evaluation_sequence_ids
exactly one assay_results row for every assay_id × evaluation_sequence_id pair
no missing or duplicate matrix pair
no unknown assay_id in assay_results
no unknown assay_id in oligo_matches
no unknown assay_id in variations
no unknown assay_id in proposals
proposal role is FORWARD, PROBE, or REVERSE
at most one proposal for each assay_id × role
```

For `SpecificityResult.status == COMPLETE`:

```text
assay_count equals Primer3 assay count
unique dataset_names
exactly one retention row for every dataset × assay × role
no missing or duplicate retention key
retention role is FORWARD, PROBE, or REVERSE
retention counts are non-negative integers and retained_hit_count <= total_hit_count
no unknown assay_id in retention
no unknown assay_id in hits
no unknown assay_id in amplicons
```

A malformed `COMPLETE` result raises `RankingError`; it never degrades to REVIEW.

- [ ] **Step 5: Generate reasons and finalize class**

Use severity priority:

```python
_SEVERITY_ORDER = {"HIGH_RISK": 0, "REVIEW": 1, "ADVISORY": 2}
```

Generate inclusivity, specificity, and proposal reasons from the approved rules. Create one aggregate `EVIDENCE_INCOMPLETE` reason with `source="ranking"` whose evidence contains a sorted tuple of unavailable component names. Raw prerequisites are known before scoring:

```text
inclusivity SKIPPED -> inclusivity, robustness unavailable
empty Evaluation Set -> inclusivity, robustness unavailable
specificity SKIPPED -> specificity unavailable
pair_penalty is None -> primer3_quality unavailable
```

Deduplicate reasons by `(code, source)`. For proposal reasons, aggregate all matching proposals into one reason before deduplication. Sort final reasons by severity priority, then `source`, then `code`.

Classify only from the reason set:

```python
if any(reason.severity == "HIGH_RISK" for reason in reasons):
    classification = "HIGH_RISK"
elif any(reason.severity == "REVIEW" for reason in reasons):
    classification = "REVIEW"
else:
    classification = "IN SILICO PASS"
```

No score function is called in this task.

- [ ] **Step 6: Add malformed-evidence and determinism regression tests**

Add concrete cases for duplicate assay ID, missing candidate region, missing/duplicate inclusivity matrix row, unknown assay in `oligo_matches`, unknown assay in `variations`, duplicate proposal `(assay_id, role)`, mismatched specificity `assay_count`, missing/duplicate retention row, unknown assay in specificity hits, and unknown assay in amplicons. Every case must raise `RankingError`.

Run equivalent valid inputs with proposal/evidence tuples in different input orders and assert the final ordered `(severity, source, code)` reason tuples are identical. This verifies deterministic reason ordering rather than relying on construction order.

- [ ] **Step 7: Verify GREEN**

```bash
python -m pytest tests/test_ranking_classification.py -q
```

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add qpcr_pipeline/ranking.py tests/ranking_fixtures.py tests/test_ranking_classification.py
git commit -m "feat: classify final assay evidence"
```

---

### Task 3: Calculate score components and apply deterministic class-first ranking

**Files:**
- Modify: `qpcr_pipeline/ranking.py`
- Create: `tests/test_ranking_scoring.py`
- Create: `tests/test_ranking_ordering.py`

**Interfaces:**

```python
ScoreStatus = Literal["COMPLETE", "INCOMPLETE"]


@dataclass(frozen=True, slots=True)
class ScoreComponents:
    inclusivity: float | None
    specificity: float | None
    conservation: float | None
    primer3_quality: float | None
    robustness: float | None


@dataclass(frozen=True, slots=True)
class RankedAssay:
    rank: int
    assay_id: str
    region_id: str
    primer3_index: int
    classification: AssayClassification
    final_score: float | None
    score_status: ScoreStatus
    components: ScoreComponents
    reasons: tuple[RankingReason, ...]
    original_compatible_count: int | None
    evaluation_sequence_count: int | None
    compatible_off_target_hit_count: int | None
    plausible_off_target_count: int | None
    detectable_off_target_count: int | None
    pair_penalty: float | None
```

Produces `rank_assays(primer_design: PrimerDesignResult, inclusivity: InclusivityResult, specificity: SpecificityResult, config: RankingConfig) -> tuple[RankedAssay, ...]`. Its first phase calls `classify_assays(...)`; component calculation and sorting happen only after classification returns.

- [ ] **Step 1: Write RED score-formula tests**

In `tests/test_ranking_scoring.py`, build a 9/10-compatible assay, two isolated off-target hits, perfect conservation, `pair_penalty=1.0`, and one accepted forward proposal with degeneracy `1 -> 2`. Assert:

```python
self.assertAlmostEqual(item.components.inclusivity, 0.9)
self.assertAlmostEqual(item.components.specificity, 0.96)
self.assertAlmostEqual(item.components.conservation, 1.0)
self.assertAlmostEqual(item.components.primer3_quality, 0.5)
self.assertAlmostEqual(item.components.robustness, (0.5 + 1.0 + 1.0) / 3.0)
self.assertAlmostEqual(item.final_score, 88.83333333333333)
self.assertEqual(item.score_status, "COMPLETE")
```

Separately assert specificity scores:

```text
0 hits, no amplicon -> 1.00
1 isolated hit -> 0.98
25 isolated hits -> 0.80 floor
plausible F/R, no detectable probe -> 0.40
detectable off-target -> 0.00
```

The `25` case must use `HitRetentionSummary.total_hit_count=25` even when `retained_hit_count=1`, proving detailed-hit truncation does not change scientific scoring.

- [ ] **Step 2: Write RED incomplete-score tests**

Assert `pair_penalty=None` gives:

```text
score_status == INCOMPLETE
final_score is None
classification == REVIEW
EVIDENCE_INCOMPLETE is present
```

Repeat with `primer3_quality` weight set to zero while the other weights still sum to one; the score remains incomplete. Use a complete inclusivity result with an empty Evaluation Set and assert inclusivity/robustness components are `None`, no division by zero occurs, and final score is incomplete.

- [ ] **Step 3: Verify RED**

```bash
python -m pytest tests/test_ranking_scoring.py -q
```

Expected: missing score models/functions or failed component assertions.

- [ ] **Step 4: Implement exact component formulas**

Inclusivity is `original_compatible_count / evaluation_sequence_count` when the denominator is nonzero, otherwise `None`.

Specificity is `None` when specificity is unavailable; otherwise use exactly `0.0` for any detectable off-target, `0.40` for any plausible amplicon without detectable probe, `max(0.80, 1.0 - 0.02 * compatible_hit_count)` for isolated hits, and `1.0` for zero compatible hits.

Conservation is:

```python
(
    item.region.mean_conservation
    + item.region.minimum_conservation
    + item.region.mean_coverage
    + (1.0 - item.region.mean_gap_frequency)
    + (1.0 - min(item.region.mean_entropy_bits / 2.0, 1.0))
) / 5.0
```

Primer3 quality is `None` for missing pair penalty, otherwise `1.0 / (1.0 + pair_penalty)`.

Robustness uses `1.0` for every role without an accepted proposal and `original_degeneracy / proposed_degeneracy` for an accepted proposal; take the arithmetic mean across forward/probe/reverse. If inclusivity is unavailable or its Evaluation Set is empty, robustness is `None`.

Validate present candidate fractions as finite values in `[0, 1]`, entropy as finite and non-negative, pair penalty as finite and non-negative, and accepted degeneracies as positive with proposed degeneracy not lower than original. Invalid present values raise `RankingError`.

- [ ] **Step 5: Compute final score only when all five components exist**

When any component is `None`, set `score_status="INCOMPLETE"` and `final_score=None`. Otherwise calculate full-precision:

```python
100.0 * (
    config.weights.inclusivity * components.inclusivity
    + config.weights.specificity * components.specificity
    + config.weights.conservation * components.conservation
    + config.weights.primer3_quality * components.primer3_quality
    + config.weights.robustness * components.robustness
)
```

Do not round model values.

- [ ] **Step 6: Verify scoring GREEN**

```bash
python -m pytest tests/test_ranking_scoring.py -q
```

Expected: PASS.

- [ ] **Step 7: Write RED deterministic ordering tests**

In `tests/test_ranking_ordering.py`, create at least one assay from each class and make the HIGH_RISK assay's numeric inputs stronger than the REVIEW assay. Assert final ordering remains PASS, REVIEW, HIGH_RISK. Add same-class tests proving COMPLETE score sorts before INCOMPLETE, then final score descending, inclusivity descending, pair penalty ascending, Primer3 index ascending, and assay ID ascending. Assert ranks are contiguous `1..N` and every Primer3 assay remains present.

Also rank a fixed assay alone and beside an additional assay; its `final_score` must be identical in both runs, proving there is no relative normalization.

- [ ] **Step 8: Implement sort key and contiguous rank assignment**

Use class priorities `PASS=0`, `REVIEW=1`, `HIGH_RISK=2`. Missing optional numeric values sort after present values. Sort by class, score status, final score descending, inclusivity descending, pair penalty ascending, Primer3 index ascending, assay ID ascending. Assign `rank` afterward with `dataclasses.replace` and 1-based contiguous indexes.

- [ ] **Step 9: Verify all core ranking tests**

```bash
python -m pytest \
  tests/test_ranking_classification.py \
  tests/test_ranking_scoring.py \
  tests/test_ranking_ordering.py -q
```

Expected: PASS.

- [ ] **Step 10: Commit**

```bash
git add qpcr_pipeline/ranking.py tests/test_ranking_scoring.py tests/test_ranking_ordering.py
git commit -m "feat: score and order final assays"
```

---

### Task 4: Render and atomically publish final ranking artifacts

**Files:**
- Create: `qpcr_pipeline/assay_report_html.py`
- Modify: `qpcr_pipeline/ranking.py`
- Create: `tests/test_assay_report_html.py`
- Create: `tests/test_ranking_artifacts.py`

**Interfaces:**

Produces `render_assay_report_html(target_name: str, primer_design: PrimerDesignResult, inclusivity: InclusivityResult, specificity: SpecificityResult, assays: tuple[RankedAssay, ...]) -> str`.

```python
@dataclass(frozen=True, slots=True)
class RankingResult:
    status: Literal["SKIPPED", "COMPLETE"]
    assays: tuple[RankedAssay, ...]
    ranking_tsv_path: Path | None
    ranking_report_path: Path
    html_report_path: Path | None
```

Produces `evaluate_ranking(primer_design: PrimerDesignResult, inclusivity: InclusivityResult, specificity: SpecificityResult, config: RankingConfig, output_dir: Path, *, target_name: str) -> RankingResult`.

- [ ] **Step 1: Write RED HTML tests**

Create `tests/test_assay_report_html.py` with one ranked assay and an accepted IUPAC proposal. Assert the returned HTML includes class, rank, total score, all five named components, original F/Probe/R sequences, product size, pair penalty, candidate-region conservation metrics, original inclusivity count/fraction, proposal sequence/status as contextual evidence, and specificity counts by dataset.

Use `target_name='<img src=x onerror="boom">'`; assert the literal tag is absent and `&lt;img` is present. Assert the document contains no `http://`, no `https://`, and no `<script`. With zero assays, assert visible text `No assay candidates`.

- [ ] **Step 2: Verify HTML RED**

```bash
python -m pytest tests/test_assay_report_html.py -q
```

Expected: import failure because the renderer does not exist.

- [ ] **Step 3: Implement a static renderer with escaping**

Use stdlib `html.escape` for every dynamic string and a static CSS-only document. The summary table columns are rank, assay, class, score, reasons. One `<details>` block per assay contains four panels: assay oligos/design, conservation, inclusivity/degeneracy, specificity/score. Dynamic rows follow final assay order; oligo roles always use `FORWARD`, `PROBE`, `REVERSE`; proposal rows sort by role. Do not use a script tag, CDN, remote font, remote stylesheet, or network fetch.

- [ ] **Step 4: Verify HTML GREEN**

```bash
python -m pytest tests/test_assay_report_html.py -q
```

Expected: PASS.

- [ ] **Step 5: Write RED stage-publication tests**

Create `tests/test_ranking_artifacts.py` and assert disabled behavior with deliberately invalid `object()` upstream values. Before the call, create stale `ranking/assay_ranking.tsv` and root `report.html` containing `conservation report`. After calling disabled ranking, assert:

```text
status == SKIPPED
stale ranking TSV removed
ranking/ranking_report.json exists with SKIPPED
root report.html still contains conservation report
upstream objects were never inspected
```

For enabled failure, pre-create stale TSV/JSON/root HTML, pass `PrimerDesignResult.status=SKIPPED`, assert `RankingError`, then assert all three stale ranking-owned artifacts were removed.

For enabled success, assert exact paths:

```text
ranking/assay_ranking.tsv
ranking/ranking_report.json
report.html
```

TSV must contain rank, assay ID, region ID, class, score status, final score, all five components, inclusivity counts/fraction, compatible-hit count, plausible/detectable counts, pair penalty, and semicolon-separated reason codes.

JSON must contain `schema_version: 1`, effective config, class counts, score-state counts, ordered assays, full-precision components, structured reasons/evidence, and relative artifact paths.

A zero-assay complete primer result must publish header-only TSV, zero-count JSON, and empty-state HTML.

- [ ] **Step 6: Verify publication RED**

```bash
python -m pytest tests/test_ranking_artifacts.py -q
```

Expected: `RankingResult` or `evaluate_ranking` missing.

- [ ] **Step 7: Implement stage publication order and atomic writes**

Artifact paths are `ranking/assay_ranking.tsv`, `ranking/ranking_report.json`, and root `report.html`. Validate `RankingConfig` before any cleanup.

Disabled path: remove ranking TSV only, leave root HTML untouched, atomically replace ranking JSON with SKIPPED report, return no TSV/HTML path.

Enabled path: remove stale TSV/ranking JSON/root HTML, rank assays, render all three final text payloads in memory, atomically write TSV, atomically write HTML, atomically write JSON report last, then return COMPLETE result.

Use a UUID-named sibling temporary file opened with mode `x`, UTF-8, newline `\n`, then `Path.replace()`. Always remove a leftover temporary file in `finally`.

- [ ] **Step 8: Add determinism tests**

Run identical complete ranking inputs into two clean output directories and assert TSV and JSON bytes are identical. Ensure JSON stores relative artifact paths so absolute temporary directory names do not enter output. Assert scoring with a truncated retained-hit list still uses retention `total_hit_count`.

- [ ] **Step 9: Verify artifact + renderer suites**

```bash
python -m pytest \
  tests/test_assay_report_html.py \
  tests/test_ranking_artifacts.py \
  tests/test_ranking_classification.py \
  tests/test_ranking_scoring.py \
  tests/test_ranking_ordering.py -q
```

Expected: PASS.

- [ ] **Step 10: Commit**

```bash
git add qpcr_pipeline/assay_report_html.py qpcr_pipeline/ranking.py \
  tests/test_assay_report_html.py tests/test_ranking_artifacts.py
git commit -m "feat: publish final assay ranking"
```

---

### Task 5: Integrate ranking into the pipeline and document user-facing semantics

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Create: `tests/test_pipeline_ranking.py`
- Modify: `README.md`

**Interfaces:**

`run_pipeline()` calls `evaluate_ranking()` after specificity and before final QC/run-summary publication. `qc_report.json["ranking"]` contains exactly:

```text
status
assay_count
in_silico_pass_count
review_count
high_risk_count
complete_score_count
incomplete_score_count
top_recommended_assay_id
```

- [ ] **Step 1: Write RED pipeline tests**

Default pipeline config must produce:

```python
{
    "status": "SKIPPED",
    "assay_count": 0,
    "in_silico_pass_count": 0,
    "review_count": 0,
    "high_risk_count": 0,
    "complete_score_count": 0,
    "incomplete_score_count": 0,
    "top_recommended_assay_id": None,
}
```

For enabled ranking, patch clustering/alignment/conservation/primer design/inclusivity/specificity to return valid public result objects from `tests/ranking_fixtures.py`; leave ranking itself unpatched. Assert one PASS assay produces `top_recommended_assay_id == "a1"` and the three ranking artifacts exist. Add a REVIEW-only case and assert `top_recommended_assay_id is None`.

Use mock side effects or call assertions to prove `evaluate_specificity` completes before ranking is invoked.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_pipeline_ranking.py -q
```

Expected: ranking summary missing or ranking not invoked.

- [ ] **Step 3: Integrate ranking and QC summary**

Import `evaluate_ranking`, call it with primer design, inclusivity, specificity, ranking config, output directory, and target name. Add counts from `ranking.assays`. `top_recommended_assay_id` is the first final ordered assay whose classification is `IN SILICO PASS`; return `None` when there is no PASS.

- [ ] **Step 4: Verify pipeline and prior specificity regressions**

```bash
python -m pytest \
  tests/test_pipeline_ranking.py \
  tests/test_pipeline_specificity.py \
  tests/test_specificity_artifacts.py -q
```

Expected: PASS.

- [ ] **Step 5: Update README**

Add `## Classificação e ranking final dos assays` after specificity with the exact default YAML, inclusivity boundaries, specificity class rules, five score components and weights, class-before-score ordering, incomplete-score behavior, reason-code auditability, and artifact paths.

State explicitly that accepted IUPAC proposals remain contextual and are not independently re-tested for specificity. State that `IN SILICO PASS` is not experimental validation.

Add this report-ownership rule in Portuguese:

```text
Com ranking desabilitado, o estágio não altera report.html; portanto o relatório
publicado pela conservação continua disponível. Com ranking habilitado, o ranking
assume a propriedade de report.html e substitui o relatório de conservação pelo
relatório final consolidado dos assays.
```

Also adjust the conservation section to mention that enabled ranking later replaces root `report.html` in the same run.

- [ ] **Step 6: Commit**

```bash
git add qpcr_pipeline/pipeline.py tests/test_pipeline_ranking.py README.md
git commit -m "feat: integrate final assay ranking"
```

---

### Task 6: Final regression, review, and completion evidence

**Files:**
- Review every file changed by Tasks 1-5.
- Keep `.circleci/config.yml` unchanged.

**Interfaces:**
- Produces a reviewed feature branch ready for `superpowers:verification-before-completion` and `superpowers:finishing-a-development-branch`.

- [ ] **Step 1: Run the full normal suite**

```bash
python -m pytest -q
```

Expected: PASS.

- [ ] **Step 2: Run diff sanity checks**

```bash
git diff --check develop...HEAD
git diff --stat develop...HEAD
git status --short
```

Expected: `git diff --check` exits zero and there are no uncommitted issue-#10 implementation files.

- [ ] **Step 3: Perform a direct diff review against the issue and spec**

Inspect all issue-#10 changes and explicitly verify:

```text
classification occurs before scoring
IN SILICO PASS / REVIEW / HIGH_RISK rules match approved thresholds
HIGH_RISK cannot be outranked by score
five components are named and decomposable
all Primer3 assays remain in output
reason codes are structured and deterministic
COMPLETE evidence matrices fail on inconsistency
SKIPPED evidence yields REVIEW rather than fake negative evidence
root final report contains F/Probe/R, degeneracy, inclusivity, specificity, class and ranking
disabled ranking preserves a conservation report.html
no network or new dependency was introduced
```

No independent subagent reviewer is available in this environment, so record this as a direct diff review and do not claim an independent review.

- [ ] **Step 4: Re-run focused issue-#10 tests after review**

```bash
python -m pytest \
  tests/test_ranking_config.py \
  tests/test_ranking_classification.py \
  tests/test_ranking_scoring.py \
  tests/test_ranking_ordering.py \
  tests/test_ranking_artifacts.py \
  tests/test_assay_report_html.py \
  tests/test_pipeline_ranking.py -q
```

Expected: PASS from a fresh post-review run.

- [ ] **Step 5: Fix any review defect through RED → GREEN only**

For every discovered defect, first add one regression test that fails for the defect, run that focused test to confirm RED, implement the minimum fix, rerun to GREEN, then rerun the complete issue-#10 focused suite. If no defect is found, create no review-only commit.

- [ ] **Step 6: Finish the branch with verification evidence**

Invoke `superpowers:verification-before-completion`, then `superpowers:finishing-a-development-branch`. Issue #10 is not complete until the merged `develop` SHA receives a successful normal CircleCI run. After that run is green, mark the issue acceptance checkboxes complete, comment with the merged SHA and verification result, and close issue #10 as completed.
