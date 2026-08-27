# Issue #10 Assay Ranking and Classification Design

## Status and goal

Approved in conversation on 2026-08-27. This document defines issue #10, **Classificar, ranquear e explicar os assays finais**.

The stage transforms Primer3 assays plus conservation, inclusivity, degeneracy, and specificity evidence into an explainable final classification and deterministic ranking. Classification is evaluated before score, so a `HIGH_RISK` assay can never outrank `REVIEW` or `IN SILICO PASS` because of a high numeric score.

All Primer3 assays remain visible. `HIGH_RISK` means rejected for recommendation, not deleted.

## Scope and scientific boundaries

Issue #10 adds:

- `IN SILICO PASS`, `REVIEW`, and `HIGH_RISK` classes;
- structured reason codes with source and evidence;
- an absolute, decomposable 0-100 score;
- deterministic class-first ordering;
- `ranking/assay_ranking.tsv` and `ranking/ranking_report.json`;
- a root `report.html` consolidating assay design, conservation, inclusivity, degeneracy, specificity, classification, reasons, and score;
- a compact ranking summary in `qc_report.json`.

It does not:

- turn an accepted IUPAC proposal into a new assay;
- recompute specificity for proposed degenerate oligos;
- allow score to override class;
- hide rejected candidates;
- perform network access;
- add external report dependencies.

Only the original Primer3 assay is ranked. IUPAC proposals from #8 are contextual robustness evidence and never silently replace the original oligos.

## Architecture

Ranking is a new stage after specificity:

```text
Primer design / candidate regions
          |
          +------------------+
          v                  v
     Inclusivity        Specificity
          |                  |
          +--------+---------+
                   v
                Ranking
                   |
        reasons -> class -> score -> deterministic order
                   |
               TSV / JSON / HTML
```

The stage consumes public result objects and joins them by `assay_id`; it does not parse upstream artifact files.

Recommended files:

- `qpcr_pipeline/ranking.py`: typed models, evidence aggregation, reasons, classification, scoring, ordering, artifacts;
- `qpcr_pipeline/assay_report_html.py`: final self-contained assay report;
- `qpcr_pipeline/config.py`: ranking configuration;
- `qpcr_pipeline/pipeline.py`: orchestration and QC summary.

The existing `qpcr_pipeline/report_html.py` remains specific to conservation.

Processing order is fixed:

```text
1. collect and validate evidence
2. generate reason codes
3. determine classification
4. calculate five score components
5. calculate final score
6. order by class, score state, score, and deterministic tie breakers
7. publish artifacts
```

## Configuration

Ranking is disabled by default for backward compatibility.

```yaml
ranking:
  enabled: false
  min_inclusivity_for_pass: 1.0
  min_inclusivity_before_high_risk: 0.90
  weights:
    inclusivity: 0.35
    specificity: 0.25
    conservation: 0.20
    primer3_quality: 0.10
    robustness: 0.10
```

Validation:

- `enabled` is boolean;
- thresholds are finite numbers in `[0, 1]`;
- `min_inclusivity_before_high_risk <= min_inclusivity_for_pass`;
- weights are finite and non-negative;
- weights sum to `1.0` within numerical tolerance;
- enabled ranking requires enabled primer design;
- inclusivity and specificity are not enabling prerequisites because `SKIPPED` evidence is allowed.

The MVP makes inclusivity thresholds configurable. Specificity severity is a minimum safety contract: detectable off-targets cannot be configured into PASS. Future options may make specificity handling more conservative, not weaker.

## Public contract

Use immutable typed models consistent with existing stages.

```text
RankingReason
- code
- severity: HIGH_RISK | REVIEW | ADVISORY
- source
- message
- evidence

ScoreComponents
- inclusivity: float | None
- specificity: float | None
- conservation: float | None
- primer3_quality: float | None
- robustness: float | None

RankedAssay
- rank
- assay_id
- classification: IN SILICO PASS | REVIEW | HIGH_RISK
- final_score: float | None
- score_status: COMPLETE | INCOMPLETE
- components
- reasons

RankingResult
- status: SKIPPED | COMPLETE
- assays
- ranking_tsv_path
- ranking_report_path
- html_report_path
```

## Evidence integrity

For every Primer3 assay, ranking resolves the candidate region, inclusivity rows/proposals, specificity retention totals, and plausible amplicons.

The following are internal inconsistencies and fail explicitly rather than becoming `REVIEW`:

- duplicate Primer3 `assay_id`;
- duplicate candidate `region_id`;
- assay references a missing candidate region;
- a `COMPLETE` inclusivity result references an unknown assay;
- a `COMPLETE` specificity result references an unknown assay;
- duplicate degeneracy proposal for the same `(assay_id, role)`;
- non-finite or structurally invalid ranking metrics.

A `COMPLETE` inclusivity result must contain exactly one assay evaluation for every `(assay_id, evaluation_sequence_id)` pair. Missing or duplicate rows are invalid upstream evidence. If its Evaluation Set is empty, the inclusivity fraction is undefined rather than zero.

A `COMPLETE` specificity result must be structurally consistent with the Primer3 assays. Its `assay_count` must equal the Primer3 assay count, and retention summaries must contain exactly one row for every `(dataset, assay_id, role)` for `FORWARD`, `PROBE`, and `REVERSE`. Unknown, missing, or duplicate retention rows fail explicitly. Amplicon rows may legitimately be empty.

`SKIPPED` is different from inconsistent `COMPLETE`: skipped evidence is allowed and produces an incomplete score.

## Reasons and classification

Reason severities:

```text
HIGH_RISK -> forces HIGH_RISK
REVIEW    -> prevents IN SILICO PASS
ADVISORY  -> explanation only
```

Initial stable codes:

```text
HIGH_RISK
- INCLUSIVITY_BELOW_MINIMUM
- DETECTABLE_OFF_TARGET

REVIEW
- INCLUSIVITY_BELOW_PASS
- PLAUSIBLE_OFF_TARGET_AMPLICON
- EVIDENCE_INCOMPLETE
- INCLUSIVITY_EVIDENCE_MISSING
- SPECIFICITY_EVIDENCE_MISSING

ADVISORY
- ISOLATED_OFF_TARGET_HITS
- IUPAC_PROPOSAL_ACCEPTED
- IUPAC_PROPOSAL_REJECTED
```

Each reason carries `code`, `severity`, `source`, `message`, and structured `evidence`.

Reasons are deduplicated by stable code/source identity and serialized in deterministic order: severity priority (`HIGH_RISK`, `REVIEW`, `ADVISORY`), then source, then code. This prevents duplicate `EVIDENCE_INCOMPLETE` entries or output-order drift.

Classification:

```text
if any HIGH_RISK reason:
    HIGH_RISK
else if any REVIEW reason:
    REVIEW
else:
    IN SILICO PASS
```

### Inclusivity

Using only original assay compatibility:

```text
fraction >= min_inclusivity_for_pass
    -> no downgrade

min_inclusivity_before_high_risk <= fraction < min_inclusivity_for_pass
    -> REVIEW / INCLUSIVITY_BELOW_PASS

fraction < min_inclusivity_before_high_risk
    -> HIGH_RISK / INCLUSIVITY_BELOW_MINIMUM
```

Defaults: PASS requires `1.0`; below `0.90` is HIGH_RISK.

### Specificity

Across all configured off-target datasets:

```text
any detectable_off_target
    -> HIGH_RISK / DETECTABLE_OFF_TARGET

else any plausible F/R amplicon without compatible probe
    -> REVIEW / PLAUSIBLE_OFF_TARGET_AMPLICON

else compatible oligo hits but no plausible amplicon
    -> ADVISORY / ISOLATED_OFF_TARGET_HITS

else
    -> no specificity reason
```

Use one highest-risk specificity classification reason per assay. Detailed counts by dataset remain in evidence and reports.

### Degeneracy proposals

- `ACCEPTED` -> `ADVISORY / IUPAC_PROPOSAL_ACCEPTED`;
- `REJECTED` -> `ADVISORY / IUPAC_PROPOSAL_REJECTED`;
- `UNCHANGED` -> no reason required.

A proposal alone never changes class.

### Missing evidence

If inclusivity or specificity is `SKIPPED`, add `EVIDENCE_INCOMPLETE` plus the corresponding source-specific missing-evidence code. This forces at least `REVIEW` but never creates `HIGH_RISK` by itself.

If another available source has HIGH_RISK evidence, HIGH_RISK still wins.

## Score model

The score is absolute and deterministic, not normalized against the other assays in the run.

```text
final_score = 100 * (
    weight_inclusivity     * inclusivity
  + weight_specificity     * specificity
  + weight_conservation    * conservation
  + weight_primer3_quality * primer3_quality
  + weight_robustness      * robustness
)
```

Default weights are `0.35 / 0.25 / 0.20 / 0.10 / 0.10`.

Components are calculated with full precision in `[0, 1]`; presentation rounds to two decimals only. Score never influences class.

### Inclusivity component

```text
inclusivity = original_compatible_count / evaluation_sequence_count
```

If inclusivity is `SKIPPED` or the Evaluation Set is empty, the component is `None`.

### Specificity component

Use complete #9 evidence, never the truncated detailed hit list. When there is no plausible amplicon, compatible hit count is the sum of `HitRetentionSummary.total_hit_count` for the assay.

```text
1.00                              -> zero compatible off-target hits
max(0.80, 1 - 0.02 * hit_count)  -> hits but no plausible amplicon
0.40                              -> plausible F/R amplicon, no detectable probe
0.00                              -> any detectable_off_target
```

If specificity is `SKIPPED`, the component is `None`.

### Conservation component

Resolve the assay's candidate region:

```text
conservation = mean(
    mean_conservation,
    minimum_conservation,
    mean_coverage,
    1 - mean_gap_frequency,
    1 - min(mean_entropy_bits / 2, 1)
)
```

Two bits is the theoretical maximum Shannon entropy for four equiprobable canonical bases. Fraction metrics must be finite and in `[0, 1]`; entropy must be finite and non-negative. Invalid upstream values fail instead of being silently clamped, except for the explicit entropy normalization cap.

### Primer3 quality component

```text
primer3_quality = 1 / (1 + pair_penalty)
```

A present `pair_penalty` must be finite and non-negative. `None` makes this component incomplete.

### Robustness component

With complete inclusivity evidence:

```text
no accepted proposal for role -> role_robustness = 1.0
accepted proposal             -> original_degeneracy / proposed_degeneracy
robustness = mean(FORWARD, PROBE, REVERSE)
```

Accepted proposal degeneracies must be positive and `proposed_degeneracy >= original_degeneracy`.

If inclusivity is `SKIPPED` or its Evaluation Set is empty, robustness is `None`.

## Incomplete score behavior

Missing evidence is never replaced by zero.

If any of the five components is unavailable:

```text
score_status = INCOMPLETE
final_score = null
```

This remains true even if the configured weight of the missing component is zero. `score_status=COMPLETE` means all five evidence components were actually computable; weights cannot hide missing evidence.

The assay receives `REVIEW / EVIDENCE_INCOMPLETE` unless an available HIGH_RISK reason already wins classification.

## Deterministic ordering

Sort by:

```text
1. class: IN SILICO PASS, REVIEW, HIGH_RISK
2. score status: COMPLETE before INCOMPLETE
3. final_score descending when present
4. inclusivity component descending when present
5. pair_penalty ascending when present
6. primer3_index ascending
7. assay_id ascending
```

Missing numeric values sort after present values. Ranks are contiguous from 1 after sorting.

Therefore a `HIGH_RISK 99.8` remains behind every REVIEW/PASS assay.

## Artifacts

Enabled and complete ranking publishes:

```text
ranking/
├── assay_ranking.tsv
└── ranking_report.json

report.html
```

### `assay_ranking.tsv`

One row per Primer3 assay with at least:

- rank, assay_id, region_id, class, score status, final score;
- all five components;
- inclusivity count/fraction when available;
- compatible off-target hit count;
- plausible and detectable off-target counts;
- pair penalty;
- semicolon-separated reason codes.

TSV is presentation/audit output, not the source of truth for downstream computation.

### `ranking_report.json`

Contains schema version, effective config, class and score-state counts, ordered assays, full score components, structured reasons/evidence, and artifact paths.

### `report.html`

The root report is deterministic, self-contained, and offline. It uses safe HTML escaping for all dynamic values and contains no CDN, remote fonts, external scripts, or network fetches.

For each assay it shows:

- rank, class, score, reason codes;
- forward/probe/reverse sequences, coordinates, Tm, GC, and available penalties;
- product size and pair penalty;
- region conservation metrics;
- original inclusivity count/fraction;
- IUPAC proposals as contextual evidence only;
- specificity hit, plausible amplicon, and detectable counts by dataset;
- all five score components.

A static ordered table plus `<details>` sections is sufficient; no charting library is required.

## Disabled, empty, and failure behavior

### Disabled

`ranking.enabled: false`:

- returns `SKIPPED`;
- does not inspect upstream scientific result integrity;
- removes stale `ranking/assay_ranking.tsv` and root `report.html`;
- publishes only `ranking/ranking_report.json` with `SKIPPED` status.

### Enabled execution

After ranking config itself is validated, remove stale `ranking/assay_ranking.tsv`, `ranking/ranking_report.json`, and root `report.html` before validating/aggregating enabled upstream evidence. Therefore a failed current ranking run cannot leave a previous successful report looking current.

All final artifact contents are computed and validated before publication. Each text file is written through a temporary file and atomically replaced, following existing project patterns.

### Zero assays

Enabled ranking with `PrimerDesignResult.status == COMPLETE` and zero assays returns `COMPLETE`, publishes a header-only TSV, zero-assay JSON, and empty-state HTML. No evidence rows are required for nonexistent assays.

### Invalid primer design

Enabled ranking with a non-`COMPLETE` PrimerDesignResult fails explicitly.

## Pipeline integration

`run_pipeline()` calls ranking after specificity and before final summary publication.

`qc_report.json` receives only:

```text
ranking.status
ranking.assay_count
ranking.in_silico_pass_count
ranking.review_count
ranking.high_risk_count
ranking.complete_score_count
ranking.incomplete_score_count
ranking.top_recommended_assay_id
```

`top_recommended_assay_id` is the first `IN SILICO PASS` assay. If there is no PASS, it is `null`; REVIEW is not silently promoted to recommended.

## Test strategy

Implementation follows RED -> GREEN -> REFACTOR TDD.

Suggested files:

```text
tests/
├── test_ranking_config.py
├── test_ranking_classification.py
├── test_ranking_scoring.py
├── test_ranking_ordering.py
├── test_ranking_artifacts.py
├── test_assay_report_html.py
└── test_pipeline_ranking.py
```

Required coverage:

- default disabled behavior and stale cleanup;
- config thresholds and weight validation;
- inclusivity at 100%, 90-<100%, and below 90%;
- configurable inclusivity thresholds;
- detectable off-target HIGH_RISK;
- plausible F/R without probe REVIEW;
- isolated-hit advisory and score penalty;
- specificity scoring from full retention totals despite detailed-hit truncation;
- accepted/rejected IUPAC proposal advisories and robustness scoring;
- skipped evidence -> REVIEW/EVIDENCE_INCOMPLETE;
- HIGH_RISK evidence wins even with another missing source;
- empty Evaluation Set without division by zero;
- missing pair penalty -> incomplete score;
- a zero-weight missing component still yields incomplete score;
- all five component formulas;
- no relative normalization between assays;
- class-first ordering and deterministic tie breaks;
- HIGH_RISK never outranks REVIEW/PASS by score;
- complete score before incomplete score within a class;
- all Primer3 assays preserved;
- reason deduplication and deterministic ordering;
- unknown/missing/duplicate COMPLETE evidence fails explicitly;
- specificity retention matrix integrity;
- zero-assay result;
- deterministic TSV/JSON;
- HTML escaping and no remote resources;
- QC summary and top PASS assay;
- full pipeline path `primer_design -> inclusivity -> specificity -> ranking`.

Everything runs offline in normal `pytest` on `develop`. No new CircleCI dependency, integration executable, cache, schedule, parallelism, or network test is required.

## Completion criteria

Issue #10 closes only when:

1. acceptance criteria are implemented and tested;
2. classification is demonstrably computed before score ordering;
3. every Primer3 assay remains auditable;
4. score components and structured reasons appear in JSON and final HTML;
5. README documents class/score semantics and limitations;
6. the normal test suite passes;
7. the merged `develop` commit is green in CircleCI;
8. review has no unresolved critical or important finding.
