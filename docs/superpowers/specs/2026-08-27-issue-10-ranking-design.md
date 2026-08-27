# Issue #10 Assay Ranking and Classification Design

## Status

Approved in conversation on 2026-08-27. This document defines the design for issue #10, "Classificar, ranquear e explicar os assays finais".

## Goal

Transform Primer3 assay candidates plus inclusivity, degeneracy, conservation, and specificity evidence into a transparent final classification and deterministic ranking while preserving every candidate, including those not recommended.

The ranking stage must explain why an assay received its class and score. Classification always happens before quantitative ranking, so a numerically strong HIGH_RISK assay can never outrank REVIEW or IN SILICO PASS.

## Scope

Issue #10 adds:

- explicit assay classes: `IN SILICO PASS`, `REVIEW`, and `HIGH_RISK`;
- structured reason codes with evidence;
- a decomposable 0-100 score used only within the same class;
- deterministic ordering;
- a final offline `report.html` consolidating assay design, conservation, inclusivity, degeneracy, specificity, classification, and score;
- auditable TSV and JSON artifacts;
- a small ranking summary in `qc_report.json`.

Issue #10 does not:

- promote an accepted IUPAC proposal into a new assay candidate;
- recompute specificity for a proposed degenerate oligo set;
- allow quantitative score to override class;
- hide or delete rejected/high-risk candidates;
- perform network access;
- add external reporting or visualization dependencies.

## Scientific contract

### Original assay remains the ranked unit

Only assays produced by Primer3 are ranked. An accepted IUPAC proposal from inclusivity is contextual evidence about robustness. It never silently replaces the original oligo sequence and never becomes a separate ranked candidate in this issue.

### Evidence absence is not negative evidence

If inclusivity or specificity was `SKIPPED`, ranking still runs. The affected assay receives `REVIEW` with `EVIDENCE_INCOMPLETE` plus a source-specific reason. Missing evidence is never interpreted as a biological failure and never produces `HIGH_RISK` by itself.

If another available source provides real HIGH_RISK evidence, HIGH_RISK still wins even when some other evidence is missing.

### High-risk candidates remain visible

`HIGH_RISK` means rejected for recommendation, not deleted. Every Primer3 assay remains in the TSV, JSON, and HTML report with its reasons and score state.

## Architecture

Add a dedicated ranking stage after specificity:

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
        +----------+----------+
        v          v          v
 classification   score   reason codes
        |          |          |
        +----------+----------+
                   v
          deterministic order
                   |
          TSV + JSON + HTML
```

The stage consumes public result objects and joins evidence by `assay_id`. It does not parse upstream TSV/JSON artifacts.

Recommended modules:

- `qpcr_pipeline/ranking.py`: public ranking models, evidence aggregation, classification, scoring, ordering, artifact publication;
- `qpcr_pipeline/assay_report_html.py`: deterministic, self-contained HTML rendering for the final assay report;
- `qpcr_pipeline/config.py`: `RankingConfig` and `RankingWeights`;
- `qpcr_pipeline/pipeline.py`: stage orchestration and `qc_report.json` summary.

The existing `qpcr_pipeline/report_html.py` remains specific to conservation and is not generalized.

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

Validation rules:

- `enabled` is boolean;
- inclusivity thresholds are finite numbers in `[0, 1]`;
- `min_inclusivity_before_high_risk <= min_inclusivity_for_pass`;
- all weights are finite and non-negative;
- weights sum to `1.0` within a small numerical tolerance;
- ranking enabled requires primer design enabled;
- inclusivity and specificity are optional dependencies for execution, not enabling prerequisites.

The configurable classification behavior in the MVP is the inclusivity threshold pair. The specificity safety mapping below is an explicit minimum safety contract and cannot be weakened by configuration in this issue. Future configuration may make these rules more conservative but must not allow a detectable off-target to become PASS.

## Public data contract

Use typed, immutable models analogous to existing pipeline stages.

```text
RankingReason
- code: str
- severity: HIGH_RISK | REVIEW | ADVISORY
- source: str
- message: str
- evidence: structured mapping

ScoreComponents
- inclusivity: float | None
- specificity: float | None
- conservation: float | None
- primer3_quality: float | None
- robustness: float | None

RankedAssay
- rank: int
- assay_id: str
- classification: IN SILICO PASS | REVIEW | HIGH_RISK
- final_score: float | None
- score_status: COMPLETE | INCOMPLETE
- components: ScoreComponents
- reasons: tuple[RankingReason, ...]

RankingResult
- status: SKIPPED | COMPLETE
- assays: tuple[RankedAssay, ...]
- ranking_tsv_path: Path | None
- ranking_report_path: Path
- html_report_path: Path | None
```

The ranking result may additionally retain small typed evidence summaries needed by the HTML renderer, but it must not duplicate complete upstream scientific datasets unnecessarily.

## Evidence aggregation and integrity checks

For each Primer3 assay:

1. resolve its `region_id` against `PrimerDesignResult.candidates`;
2. collect inclusivity rows for the same `assay_id`;
3. collect degeneracy proposals for the same `assay_id`;
4. collect specificity retention and plausible amplicons for the same `assay_id`;
5. derive reasons, classification, and score components.

Internal inconsistencies fail explicitly rather than being converted into `REVIEW`:

- duplicate Primer3 `assay_id`;
- duplicate candidate `region_id`;
- assay references a missing region;
- a `COMPLETE` inclusivity result contains unknown assay IDs;
- a `COMPLETE` specificity result contains unknown assay IDs;
- duplicate degeneracy proposal for the same assay/oligo role;
- non-finite or structurally invalid metrics needed for ranking.

A `COMPLETE` inclusivity result is expected to contain exactly one assay evaluation per Evaluation Set sequence for each Primer3 assay. Missing or duplicate rows under a `COMPLETE` status are treated as inconsistent upstream evidence and fail explicitly. This is different from the stage being `SKIPPED`, which is allowed and yields incomplete evidence.

If the Evaluation Set is empty, inclusivity fraction is undefined. Ranking does not divide by zero; the inclusivity and robustness components become incomplete and the assay receives `REVIEW / EVIDENCE_INCOMPLETE`.

## Structured reasons and classification

### Reason severities

```text
HIGH_RISK -> forces HIGH_RISK
REVIEW    -> prevents IN SILICO PASS
ADVISORY  -> explanatory only
```

Initial reason code vocabulary:

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

Each reason includes source, a stable human-readable message, and structured evidence. Example:

```json
{
  "code": "DETECTABLE_OFF_TARGET",
  "severity": "HIGH_RISK",
  "source": "specificity",
  "message": "Detectable off-target found",
  "evidence": {
    "dataset_count": 1,
    "amplicon_count": 2
  }
}
```

Classification precedence:

```text
if any HIGH_RISK reason:
    HIGH_RISK
else if any REVIEW reason:
    REVIEW
else:
    IN SILICO PASS
```

### Inclusivity rules

Using the original assay only:

```text
fraction >= min_inclusivity_for_pass
    -> no class downgrade

min_inclusivity_before_high_risk <= fraction < min_inclusivity_for_pass
    -> REVIEW / INCLUSIVITY_BELOW_PASS

fraction < min_inclusivity_before_high_risk
    -> HIGH_RISK / INCLUSIVITY_BELOW_MINIMUM
```

Defaults are `1.0` for PASS and `0.90` for the HIGH_RISK floor.

### Specificity rules

Per assay across all configured off-target datasets:

```text
any detectable_off_target
    -> HIGH_RISK / DETECTABLE_OFF_TARGET

else any plausible F/R amplicon without compatible probe
    -> REVIEW / PLAUSIBLE_OFF_TARGET_AMPLICON

else compatible oligo hits only, with no plausible amplicon
    -> ADVISORY / ISOLATED_OFF_TARGET_HITS

else
    -> no specificity reason
```

A detectable off-target is the #9 contract: compatible inward-facing F/R geometry plus a compatible probe inside the plausible amplicon.

### Degeneracy proposal rules

An accepted proposal does not downgrade the class by itself.

- `ACCEPTED` -> `ADVISORY / IUPAC_PROPOSAL_ACCEPTED`;
- `REJECTED` -> `ADVISORY / IUPAC_PROPOSAL_REJECTED`;
- `UNCHANGED` -> no reason required.

The original assay remains the only ranked candidate.

## Score model

Score is absolute, deterministic, and decomposable. It is not normalized against other assays in the same run.

```text
final_score = 100 * (
    0.35 * inclusivity
  + 0.25 * specificity
  + 0.20 * conservation
  + 0.10 * primer3_quality
  + 0.10 * robustness
)
```

Weights come from configuration. Components use full internal precision and are in `[0, 1]`. Presentation rounds scores to two decimal places only.

The score never determines class and is used only to order assays within the same class.

### Inclusivity component

```text
inclusivity = original_compatible_count / Evaluation_Set_count
```

Use only the original assay compatibility from #8, not proposed IUPAC compatibility.

If inclusivity is `SKIPPED` or the Evaluation Set is empty, this component is `None`.

### Specificity component

Use complete #9 evidence, not the truncated detailed hit artifact.

When no plausible amplicon exists, compatible hit count is taken from the complete `HitRetentionSummary.total_hit_count` values for the assay, so `max_hits_per_oligo_per_dataset` cannot artificially improve the ranking.

```text
1.00                              -> zero compatible off-target hits
max(0.80, 1 - 0.02 * hit_count)  -> compatible hits but no plausible amplicon
0.40                              -> plausible F/R amplicon, no detectable probe
0.00                              -> any detectable_off_target
```

If specificity is `SKIPPED`, this component is `None`.

### Conservation component

Resolve the assay's candidate region and calculate:

```text
conservation = mean(
    mean_conservation,
    minimum_conservation,
    mean_coverage,
    1 - mean_gap_frequency,
    1 - min(mean_entropy_bits / 2, 1)
)
```

Two bits is the theoretical maximum Shannon entropy for four equiprobable canonical bases. Metrics must be finite and within their scientifically valid ranges; invalid upstream metrics fail rather than being silently clamped, except for the explicit entropy normalization cap above.

### Primer3 quality component

```text
primer3_quality = 1 / (1 + pair_penalty)
```

`pair_penalty` must be finite and non-negative when present. If it is `None`, this component is incomplete rather than treated as zero.

### Robustness component

When inclusivity evidence is complete and no accepted proposal exists:

```text
robustness = 1.0
```

For each role (`FORWARD`, `PROBE`, `REVERSE`):

```text
no accepted proposal -> role_robustness = 1.0
accepted proposal    -> role_robustness = original_degeneracy / proposed_degeneracy
```

Then:

```text
robustness = mean(FORWARD, PROBE, REVERSE)
```

Accepted proposal degeneracies must be positive and `proposed_degeneracy >= original_degeneracy`; structurally impossible values fail explicitly.

If inclusivity is `SKIPPED` or the Evaluation Set is empty, robustness is `None`.

## Incomplete score behavior

Do not substitute zero for missing evidence.

If any component is unavailable:

```text
score_status = INCOMPLETE
final_score = null
```

The assay receives `REVIEW / EVIDENCE_INCOMPLETE` unless another available reason already forces HIGH_RISK. HIGH_RISK still wins classification precedence.

Within a class, complete-score assays are ordered before incomplete-score assays.

## Deterministic ordering

Sort by:

```text
1. classification priority:
   IN SILICO PASS
   REVIEW
   HIGH_RISK

2. score status:
   COMPLETE before INCOMPLETE

3. final_score descending when present
4. inclusivity component descending when present
5. pair_penalty ascending when present
6. primer3_index ascending
7. assay_id ascending
```

For tie breaking, missing numeric values sort after present values. Rank values are assigned after this ordering and are contiguous from 1.

This guarantees that no HIGH_RISK assay outranks REVIEW or PASS regardless of score.

## Artifacts

When ranking is enabled and completes:

```text
ranking/
├── assay_ranking.tsv
└── ranking_report.json

report.html
```

### `assay_ranking.tsv`

One row per Primer3 assay. Minimum fields:

- rank;
- assay_id;
- region_id;
- classification;
- score_status;
- final_score;
- all five score components;
- inclusivity numerator/denominator/fraction when available;
- compatible off-target hit count;
- plausible amplicon count;
- detectable off-target count;
- pair_penalty;
- semicolon-separated stable reason codes.

TSV uses presentation rounding where applicable but must not be the source of truth for downstream computation.

### `ranking_report.json`

Contains:

- schema version;
- effective ranking configuration;
- class counts;
- complete/incomplete score counts;
- ordered assays;
- full score components;
- structured reasons and evidence;
- artifact paths.

### `report.html`

The final run report is self-contained, deterministic, and offline. It uses HTML escaping for all dynamic text and no CDN, remote font, external JavaScript, or network fetch.

For each ranked assay it shows:

- rank, class, score, and reason codes;
- forward/probe/reverse sequence, reference coordinates, Tm, GC, and available penalties;
- product size and pair penalty;
- candidate-region conservation metrics;
- original inclusivity count and fraction;
- IUPAC proposals as contextual evidence only;
- off-target compatible hits, plausible amplicons, and detectable off-targets by dataset;
- five named score components.

A static ordered table plus `<details>` sections is sufficient for the MVP. No charting library is required.

## Disabled and empty behavior

### Ranking disabled

`ranking.enabled: false`:

- status `SKIPPED`;
- does not inspect or validate upstream scientific results;
- publishes only `ranking/ranking_report.json` with `SKIPPED` status;
- removes stale `ranking/assay_ranking.tsv` and root `report.html` from a prior run.

### Zero Primer3 assays

Ranking enabled with `PrimerDesignResult.status == COMPLETE` and zero assays:

- status `COMPLETE`;
- publishes header-only `assay_ranking.tsv`;
- publishes `ranking_report.json` with zero assays;
- publishes an empty-state `report.html`;
- does not require inclusivity/specificity evidence rows for nonexistent assays.

### Invalid primer design state

Ranking enabled with a non-`COMPLETE` `PrimerDesignResult` fails explicitly.

## Pipeline integration

`run_pipeline()` calls ranking after specificity and before final run summary publication.

`qc_report.json` receives only a compact summary:

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

`top_recommended_assay_id` means the first `IN SILICO PASS` assay after deterministic ordering. If there is no PASS, it is `null`; REVIEW is not silently promoted to recommended.

## Error and publication behavior

Errors are contextualized with assay/source identifiers where possible.

The stage computes and validates all ranking results before publishing final artifacts. Text artifacts are written through temporary files and atomically replaced using the same project pattern as existing stages. A failed ranking must not present stale output as the result of the current run.

No source sequence database is queried and no network access is introduced.

## Test strategy

Implementation follows RED -> GREEN -> REFACTOR TDD.

Suggested test files:

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

Required coverage includes:

- ranking disabled defaults and stale cleanup;
- config validation, thresholds, and weight sum;
- 100% inclusivity PASS eligibility;
- inclusivity in `[0.90, 1.0)` -> REVIEW;
- inclusivity below `0.90` -> HIGH_RISK;
- configurable inclusivity thresholds;
- detectable off-target -> HIGH_RISK;
- plausible F/R without probe -> REVIEW;
- isolated hits -> ADVISORY and specificity score penalty;
- specificity hit scoring uses full retention totals, not truncated published hits;
- accepted/rejected IUPAC proposal advisories;
- proposal robustness calculation;
- skipped inclusivity/specificity -> REVIEW/EVIDENCE_INCOMPLETE;
- HIGH_RISK evidence still wins when another source is missing;
- empty Evaluation Set -> incomplete score without division by zero;
- missing pair penalty -> incomplete score;
- all five component formulas;
- no relative normalization against other candidates;
- class-first ordering;
- HIGH_RISK cannot outrank REVIEW or PASS by score;
- complete score before incomplete score within a class;
- deterministic tie breaks;
- all Primer3 assays preserved;
- unknown/duplicate assay evidence fails explicitly;
- missing candidate region fails explicitly;
- zero-assay empty result;
- TSV/JSON deterministic content;
- safe HTML escaping for sequences, IDs, labels, messages, and evidence text;
- offline HTML with no remote resources;
- `qc_report.json` summary and top PASS assay;
- full pipeline integration `primer_design -> inclusivity -> specificity -> ranking`.

All tests run offline in the normal `pytest` suite on `develop`. No new CircleCI dependency, integration executable, cache, schedule, parallelism, or network test is needed.

## Completion criteria

Issue #10 is complete only when:

1. all issue acceptance criteria are covered by implementation and tests;
2. classification is demonstrably evaluated before ordering by score;
3. all Primer3 assays remain auditable, including HIGH_RISK candidates;
4. score components and reason evidence are present in JSON and report HTML;
5. the normal test suite passes;
6. the merged `develop` commit is green in CircleCI;
7. README documents classification, score semantics, and the non-diagnostic/non-experimental nature of the result;
8. review finds no unresolved critical or important issue.
