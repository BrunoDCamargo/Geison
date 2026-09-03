# Contrastive Conservation + Guided Colab UX

Date: 2026-09-03
Status: Design approved in chat; implementation not started

## 1. Purpose

Extend Geison with an explicit scientific stage between target conservation and assay design that answers a question the current pipeline does not model directly:

> Which regions are stable within the target diversity and, at the same time, differentiated from the approved challenge panel?

The new stage is named `contrastive_conservation` in the core and presented to notebook users as **Target vs non-target contrast** / **Contraste alvo x nao alvo**.

The feature must remain CLI-first and reproducible. The notebook is a guided user interface over Geison artifacts and commands, not a second implementation of scientific logic.

The primary demonstration scenario is a **West Nile-like broad-detection workflow using synthetic data**. The story and panel structure may refer to West Nile and differential arbovirus categories, but the bundled demo must not contain organism-specific assay sequences, tuned biological cutoffs, or an embedded recipe for a real diagnostic assay.

## 2. Problem in the current architecture

The current stage order is:

`panel -> input -> qc -> clustering -> alignment -> conservation -> primer_design -> inclusivity -> specificity -> ranking`

`primer_design` currently derives candidate regions directly from `ConservationResult`. Candidate eligibility and ordering therefore answer only target-side questions such as conservation, coverage, gaps, entropy and usable fraction.

The later `specificity` stage works at a different level. It evaluates already-designed forward primer, probe and reverse primer candidates against off-target datasets and reports oligo hits and plausible off-target amplicons.

These are different questions:

- **Conservation:** is a region stable within the target set?
- **Contrastive conservation:** is a target-stable region differentiated from the approved challenge datasets?
- **Assay specificity:** after oligos exist, do those concrete oligos create concerning off-target evidence?

Geison needs all three. The new stage does not replace `specificity`.

## 3. Goals

1. Add `contrastive_conservation` as a first-class resumable pipeline stage.
2. Reuse the approved Panel Model as the semantic definition of target groups, non-targets, dataset roles and criticality.
3. Reuse the same resolved CHALLENGE datasets later used by assay specificity rather than creating a second off-target registry.
4. Preserve per-dataset evidence instead of reducing the entire challenge panel to a single average.
5. Produce consolidated candidate regions suitable for `primer_design`, not dozens of overlapping windows exposed as separate user choices.
6. Preserve existing behavior when the new stage is disabled.
7. Deliver a guided Colab notebook with progressive disclosure, human-readable results and advanced audit details.
8. Keep all scientific computation in the Geison package/CLI; the notebook only configures, invokes and renders.

## 4. Non-goals

This subproject does not:

- automatically construct a biological challenge panel;
- decide which real organisms must be included for a particular diagnostic use case;
- ship organism-specific diagnostic thresholds;
- download or curate a real West Nile assay dataset as part of the demo;
- replace the existing `specificity` stage;
- implement a web application;
- replace Panel approval with automatic decisions;
- change the meaning of `IN SILICO PASS` or claim experimental validation.

## 5. Pipeline architecture

The new stage order is:

`panel -> input -> qc -> clustering -> alignment -> conservation -> contrastive_conservation -> primer_design -> inclusivity -> specificity -> ranking`

Dependencies:

- `contrastive_conservation` depends on `conservation`.
- `primer_design` depends on `contrastive_conservation` when contrastive analysis is enabled.
- In disabled/legacy mode, `primer_design` retains the current target-conservation selection behavior.
- `specificity` continues to depend on `primer_design` and uses the same CHALLENGE datasets independently at the oligo/amplicon level.

The contrastive stage receives target-side evidence from `ConservationResult` and challenge-side sequence datasets resolved from existing off-target configuration.

### 5.1 Backward compatibility

`contrastive_conservation.enabled` defaults to `false`.

When disabled:

- the stage publishes a `SKIPPED` report;
- existing configurations remain valid;
- `primer_design` uses the existing conservation-only candidate selection path;
- existing CLI commands and manifests remain readable;
- older runs without this checkpoint remain valid historical runs, not corrupt runs.

When enabled:

- an approved frozen panel is required;
- at least one approved non-target with `CHALLENGE` role must resolve to a dataset;
- target-side conservation must be complete;
- each dataset used for contrast must be fingerprinted for reproducibility.

## 6. Panel and challenge-dataset contract

The Panel Model remains the semantic source of truth. It already carries:

- target identity and target mode;
- target groups and `required` status;
- dataset roles `DESIGN` and `CHALLENGE`;
- non-target criticality `CRITICAL`, `IMPORTANT` and `BACKGROUND`;
- reasons and provenance.

The existing `off_targets` configuration remains the physical dataset source used by both contrastive analysis and final specificity.

### 6.1 Name mapping

For v1, challenge datasets are mapped by exact normalized name between:

- approved `PanelNonTarget.name` entries that contain the `CHALLENGE` role; and
- configured `off_targets[].name` entries.

Rules:

1. Every enabled CHALLENGE non-target selected for contrast must resolve to exactly one off-target dataset.
2. Duplicate normalized names are invalid.
3. A configured off-target not represented in the approved panel may still be used by final specificity only if explicitly allowed by existing configuration rules, but it is not silently promoted into the contrastive panel.
4. Dataset provenance and hashes are recorded in contrastive artifacts.
5. Criticality comes from the approved panel, never from the off-target file name or notebook UI.

This avoids two independent registries for the same scientific evidence.

## 7. Contrastive evidence model

The stage operates on target-conserved windows in reference coordinates and evaluates how strongly each window is represented in each challenge dataset.

The production contract intentionally separates **measurement** from **policy**.

### 7.1 Measurements retained per window

For each target window Geison records at least:

- reference start/end;
- target mean conservation;
- target minimum conservation;
- target mean coverage;
- target mean gap frequency;
- target mean entropy;
- per-challenge-dataset similarity evidence;
- dataset name;
- dataset criticality;
- challenge sequence count;
- dataset source fingerprint;
- worst non-target dataset;
- worst non-target criticality;
- worst observed similarity;
- an ordering metric derived from target stability plus challenge differentiation.

The exact low-level similarity engine is an implementation detail of this subproject and must be testable behind a focused interface. The stage output must not discard per-dataset values after producing an aggregate ordering metric.

### 7.2 Criticality semantics

Criticality affects interpretation and ordering, not raw measurements.

A region with concerning evidence in a `CRITICAL` dataset must not look equivalent in the UI to a region whose only concerning evidence occurs in `BACKGROUND`.

The report must therefore expose at least:

- worst evidence across all challenge datasets;
- worst evidence among `CRITICAL` datasets;
- worst evidence among `IMPORTANT` datasets;
- supporting per-dataset rows.

The first implementation must not hide these categories inside a single opaque score.

### 7.3 No hidden biological defaults

The core must not silently embed organism-specific cutoffs.

The stage may produce continuous evidence and deterministic ranking without claiming a universal biological PASS threshold. If policy thresholds are supported, they must be explicit configuration values and recorded in the run manifest/report.

The guided synthetic demo may use illustrative values solely to make the architecture visible; those values must be labeled as synthetic-demo parameters and must not become production defaults.

## 8. Window consolidation into candidate regions

Raw sliding windows are an analysis primitive, not the primary user object.

The stage consolidates overlapping high-ranking windows into candidate regions before handing them to assay design.

Example conceptually:

`721-800, 731-810, 741-820, 751-830 -> one candidate region`

The consolidation contract must be deterministic and auditable:

- every candidate region lists its contributing windows;
- reference coordinates remain 1-based inclusive in public artifacts;
- candidate regions have stable deterministic IDs such as `contrast-region-001`;
- ordering is deterministic;
- overlap consolidation parameters are recorded in configuration/report;
- raw windows remain available in a separate artifact.

No candidate-region consolidation may delete the underlying evidence needed to explain why a region was ranked.

## 9. Core result type

Introduce a typed `ContrastiveConservationResult` with, conceptually:

- `status`: `SKIPPED` or `COMPLETE`;
- `reference_id`;
- window evidence rows;
- per-dataset evidence rows;
- consolidated candidate regions;
- paths to published artifacts;
- challenge dataset names and provenance summary.

The result must be serializable through the checkpoint codec without depending on notebook state.

## 10. Artifacts

When enabled, publish under `contrastive_conservation/`:

1. `window_metrics.tsv`
   - one row per target window with aggregate contrast fields;
2. `dataset_metrics.tsv`
   - one row per window x challenge dataset;
3. `candidate_regions.tsv`
   - consolidated candidate regions and their traceable metrics;
4. `contrastive_conservation_report.json`
   - configuration, panel identity, dataset provenance, counts, summaries and artifact paths;
5. `report.html`
   - stage-specific human-readable contrast report under `contrastive_conservation/report.html`, not the root final report.

The root `report.html` ownership remains unchanged:

- conservation may publish it earlier;
- ranking owns/replaces the root final report when ranking is enabled.

The notebook may render JSON/TSV directly, but must not require notebook-only artifacts for reproducibility.

## 11. Checkpoint and invalidation rules

Add a checkpoint codec for `contrastive_conservation`.

Its request fingerprint includes:

- the completed conservation dependency fingerprint;
- Geison version;
- contrastive configuration;
- approved panel manifest SHA-256;
- resolved CHALLENGE dataset names;
- challenge dataset record hashes/manifests;
- criticality and role information relevant to interpretation.

Expected invalidation behavior:

- changing only a Primer3 parameter does not invalidate contrastive conservation;
- changing target input/alignment/conservation invalidates contrastive and all descendants;
- changing a challenge dataset invalidates contrastive, primer design, inclusivity where dependent, specificity and ranking, but does not rerun target acquisition/alignment/conservation;
- changing non-target criticality or CHALLENGE membership invalidates contrastive and descendants, but not target conservation;
- changing only final specificity tolerances must not invalidate contrastive conservation;
- corruption or removal of declared contrastive artifacts makes the checkpoint non-reusable.

## 12. Primer-design integration

Current primer design selects regions directly from conservation windows. The new design introduces two paths.

### 12.1 Legacy/disabled path

`contrastive_conservation.enabled = false`

`ConservationResult -> existing candidate selection -> Primer3`

Behavior remains unchanged.

### 12.2 Contrastive path

`contrastive_conservation.enabled = true`

`ConservationResult + ContrastiveConservationResult -> contrast-ranked candidate regions -> Primer3`

`primer_design` still receives target consensus/reference information from conservation because Primer3 needs the target sequence context. Candidate-region choice, however, is constrained/ordered by the contrastive result rather than recomputed independently from all conservation windows.

The primer-design report must record whether each candidate came from:

- `CONSERVATION_ONLY`; or
- `CONTRASTIVE_CONSERVATION`.

Each assay keeps a `region_id` that traces back to the contrastive candidate region and its contributing windows.

## 13. Relationship with final specificity

Final assay specificity remains independent confirmation.

The contrastive stage asks a region-level question before oligo design. Final specificity asks an assay-level question after concrete oligos exist.

The same CHALLENGE datasets are reused, but the evidence types differ:

- contrastive: regional differentiation evidence;
- specificity: forward/probe/reverse hits, orientation, mismatches and plausible amplicon geometry.

Ranking may consume both evidence streams, but it must not infer that strong regional contrast guarantees assay specificity.

## 14. Guided Colab experience

The notebook becomes a workbench-style interface over the CLI.

File target for the new experience:

`notebooks/geison_guided_colab.ipynb`

The current `geison_colab.ipynb` remains as the low-level operational notebook until the guided notebook reaches acceptance criteria. The guided notebook must not break or replace it prematurely.

### 14.1 Design principles

- progressive disclosure;
- plain-language stage names;
- visible scientific provenance;
- no scientific functions implemented in notebook cells;
- clear separation between demo and real-project inputs;
- user decisions shown before execution;
- technical details available but collapsed/secondary;
- deterministic reruns;
- compatible with Colab `Run all` except explicit human review gates already required by Panel approval.

### 14.2 Colab-native forms

Use Colab form metadata (`#@title`, `#@param`) for common inputs rather than asking users to edit YAML directly.

Generated YAML is displayed in an **Advanced / generated config** section and saved as an artifact. Users can still edit YAML in advanced mode.

No custom web frontend is introduced in this subproject.

### 14.3 User journey

#### Screen/section 0: Welcome

Title: `Geison - Assay discovery workbench`

Show:

- project purpose;
- current Geison commit/version;
- `Demo (synthetic)` vs `Project` mode;
- statement that in-silico evidence does not replace experimental validation.

The bundled demo is labeled:

`West Nile-like broad-detection scenario (synthetic sequences)`

It demonstrates target diversity and arbovirus-like challenge categories without bundling real assay sequences or tuned organism-specific parameters.

#### Section 1: Environment

One setup cell:

- clone/update `main`;
- install Geison;
- install required external tools;
- run `qpcr-pipeline doctor`;
- render a compact status card.

Normal users see `Ready`/`Action needed`. Raw versions and paths are under Advanced.

#### Section 2: Project and panel

User-facing form fields:

- project name;
- target display name;
- target mode (`broad detection` / `subtype specific`);
- diagnostic context text;
- target groups;
- challenge panel entries;
- challenge criticality.

For the synthetic West Nile-like demo, the notebook presents example challenge categories such as:

- related flavivirus-like synthetic groups marked `CRITICAL`;
- clinically/epidemiologically overlapping arbovirus-like synthetic groups marked `IMPORTANT`;
- optional background groups.

The UI may display familiar organism labels as narrative examples, but the bundled data remain synthetic.

The notebook renders a panel review card before approval.

Existing `panel proposal -> ACTION_REQUIRED -> panel approve -> frozen_manifest -> resume` behavior remains the scientific approval mechanism.

#### Section 3: Data readiness

Show compact cards:

- target sequences available;
- target groups represented;
- challenge datasets resolved;
- critical/important/background coverage;
- warnings for missing mappings.

Do not proceed with enabled contrastive conservation if required CHALLENGE mappings are unresolved.

#### Section 4: Target conservation

Render:

- sequence count;
- target-group representation;
- reference ID;
- number of conservation windows;
- interactive or static genome conservation view;
- plain-language summary: `Where is the target stable?`

The user can expand raw conservation artifacts.

#### Section 5: Target vs non-target contrast

This is the central notebook view.

Render two complementary visuals:

1. **Quadrant plot**
   - Y: target conservation;
   - X: challenge similarity/evidence;
   - highlight contrast-ranked regions;
   - plain-language quadrant labels.

2. **Reference track**
   - target conservation across reference;
   - worst challenge similarity across reference;
   - candidate-region overlays;
   - click/selection support if feasible without custom frontend.

Below the plots, show candidate-region cards/table with:

- region ID and coordinates;
- target stability summary;
- worst challenge dataset;
- worst challenge criticality;
- contrast interpretation;
- evidence completeness;
- expandable per-dataset details.

Users primarily see consolidated regions, not raw overlapping windows.

#### Section 6: Assay design

Show how many regions are passed to Primer3 and how many assay candidates are returned.

Normal view:

- assay ID;
- source region;
- design status;
- evidence completeness.

Advanced view:

- oligo coordinates and technical Primer3 artifacts already produced by Geison.

The notebook itself does not compute or modify oligo sequences.

#### Section 7: Inclusivity and specificity

Render challenge results by assay and dataset.

Normal view uses statuses such as:

- `No concerning evidence found`;
- `Review`;
- `High-risk evidence`;
- `Not evaluated`.

Expandable detail links back to Geison specificity artifacts. Raw hits/amplicons remain available under Advanced.

#### Section 8: Final ranking

Present a short ranked list with explanations, not only a numeric score.

For each assay show reasons such as:

- target representation evidence complete;
- high target conservation;
- strong regional contrast;
- no critical specificity signal;
- missing evidence requiring review.

The final view also shows run status (`COMPLETED` / `PARTIAL`) and missing evidence from `run_manifest.json`.

#### Section 9: Reproducibility package

One cell lists/downloads or packages:

- generated config;
- approved panel manifest;
- run manifest/log;
- stage reports;
- candidate-region evidence;
- assay/ranking artifacts;
- commit/version information.

The notebook does not hide the underlying files.

## 15. User-facing language

Core/internal -> Notebook label:

- `panel` -> `Panel`
- `conservation` -> `Target conservation`
- `contrastive_conservation` -> `Target vs non-target contrast`
- `primer_design` -> `Assay design`
- `inclusivity` -> `Target coverage`
- `specificity` -> `Specificity`
- `ranking` -> `Final candidates`

Avoid unexplained jargon in the primary view. Terms such as entropy, mismatch position, amplicon geometry, checkpoint fingerprint and manifest SHA are shown in Advanced details with short definitions.

## 16. Error and action-required UX

The notebook maps CLI/runtime outcomes to explicit cards:

- `ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED` -> `Review and approve the panel before scientific execution.`
- missing challenge mapping -> `Challenge dataset missing for an approved panel entry.`
- target mismatch -> `The approved panel belongs to a different target.`
- missing external tool -> `Environment setup incomplete.`
- `PARTIAL` -> show exact `missing_evidence` from the run manifest;
- `FAILED` -> show failed stage and sanitized diagnostic message;
- stale checkpoint -> explain which stages will rerun before execution.

The notebook must never convert a failed or partial run into a visually successful result.

## 17. Reporting and explanation requirements

Every candidate region must be explainable without reading source code.

A region detail view must answer:

1. Why was this region considered stable in the target?
2. Which challenge dataset was most similar?
3. What was that dataset's criticality?
4. Which other challenge datasets were evaluated?
5. Which raw windows contributed to this region?
6. Which assay candidates were designed from it?
7. Which final specificity results relate to those assays?

No opaque score may be the only explanation for candidate ordering.

## 18. Testing strategy

### 18.1 Unit tests

Cover:

- config parsing/validation;
- panel-to-off-target mapping;
- per-dataset metric preservation;
- criticality summaries;
- deterministic ordering;
- window consolidation;
- empty target/challenge inputs;
- duplicate dataset names;
- missing CHALLENGE dataset mapping;
- disabled backward-compatible path;
- serialization/checkpoint codec.

### 18.2 Pipeline tests

Cover:

- stage order and dependencies;
- contrastive checkpoint creation;
- resume/reuse;
- invalidation on challenge dataset changes;
- invalidation on panel criticality/role changes;
- no invalidation of target conservation when only challenge evidence changes;
- downstream invalidation of primer design/specificity/ranking;
- legacy config behavior unchanged.

### 18.3 Synthetic scientific fixtures

Keep a deterministic synthetic fixture with:

- a target-shared conserved region that should rank poorly for contrast;
- a target-stable/challenge-divergent region that should rank highly;
- overlapping windows that consolidate into one region;
- at least two challenge criticality classes.

Tests assert relational behavior, not organism-specific biological thresholds.

### 18.4 Notebook tests

Static tests verify that the guided notebook:

- calls the Geison CLI rather than importing scientific internals;
- contains no scientific helper functions;
- exposes demo/project mode;
- performs Panel approval;
- renders target conservation and contrast sections;
- surfaces `run_manifest.json` status/missing evidence;
- exposes generated config and artifact paths;
- uses synthetic data for the bundled West Nile-like demo.

A lightweight execution/smoke path verifies notebook-supporting CLI commands in CI without requiring a browser.

## 19. Migration and rollout

Implementation should be staged:

1. core result/config/artifact model;
2. contrastive computation and deterministic synthetic tests;
3. checkpoint/resume integration;
4. primer-design integration with backward-compatible fallback;
5. report rendering;
6. guided Colab notebook;
7. documentation and CI smoke coverage.

The existing official Colab notebook remains available throughout rollout.

The synthetic spike notebook remains under `notebooks/spikes/` as design evidence and is not promoted to production UI.

## 20. Acceptance criteria

The subproject is accepted when all of the following are true:

1. `contrastive_conservation` is a typed first-class stage in the execution graph.
2. Existing configs with the stage omitted/disabled preserve previous behavior.
3. Enabled contrastive runs require and record approved panel provenance.
4. CHALLENGE non-targets map deterministically to the same datasets available to specificity.
5. Per-dataset evidence and criticality are preserved in artifacts.
6. Overlapping windows are consolidated into deterministic candidate regions.
7. Primer design can consume contrastive candidate regions while retaining target consensus context.
8. Final specificity remains independently evaluated after oligo design.
9. Resume/checkpoint invalidation follows the rules in section 11.
10. A deterministic synthetic fixture distinguishes a shared-conserved region from a target-stable/challenge-divergent region.
11. The guided notebook contains no duplicated scientific implementation.
12. A normal notebook user can complete setup, panel review, run, contrast review and result inspection without manually writing YAML.
13. Advanced users can inspect generated YAML, manifests, tool versions and raw artifacts.
14. `PARTIAL`, `FAILED` and `ACTION_REQUIRED` states are never presented as completed success.
15. Full regression and focused tests pass in CircleCI.
16. The bundled West Nile-like demonstration remains synthetic and does not embed organism-specific assay sequences or hidden biological cutoffs.

## 21. Resulting conceptual model

The final architecture is:

`Target diversity -> Conservation -> Target/non-target contrast -> Candidate regions -> Primer3 -> Inclusivity -> Assay specificity -> Ranking`

with the approved Panel Model governing target groups, challenge membership, criticality and provenance across the entire run.

This makes the scientific questions explicit instead of asking final assay specificity to compensate for a missing pre-design discrimination stage.