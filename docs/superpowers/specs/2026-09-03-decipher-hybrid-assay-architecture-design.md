# DECIPHER hybrid assay architecture for Geison

## Goal

Evolve Geison from a target-only conservation-first qPCR assay pipeline into a reproducible diagnostic assay-design system that explicitly optimizes target inclusivity and non-target exclusivity.

The new architecture keeps Python as the pipeline orchestrator, adds DECIPHER/R as an optional specialized primer-design engine, separates primer-pair design from hydrolysis-probe design, introduces explicit target/non-target panel manifests, and validates final assays against independent challenge data.

The scientific goal is not merely to find regions conserved within a target. It is to identify assays that remain inclusive across the intended target diversity while avoiding clinically, phylogenetically, and matrix-relevant non-targets.

## Scope

This design covers four related subsystems:

1. **Panel model and manifests**: target/non-target definitions, criticality, design/challenge roles, approval, freezing, and provenance.
2. **DECIPHER engine adapter**: Python-to-R boundary for `TileSeqs`, `DesignPrimers`, and `CalculateEfficiencyPCR` while preserving the existing classic Primer3 path.
3. **qPCR assay design**: primer-pair design followed by separate hydrolysis-probe generation and assay construction.
4. **Panel intelligence and validation**: assisted panel construction, arbovirus-focused clinical knowledge, independent challenge evaluation, hard gates, ranking, and reporting.

The implementation must be incremental. Existing acquisition, QC, MAFFT alignment, conservation, checkpoints, provenance, diagnostics, and reporting infrastructure should be reused where practical rather than replaced wholesale.

## Non-goals for the first implementation cycle

The first DECIPHER MVP does **not**:

- replace MAFFT with DECIPHER `AlignSeqs`;
- replace Geison's current conservation calculations;
- use DECIPHER `DesignProbes` for hydrolysis-probe design, because that function addresses a different probe-design use case;
- make AI or web research a mandatory dependency of the scientific core;
- silently infer a clinically authoritative panel without human review;
- define universal numeric thresholds for coverage, efficiency, or cross-amplification;
- claim clinical validation from in silico evidence alone;
- require the MVP to produce a successful assay for every target.

A scientifically valid `NO_VALID_ASSAY` outcome is an acceptable successful execution.

## Architectural choice

Geison remains a Python application and pipeline orchestrator.

DECIPHER is integrated as an external R-based scientific engine behind a stable adapter boundary. Python prepares deterministic input artifacts, invokes an `Rscript` entrypoint without a shell, and parses structured outputs. Python does not depend on DECIPHER's internal R object model.

The first DECIPHER integration is deliberately limited to:

- `TileSeqs`;
- `DesignPrimers`;
- `CalculateEfficiencyPCR`.

MAFFT and Geison's conservation stage remain unchanged for the initial MVP.

The DECIPHER engine is optional initially. The classic Primer3 path remains available so the same frozen dataset can be compared across engines.

## High-level data flow

```text
Target definition
      |
      v
Target + non-target panel proposal
      |
      v
Human review / frozen panel manifest
      |
      v
Sequence acquisition + QC
      |
      v
Design/challenge partition
      |
      v
MAFFT + Geison conservation on TARGET DESIGN data
      |
      v
Primer-pair design
  |                 |
  | classic         | DECIPHER
  | Primer3         | TileSeqs -> DesignPrimers
  |                 | -> CalculateEfficiencyPCR
  +--------+--------+
           |
           v
Common PrimerPairCandidate contract
           |
           v
Hydrolysis-probe design
           |
           v
F + Probe + R assays
           |
           v
Design-set inclusivity/specificity
           |
           v
Independent challenge
           |
           v
Hard gates + ranking + report
```

## Target modes

Geison supports two explicit biological intents.

### `broad_detection`

This is the default. Geison attempts to design an assay inclusive across the relevant diversity of the named target.

Examples of relevant strata can include lineage/genotype, geography, time, host, and other metadata when available. The selection strategy should seek representation rather than merely accepting the first N NCBI results.

### `subtype_specific`

The named subtype/lineage becomes the target. Other subtypes of the same species can become non-targets when discrimination from them is part of the user's intent.

This allows the same architecture to support broad detection and discriminatory/genotyping use cases without changing the core pipeline.

## Diagnostic context

Panel construction can use contextual fields such as:

```yaml
diagnostic_context:
  syndrome: arboviral febrile disease
  geography: Brazil
  sample_type: human serum
  vector: mosquito
```

The default behavior is hybrid inference:

- Geison infers context when the target makes the inference sufficiently clear;
- ambiguous context is presented for review rather than silently treated as authoritative;
- all inferred fields are editable;
- final approved context is preserved in the frozen manifest.

## Panel model

The panel is a first-class scientific artifact rather than an incidental list of FASTA files.

Each organism/group entry records at least:

- canonical identifier, preferably NCBI TaxID when applicable;
- accepted display name;
- target/non-target role;
- design/challenge role;
- target required-group status when relevant;
- non-target criticality;
- reason for inclusion;
- source(s) that proposed it;
- approval state;
- sequence-selection provenance.

Conceptual example:

```yaml
target:
  name: West Nile virus
  mode: broad_detection
  groups:
    - name: lineage_1
      required: true
    - name: lineage_2
      required: true

non_targets:
  - organism: Usutu virus
    taxid: 64286
    criticality: CRITICAL
    reasons:
      - phylogenetic_neighbor
      - diagnostic_cross_reactivity_risk
    proposed_by:
      - taxonomy
      - geison_clinical_knowledge_base
    approved_by_user: true
```

## Non-target criticality

Three levels are supported.

### `CRITICAL`

A single clinically unacceptable predicted cross-reaction can reject an assay. Risk against one CRITICAL organism must never be hidden by averaging across many safe organisms.

Examples can include especially close phylogenetic neighbors or critical diagnostic differentials.

### `IMPORTANT`

Evidence contributes a strong penalty or warning but is not automatically a hard fail solely because a configurable intermediate-risk metric is observed.

### `BACKGROUND`

Provides additional protection against host, vector, matrix, environmental, or other background sequences. It contributes evidence and ranking but normally has less weight than CRITICAL and IMPORTANT organisms.

Numeric cutoffs remain configurable and must not be presented as universal biological constants.

## Region-level versus assay-level specificity

Geison distinguishes two concepts.

### Region-level exclusivity

Before or during primer-pair generation, Geison may use non-target evidence to prioritize target regions capable of supporting discriminating oligos.

Region-level similarity is not by itself final proof of assay cross-reactivity. A broadly similar region can still contain discriminating primer positions.

Therefore region-level exclusivity is primarily a prioritization/risk signal. Early hard rejection is reserved for cases where the configured evidence indicates that a region does not contain realistic discriminatory design space.

### Assay-level specificity

The final scientific decision is based on the actual assay geometry and oligos.

The complete forward primer, hydrolysis probe, and reverse primer are evaluated against non-target sequences. A detectable off-target for a CRITICAL organism is a hard fail regardless of aggregate score.

This preserves and extends the existing Geison concept of distinguishing isolated oligo hits, plausible amplicons, and detectable off-targets.

## DECIPHER design panel versus full specificity panel

DECIPHER's primer-design logic is used only where its aligned-group model is biologically appropriate.

The complete approved non-target panel is therefore divided conceptually into two scopes.

### DECIPHER design non-targets

This subset contains homologous/alignment-compatible non-target groups for which target-versus-non-target primer discrimination can be modeled meaningfully during DECIPHER design.

For a West Nile virus project, examples may include sufficiently homologous flaviviral neighbors such as Usutu virus or Japanese encephalitis virus when the selected datasets support a biologically meaningful alignment.

Membership is based on homology/alignment suitability, not merely on CRITICAL/IMPORTANT/BACKGROUND criticality.

### Full specificity panel

The final assay is evaluated against **all** approved non-target groups, including phylogenetically distant but clinically relevant organisms and matrix/background sequences.

A West Nile project may therefore include Dengue, Zika, Chikungunya, Mayaro, Oropouche, human, vector, and other approved groups in full specificity even when they are not appropriate members of the DECIPHER design alignment.

This avoids forcing biologically dissimilar genomes into one DECIPHER design alignment while still preserving comprehensive diagnostic specificity testing.

## Automatic panel construction

Panel automation has deterministic and assistive layers.

### Deterministic core

The core can propose a useful panel without AI or live web research by combining:

1. structured taxonomy/phylogenetic relationships;
2. a versioned Geison clinical-differential knowledge base;
3. diagnostic-context rules;
4. sample host/vector/background rules.

### Assistive discovery

Literature, web research, or AI-assisted discovery may suggest additional non-targets, but suggestions never enter the final scientific dataset without explicit provenance and user approval.

The first clinical knowledge-base scope is arboviruses, while the data model remains pathogen-agnostic.

The arbovirus seed can include relevant entities such as West Nile, Usutu, Japanese encephalitis virus, Dengue, Zika, Chikungunya, Mayaro, Oropouche, and other justified entries. Inclusion in a specific run still depends on target/context and approval; the knowledge base is not a universal fixed panel.

## Human approval gate

A newly proposed panel cannot silently enter assay design.

The first run for a new proposal produces an artifact such as:

```text
panel_proposal.yaml
status: NEEDS_PANEL_APPROVAL
```

The pipeline stops in an action-required state rather than prompting interactively inside `qpcr-pipeline run`.

After review, an approved immutable manifest is produced. Reproducible reruns can point to the frozen manifest and proceed without another interactive approval.

This is intentionally friendly to CLI automation, Colab, CI, and checkpoint/resume workflows.

## Sequence representativeness

`max_records: N` is not a sufficient scientific representation strategy by itself.

For target and important non-target groups, Geison should be able to construct a broader candidate pool, perform QC and redundancy reduction, inspect available metadata, and select a representative subset.

Useful selection dimensions include:

- lineage/genotype;
- geography;
- collection time;
- host;
- source metadata availability;
- sequence completeness/quality;
- redundancy.

The objective is representation, not mechanically equal counts in every category.

If available metadata cannot support a requested stratification, Geison records the evidence limitation rather than inventing balance.

## Design versus independent challenge

Data used to optimize an assay must be separated from data used to challenge it.

Both target and non-target pools therefore support:

- **design/discovery data**, visible to conservation and assay design;
- **challenge/holdout data**, excluded from design and used only in independent evaluation.

Target challenge data estimate out-of-design inclusivity. Non-target challenge data test exclusivity against independent examples and potentially broader diversity.

The system must enforce that challenge sequences do not participate in conservation, DECIPHER primer design, classic primer design, or probe optimization.

When insufficient sequence diversity exists to create a meaningful independent holdout, Geison proceeds only if the configured workflow permits it and records an explicit `INSUFFICIENT_EVIDENCE` limitation for that group.

## DECIPHER adapter boundary

Python owns orchestration and durable artifacts. R owns the DECIPHER calls.

A deterministic workspace contains structured inputs such as:

```text
aligned_sequences.fasta
groups.tsv
decipher_config.json
```

An R entrypoint consumes those inputs and emits structured outputs such as JSON/TSV. The adapter records stdout/stderr safely for diagnostics but does not use fragile human-readable console output as the primary contract.

The adapter must record at least:

- R version;
- DECIPHER version;
- relevant supporting-tool versions, including OligoArrayAux when used;
- exact effective engine parameters;
- input artifact hashes;
- selected engine;
- output artifact hashes.

R/DECIPHER failures are domain-specific execution failures with actionable sanitized diagnostics.

## Engine selection

Configuration concept:

```yaml
primer_pair_design:
  engine: auto
```

Supported values:

### `auto`

Use DECIPHER when the required DECIPHER environment is available and the run is compatible with the DECIPHER path. Otherwise use the classic engine and record a warning and explicit fallback provenance.

### `decipher`

DECIPHER is mandatory. Missing R/DECIPHER/OligoArrayAux or an engine failure is a technical failure such as `DECIPHER_UNAVAILABLE` or `DECIPHER_EXECUTION_FAILED`. There is no silent classic fallback.

### `classic`

Use the classic Primer3-derived path explicitly.

The selected engine must always be visible in reports and run provenance.

## Common primer-pair contract

Downstream stages must not need to know whether a primer pair came from DECIPHER or the classic engine.

Both engines normalize into a common immutable model conceptually equivalent to:

```text
PrimerPairCandidate
  id
  forward
  reverse
  amplicon coordinates/sequence reference
  target coverage evidence
  target efficiency evidence
  non-target evidence
  source engine
  source engine metrics
  provenance
```

Engine-specific metrics may be preserved in an extensible structured field, but downstream probe design, inclusivity, specificity, ranking, checkpointing, and reporting consume the stable shared fields.

## Primer-pair design responsibility

The DECIPHER engine receives target design data plus the DECIPHER-compatible non-target design subset.

Its role is to generate or rank forward/reverse primer candidates using target coverage and non-target discrimination evidence. `TileSeqs`, `DesignPrimers`, and `CalculateEfficiencyPCR` are the first integration targets.

The classic engine remains available for comparison and compatibility.

The first MVP must make side-by-side comparison possible on the same frozen panel and datasets.

## Hydrolysis-probe design

Probe design becomes a separate Geison stage.

DECIPHER primer-pair selection does not imply use of DECIPHER `DesignProbes`. Geison designs qPCR hydrolysis-probe candidates inside each surviving amplicon using a dedicated probe-design boundary.

For each primer pair, the probe stage should be able to consider multiple probe candidates rather than asking one upstream command to return an all-or-nothing forward/probe/reverse triple.

Each probe candidate is evaluated for configurable properties including:

- size;
- melting-temperature evidence;
- GC composition;
- target conservation/coverage;
- variation within the probe site;
- non-target matches and mismatch patterns;
- compatibility with the primer-defined amplicon.

If one primer pair has no acceptable probe, only that pair is rejected. Other primer pairs remain eligible.

This specifically avoids the current failure mode where otherwise usable primer options can be hidden because Primer3's internal-oligo constraints produce zero complete triples.

## Final assay model

A final assay is the explicit trio:

```text
Forward primer + Hydrolysis probe + Reverse primer
```

Scientific evaluation is performed on the trio, not merely on individual oligos.

For a CRITICAL non-target, the default classification concept is:

1. no plausible primer amplicon: `PASS`;
2. plausible primer amplicon but probe does not support detection: `HIGH_RISK` or strong penalty, subject to configurable policy;
3. plausible primer amplicon plus compatible/detectable probe: `HARD_FAIL`.

The final hard fail is based on detectable off-target evidence, not merely on the existence of any local sequence match.

## Inclusivity policy

For `broad_detection`, each configured required target group must meet the applicable inclusivity criteria. High average target coverage cannot hide systematic failure of a required lineage/genotype.

For `subtype_specific`, intended target subtype coverage is evaluated while relevant non-target subtypes can be assessed for discriminatory exclusion.

Design-set inclusivity and challenge-set inclusivity are reported separately.

## Specificity policy

Specificity is evaluated individually by non-target group and then summarized.

A global average cannot override a CRITICAL hard fail.

`IMPORTANT` and `BACKGROUND` evidence contributes warnings, penalties, and rank signals according to configured policy. Thresholds are auditable configuration, not hidden constants.

Both design non-target data and independent challenge non-target data must remain distinguishable in outputs.

## Scientific outcomes versus execution status

Technical execution and scientific result are separate dimensions.

### Execution status

Possible top-level concepts include:

- `COMPLETED`;
- `ACTION_REQUIRED`;
- `FAILED`.

`ACTION_REQUIRED` includes cases such as first-time panel approval.

`FAILED` is reserved for technical failures such as invalid input, corrupted artifacts, missing mandatory tools, malformed engine output, or subprocess failure.

### Scientific outcome

A technically completed run can produce:

- `RECOMMENDED`;
- `CANDIDATES_WITH_WARNINGS`;
- `NO_VALID_ASSAY`;
- `INSUFFICIENT_EVIDENCE`.

For example, a run in which every candidate is rejected due to CRITICAL detectable off-targets is `COMPLETED / NO_VALID_ASSAY`, not a pipeline failure.

This replaces the assumption that zero assays must always mean an incomplete technical run.

## Minimum structural requirements for `RECOMMENDED`

Exact numeric thresholds remain configurable, but a recommended assay structurally requires:

- approved/frozen panel evidence;
- required target design groups passing configured inclusivity criteria;
- target challenge evidence passing configured policy when challenge data are available and required;
- an acceptable hydrolysis probe;
- complete specificity evaluation required by the configured run;
- no CRITICAL detectable off-target hard fail;
- sufficient evidence for final ranking/recommendation.

If independent challenge evidence is unavailable, the system must not silently describe the result as having independent validation.

## Checkpoints and invalidation

The existing deterministic checkpoint/resume system remains the foundation.

New stages and manifests must participate in dependency hashes so scientifically relevant changes invalidate descendants.

Examples:

- changing an approved panel invalidates affected sequence selection and all dependent design/evaluation stages;
- changing design/challenge membership invalidates dependent stages;
- changing `primer_pair_design.engine` invalidates primer-pair design and all downstream stages;
- changing DECIPHER parameters or tool identity invalidates DECIPHER-derived primer-pair checkpoints;
- changing probe constraints invalidates probe design and downstream assay evaluation, but need not invalidate upstream conservation;
- an identical frozen manifest and effective configuration can reuse valid checkpoints.

## Proposed stage evolution

The current linear pipeline is evolved conceptually from:

```text
input -> qc -> clustering -> alignment -> conservation
      -> primer_design -> inclusivity -> specificity -> ranking
```

toward a staged model containing explicit boundaries such as:

```text
panel
  -> input/acquisition
  -> qc
  -> dataset_partition
  -> clustering
  -> alignment
  -> conservation
  -> primer_pair_design
  -> probe_design
  -> assay_evaluation
  -> challenge
  -> ranking
```

The exact internal stage decomposition can be refined in implementation plans, but the following invariants are binding:

- panel approval precedes scientific design;
- challenge data are excluded from design stages;
- primer-pair design and probe design are separate;
- engine selection applies to primer-pair design;
- final specificity evaluates complete assays against the full approved non-target scope;
- checkpoints respect these scientific boundaries.

## Implementation decomposition

This architecture is intentionally too large for one implementation PR. It is decomposed into four subprojects, delivered in order.

### 1. Panel model and manifests

Deliver first:

- target modes;
- non-target criticality;
- design/challenge roles;
- panel proposal/approval/freeze contract;
- manual West Nile fixtures;
- provenance and checkpoint integration.

This subproject does not yet need automatic panel intelligence.

### 2. DECIPHER engine and common primer-pair contract

Deliver next:

- DECIPHER environment detection in `doctor`;
- Python/R adapter;
- structured input/output contract;
- `TileSeqs`, `DesignPrimers`, `CalculateEfficiencyPCR` integration;
- `auto` / `decipher` / `classic` engine configuration;
- normalization into `PrimerPairCandidate`;
- classic-versus-DECIPHER comparison fixtures.

### 3. Separate hydrolysis-probe and assay design

Deliver next:

- probe-design stage;
- multiple probe candidates per primer pair;
- final F + Probe + R assay model;
- assay-level inclusivity/specificity hard gates;
- integration with existing off-target geometry where reusable.

### 4. Panel intelligence and independent validation

Deliver last:

- taxonomy-based panel proposal;
- versioned arbovirus clinical knowledge base;
- diagnostic-context inference;
- optional literature/research suggestions with provenance;
- representative sequence selection;
- design/challenge partitioning;
- final ranking/report changes;
- explicit scientific outcomes.

Each subproject receives its own implementation plan and can be validated before proceeding to the next.

## Error handling and observability

Failures at external-tool boundaries must remain sanitized, deterministic, and actionable.

The run record should distinguish at least:

- panel approval required;
- panel manifest invalid;
- DECIPHER unavailable;
- DECIPHER execution failed;
- DECIPHER output invalid;
- no eligible primer pairs;
- no acceptable probe for a primer pair;
- no valid assays after hard gates;
- insufficient challenge evidence.

The last three are scientific/evidence outcomes when the pipeline itself executed correctly, not necessarily technical exceptions.

Tool versions, selected engine, fallback decisions, frozen panel identity, dataset hashes, design/challenge membership, and effective thresholds must be preserved in run provenance.

## Testing strategy

### Panel tests

Verify:

- target/non-target roles;
- CRITICAL/IMPORTANT/BACKGROUND preservation;
- canonical identifiers and provenance;
- approval/freeze behavior;
- deterministic manifest serialization;
- rejection of invalid or contradictory memberships.

### Design/challenge separation tests

Verify that challenge sequences cannot enter:

- conservation;
- DECIPHER design;
- classic primer-pair design;
- probe optimization.

### DECIPHER adapter unit tests

Use fake/injected process boundaries or frozen structured outputs to verify:

- deterministic workspace generation;
- correct group membership serialization;
- output parsing;
- engine metadata capture;
- domain-specific error mapping;
- no dependency on human-readable R console text.

### Real DECIPHER integration tests

Keep separate from ordinary fast unit tests. Use a small frozen fixture and a locally installed R/DECIPHER environment. Do not require live NCBI/network access.

### Engine contract tests

Verify classic and DECIPHER paths both normalize into a compatible `PrimerPairCandidate` interface and preserve engine-specific provenance.

### Probe tests

Verify:

- multiple probes can be explored per primer pair;
- a failed probe does not eliminate unrelated pairs;
- a pair with no valid probe is rejected cleanly;
- probe provenance and metrics remain available.

### Hard-gate tests

Verify:

- CRITICAL detectable off-target rejects an assay;
- one CRITICAL failure cannot be hidden by safe averages;
- IMPORTANT risk is not automatically promoted to a CRITICAL hard fail;
- primer-only off-target risk and detectable F+Probe+R off-target remain distinguishable.

### Checkpoint tests

Verify invalidation when changing:

- panel manifest;
- design/challenge membership;
- engine selection;
- DECIPHER version/parameters;
- probe constraints;
- specificity policy.

Verify identical reruns reuse valid checkpoints.

### End-to-end MVP fixture

Use a small, frozen West Nile-centered fixture containing target diversity plus selected non-targets. The test must be deterministic and offline.

The end-to-end test succeeds if the pipeline produces an explainable scientific outcome, even when that outcome is `NO_VALID_ASSAY`.

## MVP acceptance criteria

The DECIPHER hybrid MVP is accepted when Geison can, using frozen inputs:

1. load an approved target/non-target panel manifest;
2. distinguish design from challenge data;
3. calculate target conservation using the existing path;
4. select primer pairs through the DECIPHER engine;
5. normalize those pairs through a stable Python contract;
6. generate multiple hydrolysis-probe candidates per viable primer pair;
7. evaluate complete assays against design data;
8. challenge survivors against independent target and non-target data;
9. apply CRITICAL hard gates without hiding single-organism risk in averages;
10. report selected engine, versions, panel provenance, design/challenge evidence, candidate rejection reasons, and scientific outcome;
11. reproduce the same result from the same frozen manifest, data, configuration, and tool versions.

The MVP does not require a `RECOMMENDED` assay. A deterministic `NO_VALID_ASSAY` result with complete rejection evidence satisfies the technical architecture.

## Scientific interpretation boundary

All outputs remain in silico design evidence. `RECOMMENDED` means recommended by the configured computational workflow, not clinically validated, analytically validated, or approved for diagnostic use.

Reports must preserve this distinction.

## Relationship to existing Geison design documents

This architecture supersedes only the assumptions that conflict with the new primer-pair/probe separation, explicit panel model, independent challenge model, and multi-engine design path.

Existing designs for NCBI acquisition, CD-HIT, MAFFT, conservation, inclusivity, specificity, ranking, checkpoints, and reproducibility remain useful implementation foundations unless a later subproject spec explicitly changes their contract.

In particular, the current Primer3 assay design remains the basis of the `classic` path during migration. It is not deleted as part of the initial DECIPHER integration.

## Approved architectural decisions

The following decisions are binding for implementation planning:

- Python remains the Geison orchestrator.
- DECIPHER/R is an optional specialized engine, recommended/default through `auto` when available and compatible.
- The DECIPHER MVP covers `TileSeqs`, `DesignPrimers`, and `CalculateEfficiencyPCR` only.
- MAFFT and Geison conservation remain in place initially.
- `broad_detection` is the default target mode; `subtype_specific` is explicitly supported.
- The panel uses CRITICAL, IMPORTANT, and BACKGROUND non-target criticality.
- A new panel requires human approval before design; frozen manifests rerun without interactive approval.
- Target and non-target data are separated into design and independent challenge subsets.
- The automatic panel architecture is generic, while the first clinical knowledge-base content focuses on arboviruses.
- DECIPHER design non-targets are limited to biologically meaningful homologous/alignment-compatible groups.
- Final specificity always uses the full approved non-target panel.
- DECIPHER designs primer pairs; Geison separately designs hydrolysis probes.
- The complete F + Probe + R assay is the unit of final specificity decision.
- A detectable CRITICAL off-target is a hard fail.
- Execution status and scientific outcome are separate.
- `NO_VALID_ASSAY` can be a successful completed run.
- Silent engine fallback is prohibited when `engine: decipher`; `engine: auto` may fall back only with explicit provenance and warning.
- The architecture is implemented as four sequential subprojects rather than one monolithic change.
