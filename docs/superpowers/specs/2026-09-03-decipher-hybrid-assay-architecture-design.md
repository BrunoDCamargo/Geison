# DECIPHER hybrid assay architecture for Geison

## Goal

Evolve Geison from a target-only conservation-first qPCR assay pipeline into a reproducible diagnostic assay-design system that explicitly optimizes **target inclusivity** and **non-target exclusivity**.

Geison remains the Python orchestrator. DECIPHER/R becomes an optional specialized primer-pair engine. Primer-pair design is separated from hydrolysis-probe design. Target/non-target panels become explicit, approved, frozen scientific artifacts, and final assays are challenged against data that did not participate in design.

The scientific question becomes:

> Which assay remains inclusive across the intended target diversity while avoiding clinically, phylogenetically, and matrix-relevant non-targets?

## Scope

The architecture is split into four sequential subprojects:

1. **Panel model and manifests**: target/non-target roles, criticality, design/challenge membership, approval, freezing, and provenance.
2. **DECIPHER engine**: Python/R adapter for `TileSeqs`, `DesignPrimers`, and `CalculateEfficiencyPCR`, while preserving the classic Primer3 path.
3. **qPCR assay design**: common primer-pair contract, separate hydrolysis-probe generation, and complete F + Probe + R assays.
4. **Panel intelligence and validation**: assisted panel construction, arbovirus knowledge base, representative datasets, independent challenge, hard gates, ranking, and reporting.

Existing NCBI acquisition, QC, CD-HIT, MAFFT, conservation, checkpointing, provenance, diagnostics, and reporting infrastructure should be reused where practical.

## Non-goals for the first DECIPHER MVP

The first integration does not:

- replace Geison's target alignment/conservation path with DECIPHER `AlignSeqs`;
- replace Geison's conservation calculations;
- use DECIPHER `DesignProbes` for qPCR hydrolysis probes;
- make AI or web research mandatory for the scientific core;
- silently create an authoritative clinical panel without human approval;
- define universal hidden numeric thresholds;
- claim clinical or analytical validation from in silico evidence;
- require every run to produce a valid assay.

`COMPLETED / NO_VALID_ASSAY` is a legitimate successful scientific result.

## Core architecture

```text
Target definition
      |
      v
Panel proposal
(target + non-target + context)
      |
      v
Human review / frozen manifest
      |
      v
Acquisition + QC + representative selection
      |
      v
Design / challenge partition
      |
      +------------------------------+
      |                              |
      v                              v
TARGET DESIGN                  CHALLENGE DATA
      |                              |
      v                              |
MAFFT target alignment               |
      |                              |
Geison conservation                  |
      |                              |
      v                              |
Primer-pair design                   |
  |              |                   |
  | classic      | DECIPHER          |
  |              |                   |
  |              + target design     |
  |              + homologous        |
  |                non-target design |
  |              -> separate MAFFT   |
  |                 design alignment |
  |              -> TileSeqs         |
  |              -> DesignPrimers    |
  |              -> EfficiencyPCR    |
  +-------+------+                   |
          |                          |
          v                          |
Common PrimerPairCandidate           |
          |                          |
          v                          |
Hydrolysis-probe design              |
          |                          |
          v                          |
Complete F + Probe + R assays        |
          |                          |
          v                          |
Design-set evaluation                |
          |                          |
          +------------+-------------+
                       |
                       v
Independent challenge
                       |
                       v
Hard gates + ranking + report
```

### Two distinct alignments

The existing **target alignment** remains responsible for target conservation and is not replaced.

DECIPHER additionally needs an alignment containing the **target design sequences plus the homologous/alignment-compatible non-target design sequences** used for discrimination. The DECIPHER adapter therefore prepares a separate deterministic **DECIPHER design alignment**, initially using MAFFT.

This preserves the decision not to migrate Geison's main alignment/conservation engine while satisfying DECIPHER's aligned-group input model.

The two alignments have separate artifacts, hashes, provenance, and checkpoint identities.

## Target modes

### `broad_detection`

Default mode. Geison attempts to cover the relevant diversity of the target, such as lineages/genotypes, geography, time, and host when suitable metadata exist.

Selection seeks representation rather than simply taking the first N NCBI results.

### `subtype_specific`

A requested subtype/lineage becomes the target. Other subtypes can become non-targets when discrimination among them is part of the design intent.

## Diagnostic context

Panel construction can use fields such as:

```yaml
diagnostic_context:
  syndrome: arboviral febrile disease
  geography: Brazil
  sample_type: human serum
  vector: mosquito
```

The behavior is hybrid:

- infer context when sufficiently clear;
- expose inferred values for review;
- ask only where ambiguity matters;
- preserve approved context in the frozen manifest.

## Panel as a first-class artifact

Each group records at least:

- canonical identifier, preferably NCBI TaxID when applicable;
- display name;
- target/non-target role;
- design/challenge role;
- required target-group status when applicable;
- non-target criticality;
- reasons for inclusion;
- proposal source(s);
- approval status;
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

### `CRITICAL`

A single unacceptable predicted cross-reaction can reject an assay. One CRITICAL failure cannot be hidden by a good average across many other organisms.

### `IMPORTANT`

Produces strong penalty/warning evidence but does not automatically behave like a CRITICAL hard fail.

### `BACKGROUND`

Adds protection against host, vector, matrix, environmental, or other background sequences with normally lower weighting.

Numeric thresholds remain explicit configuration and must not be presented as universal biological constants.

## Region-level versus assay-level specificity

### Region-level exclusivity

Non-target evidence can prioritize regions with useful discriminatory design space. Similarity alone is not final proof of assay cross-reactivity because discriminatory oligo positions may still exist inside a broadly homologous region.

Region-level evidence is therefore primarily prioritization/risk information. Early rejection is reserved for configured evidence showing that realistic discriminatory design space is absent.

### Assay-level specificity

The final decision is made on the actual **Forward + Probe + Reverse** assay and its geometry against non-target sequences.

A detectable off-target involving a CRITICAL non-target is a hard fail regardless of aggregate score.

## DECIPHER design panel versus full specificity panel

The complete approved non-target universe is not forced into one DECIPHER alignment.

### DECIPHER design non-targets

Only homologous/alignment-compatible groups suitable for meaningful target-versus-non-target discrimination during DECIPHER primer design are included.

For West Nile virus, sufficiently homologous flaviviral neighbors such as Usutu virus or Japanese encephalitis virus may be appropriate when the selected datasets support meaningful alignment.

Membership depends on homology/alignment suitability, not merely CRITICAL/IMPORTANT/BACKGROUND status.

The adapter constructs the separate MAFFT **DECIPHER design alignment** from:

```text
TARGET DESIGN
+
DECIPHER-COMPATIBLE NON-TARGET DESIGN
```

### Full specificity panel

Every approved non-target can participate in final assay specificity, including distant but clinically relevant organisms and host/vector/background sequences.

A West Nile run may therefore test Usutu, JEV, Dengue, Zika, Chikungunya, Mayaro, Oropouche, human, vector, and other approved groups even though only a homologous subset participates in DECIPHER primer design.

## Automatic panel construction

### Deterministic core

A useful proposal should be possible without AI by combining:

1. structured taxonomy/phylogenetic relationships;
2. a versioned Geison clinical-differential knowledge base;
3. diagnostic-context rules;
4. host/vector/background rules.

### Assistive discovery

Literature, web research, or AI can suggest additional non-targets, but suggestions require explicit provenance and human approval before entering the frozen panel.

The first knowledge-base content focuses on arboviruses while the architecture remains pathogen-agnostic.

The seed knowledge can include relevant entities such as West Nile, Usutu, JEV, Dengue, Zika, Chikungunya, Mayaro, and Oropouche. Presence in the knowledge base does not mean automatic inclusion in every run.

## Human approval gate

A newly generated proposal produces an action-required artifact such as:

```text
panel_proposal.yaml
status: NEEDS_PANEL_APPROVAL
```

`qpcr-pipeline run` does not hide interactive `[y/n]` questions inside execution.

After review, Geison produces an approved immutable manifest. Reruns can reference that frozen manifest and proceed automatically.

## Representative sequence selection

`max_records: N` alone is not a scientific representation strategy.

For target and relevant non-target groups, Geison should be able to:

1. acquire a broader candidate pool;
2. run QC;
3. inspect available metadata;
4. reduce redundancy;
5. select representative sequences.

Useful strata can include lineage/genotype, geography, collection time, host, sequence quality/completeness, and metadata availability.

The goal is representation, not artificial equal counts. Missing metadata or limited diversity produces an explicit evidence limitation rather than fabricated balance.

## Design versus independent challenge

Both target and non-target pools support:

- **design/discovery data**, visible to optimization;
- **challenge/holdout data**, never visible to design.

Challenge sequences must not participate in:

- conservation;
- the DECIPHER design alignment;
- DECIPHER primer design;
- classic primer-pair design;
- probe optimization.

Target challenge evaluates out-of-design inclusivity. Non-target challenge evaluates exclusivity against independent examples and potentially broader diversity.

If a meaningful holdout cannot be created, the limitation is explicitly recorded. The system must not claim independent validation where no independent challenge evidence exists.

## DECIPHER adapter boundary

Python owns orchestration and durable artifacts. R owns DECIPHER calls.

A deterministic adapter workspace contains inputs such as:

```text
decipher_design_alignment.fasta
groups.tsv
decipher_config.json
```

The design alignment is created from target design plus DECIPHER-compatible non-target design sequences using the configured aligner, initially MAFFT.

An `Rscript` entrypoint consumes structured inputs and emits structured JSON/TSV outputs. Human-readable console text is diagnostic only, not the primary API contract.

Provenance records at least:

- MAFFT version used for the DECIPHER design alignment;
- R version;
- DECIPHER version;
- supporting-tool versions such as OligoArrayAux when used;
- effective engine parameters;
- input/output hashes;
- selected engine.

## Engine selection

```yaml
primer_pair_design:
  engine: auto
```

### `auto`

Use DECIPHER when the required environment is available and the run supports the DECIPHER path. Otherwise use classic and record an explicit warning and fallback provenance.

### `decipher`

DECIPHER is mandatory. Missing dependencies or engine failure produces a technical failure such as `DECIPHER_UNAVAILABLE` or `DECIPHER_EXECUTION_FAILED`. No silent fallback is allowed.

### `classic`

Use the classic Primer3-derived path explicitly.

The report and run provenance always show both requested and selected engine.

## Common primer-pair contract

Classic and DECIPHER outputs normalize into the same downstream model:

```text
PrimerPairCandidate
  id
  forward
  reverse
  amplicon coordinates/reference
  target coverage evidence
  target efficiency evidence
  non-target evidence
  source engine
  engine-specific metrics
  provenance
```

Probe design, inclusivity, specificity, ranking, checkpointing, and reporting consume this stable contract rather than engine internals.

## DECIPHER primer-pair responsibility

The DECIPHER path initially integrates:

- `TileSeqs`;
- `DesignPrimers`;
- `CalculateEfficiencyPCR`.

It consumes the separate DECIPHER design alignment and produces forward/reverse pairs with target coverage/efficiency and non-target discrimination evidence.

The first MVP must allow side-by-side classic-versus-DECIPHER comparison on the same frozen target/panel datasets.

## Separate hydrolysis-probe design

Probe design becomes its own Geison stage.

For every viable primer pair, Geison explores multiple hydrolysis-probe candidates inside the amplicon. Configurable evidence can include:

- length;
- Tm;
- GC composition;
- target conservation/coverage;
- variation within the probe site;
- non-target matches/mismatch patterns;
- compatibility with the primer-defined amplicon.

If a primer pair has no acceptable probe, that pair is rejected and other pairs continue.

This avoids the current all-or-nothing behavior where otherwise usable primers can disappear because Primer3 returns zero internal oligos under the configured probe constraints.

## Final assay and specificity policy

A final assay is explicitly:

```text
Forward primer + Hydrolysis probe + Reverse primer
```

For a CRITICAL non-target, the conceptual classification is:

1. no plausible primer amplicon: `PASS`;
2. plausible primer amplicon but probe does not support detection: `HIGH_RISK` / strong penalty according to policy;
3. plausible amplicon plus compatible/detectable probe: `HARD_FAIL`.

The definitive hard fail is therefore based on detectable complete-assay off-target evidence, not merely any local alignment hit.

## Inclusivity policy

In `broad_detection`, each required target group must meet the configured inclusivity policy. High average coverage cannot hide failure of a required lineage/genotype.

In `subtype_specific`, the intended subtype is evaluated for inclusion while selected sibling subtypes can be evaluated for exclusion.

Design-set and challenge-set inclusivity are reported separately.

## Execution status versus scientific outcome

### Execution status

Top-level concepts:

- `COMPLETED`;
- `ACTION_REQUIRED`;
- `FAILED`.

`ACTION_REQUIRED` includes panel approval. `FAILED` is reserved for technical problems such as invalid inputs, corrupted artifacts, mandatory tool absence, malformed engine output, or subprocess failure.

### Scientific outcome

A completed run can produce:

- `RECOMMENDED`;
- `CANDIDATES_WITH_WARNINGS`;
- `NO_VALID_ASSAY`;
- `INSUFFICIENT_EVIDENCE`.

A run in which all candidates are rejected by CRITICAL off-target evidence is `COMPLETED / NO_VALID_ASSAY`, not a broken pipeline.

## Structural requirements for `RECOMMENDED`

Exact numeric thresholds remain configurable. Structurally, recommendation requires:

- an approved/frozen panel;
- required target design groups passing inclusivity policy;
- required independent target challenge evidence when configured and available;
- an acceptable hydrolysis probe;
- complete required specificity evidence;
- no CRITICAL detectable off-target hard fail;
- sufficient evidence for final ranking.

## Checkpoints and invalidation

The existing deterministic checkpoint/resume system remains the foundation.

Scientifically relevant changes invalidate descendants, including:

- panel manifest changes;
- design/challenge membership changes;
- DECIPHER design non-target membership changes;
- target or DECIPHER design alignment changes;
- primer-pair engine changes;
- DECIPHER/MAFFT tool identity or parameter changes;
- probe-constraint changes;
- specificity-policy changes.

Identical frozen manifests, datasets, configurations, and tool identities can reuse valid checkpoints.

## Proposed stage evolution

Current conceptual pipeline:

```text
input -> qc -> clustering -> alignment -> conservation
      -> primer_design -> inclusivity -> specificity -> ranking
```

Target architecture:

```text
panel
  -> input/acquisition
  -> qc
  -> dataset_partition
  -> clustering
  -> target_alignment
  -> conservation
  -> decipher_design_alignment (when required)
  -> primer_pair_design
  -> probe_design
  -> assay_evaluation
  -> challenge
  -> ranking
```

Binding invariants:

- panel approval precedes design;
- challenge data never enter optimization;
- target conservation alignment and DECIPHER design alignment are separate artifacts;
- primer-pair design and probe design are separate;
- engine selection applies to primer-pair design;
- final specificity uses the full approved non-target scope;
- checkpoints respect all these boundaries.

## Implementation decomposition

### 1. Panel model and manifests

Deliver:

- target modes;
- criticality;
- design/challenge roles;
- panel proposal/approval/freeze contract;
- manual West Nile fixtures;
- provenance/checkpoint integration.

Automatic panel intelligence is not required yet.

### 2. DECIPHER engine and common primer-pair contract

Deliver:

- DECIPHER detection in `doctor`;
- separate MAFFT DECIPHER-design-alignment preparation;
- Python/R adapter;
- structured I/O;
- `TileSeqs`, `DesignPrimers`, `CalculateEfficiencyPCR`;
- `auto` / `decipher` / `classic`;
- normalization into `PrimerPairCandidate`;
- classic-versus-DECIPHER comparison fixtures.

### 3. Separate probe and complete-assay design

Deliver:

- probe-design stage;
- multiple probes per primer pair;
- F + Probe + R assay model;
- assay-level inclusivity/specificity gates;
- reuse of existing off-target geometry where appropriate.

### 4. Panel intelligence and independent validation

Deliver:

- taxonomy-based proposals;
- versioned arbovirus clinical knowledge base;
- context inference;
- optional research suggestions with provenance;
- representative sequence selection;
- automatic design/challenge partitioning;
- final ranking/report changes;
- explicit scientific outcomes.

Each subproject receives its own implementation plan and validation before the next begins.

## Error handling and observability

The run record distinguishes at least:

- panel approval required;
- invalid panel manifest;
- DECIPHER unavailable;
- DECIPHER design-alignment failure;
- DECIPHER execution failure;
- invalid DECIPHER output;
- no eligible primer pairs;
- no acceptable probe for a pair;
- no valid assays after hard gates;
- insufficient challenge evidence.

The last three can be scientific/evidence outcomes rather than technical exceptions.

Run provenance preserves selected engine, fallback reason, tool versions, panel identity, dataset hashes, design/challenge membership, DECIPHER design-panel membership, effective thresholds, and alignment identities.

## Testing strategy

### Panel tests

Verify roles, criticality, canonical identifiers, provenance, approval/freeze behavior, deterministic serialization, and invalid-membership rejection.

### Design/challenge isolation

Verify challenge sequences cannot enter conservation, either design alignment, either primer-pair engine, or probe optimization.

### DECIPHER design-alignment tests

Verify:

- only target design + DECIPHER-compatible non-target design members enter the alignment;
- distant/full-specificity-only groups do not enter it;
- deterministic MAFFT artifacts and hashes;
- alignment/tool changes invalidate DECIPHER design descendants.

### DECIPHER adapter unit tests

Use fake/injected process boundaries or frozen structured output to verify deterministic workspace generation, group serialization, output parsing, version capture, error mapping, and independence from human-readable R console output.

### Real DECIPHER integration tests

Keep separate from ordinary fast unit tests. Use a small frozen fixture and local R/DECIPHER installation without live NCBI/network access.

### Engine-contract tests

Verify classic and DECIPHER both normalize to compatible `PrimerPairCandidate` objects while preserving engine-specific provenance.

### Probe tests

Verify multiple probes per pair, local pair rejection when no probe succeeds, and preservation of probe metrics/provenance.

### Hard-gate tests

Verify CRITICAL detectable off-target rejection, no averaging away a single CRITICAL failure, IMPORTANT not automatically becoming CRITICAL, and primer-only risk remaining distinct from detectable F+Probe+R off-target.

### Checkpoint tests

Verify invalidation for panel, partition, design-panel membership, engine, alignment/tool versions, DECIPHER parameters, probe constraints, and specificity policy. Verify identical reruns reuse checkpoints.

### End-to-end MVP fixture

Use a small frozen West Nile-centered fixture with target diversity and selected non-targets. The run is deterministic and offline.

Success means an explainable scientific outcome, including `NO_VALID_ASSAY` when appropriate.

## MVP acceptance criteria

Using frozen inputs, Geison can:

1. load an approved target/non-target manifest;
2. distinguish design from challenge data;
3. calculate target conservation using the existing path;
4. create a separate DECIPHER design alignment from target design + homologous non-target design data;
5. select primer pairs through DECIPHER;
6. normalize them through a common Python contract;
7. generate multiple hydrolysis-probe candidates per viable pair;
8. evaluate complete assays against design data;
9. challenge survivors against independent target and non-target data;
10. apply CRITICAL hard gates without hiding single-organism risk in averages;
11. report engine, versions, panel provenance, alignments, design/challenge evidence, rejection reasons, and scientific outcome;
12. reproduce the same result from the same frozen manifest, datasets, configuration, and tool versions.

The MVP does not require a `RECOMMENDED` assay.

## Scientific interpretation boundary

All results remain **in silico design evidence**. `RECOMMENDED` means recommended by the configured computational workflow, not clinically validated, analytically validated, or approved for diagnostic use.

Reports must preserve this distinction.

## Relationship to existing Geison designs

This architecture changes only assumptions that conflict with:

- explicit panel modeling;
- design/challenge separation;
- DECIPHER multi-engine primer-pair design;
- separate probe design;
- full-assay specificity;
- scientific-outcome semantics.

Existing NCBI, CD-HIT, MAFFT, conservation, inclusivity, specificity, ranking, checkpoint, and reproducibility designs remain implementation foundations unless a later subproject spec explicitly changes them.

The current Primer3 assay design remains the migration basis for the `classic` path and is not deleted by the initial DECIPHER work.

## Approved architectural decisions

- Python remains the Geison orchestrator.
- DECIPHER/R is optional and selected through `auto` when available/compatible.
- The DECIPHER MVP integrates `TileSeqs`, `DesignPrimers`, and `CalculateEfficiencyPCR` only.
- Existing target MAFFT alignment and Geison conservation remain in place initially.
- DECIPHER uses a separate MAFFT design alignment containing target design plus homologous/alignment-compatible non-target design data.
- `broad_detection` is the default; `subtype_specific` is supported.
- Non-target criticality is CRITICAL / IMPORTANT / BACKGROUND.
- New panels require human approval; frozen manifests rerun automatically.
- Target and non-target data are separated into design and independent challenge subsets.
- The first clinical knowledge-base content focuses on arboviruses while the model remains generic.
- DECIPHER receives only biologically meaningful alignment-compatible non-targets.
- Final specificity uses the full approved non-target panel.
- DECIPHER designs primer pairs; Geison separately designs hydrolysis probes.
- F + Probe + R is the unit of final specificity decision.
- A detectable CRITICAL off-target is a hard fail.
- Execution status and scientific outcome are separate.
- `NO_VALID_ASSAY` can be a successful completed run.
- `engine: decipher` never silently falls back; `engine: auto` may fall back only with explicit warning/provenance.
- Implementation proceeds as four sequential subprojects, not one monolithic change.
