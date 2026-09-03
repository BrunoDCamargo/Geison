# Contrastive Conservation Checkpoint Clarification

Date: 2026-09-03
Status: Consistency clarification for the approved contrastive-conservation design

This note resolves one dependency contradiction discovered while writing the implementation plan for `2026-09-03-contrastive-conservation-guided-colab-design.md`.

## Problem

The design describes the visible stage order as:

`panel -> input -> qc -> clustering -> alignment -> conservation -> contrastive_conservation -> primer_design -> inclusivity -> specificity -> ranking`

but also requires this invalidation behavior:

- changing non-target criticality or `CHALLENGE` membership invalidates contrastive conservation and its descendants;
- target acquisition, QC, clustering, alignment and conservation remain reusable when only challenge-panel semantics change.

If the `input` checkpoint directly depends on the `panel` checkpoint fingerprint, any approved-panel change invalidates the entire target-side chain. That contradicts the explicit invalidation requirement.

## Resolution

Execution order and checkpoint dependency are separated.

The panel approval preflight still occurs before any scientific execution. A proposal still returns `ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED` and stops before target input acquisition.

After an approved frozen panel exists, the checkpoint graph is:

```text
panel ───────────────────────────────┐
                                     v
input -> qc -> clustering -> alignment -> conservation -> contrastive_conservation -> primer_design -> inclusivity -> specificity -> ranking
```

with these exact dependency rules:

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

This preserves all approved behavior:

- panel approval remains mandatory before execution when contrastive analysis is enabled;
- panel identity/criticality/CHALLENGE membership participates in the contrastive fingerprint;
- target-side acquisition and conservation do not depend on challenge-panel semantics;
- changing the target input still invalidates the entire target-side chain and then contrastive descendants;
- changing only the panel challenge semantics invalidates `panel`, `contrastive_conservation`, and descendants while reusing `input` through `conservation`;
- changing only a challenge dataset invalidates contrastive descendants while reusing the target-side chain.

## Proposal-to-approved configuration transition

A proposal is not an approved frozen panel. Therefore a configuration containing `panel.proposal` must not enable `contrastive_conservation` yet.

The guided flow uses two generated configurations:

1. `config-proposal.yaml`
   - contains `panel.proposal`;
   - keeps `contrastive_conservation.enabled: false`;
   - exists only to produce the review artifact and `ACTION_REQUIRED` gate.
2. `config-approved.yaml`
   - replaces the proposal with `panel.frozen_manifest`;
   - enables `contrastive_conservation` and the downstream scientific stages selected by the user;
   - resumes the run after explicit approval.

This is a state transition, not a hidden configuration mutation. Both generated files are preserved for auditability.

## Result-type provenance clarification

`ContrastiveConservationResult` should preserve typed challenge-dataset provenance, not only dataset names. The result therefore carries a tuple of summaries containing at least:

- dataset name;
- criticality;
- source type;
- record SHA-256;
- sequence count.

The report JSON may include richer source-manifest details, but the checkpointed typed result must retain enough provenance to explain and render a resumed run without re-reading notebook state.
