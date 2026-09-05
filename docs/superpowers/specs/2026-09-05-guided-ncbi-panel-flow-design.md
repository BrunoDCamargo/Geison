# Guided NCBI Panel Workflow Design

## Goal

Make Geison's default Guided Colab workflow match the intended product experience: the researcher chooses a supported biological target, Geison prepares an explicit panel proposal, then — only after approval — Geison itself acquires target and challenge sequences from NCBI and runs the existing assay-discovery pipeline. Manual FASTA download/upload and path management disappear from the default path.

The current local FASTA/GenBank path remains available as an advanced workflow.

## Existing capability to reuse

Geison already has a mature NCBI acquisition subsystem for target inputs:

- `NcbiInputConfig` query/accession/frozen modes;
- `BioEntrezClient`;
- bounded retry/backoff for HTTP 429/5xx and network failures;
- batched acquisition;
- resumable/frozen `dataset_manifest.json` + `records.gb`;
- accession/version and SHA-256 provenance.

The Guided Colab currently bypasses this and constructs local FASTA configuration. Off-targets currently accept local FASTA or an already-frozen NCBI dataset, so the cleanest low-risk implementation is to make Guided mode materialize challenge datasets through the same existing NCBI service and then feed the ordinary pipeline its already-supported frozen off-target sources.

This restores the intended experience without changing contrast, specificity, ranking, or checkpoint semantics.

## User experience

The notebook exposes three modes:

1. **Guided (NCBI)** — default. User selects/types a supported target and provides the NCBI contact email required by the existing E-utilities adapter. No FASTA paths or uploads are requested.
2. **Demo (synthetic)** — deterministic synthetic walkthrough.
3. **Advanced (local sequences)** — current bring-your-own-FASTA workflow.

Guided flow:

```text
Target selection
      |
      v
Geison guided preset -> proposal config
      |
      v
PANEL_APPROVAL_REQUIRED
      |
      v
User reviews target + challenge organisms + criticality
      |
      v
APROVAR
      |
      v
Geison guided finalizer
  - validates approved panel
  - freezes challenge datasets from NCBI
  - writes approved pipeline config
      |
      v
qpcr-pipeline run
  - target NCBI acquisition through existing input stage
  - QC/alignment/conservation/contrast
  - primer/probe design/inclusivity/specificity/ranking/report
```

The only scientific human gate remains panel approval. NCBI contact identity is operational configuration, not a scientific decision.

## Supported guided knowledge

The architecture is generic, but the first versioned knowledge seed is deliberately small. It supports **West Nile virus** first because the current acceptance exercise uses WNV and the repository already has a WNV panel fixture.

Initial WNV proposal:

- Target: West Nile virus, `broad_detection`.
- CRITICAL challenge: Usutu virus — `phylogenetic_neighbor`.
- CRITICAL challenge: Japanese encephalitis virus — `phylogenetic_neighbor`.
- IMPORTANT challenge: Dengue virus — `clinical_differential`.

Every proposal records `proposed_by: geison_guided_knowledge` and a knowledge version. Unsupported targets fail clearly with the supported target list; Geison must not silently invent an authoritative panel.

## NCBI query policy

Guided datasets use the existing `NcbiInputConfig` contract. Initial operational defaults for the real WNV smoke workflow are:

- target query: `"West Nile virus"[Organism] AND complete genome[Title]`, `max_records: 50`;
- each challenge: organism-specific `complete genome[Title]`, `max_records: 20`;
- `batch_size: 20`;
- `retries: 5`.

These are bounded smoke/workbench defaults, not claims that the first N results form a scientifically representative population. The guided metadata records that limitation. Future representative-selection work can replace this policy without changing the notebook UX.

## NCBI contact identity

The current production adapter requires `NCBI_EMAIL`; this remains unchanged. The notebook collects the value once and exports `NCBI_EMAIL` before any live acquisition. `NCBI_API_KEY` stays optional via environment.

No email/API key is written to YAML, panel manifests, frozen dataset manifests, reports, logs, or evidence bundles.

## Guided module

Add `qpcr_pipeline.guided` as a deterministic, network-free knowledge/config layer.

It owns:

- versioned supported-target presets;
- target NCBI query/config;
- challenge NCBI query/config;
- panel proposal metadata;
- operational limits and explicit limitations.

Public interfaces:

```python
supported_guided_targets() -> tuple[str, ...]
build_guided_proposal_config(target_name: str) -> dict[str, object]
finalize_guided_project(
    target_name: str,
    approved_panel_path: str | Path,
    workspace: str | Path,
    *,
    ncbi_client: NcbiClient | None = None,
) -> Path
```

`build_guided_proposal_config()` makes no network call. It creates the normal pipeline proposal config with `input.ncbi` and `contrastive_conservation.enabled: false`.

`finalize_guided_project()` is called only after panel approval. It:

1. loads and validates the approved panel manifest;
2. verifies the approved target matches the guided preset;
3. maps each approved CHALLENGE non-target to an explicit preset NCBI query;
4. calls the existing `acquire_ncbi_dataset()` for each challenge under `<workspace>/guided_challenges/<NNN-slug>`;
5. writes `<workspace>/config-approved.yaml` using ordinary `off_targets[].frozen_dataset` entries and `input.ncbi` for the target;
6. writes a small `guided_acquisition_manifest.json` containing knowledge version, approved-panel SHA-256, dataset names, frozen paths, record counts and hashes, but no credentials.

If an approved challenge is not present in the preset, finalization fails rather than substituting or silently skipping it.

## CLI integration

Extend the existing CLI with two non-scientific orchestration commands:

```text
qpcr-pipeline guided prepare --target "West Nile virus" --workspace <dir>
qpcr-pipeline guided finalize --target "West Nile virus" --workspace <dir> --approved-panel <path>
```

`guided prepare` writes `config-proposal.yaml` and prints its path.

`guided finalize` performs approved challenge acquisition and writes `config-approved.yaml` plus `guided_acquisition_manifest.json`.

The existing commands remain responsible for scientific execution:

```text
qpcr-pipeline run config-proposal.yaml --outdir ...
qpcr-pipeline panel approve panel_proposal.yaml --output approved_panel.json
qpcr-pipeline run config-approved.yaml --outdir ... --resume
```

The notebook never implements NCBI HTTP calls itself.

## Notebook behavior

The first Python cell explicitly inserts `/content/Geison` into the active kernel `sys.path`. This closes the observed Colab editable-install/path-hook issue and makes the final report/evidence-bundle helper import reliable.

Default Guided (NCBI) mode shows only:

- target;
- NCBI contact email;
- workspace.

It does not show target FASTA, challenge names, challenge FASTAs, or upload instructions.

The proposal and approval gate remain visible. After `APROVAR`, the notebook calls `qpcr-pipeline guided finalize`, displays the generated approved config, then resumes the normal pipeline.

Advanced local mode retains the current manual fields.

## Provenance and reproducibility

Target NCBI provenance remains the existing pipeline-owned `ncbi_dataset_manifest.json`.

Challenge provenance is preserved by the existing NCBI frozen dataset manifests created under `guided_challenges/` plus `guided_acquisition_manifest.json` linking them to the approved panel and knowledge version. The final pipeline uses those frozen paths, so contrast and specificity consume the exact same challenge files.

No challenge network lookup occurs during contrast or specificity.

## Error behavior

- Missing `NCBI_EMAIL`: actionable failure before approved live acquisition.
- Unsupported guided target: fail before proposal with supported targets.
- NCBI 429/5xx/network failures: existing bounded retry/backoff applies.
- Empty/inconsistent NCBI dataset: existing acquisition validation fails explicitly.
- Approved challenge absent from preset: finalization fails; no silent omission.
- Panel approval remains mandatory before challenge acquisition and scientific execution.

## Testing

Offline tests cover:

- deterministic WNV proposal generation;
- unsupported target behavior;
- approved-panel mismatch/rejection;
- challenge acquisition through a deterministic fake `NcbiClient`;
- final approved YAML uses `input.ncbi` for target and `off_targets[].frozen_dataset` for challenges;
- acquisition manifest contains dataset hashes/counts and no credentials;
- CLI `guided prepare`/`guided finalize` contracts;
- notebook contract: Guided (NCBI) default, no manual FASTA requirement in that path, NCBI email exported, kernel source path inserted, existing approval/report sections preserved.

Live NCBI testing stays opt-in under `network_tests/` and is never required by the normal suite.

## Out of scope

- universal panel intelligence for arbitrary pathogens;
- literature/AI-based panel suggestions;
- metadata-stratified representative sampling;
- DECIPHER integration;
- report dark/light theme redesign;
- analytical or clinical validation claims.

This iteration is specifically about restoring the intended low-intrusion guided acquisition workflow with the NCBI infrastructure Geison already has.