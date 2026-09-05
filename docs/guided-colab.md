# Guided Google Colab workbench

`notebooks/geison_guided_colab.ipynb` is the researcher-facing Colab flow for Geison. Scientific computation remains in the installed `qpcr-pipeline` CLI; notebook cells collect the small amount of user input needed, invoke official Geison commands, and render published YAML, TSV, JSON, and HTML artifacts.

The default experience is now **Guided (NCBI)**. A researcher selects a supported target, reviews the proposed scientific panel, approves it explicitly, and Geison acquires the target and challenge sequence datasets itself. Manual FASTA download/upload is not part of normal guided use.

The existing `notebooks/geison_colab.ipynb` remains available for lower-level operational validation.

## Modes

### Guided (NCBI) — recommended

The default mode starts from a target name. The first supported guided knowledge preset is **West Nile virus**.

The user supplies only:

- target name;
- workspace;
- an NCBI contact email, exposed to the runtime as `NCBI_EMAIL` before live acquisition.

Geison then creates the proposal configuration through:

```text
qpcr-pipeline guided prepare
```

The WNV preset uses the existing Geison NCBI acquisition subsystem. It proposes a reviewable challenge panel containing Usutu virus (`CRITICAL`), Japanese encephalitis virus (`CRITICAL`), and Dengue virus (`IMPORTANT`). This is a versioned guided knowledge seed, not an authoritative clinical panel.

The target is configured through `input.ncbi`; no target FASTA path is required. Before approval, challenge entries use only future frozen-dataset locations, so the panel proposal can be reviewed without downloading biological sequences.

After explicit approval, the notebook runs:

```text
qpcr-pipeline guided finalize
```

The finalizer uses Geison's existing NCBI acquisition code to download and freeze the approved challenge datasets under the workspace. It then writes a normal `config-approved.yaml` whose off-targets point to those frozen datasets. Contrast and specificity therefore consume the same immutable challenge evidence.

The current Guided workbench deliberately caps NCBI records to keep the Colab run bounded. Those caps are operational smoke/workbench limits and must not be interpreted as a representative-sampling strategy. Representative selection across lineage, geography, time, host, quality, and metadata remains separate scientific work.

### Demo (synthetic)

The notebook runs `examples/guided_demo/generate_demo_data.py` and creates a deterministic West Nile-like synthetic scenario. The sequences exist to exercise workflow behavior and reproducibility. They do not represent a validated diagnostic assay, a biological reference panel, or experimental evidence.

The demo contains target diversity plus `CRITICAL` and `IMPORTANT` synthetic challenge datasets. It includes target-stable regions with different challenge similarity so Geison can exercise contrastive conservation and downstream assay design reproducibly.

### Advanced (local sequences)

Advanced mode preserves bring-your-own-data behavior. It accepts one target FASTA and up to five challenge FASTA files, each with a name and criticality. The challenge name must match its panel entry because Geison resolves approved `CHALLENGE` entries against `off_targets` by normalized name.

Advanced mode validates local paths before panel approval so a missing file is surfaced before the researcher commits to the panel.

## Approval gate

Every assay-discovery mode uses the same human scientific gate. The proposal configuration contains `panel.proposal` and keeps `contrastive_conservation.enabled: false`.

The first run must stop with:

```text
ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED
```

In Guided (NCBI), this stop occurs before live target or challenge sequence acquisition. The researcher reviews target groups, challenge organisms, criticalities, and rationale, then types exactly:

```text
APROVAR
```

Approval runs `qpcr-pipeline panel approve` and creates the immutable `approved_panel.json`.

For Guided (NCBI), Geison then freezes approved challenge datasets through `guided finalize`, creates `config-approved.yaml`, and resumes the ordinary pipeline with `--resume`. Target acquisition happens through the existing `input.ncbi` stage. Demo and Advanced local modes continue to create the approved configuration without live challenge acquisition.

## NCBI acquisition behavior

Guided mode does not contain a second HTTP/download implementation in the notebook. It delegates acquisition to Geison's existing NCBI subsystem, which provides batched downloads, bounded retry/backoff for transient failures such as HTTP 429 and 5xx responses, resumable partial datasets, and frozen manifests with accession/version and SHA-256 provenance.

`NCBI_EMAIL` is required by live acquisition. `NCBI_API_KEY` remains optional through the environment. Neither value is written into configuration, panel manifests, reports, evidence bundles, or frozen-dataset manifests.

Guided challenge provenance is additionally linked in:

```text
guided_acquisition_manifest.json
```

This manifest records the guided knowledge version, approved-panel hash, challenge dataset paths, record counts, and data/manifest hashes.

## Evidence views

The workbench separates the scientific questions instead of presenting one aggregate score:

1. **Target conservation** reads target-side conservation artifacts.
2. **Target vs non-target contrast** reads `contrastive_conservation/window_metrics.tsv`, `dataset_metrics.tsv`, and `candidate_regions.tsv`.
3. **Assay design** reads Primer3 outputs and the recorded candidate source.
4. **Target coverage** reads inclusivity evidence.
5. **Specificity** reads assay-level off-target evidence and plausible amplicons.
6. **Final candidates** reads the published ranking artifacts.

Contrastive conservation prioritizes regions before oligo design. For candidates sourced from `CONTRASTIVE_CONSERVATION`, the broader candidate region gives Primer3 design space while the recorded contrast anchor is the discriminant interval that the final amplicon must contain. Final `specificity` remains an independent assay-level check after concrete primers and probes exist. Strong region-level contrast does not establish oligo specificity.

## Run states and scientific outcome

The notebook uses `run_manifest.json` as the authority for technical run state:

- `ACTION_REQUIRED`: scientific panel approval is still required.
- `PARTIAL`: inspect `missing_evidence`; do not present the run as complete.
- `FAILED`: inspect the failed stage and recorded diagnostic.
- `COMPLETED`: all evidence required by the run-completeness rules was recorded.

`COMPLETED` is not the same as an acceptable assay. A completed run can still contain only `REVIEW` or `HIGH_RISK` candidates. Final classification and reason codes remain the scientific interpretation of the recorded evidence.

## Researcher report

Assay-discovery runs publish a standalone English Researcher report at:

```text
output/report.html
```

The report is a static, self-contained view over Geison's published evidence: approved panel, conservation, target-versus-non-target contrast, assay design, target coverage, specificity, ranking, limitations, and reproducibility information.

The report separates technical completion from scientific outcome and does not turn an in-silico result into experimental or diagnostic validation.

Section **11. Researcher report** provides:

- **View report** — preview the published `report.html`;
- **Download report.html** — save the standalone report;
- **Download evidence bundle.zip** — save the report and allowlisted evidence together.

The notebook does not recalculate scientific metrics when producing these convenience outputs.

### Evidence bundle

The evidence bundle can contain, when present:

- `report.html`;
- `run_manifest.json`, `run_summary.json`, and `qc_report.json`;
- approved panel/config artifacts;
- conservation and contrastive-conservation evidence;
- primer-design, inclusivity, specificity, and ranking evidence;
- checkpoint manifests required for audit.

It uses an explicit allowlist rather than recursively packaging unrelated workspace files.

## Advanced evidence and reproducibility

The final Advanced section exposes generated configs, approved panel, Geison commit, environment/tool information, optional guided acquisition manifest, checkpoint manifests, `run_manifest.json`, and raw artifact paths.

Important contrastive artifacts live under:

```text
output/contrastive_conservation/window_metrics.tsv
output/contrastive_conservation/dataset_metrics.tsv
output/contrastive_conservation/candidate_regions.tsv
output/contrastive_conservation/contrastive_conservation_report.json
output/contrastive_conservation/report.html
```

Final assay evidence also includes the `primer_design`, `inclusivity`, `specificity`, and `ranking` directories.

For a clean restart, use a fresh Colab runtime/workspace and rerun from the top. In Guided (NCBI), there is no need to preserve or re-upload target/challenge FASTAs because Geison owns acquisition. For an offline reproducible resume, preserve `config-approved.yaml`, `approved_panel.json`, the guided frozen challenge directories when applicable, and the complete `output/` directory including `.checkpoints`.

## Scientific boundary

Guided (NCBI), Demo (synthetic), Advanced local sequences, and the Researcher report support in-silico assay discovery and evidence review. They do not replace experimental validation, wet-lab optimization, clinical validation, representative population sampling, or the biological judgment required to define a diagnostic panel. `IN SILICO PASS` is not experimental or diagnostic validation.
