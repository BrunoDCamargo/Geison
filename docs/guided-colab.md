# Guided Google Colab workbench

`notebooks/geison_guided_colab.ipynb` is the researcher-facing Colab flow for Geison. It keeps scientific computation in the installed `qpcr-pipeline` CLI and uses notebook cells to collect inputs, invoke commands, and render published YAML, TSV, JSON, and HTML artifacts.

The existing `notebooks/geison_colab.ipynb` remains available for lower-level operational validation.

## Modes

### Demo (synthetic)

The notebook runs `examples/guided_demo/generate_demo_data.py` and creates a deterministic West Nile-like synthetic scenario. The sequences exist to exercise workflow behavior and reproducibility. They do not represent a validated diagnostic assay, a biological reference panel, or experimental evidence.

The demo contains target diversity plus `CRITICAL` and `IMPORTANT` synthetic challenge datasets. It includes one target-stable region shared with challenges and another target-stable region deliberately altered in the challenges. Geison measures and ranks those regions through the normal contrastive-conservation implementation.

### Project

Project mode accepts one target FASTA and up to five challenge FASTA files. Each challenge needs a name and criticality. The challenge name must match its panel entry because Geison resolves approved `CHALLENGE` entries against `off_targets` by normalized name.

The notebook generates the YAML. Normal use does not require hand-editing the configuration.

## Approval gate

The first configuration contains `panel.proposal` and keeps `contrastive_conservation.enabled: false`.

Running it must stop with:

```text
ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED
```

No contrastive scientific execution should occur before approval. The notebook displays the generated panel proposal and requires the researcher to type exactly `APROVAR` after reviewing the target groups, challenge datasets, criticalities, and rationale.

Approval runs:

```text
qpcr-pipeline panel approve
```

The notebook then creates `approved_panel.json`, replaces the placeholder in `config-approved-template.yaml`, writes `config-approved.yaml`, enables contrastive conservation, and resumes the same output directory with `--resume`.

## Evidence views

The workbench separates the main questions instead of presenting one aggregate score:

1. **Target conservation** reads the target-side conservation artifacts.
2. **Target vs non-target contrast** reads `contrastive_conservation/window_metrics.tsv`, `dataset_metrics.tsv`, and `candidate_regions.tsv`.
3. **Assay design** reads Primer3 outputs and the recorded candidate source.
4. **Target coverage** reads inclusivity evidence.
5. **Specificity** reads assay-level off-target evidence and plausible amplicons.
6. **Final candidates** reads the published ranking artifacts.

Contrastive conservation prioritizes regions before oligo design. For candidates sourced from `CONTRASTIVE_CONSERVATION`, the broader candidate region gives Primer3 design space while the recorded contrast anchor is the discriminant interval that the final amplicon must contain. Final `specificity` remains an independent assay-level check after concrete primers and probes exist. Strong region-level contrast does not establish oligo specificity.

## Run states and scientific outcome

The notebook uses `run_manifest.json` as the authority for technical run state:

- `ACTION_REQUIRED`: review and approval are still required.
- `PARTIAL`: inspect `missing_evidence`; do not present the run as complete.
- `FAILED`: inspect the failed stage and recorded diagnostic.
- `COMPLETED`: all evidence required by the run-completeness rules was recorded.

`COMPLETED` is not the same as an acceptable assay. A completed run can still contain only `REVIEW` or `HIGH_RISK` candidates. The final classification and reason codes remain the scientific interpretation of the recorded evidence.

The presence of a chart or output file does not convert `PARTIAL` or `FAILED` into success.

## Researcher report

Assay-discovery runs publish a standalone English report at:

```text
output/report.html
```

The Researcher report is a static, self-contained, offline view over Geison's published evidence. It presents a readable summary first and then the full evidence chain: approved panel, conservation, target-versus-non-target contrast, assay design, target coverage, specificity, ranking, limitations, and reproducibility information.

The report separates technical completion from scientific outcome. It is generated for scientifically negative assay-discovery outcomes as well as successful ones. For example, a `COMPLETED` run in which every candidate is `HIGH_RISK` is presented as a completed computation with no in-silico acceptable candidate, not as a successful assay.

The HTML can be saved directly and can also be printed to PDF from a browser. PDF generation is not required to preserve the scientific record; `report.html` is the primary readable artifact.

In section **11. Researcher report**, the guided notebook provides three researcher-facing actions:

- **View report** previews the already-published `report.html`.
- **Download report.html** saves the standalone readable report.
- **Download evidence bundle.zip** saves the report together with the allowlisted run evidence and approved inputs.

The notebook does not recalculate scientific metrics when creating these outputs.

### Evidence bundle

The `evidence_bundle.zip` convenience package is assembled from already-published artifacts. It can contain, when present:

- `report.html`;
- `run_manifest.json`, `run_summary.json`, and `qc_report.json`;
- the approved panel artifact;
- conservation and contrastive-conservation evidence;
- primer-design, inclusivity, specificity, and ranking evidence;
- checkpoint manifests needed for audit;
- explicitly supplied approved input files such as `config-approved.yaml` under `inputs/`.

The evidence bundle uses an explicit allowlist. It does not recursively package unrelated workspace files and it does not recalculate scientific results.

## Advanced evidence

The final advanced section exposes the generated configs, approved panel, Geison commit, environment/tool information, checkpoint manifests, `run_manifest.json`, and raw artifact paths. Preserve these files when a run needs to be reviewed or reproduced.

Important contrastive artifacts live under:

```text
output/contrastive_conservation/window_metrics.tsv
output/contrastive_conservation/dataset_metrics.tsv
output/contrastive_conservation/candidate_regions.tsv
output/contrastive_conservation/contrastive_conservation_report.json
output/contrastive_conservation/report.html
```

Final assay evidence also includes the `primer_design`, `inclusivity`, `specificity`, and `ranking` directories.

## Clean restart and resume

For a clean restart, remove or choose a new workspace/output directory and rerun from the proposal stage.

For a reproducible resume in another Colab session, preserve at least:

- `config-approved.yaml`;
- `approved_panel.json`;
- the complete `output/` directory, including `.checkpoints`;
- the Geison commit recorded for the run.

Restoring only final tables does not restore checkpoint state.

## Scientific boundary

The guided workbench and Researcher report support in-silico assay discovery and evidence review. They do not replace experimental validation, wet-lab optimization, clinical validation, or the biological judgment needed to define a real diagnostic panel. `IN SILICO PASS` is not experimental or diagnostic validation.
