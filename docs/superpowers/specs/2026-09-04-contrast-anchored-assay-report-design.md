# Contrast-anchored assay design and researcher report

**Date:** 2026-09-04  
**Status:** Approved design, pending implementation plan  
**Scope:** Geison qPCR assay-discovery pipeline and guided Colab workbench

## 1. Problem statement

Geison currently ranks conserved target regions against approved non-target challenge datasets before Primer3 design. The contrastive stage may identify a highly discriminant target window, but the selected window is expanded to a broader candidate region so Primer3 has enough sequence context for forward primer, probe, and reverse primer design.

In the guided synthetic demo, this exposed a scientific gap: the strongest discriminant window was approximately 601-700, while the corresponding candidate region was expanded to approximately 501-800. Primer3 was allowed to design anywhere inside that expanded interval and returned assays concentrated around 714-800. Those assays remained target-inclusive but no longer traversed the evidence that justified selection of the region. Final specificity then detected plausible off-target amplification and all assays were classified HIGH_RISK.

The pipeline behaved correctly at the final specificity stage, but the design boundary between contrastive region selection and Primer3 was too weak.

A second usability problem is that the guided Colab exposes many tables and charts, but the notebook itself is not a suitable final scientific deliverable. Researchers need a readable, self-contained report that summarizes the outcome and preserves the evidence required for review and reproducibility.

## 2. Goals

This change has two goals:

1. Make contrastive evidence an explicit constraint on oligo design so that assays designed from `CONTRASTIVE_CONSERVATION` must traverse the discriminant evidence that caused their candidate region to be selected.
2. Promote the existing final HTML from a ranking-focused artifact into a complete researcher-facing study report that is always generated when enough run evidence exists, including scientifically negative outcomes.

## 3. Non-goals

This change does not:

- replace final assay-level specificity with region-level contrast;
- guarantee an `IN SILICO PASS` candidate;
- convert a completed computational run into experimental or diagnostic validation;
- introduce wet-lab, clinical, or regulatory validation logic;
- change `CONSERVATION_ONLY` candidate semantics unless needed for compatibility with the new report interface;
- add PDF generation as a primary output.

## 4. Scientific design

### 4.1 Candidate region versus contrast anchor

A contrastive candidate must preserve two distinct coordinate concepts:

- **Candidate region** (`reference_start`, `reference_end`): the broader interval in which Primer3 is allowed to design an assay.
- **Contrast anchor** (`peak_start`, `peak_end`): the discriminant target window whose contrastive evidence justified selecting that candidate region.

The existing `CandidateRegion` already carries `peak_start` and `peak_end`. For candidates sourced from `CONTRASTIVE_CONSERVATION`, these coordinates become an active design constraint instead of audit-only metadata.

### 4.2 Primer3 constraint

For `CONTRASTIVE_CONSERVATION` candidates:

- continue to send the broader candidate interval as `SEQUENCE_INCLUDED_REGION`;
- send the contrast anchor as `SEQUENCE_TARGET`, converted to Primer3's zero-based coordinate convention;
- require every returned assay product to contain the complete contrast anchor;
- independently validate returned coordinates after parsing Primer3 output so an assay entirely outside the anchor cannot be accepted even if an upstream tool behaves unexpectedly.

The post-parse validation is defense in depth. Primer3 remains responsible for design; Geison remains responsible for enforcing Geison's scientific contract.

For `CONSERVATION_ONLY` candidates, do not fabricate a contrast anchor and preserve current design behavior.

### 4.3 No silent fallback

If a `CONTRASTIVE_CONSERVATION` result is `COMPLETE`, it remains authoritative.

If a contrastive candidate cannot produce a valid assay under the anchor constraint:

- that candidate produces no accepted assay;
- Geison may continue to other contrastive candidate regions;
- Geison must not silently fall back to conservation-only regions.

A complete contrastive result with zero viable assays remains a valid scientific outcome and must be represented as such in completeness/ranking/reporting.

### 4.4 Downstream stages remain independent

After constrained Primer3 design, the normal stages remain authoritative:

1. inclusivity / target coverage;
2. assay-level specificity against approved non-target datasets;
3. final classification and ranking.

A contrast-anchored assay may still be `REVIEW` or `HIGH_RISK`. Region-level contrast is not proof of oligo specificity.

## 5. Run status and scientific outcome

Technical run state and scientific interpretation must be separated.

| Technical state | Scientific report outcome |
| --- | --- |
| `COMPLETED` with at least one `IN SILICO PASS` | `In-silico candidate(s) identified` |
| `COMPLETED` with no pass but at least one `REVIEW` | `No in-silico pass; candidate(s) require review` |
| `COMPLETED` with all ranked assays `HIGH_RISK` | `No in-silico acceptable assay candidates identified` |
| `PARTIAL` | `Inconclusive - insufficient evidence` |
| `FAILED` | No conclusive scientific outcome; present execution failure only |

`COMPLETED` means the required computational evidence was produced. It must never be presented as synonymous with scientific acceptance.

The report must also support zero-assay runs. If the pipeline cannot satisfy scientific completeness because there are no assays, the report must state that the run is incomplete/inconclusive rather than presenting a negative assay conclusion.

## 6. Researcher-facing report

### 6.1 Output

The primary final artifact remains:

```text
output/report.html
```

It must be:

- single-file;
- static;
- self-contained;
- readable offline;
- escaped safely;
- suitable for browser printing to PDF;
- generated from published Geison artifacts rather than notebook-side scientific recalculation.

The report language is English.

### 6.2 Information hierarchy

The report has two layers:

1. **Summary** for rapid interpretation.
2. **Full evidence** for scientific review and reproducibility.

### 6.3 Report sections

#### A. Run summary

Show:

- target name;
- technical run status;
- explicit scientific outcome;
- target sequence count;
- challenge dataset count;
- contrastive candidate-region count;
- designed assay count;
- counts of `IN SILICO PASS`, `REVIEW`, and `HIGH_RISK` assays.

The scientific outcome must be visually more prominent than the aggregate score.

#### B. Scientific outcome

Use a concise interpretation generated from existing classified evidence.

Example negative outcome:

> **No in-silico acceptable assay candidates identified**  
> The computational workflow completed successfully, but all designed assays were classified as HIGH_RISK because detectable off-target amplification was identified.

The wording must derive from recorded classifications and reason codes, not free-form notebook logic.

#### C. Approved panel and study context

Show:

- target;
- DESIGN groups;
- approved CHALLENGE datasets;
- challenge criticality;
- available diagnostic context;
- approved/frozen panel provenance.

#### D. Target conservation

Show:

- target-side dataset summary;
- conservation visualization across reference coordinates;
- principal conserved regions or windows;
- concise explanation that conservation describes agreement inside the target population and does not establish specificity.

#### E. Target versus non-target contrast

Show:

- contrast overview plot;
- reference-position track;
- selected contrastive regions;
- candidate-region coordinates;
- contrast-anchor coordinates;
- contributing windows;
- worst challenge dataset and criticality;
- worst similarity metrics;
- contrast margin.

Visually distinguish the broader candidate region from its contrast anchor.

#### F. Assay design

For every designed assay show:

- assay ID;
- candidate region and source;
- contrast anchor when applicable;
- forward primer sequence and coordinates;
- probe sequence and coordinates;
- reverse primer sequence and coordinates;
- Tm;
- GC%;
- oligo penalties;
- pair penalty;
- product size;
- evidence that the final amplicon contains the required contrast anchor for contrastive candidates.

#### G. Target coverage / inclusivity

Show:

- evaluation-sequence count;
- original-compatible count;
- inclusivity fraction;
- per-assay target coverage summary;
- IUPAC/degeneracy proposals and their status when present.

#### H. Specificity

For each assay and challenge dataset show:

- compatible oligo hits;
- plausible off-target amplicons;
- detectable off-targets;
- the evidence contributing to the final classification.

#### I. Final candidates

Show:

- rank;
- classification;
- final score and score status;
- score components;
- reason codes and explanatory messages;
- recommended candidate only when at least one assay is `IN SILICO PASS`.

Do not label `REVIEW` or `HIGH_RISK` assays as recommended.

#### J. Interpretation and limitations

Always state that:

- results are in silico;
- region-level contrast is not equivalent to oligo specificity;
- `IN SILICO PASS` does not establish experimental or diagnostic validity;
- candidate assays require independent scientific review and experimental validation.

#### K. Reproducibility

Show, where recorded:

- Geison version and commit;
- run identifier and timestamps;
- effective configuration;
- approved panel hash/provenance;
- environment/tool versions;
- reference identifier/mode;
- artifact locations used to build the report.

Do not expose secrets or long sequence payloads from diagnostic logs.

## 7. Evidence bundle

The guided Colab must offer an evidence bundle in addition to the standalone HTML report.

The bundle is a ZIP created from already-published run artifacts and must include, when present:

- `report.html`;
- `run_manifest.json`;
- approved panel manifest;
- effective approved configuration;
- contrastive-conservation TSV/JSON/HTML artifacts;
- primer-design artifacts;
- inclusivity artifacts;
- specificity artifacts;
- ranking TSV/JSON artifacts;
- checkpoint/run metadata needed for audit.

The ZIP is a convenience packaging layer. It must not recalculate scientific results.

## 8. Guided Colab changes

Add a final researcher-facing section after the scientific result views:

```text
11. Researcher report
```

The section must provide three clear actions:

- **View report**
- **Download report.html**
- **Download evidence bundle.zip**

The notebook may preview or link the generated HTML, but report content must be generated by Geison itself.

The notebook remains a guided execution and inspection interface, not the system of record for scientific calculations.

## 9. Error behavior

### 9.1 Reporting failures

A report-generation failure must not silently overwrite or reinterpret scientific stage state.

If all scientific stages completed but final report rendering fails:

- preserve the underlying scientific artifacts;
- expose a clear report-generation error;
- do not claim that the researcher report was generated;
- do not change assay classification merely because presentation failed.

### 9.2 Partial and failed runs

When enough metadata exists, the report renderer may produce a diagnostic report for `PARTIAL` or `FAILED` states, but it must clearly label the outcome as inconclusive/execution failure and must not fabricate missing scientific sections.

Missing sections should say that evidence is unavailable and identify the relevant run state or missing evidence code.

## 10. Acceptance criteria

### 10.1 Contrast-anchor regression

Add a regression reproducing the defect observed in the guided demo:

- candidate region: approximately 501-800;
- contrast anchor: approximately 601-700;
- a Primer3 result whose entire amplicon is downstream of the anchor, for example approximately 714-800.

Expected result: the assay is rejected as violating the contrast-anchor contract.

### 10.2 Primer3 input

For a contrastive candidate:

- `SEQUENCE_INCLUDED_REGION` represents the broad candidate region;
- `SEQUENCE_TARGET` represents the contrast anchor with correct Primer3 coordinates;
- output parsing/validation confirms the resulting amplicon contains the complete anchor.

For a conservation-only candidate:

- no synthetic contrast target is added;
- existing behavior remains valid.

### 10.3 No fallback regression

A complete contrastive result that produces zero anchored assays must not fall back to conservation-only candidate regions.

### 10.4 Guided synthetic end-to-end test

The synthetic demo must test more than technical completion.

The E2E test must confirm:

- panel approval gate still functions;
- contrastive candidate generation still functions;
- Primer3 candidate source is `CONTRASTIVE_CONSERVATION`;
- every accepted contrastive assay contains its contrast anchor;
- target-coverage evaluation completes;
- specificity completes;
- final ranking is produced;
- the final scientific classification is coherent with the synthetic challenge design;
- `output/report.html` exists and contains the expected scientific outcome;
- the evidence bundle can be created from the run artifacts.

The synthetic fixture should contain enough sequence structure that at least one anchor-constrained assay can exploit the intentionally altered challenge interval and reach a meaningful final specificity outcome. The test must not weaken specificity rules merely to force a PASS.

### 10.5 Report state coverage

Test report rendering for at least:

- one or more `IN SILICO PASS` candidates;
- only `REVIEW` candidates;
- only `HIGH_RISK` candidates;
- zero assays / incomplete evidence;
- `PARTIAL` run state;
- `FAILED` run state.

### 10.6 Report integrity

Verify that the report:

- is valid self-contained HTML;
- does not depend on external JavaScript/CSS/network access;
- HTML-escapes user-provided names and identifiers;
- does not recalculate scientific metrics independently of pipeline artifacts;
- does not expose secrets or long raw sequences from diagnostic/provenance data;
- includes explicit limitations and the in-silico boundary.

## 11. Expected code areas

Implementation is expected to touch focused parts of the existing architecture rather than introduce a parallel pipeline. Likely areas include:

- `qpcr_pipeline/primer3.py` for Primer3 target serialization and returned-assay validation;
- `qpcr_pipeline/primer_design/` for candidate-source-aware constraint handling;
- `qpcr_pipeline/assay_report_html.py` or a small reporting composition layer for the expanded study report;
- `qpcr_pipeline/ranking.py` only where report assembly needs existing classified evidence;
- `qpcr_pipeline/run_recording.py` / provenance readers for run-level metadata;
- `notebooks/geison_guided_colab.ipynb` for report viewing/download actions;
- unit/regression tests and the guided E2E integration test.

Do not move scientific calculation into notebook code or HTML rendering.

## 12. Design principles

- Evidence before presentation.
- Technical completion and scientific acceptance are separate concepts.
- Contrastive selection must constrain, not merely influence, oligo design.
- Final assay specificity remains independent and authoritative.
- Negative scientific results are valid reportable outcomes.
- No silent fallback from complete contrastive evidence to conservation-only design.
- Researcher-facing outputs must remain auditable and reproducible.
- The final report is a view over published evidence, not a second implementation of the science.
