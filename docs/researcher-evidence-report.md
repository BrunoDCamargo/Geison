# Researcher evidence report: initial migration

This read-only reporting slice implements part of #29, #30 and #32 under #28.
It does not change the scientific engines, search policies or candidate generation.

## Default route

`generate_researcher_report(output_dir)` now renders `report.html` through
`qpcr_pipeline.researcher_report_evidence`. The pipeline already calls this
public generator. A complete directory of published run artifacts can be rendered
again without rerunning the scientific stages:

```python
from pathlib import Path
from qpcr_pipeline.researcher_report import generate_researcher_report

report_path = generate_researcher_report(Path("path/to/published-run"))
print(report_path)
```

The writer replaces only `report.html`, using its existing atomic-write contract.
Copy a run directory first to keep the earlier HTML presentation. JSON and TSV
artifacts are not modified. An HTML attachment alone cannot recreate the original
per-dataset/window matrix or unrounded measurements.

## Presentation contract

- The run summary reports technical state, not scientific approval.
- Candidate order follows identifiers. No automatic recommendation, final score
  or quality class appears in the default presentation.
- Per-assay compatibility uses recorded numerator/denominator pairs from the
  assay-sequence matrix. When the matrix is absent, explicit legacy summary
  counts may be displayed with their source; fractions do not reconstruct counts.
- Aggregates are named assay-sequence pairs, not unique target sequences.
- Duplicate rows, conflicting counts, unknown identifiers and inconsistent
  declared populations cause a report error. Partial matrices do not silently
  borrow complete counts from the ranking artifact.
- Skipped or missing stages do not contribute stale measurements.
- Sequences, intervals, original penalties, policy codes, published observations,
  configuration and available provenance remain inspectable.

## Window views

The line view uses actual reference coordinates and recorded window centres.
Missing values break lines. The heatmap uses one column per recorded window,
in reference order, and one row per comparison dataset. Equal column widths
represent window index, not genomic distance or nucleotide resolution. Overlaps
are not accumulated as transparency. This first version does not synchronize
zoom/pan or add an oligo track; those requirements remain in #30.

The fixed [0, 1] color scale describes the recorded normalized local alignment
score, not percent identity or a probability. Nonfinite and out-of-range scores
are rejected rather than clamped. Missing cells, no recorded best match, and empty
datasets retain different states in the numeric appendix. The appendix includes
exact recorded values, counts, best sequence identifiers and orientations,
accessible without hovering or running JavaScript. No CDN dependency is added.

## Compatibility and remaining work

Legacy ranking artifacts, `scientific_outcome`, and
`researcher_report_html.render_researcher_report_html` remain unchanged for
callers that use the previous APIs directly. Alternative legacy reports can
therefore still contain classifications. This is a staged migration, not a claim
that every report or engine in the repository has changed.

#31, #33, #34 and the unresolved scientific-engine findings remain separate work.
No biological benchmark or experimental validation is implied by renderer tests.

## Verification

The full offline suite runs with `python -m pytest -q`. Integration tests run
with `python -m unittest discover -s integration_tests -q` after installing
CD-HIT, MAFFT and Primer3. They use deterministic synthetic fixtures, not
researcher data. The guided integration test checks the default evidence
presentation and verifies that regenerating HTML leaves JSON/TSV bytes unchanged.

Optional browser tests require Playwright and Chromium:

```bash
python -m pip install playwright
python -m playwright install --with-deps chromium
python -m unittest discover -s browser_tests -q
```

`GEISON_CHROMIUM_EXECUTABLE` can select an existing Chromium executable.
The tests load generated HTML into the browser DOM, check desktop and mobile
layout, and reject page errors or remote requests. No browser policy is changed.

The `Geison verification` GitHub Actions workflow runs the offline suite,
integration tests and browser tests before merge, and again on pushes to main.
Logs and the exact tested revision are stored as seven-day artifacts. CircleCI
remains configured for its existing main/develop jobs. Network-dataset tests are
not part of these offline jobs; biological benchmarks remain in #34.

The initial partial-checkout verification covered 37 tests only. A subsequent
full-source run covered 737 tests and exposed an integration assertion for the
superseded verdict UI. The migration updates that assertion to the descriptive
contract and adds regressions for known-zero versus missing counts, retained
correspondences, stale rule evidence and text identifiers in TSV input.
The current PR check results, rather than these historical counts, govern merge.

The author performed a separate self-review; no independent review is claimed.
No uploaded researcher report or its sequences were added to repository fixtures.
