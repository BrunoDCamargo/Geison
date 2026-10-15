# Researcher evidence report implementation plan

> Execute inline using superpowers:executing-plans and test-driven-development.

Goal: replace the default final report's automatic scientific verdict with inspectable evidence and a per-dataset window heatmap. This is the first slice of the approved maturation review, not an algorithm validation.

Authority: user-approved review of 2026-10-15; issues #28, #29, #30, #32.

Architecture: keep ResearcherReportData and the atomic writer. Add a standard-library-only neutral renderer and numeric visualization helper. Redirect only generate_researcher_report to the neutral renderer. Legacy rendering APIs and ranking artifacts remain available unchanged for compatibility, explicitly documented. No modification of scientific computations.

Constraints: no remote resources; no private uploaded data in git; no inferred zeros, counts, precision, or scientific outcomes; no overwrite of input artifacts; no new production dependencies; no merge into main without regression evidence.

## Task 1: Neutral report
- [ ] Write synthetic tests for no default verdict, raw per-assay counts, no ranking, missing/partial evidence, escaping, retention details, and immutability.
- [ ] Observe the tests fail before adding the renderer.
- [ ] Add qpcr_pipeline/researcher_report_evidence.py with render_researcher_report_html(data) -> str; render only recorded values and explicit derivations.
- [ ] Redirect the atomic generator in qpcr_pipeline/researcher_report.py; keep scientific_outcome as a documented legacy API.
- [ ] Verify tests, source hashes, input preservation and atomic failure behavior.

## Task 2: Window matrix and curves
- [ ] Test no-match versus absent versus unevaluated, numeric domain, duplicates, single/irregular windows, and all-dataset alignment.
- [ ] Add qpcr_pipeline/report_evidence_charts.py. Consume contrastive JSON dataset_evidence and windows, already published by baseline; do not infer per-dataset values from the maximum curve.
- [ ] Use equal-width columns explicitly labeled as windows; preserve true reference intervals in headings/tooltips. Curves use reference coordinates. Missing values break curves. Color is a fixed numeric encoding, not stacked opacity.
- [ ] Provide a focusable scroll region, readable numeric cells, and a selectable window subset through native browser find/scroll. Interactive zoom/pan and oligo tracks remain outstanding under #30.
- [ ] Run numeric and HTML tests and inspect a browser rendering.

## Task 3: Publish a reviewable slice
- [ ] Run all tests available in the isolated snapshot and request the repository CI on the PR. Do not call a partial snapshot a full checkout or its tests a full suite.
- [ ] Self-review the diff, security boundaries and default entrypoint; report absence of an independent reviewer.
- [ ] Commit through the connector on the isolated branch and open a draft PR. Keep issues open until their complete criteria and regression gate are met.

## Review focus
Unknown/invalid data must not become zero. Missing sequence IDs must not inflate a denominator. Skipped stages must not show stale measures. Legacy ranking alone must not make a design recommendation. Dynamic strings must not become markup. HTML uploaded by the researcher does not contain the original matrix and will not be reconstructed as if it did.

## Environment
The connector can read/write GitHub, but git clone from the execution environment fails at DNS resolution. The local snapshot includes the fetched report model, verified against its remote Git blob SHA; this does not establish a complete baseline suite. Preserve untouched repository files by creating a Git tree based on the known remote tree.
