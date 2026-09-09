# Interactive Contrast Explorer

## Goal

Turn the existing self-contained contrastive-conservation HTML report into a read-only interactive explorer. The explorer must make dense evidence understandable without changing any scientific result, candidate region, assay, score, or pipeline checkpoint.

## Scope

This change is limited to `qpcr_pipeline/contrastive_report_html.py` and its focused tests. It reuses the report data already embedded by `render_contrastive_html`; no new scientific computation, persistence layer, network request, or external JavaScript dependency is introduced.

## Interaction design

The report keeps two synchronized views:

1. **Contrast overview** — target conservation versus non-target similarity. By default it shows one point per target window using the worst challenge similarity, avoiding the current 3,294-point cloud. Individual challenge datasets can be overlaid with checkboxes.
2. **Reference track** — target conservation and worst non-target similarity by reference position, with candidate regions drawn as clickable bands.

Controls:

- mouse wheel or `+` / `-` buttons zoom the reference axis;
- drag pans the visible reference interval;
- `Reset view` restores the full reference;
- one checkbox per challenge dataset shows or hides that dataset;
- `Show all` / `Hide all` provide quick dataset control;
- clicking a candidate band or candidate table row selects the region and displays its rank, coordinates, anchor/contributing windows, worst challenge, similarity, contrast margin, and conservation metrics;
- hovering a plotted point shows its window, dataset, similarity, target conservation, sequence id, and orientation when available.

The selection is exploratory only. There are no edit handles and no write-back to the pipeline.

## Data flow

Python continues to serialize `targetPoints` and `challengePoints`. It also serializes a compact `candidateRegions` array and a deterministic list of dataset names. All filtering, zooming, panning, hit-testing, and redraws happen in browser-local JavaScript against those arrays.

The HTML remains one offline file with no network dependency.

## Rendering choices

Use the existing Canvas implementation and browser APIs rather than adding Plotly, D3, or another dependency. The default view prioritizes signal over density:

- target/worst-challenge evidence is visible immediately;
- per-dataset evidence is opt-in;
- candidate regions remain visible on the reference track;
- selected regions are visually distinct and synchronized with the table/detail card.

## Safety and reproducibility

- All embedded data remains escaped through the existing JSON/HTML helpers.
- The generated HTML remains deterministic for identical inputs.
- Interaction never mutates source evidence.
- No external URL, script, font, or CDN is introduced.
- Existing static tables remain available as a complete textual evidence view.

## Testing

Extend `tests/test_contrastive_report_html.py` to verify that generated HTML:

- contains zoom, reset, and dataset-filter controls;
- embeds candidate-region data;
- keeps the existing self-contained/offline and escaping guarantees;
- remains deterministic;
- does not introduce edit controls or external URLs.

The existing pipeline/integration suite remains the regression gate because no scientific computation is changed.

## Out of scope

- dragging/resizing candidate regions;
- recalculating Primer3 from the browser;
- changing panel membership or criticality;
- NCBI acquisition preview;
- turning accepted IUPAC proposals into new assay candidates;
- changing ranking thresholds or scientific algorithms.

Those are separate changes so this explorer can be reviewed and validated independently.
