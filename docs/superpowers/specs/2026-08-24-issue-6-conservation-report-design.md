# Issue #6: Conservation analysis and interactive genome report

## Goal

Calculate deterministic per-column conservation metrics from the completed
Discovery Set alignment, derive reference-coordinate consensus sequences and
sliding-window summaries, incorporate reference annotations when available, and
publish a self-contained interactive `report.html` in which conserved peaks can be
identified with zoom and hover.

## Scope and non-goals

This issue adds conservation configuration, a pure scientific analysis service,
immutable result models, TSV/FASTA/JSON artifacts, GenBank reference-annotation
extraction, a dependency-free interactive HTML report, pipeline integration, and
unit plus browser verification.

It does not select formal candidate regions, call Primer3, design assays, evaluate
inclusivity or specificity, rank final assays, add progress/checkpoints, or build
the desktop UI. The high-conservation table in the HTML is an explanatory view of
window metrics, not the candidate-region contract owned by issue #7.

The handoff document is historical context only. GitHub issue #6, this approved
design, and the repository are the binding sources for this work.

## Architecture choice

The scientific core remains pure Python and the report uses an inline HTML/CSS/
JavaScript Canvas implementation. No plotting library, CDN, browser service, or
network access is required to generate or open the report.

This keeps the artifact portable across local Windows runs, Google Colab, and the
future PySide6 desktop shell. Canvas is preferred over one SVG node per window so
larger reference sequences remain usable. The complete scientific data lives in
TSV/FASTA/JSON; the HTML embeds compact window and annotation data for display.

## Configuration

`PipelineConfig` gains a frozen, slotted configuration:

```python
@dataclass(frozen=True, slots=True)
class ConservationConfig:
    enabled: bool = False
    window_size: int = 100
    step_size: int = 10
```

YAML uses an optional top-level section:

```yaml
conservation:
  enabled: true
  window_size: 100
  step_size: 10
```

Validation is shared by YAML-loaded and directly constructed configurations:

- `enabled` is a boolean;
- `window_size` and `step_size` are integers from 1 through 1,000,000,
  excluding booleans;
- `step_size` cannot exceed `window_size`;
- only `enabled`, `window_size`, and `step_size` are accepted fields;
- enabled conservation requires enabled alignment and is rejected during pipeline
  configuration validation before the output directory is created.

Conservation defaults to disabled so existing workflows remain compatible. A
disabled run never requires a completed alignment, removes only exact stale
issue-#6 data artifacts, and atomically publishes a `SKIPPED` conservation report
last.

## Public service and immutable models

A new `qpcr_pipeline.conservation` module owns calculation, annotation extraction,
windowing, and artifact publication:

```python
analyze_conservation(
    records: tuple[LocalSequenceRecord, ...],
    alignment: AlignmentResult,
    config: ConservationConfig,
    output_dir: Path,
    *,
    target_name: str,
) -> ConservationResult
```

Public result models are frozen and slotted:

```python
@dataclass(frozen=True, slots=True)
class PositionConservation:
    alignment_position: int
    reference_position: int | None
    reference_base: str | None
    depth: int
    coverage: float
    frequency_a: float
    frequency_c: float
    frequency_g: float
    frequency_t: float
    gap_frequency: float
    major_allele_frequency: float
    entropy_bits: float
    major_consensus: str
    iupac_consensus: str

@dataclass(frozen=True, slots=True)
class WindowConservation:
    reference_start: int
    reference_end: int
    position_count: int
    mean_conservation: float
    minimum_conservation: float
    mean_coverage: float
    mean_gap_frequency: float
    mean_entropy_bits: float

@dataclass(frozen=True, slots=True)
class ReferenceAnnotation:
    feature_type: str
    start: int
    end: int
    strand: int | None
    label: str

@dataclass(frozen=True, slots=True)
class ConservationResult:
    status: Literal["SKIPPED", "COMPLETE"]
    reference_id: str | None
    positions: tuple[PositionConservation, ...]
    windows: tuple[WindowConservation, ...]
    annotations: tuple[ReferenceAnnotation, ...]
    major_consensus: str
    iupac_consensus: str
    position_metrics_path: Path | None
    window_metrics_path: Path | None
    major_consensus_path: Path | None
    iupac_consensus_path: Path | None
    html_report_path: Path | None
    report_path: Path
```

Enabled analysis requires a `COMPLETE` `AlignmentResult`. Its sequence IDs must
equal the unchanged Discovery Set exactly once and in the same order. The supplied
records must contain exactly those IDs so the selected reference record and its
metadata are unambiguous. Invalid inputs fail before artifact mutation.

## Per-column metrics

Metrics are calculated for every alignment column, including columns where the
reference has a gap. Let `N` be the Discovery Set sequence count, `g` the number
of gap characters in the column, and `d = N - g` the non-gap depth:

- `depth = d`;
- `coverage = d / N`;
- `gap_frequency = g / N`;
- A/C/G/T frequencies are normalized among the `d` non-gap observations;
- `major_allele_frequency` is the maximum A/C/G/T frequency;
- `entropy_bits = -sum(p * log2(p))` across positive A/C/G/T frequencies.

For enabled empty alignments, all position and window collections are empty and no
division occurs. An all-gap column in a non-empty alignment is rejected.

### IUPAC fractional observations

Canonical bases contribute one full observation. An IUPAC ambiguity contributes
one observation divided equally among the canonical bases it represents. Examples:

- `R` contributes 0.5 A and 0.5 G;
- `B` contributes one third each to C, G, and T;
- `N` contributes 0.25 to every canonical base.

This makes the four normalized frequencies sum to one whenever depth is positive,
without silently discarding ambiguity. Gaps never contribute to base frequencies
or entropy.

### Consensus rules

The majority consensus chooses the canonical base with maximum frequency. Ties
prefer the canonical reference base when it is among the tied bases; remaining
ties use stable `A`, `C`, `G`, `T` order.

The IUPAC consensus encodes the exact set of canonical bases with positive
frequency using the standard inverse IUPAC ambiguity table. Consequently, an `N`
observation legitimately broadens that column's IUPAC consensus to `N`.

Per-column metrics preserve reference-gap columns. The published majority and
IUPAC consensus FASTA sequences omit reference-gap columns so both sequences are
ungapped, reference-coordinate templates suitable for issue #7 and Primer3.

## Reference coordinates and windows

Only positions with a non-null `reference_position` participate in genomic
windows. Their `major_allele_frequency` is the position's conservation score;
coverage, gap frequency, and entropy remain separate dimensions so downstream
ranking does not double-count them.

For reference length `L`:

- if `L` is zero, no windows are produced;
- if `L <= window_size`, one partial window covers `1..L`;
- otherwise, full windows start at 1 and advance by `step_size` while they fit;
- if the last regular window does not end at `L`, one final full window is anchored
  at `L - window_size + 1`, unless that start already exists.

Every window records reference start/end, position count, mean and minimum
conservation, mean coverage, mean gap frequency, and mean entropy in bits. This
ensures both ends of the genome are represented without a cascade of misleading
tiny trailing windows.

## Reference annotations

The selected reference's `LocalSequenceRecord.metadata["features"]`, already
preserved from GenBank by the local/NCBI loaders, is the annotation source. FASTA
references simply produce no annotations.

For each supported local `SeqFeature`:

- `source` features are omitted from the plot because they normally span the
  entire record and do not identify a local biological region;
- each local location part becomes one `ReferenceAnnotation`;
- Biopython zero-based, end-exclusive coordinates become 1-based inclusive
  `start = int(part.start) + 1`, `end = int(part.end)`;
- coordinates are clipped to the reference bounds and invalid/empty parts are
  skipped;
- external-reference parts are skipped;
- strand is `-1`, `0`, `1`, or null;
- label preference is `gene`, `locus_tag`, `product`, then feature type;
- annotations are ordered by start, end, type, and label.

The JSON report records published and skipped annotation counts. Compound feature
parts are deliberately published as separate spans with the same label, which
preserves visible genomic intervals without inventing continuity across joins.

## Artifacts and schemas

Every run publishes:

```text
conservation/
  conservation_report.json
```

Enabled runs additionally publish:

```text
conservation/
  position_metrics.tsv
  window_metrics.tsv
  consensus_major.fasta
  consensus_iupac.fasta
report.html
```

`position_metrics.tsv` has exact columns matching the
`PositionConservation` fields. Null reference coordinates/bases are blank. Float
values are serialized with stable round-trip-safe decimal rendering.

`window_metrics.tsv` has exact columns matching `WindowConservation` fields.

Consensus FASTA IDs are `geison-major-consensus` and
`geison-iupac-consensus`. Enabled empty analyses publish empty FASTA files rather
than header-only records with empty biological sequences.

The schema-version-1 JSON report records:

- status and whether conservation was enabled;
- effective window parameters and metric definitions;
- reference ID and Discovery Set IDs;
- counts for sequences, alignment columns, reference positions, windows, and
  annotations;
- annotation published/skipped counts;
- consensus lengths;
- relative artifact paths or null when skipped.

All final files use unique temporary siblings and atomic replacement. Before
mutating any fixed-name data artifact, publication removes the exact prior
`conservation_report.json` so a failure cannot leave an old `COMPLETE` report
beside mixed artifacts. Data artifacts and `report.html` are published first;
`conservation_report.json` is published last.

Disabled publication removes only the exact position/window TSVs, two consensus
FASTA files, and root `report.html`, preserving unrelated siblings. It then
publishes the `SKIPPED` report last.

## Interactive HTML report

`report.html` is deterministic and self-contained:

- no CDN, external font, remote image, module import, or network request;
- inline CSS and JavaScript only;
- compact JSON arrays for windows and annotations embedded with closing-script
  sequences escaped;
- untrusted target, reference, and annotation text is inserted through
  `textContent`, never interpreted as markup.

The report includes:

- analysis identity, reference, Discovery count, and window parameters;
- a Canvas genome plot with conservation as the primary line, coverage as a
  secondary line, and a reference annotation lane;
- wheel zoom centered on the pointer, drag pan, and a reset-zoom control;
- hover nearest-window details with reference interval, mean/minimum conservation,
  coverage, gaps, entropy, and overlapping annotations;
- a deterministic top-ten high-conservation window table ordered by mean
  conservation descending, minimum conservation descending, coverage descending,
  entropy ascending, then reference start;
- clicking a top-window row zooms to that interval;
- an accessible empty-state message when no windows exist;
- textual instructions and a non-Canvas summary so core results remain readable.

The HTML does not define or persist candidate regions. It helps a user interpret
peaks; issue #7 owns scientific candidate selection and ranking.

## Pipeline integration

After `align_discovery`, `run_pipeline` calls `analyze_conservation` with the
Discovery records, alignment result, configuration, output directory, and target
name.

`RunSummary`, Evaluation Set, Discovery Set, clustering artifacts, and alignment
artifacts keep their current semantics. `qc_report.json` gains exactly:

```json
"conservation": {
  "status": "COMPLETE",
  "reference_id": "seq-3",
  "position_count": 1234,
  "window_count": 115
}
```

For disabled analysis the counts are zero and reference ID is null. Invalid direct
configuration and the enabled-conservation/disabled-alignment mismatch fail before
the output directory is created.

## Error handling

- enabled conservation with a non-COMPLETE alignment fails before artifacts;
- sequence membership, order, alignment length, coordinate, reference, and record
  mismatches fail explicitly;
- a non-empty alignment with an all-gap column fails explicitly;
- unsupported alignment symbols fail explicitly even if an injected result bypassed
  the alignment parser;
- malformed reference features are skipped and counted rather than failing the
  scientific calculation;
- artifact write/replace failures propagate after invalidating the prior report;
- report generation never fetches remote resources.

## Testing and verification

The standard suite remains deterministic and offline:

- configuration tests cover defaults, YAML, unknown fields, types, bounds,
  step/window relationships, direct construction, and the alignment dependency;
- metric tests cover exact canonical frequencies, every IUPAC expansion, gaps,
  coverage, major frequency, tie-breaking, entropy values 0, 1, and 2 bits,
  reference-gap columns, and all-gap rejection;
- consensus tests prove reference-coordinate omission of insertion columns and
  deterministic FASTA IDs;
- window tests cover short references, regular stepping, anchored final windows,
  means/minima, and empty analysis;
- annotation tests cover simple and compound locations, coordinate conversion,
  label preference, source/external/malformed skipping, and deterministic order;
- artifact tests cover schemas, safe HTML embedding, hardlink-safe replacement,
  disabled stale cleanup, report-last ordering, and prior-report invalidation at
  every data replacement failure point;
- pipeline tests prove only Discovery alignment data is consumed, old set/summary
  contracts remain unchanged, disabled runs publish SKIPPED without new
  dependencies, and top-level traceability is present;
- an in-app browser verification opens a generated fixture report and exercises
  page load, reset zoom, wheel zoom, hover details, and top-window navigation while
  checking for console errors.

## Acceptance mapping

- Coverage, A/C/G/T frequencies, gap frequency, major allele frequency, entropy:
  position metrics and immutable models.
- Majority and IUPAC consensus: deterministic reference-coordinate FASTA artifacts.
- Configurable sliding windows: configurable size/step and persisted mean/minimum
  conservation plus coverage metrics.
- Interpretable genomic peaks: Canvas conservation plot and ranked explanatory
  table.
- Zoom and hover: offline pointer-centered zoom, pan, reset, and detailed tooltips.
- Reference annotations: normalized GenBank `SeqFeature` spans in report JSON and
  the annotation lane.
- Known-calculation tests: explicit frequency, IUPAC, entropy, consensus, and window
  fixtures.
