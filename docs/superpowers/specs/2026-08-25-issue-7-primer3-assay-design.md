# Issue #7: Candidate regions and Primer3 hydrolysis-probe assays

## Goal

Transform the reference-coordinate conservation windows produced by issue #6
into a bounded, deterministic set of candidate regions, then ask the official
`primer3_core` executable for multiple forward-primer, hydrolysis-probe, and
reverse-primer assays per region while preserving the scientific and Primer3
metrics needed by later evaluation stages.

## Scope and non-goals

This issue adds candidate-region configuration and ranking, a safe Boulder-IO
Primer3 boundary, immutable candidate/oligo/assay results, TSV/JSON/raw artifacts,
pipeline integration, and deterministic unit plus conditional real-binary tests.

It does not evaluate inclusivity against the Discovery Set, evaluate specificity
against the Evaluation Set, assign final risk classifications, annotate assay
overlays in `report.html`, package Primer3, or build the desktop UI. Those concerns
belong to issues #8, #9, #10, #15, and #17. The existing issue-#6 HTML remains
unchanged until the consolidated report work in issue #10.

The handoff document is historical context only. GitHub issue #7, this approved
design, and the repository are the binding sources for this work.

## Architecture choice

Candidate selection remains pure Python. Primer design uses the official
`primer3_core` executable through its Boulder-IO stdin/stdout protocol, matching
the repository's existing external-runner pattern for CD-HIT and MAFFT. The runner
never invokes a shell and accepts an injectable executable name for tests and
future WSL2 packaging. No `primer3-py` or other compiled Python dependency is
added.

One Primer3 process receives all selected candidate records in deterministic
rank order. Each record uses the complete issue-#6 majority consensus as
`SEQUENCE_TEMPLATE` and limits permitted design space with
`SEQUENCE_INCLUDED_REGION`. This satisfies the template requirement while keeping
all Primer3 coordinates directly aligned to the 1-based reference coordinate
system after explicit conversion.

## Configuration

`PipelineConfig` gains frozen, slotted configuration models:

```python
@dataclass(frozen=True, slots=True)
class OligoConstraints:
    min_size: int
    opt_size: int
    max_size: int
    min_tm: float
    opt_tm: float
    max_tm: float
    min_gc_percent: float
    max_gc_percent: float

@dataclass(frozen=True, slots=True)
class PrimerDesignConfig:
    enabled: bool = False
    max_candidate_regions: int = 10
    assays_per_region: int = 5
    candidate_region_length: int = 300
    max_region_overlap_fraction: float = 0.5
    min_mean_conservation: float = 0.90
    min_minimum_conservation: float = 0.70
    min_mean_coverage: float = 0.90
    max_mean_gap_frequency: float = 0.05
    max_mean_entropy_bits: float = 0.50
    min_usable_fraction: float = 0.80
    product_size_min: int = 70
    product_size_max: int = 200
    primer: OligoConstraints = PRIMER_DEFAULTS
    probe: OligoConstraints = PROBE_DEFAULTS
```

The immutable defaults are:

- primer size `18 / 20 / 25`, Tm `58 / 60 / 62` degrees Celsius, and GC
  `40..60` percent;
- probe size `18 / 25 / 30`, Tm `68 / 70 / 72` degrees Celsius, and GC
  `30..80` percent.

YAML uses an optional `primer_design` section with nested `primer` and `probe`
mappings. Every field above is configurable. Unknown fields are rejected at every
mapping level. Directly constructed and YAML-loaded configurations share the same
validation:

- booleans are accepted only for `enabled`;
- integer counts and lengths exclude booleans and must be positive;
- fractions and GC percentages are finite and inside their natural inclusive
  ranges;
- minimum, optimum, and maximum sizes/Tm values must be ordered;
- product-size minimum cannot exceed its maximum;
- `candidate_region_length` must be at least `product_size_max`;
- enabled primer design requires enabled conservation, which already requires
  enabled alignment.

Primer design defaults to disabled so all existing workflows remain compatible.
The default values form an initial configurable qPCR hydrolysis-probe preset, not
a claim that one preset is optimal for every organism or laboratory protocol.

## Candidate-region model and selection

```python
@dataclass(frozen=True, slots=True)
class CandidateRegion:
    region_id: str
    rank: int
    reference_start: int
    reference_end: int
    peak_start: int
    peak_end: int
    position_count: int
    usable_length: int
    usable_fraction: float
    mean_conservation: float
    minimum_conservation: float
    mean_coverage: float
    mean_gap_frequency: float
    mean_entropy_bits: float
```

Selection consumes only a `COMPLETE` `ConservationResult`. Reference-coordinate
positions must be contiguous from 1 through the majority-consensus length, and
every conservation window must lie within those bounds.

For each conservation window:

1. Center a `candidate_region_length` interval on the window midpoint.
2. Shift, rather than shorten, intervals that meet either reference boundary when
   the reference is at least the configured length. A shorter reference produces
   one interval spanning the complete reference.
3. Deduplicate identical intervals created by adjacent windows.
4. Recalculate all aggregate metrics from issue-#6 per-position data across the
   expanded interval.
5. Count a position as usable only when its major-allele frequency, coverage, gap
   frequency, and entropy each meet the corresponding configured threshold.
6. Keep an interval only when its aggregate metrics and usable fraction meet all
   configured thresholds.

Eligible intervals use a transparent lexicographic ordering instead of an
arbitrary weighted score:

1. mean conservation descending;
2. minimum conservation descending;
3. mean coverage descending;
4. mean entropy ascending;
5. mean gap frequency ascending;
6. usable length descending;
7. reference start ascending;
8. reference end ascending.

A greedy pass accepts that order while rejecting a region whose overlap divided
by the shorter region's length is greater than
`max_region_overlap_fraction`. It stops at `max_candidate_regions`. Final stable
IDs and ranks are `region-001`, `region-002`, and so on. The originating
conservation window is retained as `peak_start` and `peak_end` for traceability;
when multiple windows generated an identical interval, the best-ranked window by
the same metric ordering is retained.

An empty conservation result or a result with no eligible regions is a valid
complete design result with zero candidates and zero assays.

## Primer3 request boundary

The public runner protocol and default implementation are:

```python
class Primer3Runner(Protocol):
    def run(self, input_text: str) -> str: ...

class SubprocessPrimer3Runner:
    def __init__(self, executable: str = "primer3_core") -> None: ...
    def run(self, input_text: str) -> str: ...
```

The subprocess command is fixed to `primer3_core --strict_tags --io_version=4`
apart from the injected executable. Input is sent through stdin; stdout and stderr
are captured; `shell=False` is mandatory. A missing executable, nonzero exit,
invalid UTF-8, or empty output raises a domain-specific `PrimerDesignError` whose
message contains actionable context but never includes the full consensus.

Each candidate becomes one Boulder-IO record containing at least:

- `SEQUENCE_ID`, `SEQUENCE_TEMPLATE`, and zero-based
  `SEQUENCE_INCLUDED_REGION`;
- `PRIMER_TASK=generic`, `PRIMER_FIRST_BASE_INDEX=0`, and
  `PRIMER_EXPLAIN_FLAG=1`;
- left, right, and internal oligo pick flags set to `1`;
- `PRIMER_NUM_RETURN` and `PRIMER_PRODUCT_SIZE_RANGE`;
- primer size, Tm, and GC limits through `PRIMER_*` tags;
- hydrolysis-probe limits through the corresponding `PRIMER_INTERNAL_*` tags.

The complete majority consensus is repeated in every record as required by
Primer3. No `SEQUENCE_TARGET` is used: the included region is the candidate
contract, and requiring primers to flank the complete conservation peak would
make ordinary qPCR product sizes impossible when a peak window is wide.

The parser requires exactly one output record per requested candidate and matches
records by `SEQUENCE_ID`, independent of output order. Primer3 global errors or
warnings and left/right/internal/pair explanation fields are preserved. A record
with no returned pair is valid. Malformed coordinates, missing members of a
partially returned assay, duplicate IDs, out-of-range oligos, inconsistent lengths,
or product sizes that disagree with the reported primer coordinates fail the
entire enabled run before final publication.

## Assay models and coordinates

```python
@dataclass(frozen=True, slots=True)
class DesignedOligo:
    sequence: str
    reference_start: int
    reference_end: int
    length: int
    tm: float
    gc_percent: float
    penalty: float | None
    metrics: tuple[tuple[str, str], ...]

@dataclass(frozen=True, slots=True)
class AssayCandidate:
    assay_id: str
    region_id: str
    primer3_index: int
    forward_primer: DesignedOligo
    probe: DesignedOligo
    reverse_primer: DesignedOligo
    product_size: int
    pair_penalty: float | None
    metrics: tuple[tuple[str, str], ...]

@dataclass(frozen=True, slots=True)
class PrimerDesignResult:
    status: Literal["SKIPPED", "COMPLETE"]
    reference_id: str | None
    candidates: tuple[CandidateRegion, ...]
    assays: tuple[AssayCandidate, ...]
    candidate_regions_path: Path | None
    assays_path: Path | None
    primer3_input_path: Path | None
    primer3_output_path: Path | None
    report_path: Path
```

Primer3 uses zero-based positions. Left-primer and internal-oligo locations are
converted to 1-based inclusive `[start + 1, start + length]`. Primer3 reports the
right primer's zero-based 3-prime position, so its reference interval is
`[position - length + 2, position + 1]`. Every stored oligo sequence remains in
the 5-prime-to-3-prime synthesis orientation returned by Primer3.

Assays retain Primer3's numeric index and receive stable IDs such as
`region-001-assay-001`. Known Tm, GC, size, product size, oligo penalties, and pair
penalty have typed fields. All additional scalar per-oligo and per-pair Primer3
fields are retained as sorted key/value tuples so later issues do not lose
version-specific metrics.

## Public service and artifacts

The new `qpcr_pipeline.primer_design` module exposes:

```python
design_primers(
    conservation: ConservationResult,
    config: PrimerDesignConfig,
    output_dir: Path,
    *,
    runner: Primer3Runner | None = None,
) -> PrimerDesignResult
```

Every run publishes:

```text
primer_design/
  primer_design_report.json
```

Enabled runs with at least one candidate additionally publish:

```text
primer_design/
  candidate_regions.tsv
  assays.tsv
  primer3_input.txt
  primer3_output.txt
```

An enabled run with no candidates publishes a header-only candidate TSV and assay
TSV, but does not start Primer3 and publishes no raw Primer3 files. An enabled run
with candidates but no returned assays publishes both raw files and a header-only
assay TSV.

Candidate TSV columns match the `CandidateRegion` fields. The assay TSV flattens
the three oligos into stable forward/probe/reverse sequence, coordinate, length,
Tm, GC, and penalty columns plus assay ID, region ID, Primer3 index, product size,
and pair penalty. Version-specific metrics and explanation text live in the JSON
report rather than producing an unstable TSV schema.

The schema-version-1 JSON report records status, enabled state, effective
configuration, reference ID, candidate and assay counts, complete candidate data,
complete typed assay data, all preserved Primer3 metrics/explanations, and relative
artifact paths. Raw input/output artifacts make the exact external exchange
auditable but are never copied into the JSON report or top-level QC summary.

All files use unique temporary siblings and atomic replacement. Before mutating a
fixed-name issue-#7 data artifact, publication removes the exact previous
`primer_design_report.json`. Data artifacts are published first and the JSON
report last. Disabled publication removes only the four exact issue-#7 data/raw
artifacts, preserves unrelated siblings and all conservation artifacts, then
publishes `SKIPPED` last.

## Pipeline integration

After `analyze_conservation`, `run_pipeline` calls `design_primers` with the
conservation result, primer-design configuration, output directory, and an
optional injected runner used by tests.

Existing sequence-set, clustering, alignment, conservation, and run-summary
semantics remain unchanged. `qc_report.json` gains exactly:

```json
"primer_design": {
  "status": "COMPLETE",
  "reference_id": "seq-3",
  "candidate_region_count": 4,
  "assay_count": 18
}
```

For disabled analysis, counts are zero and reference ID is null. Invalid direct
configuration and the enabled-primer-design/disabled-conservation mismatch fail
before the output directory is created.

## Error handling

- enabled primer design with a non-`COMPLETE` conservation result fails before
  issue-#7 artifact mutation;
- malformed or noncontiguous conservation coordinates and inconsistent majority
  consensus lengths fail explicitly;
- a missing or failed Primer3 executable, malformed Boulder-IO, duplicate or
  missing records, invalid assay members, and impossible coordinates fail
  explicitly;
- a valid empty candidate set or a valid Primer3 response with zero assays is
  complete, not an error;
- all publication failures propagate after invalidating the prior issue-#7 report;
- stderr excerpts in errors are bounded, while full successful exchanges remain
  available only in raw artifacts.

## Testing and verification

The standard offline suite follows test-driven development:

- configuration tests cover defaults, nested YAML, unknown fields, exact types,
  finite bounds, ordered ranges, region/product relationships, direct
  construction, and the conservation dependency;
- pure selection tests cover boundary shifting, short references, deduplication,
  aggregate metrics, usable length, every eligibility gate, lexicographic ties,
  overlap suppression, stable IDs, caps, empty data, and invalid coordinates;
- Boulder writer tests assert exact deterministic records, complete majority
  templates, included-region conversion, probe tags, record terminators, and
  escaping/rejection of unsafe identifiers;
- parser tests cover multiple records and assays, out-of-order records, zero-pair
  explanations, all coordinate conversions, additional metric preservation, and
  each malformed-output failure;
- runner tests use a temporary fake executable to exercise stdin, arguments,
  stdout, stderr, missing executable, and nonzero exit without a shell;
- artifact tests cover stable schemas, raw exchanges, header-only empty results,
  disabled stale cleanup, hardlink-safe replacement, report-last ordering, and
  prior-report invalidation at each replacement failure point;
- pipeline tests prove majority consensus reaches Primer3, multiple assays are
  summarized, disabled runs do not invoke the runner, and all earlier contracts
  remain unchanged;
- a conditional integration test sends a fixed synthetic majority-consensus
  fixture through a real `primer3_core` when it is installed. It validates
  structural scientific invariants rather than version-sensitive penalty values
  and skips with an explicit reason when the binary is unavailable.

The current Windows development environment does not expose `primer3_core`, so
the fake-runner suite is mandatory everywhere and the real-binary integration test
is expected to skip locally until the packaged WSL2/backend toolchain exists.

## Acceptance mapping

- Candidate ranking using conservation, coverage, entropy, gaps, and usable
  length: deterministic region aggregation, eligibility, and lexicographic order.
- Configurable maximum regions: `max_candidate_regions` and greedy overlap cap.
- Majority consensus supplied to Primer3: full `SEQUENCE_TEMPLATE` in every
  Boulder record.
- Configurable hydrolysis-probe qPCR preset: nested primer/probe constraints and
  product-size range.
- Multiple assays per region: `assays_per_region` maps to `PRIMER_NUM_RETURN`.
- Persisted Tm, GC, oligo sizes, amplicon size, and Primer3 metrics: typed models,
  stable TSV fields, complete JSON metrics, and raw exchanges.
- Primer3 integration fixture: conditional real-binary test plus deterministic
  offline runner/parser coverage.
