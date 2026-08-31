# Issue #5: MAFFT alignment and genomic coordinate design

## Goal

Align the Discovery Set with MAFFT, deterministically select or accept a reference
sequence, normalize reverse-complemented inputs, and publish an auditable mapping
between 1-based alignment and reference coordinates without changing Discovery or
Evaluation Set membership.

## Scope and non-goals

This issue adds alignment configuration, a safe MAFFT process boundary, strict
aligned-FASTA validation, deterministic reference selection, orientation records,
coordinate artifacts, and pipeline integration.

It does not add an external reference FASTA, reference accession acquisition,
primer/probe discovery, consensus generation, variant statistics, visualization,
checkpoint/resume behavior, or dependency diagnostics. An explicit reference is an
existing Discovery Set sequence ID. Later issues may add external canonical
references without changing the coordinate model defined here.

Alignment is disabled by default so existing local and NCBI workflows remain usable
without MAFFT. The pipeline still publishes an alignment report with status
`SKIPPED` and removes exact stale alignment artifacts from a previous enabled run.

## Configuration

`PipelineConfig` gains a frozen, slotted configuration:

```python
@dataclass(frozen=True, slots=True)
class AlignmentConfig:
    enabled: bool = False
    threads: int = 1
    reference_id: str | None = None
```

YAML uses an optional top-level section:

```yaml
alignment:
  enabled: true
  threads: 4
  reference_id: NC_000000.1
```

Validation is shared by YAML-loaded and directly constructed configurations:

- `enabled` is a boolean;
- `threads` is an integer from 1 through 256, excluding booleans;
- `reference_id` is `null` or a non-blank string;
- only `enabled`, `threads`, and `reference_id` are accepted fields.

The initial strategy is intentionally fixed rather than configurable. Enabled
multi-sequence runs use MAFFT `--auto`, nucleotide mode, accurate direction
adjustment, input ordering, the configured thread count, and deterministic iterative
refinement threading.

## Alignment boundary

A new `qpcr_pipeline.alignment` module owns reference selection, MAFFT execution,
output parsing, coordinate construction, and artifact publication. It defines an
injectable runner protocol so standard tests remain offline:

```python
class MafftRunner(Protocol):
    def run(
        self,
        input_path: Path,
        output_path: Path,
        config: AlignmentConfig,
    ) -> None: ...

align_discovery(
    records: tuple[LocalSequenceRecord, ...],
    discovery_set: DiscoverySet,
    config: AlignmentConfig,
    output_dir: Path,
    *,
    runner: MafftRunner | None = None,
) -> AlignmentResult
```

The production runner resolves `mafft` from `PATH` and invokes a structured argument
list equivalent to:

```text
mafft --auto --nuc --inputorder --adjustdirectionaccurately
  --thread <threads> --threadit 0 --quiet input.fasta
```

The input path is one argument and no shell command string is evaluated. Standard
output is the aligned FASTA. Standard error is captured only for bounded diagnostics.
The runner requires exit code zero and a non-empty aligned output file.

`--threadit 0` disables multithreading only for the iterative refinement stage,
which avoids the nondeterminism MAFFT documents for multithreaded iterative
refinement. The all-to-all and progressive stages still use the configured threads.

## Immutable result models

The public module models are frozen and slotted:

```python
@dataclass(frozen=True, slots=True)
class AlignedSequence:
    sequence_id: str
    aligned_sequence: str
    orientation: Literal["forward", "reverse_complemented"]

@dataclass(frozen=True, slots=True)
class AlignmentCoordinate:
    alignment_position: int
    reference_position: int | None
    reference_base: str | None

@dataclass(frozen=True, slots=True)
class AlignmentResult:
    status: Literal["SKIPPED", "COMPLETE"]
    discovery_set: DiscoverySet
    reference_id: str | None
    reference_mode: Literal["explicit", "automatic"] | None
    sequences: tuple[AlignedSequence, ...]
    coordinates: tuple[AlignmentCoordinate, ...]
    alignment_fasta_path: Path | None
    coordinate_map_path: Path | None
    report_path: Path
```

The Discovery Set object in the result is the unchanged input contract.

## Reference selection

An explicit `reference_id` must occur exactly once in the Discovery Set. It is not
silently taken from the Evaluation Set and clustering is never changed to force it
into Discovery.

Without an explicit reference, the service selects one deterministically using this
ascending key:

1. ambiguous nucleotide fraction, where every base outside `A`, `C`, `G`, and `T`
   is ambiguous;
2. negative ungapped sequence length, so longer sequences win ties;
3. original Discovery Set position.

Empty sequences are rejected before selection. The chosen reference is placed first
in MAFFT's temporary input so direction normalization is anchored to it. All other
sequences retain Discovery Set order. Stable internal IDs are derived from original
Discovery positions, not temporary input positions.

## Direction normalization and identity validation

Temporary input uses stable IDs such as `geison-00000000`. MAFFT marks a sequence it
reverse-complements by prepending `_R_` to its title. The strict parser accepts each
expected internal ID exactly once, with either no prefix or one `_R_` prefix.

The reference must remain forward. Every aligned sequence must:

- have the same positive alignment length;
- contain only uppercase IUPAC DNA symbols and `-` gaps;
- contain at least one non-gap base;
- equal its original sequence after gap removal when forward;
- equal the IUPAC reverse complement of its original sequence after gap removal when
  marked reverse.

Unknown, missing, duplicate, multiply prefixed, mutated, or inconsistently sized
records fail before final artifacts are published. Final aligned FASTA records use
original sequence IDs in original Discovery Set order. `_R_` and internal IDs do not
appear in final artifacts.

## Coordinate model

Coordinates are 1-based. Iterate across the aligned reference from left to right:

- `alignment_position` is the alignment column, starting at 1;
- a non-gap reference base increments `reference_position` and records that value;
- a reference gap records `reference_position` and `reference_base` as empty/null;
- `reference_base` otherwise records the aligned reference base.

The final non-gap reference count must equal the original reference sequence length.
This column mapping is sufficient for downstream stages to convert an alignment
interval into reference coordinates without guessing around insertions.

## Empty, singleton, and disabled behavior

- Disabled alignment never constructs or calls a runner. It removes only the exact
  stale aligned FASTA and coordinate TSV, then atomically publishes a `SKIPPED`
  report last.
- Enabled empty Discovery publishes an empty aligned FASTA, a header-only coordinate
  TSV, and a `COMPLETE` report with no reference, without invoking MAFFT.
- Enabled singleton Discovery selects or validates that record as reference, publishes
  an identity alignment and coordinate map, and does not invoke MAFFT.
- Enabled Discovery with two or more records lazily constructs the production runner
  only when no runner was injected.

## Artifacts

Every run publishes:

```text
alignment/
  alignment_report.json
```

Enabled runs additionally publish:

```text
alignment/
  discovery_alignment.fasta
  coordinate_map.tsv
```

`coordinate_map.tsv` has the exact header:

```text
alignment_position\treference_position\treference_base
```

Reference gaps use empty fields for the latter two columns.

The schema-version-1 JSON report contains:

- status and whether alignment was enabled;
- tool name and effective parameters;
- exact ordered Discovery Set IDs;
- reference mode, ID, and automatic-selection rule when applicable;
- ordered sequence orientation records;
- Discovery count, alignment length, and reverse-complement count;
- relative artifact filenames or `null` when skipped.

All final files use unique temporary siblings plus atomic replacement. FASTA and TSV
are published only after complete validation. The report is published last. A failed
run does not publish a new `COMPLETE` report.

## Pipeline integration

After clustering, `run_pipeline` selects records for the unchanged Discovery Set by
iterating approved records in their original order. It does not build a dictionary
that could hide duplicate IDs. It then calls `align_discovery` and accepts an
injectable `mafft_runner`:

```python
run_pipeline(
    config,
    outdir,
    *,
    ncbi_client=None,
    cdhit_runner=None,
    mafft_runner=None,
)
```

`RunSummary`, Evaluation Set fields, Discovery Set fields, and cluster artifacts keep
their existing semantics. `qc_report.json` gains an `alignment` object containing
status, reference ID, and reference mode so the run's major stages remain traceable
from the top-level report.

Invalid direct alignment configuration and a missing explicit Discovery reference
fail before MAFFT execution. Configuration validation still occurs before output
directory creation.

## Error handling

- Missing `mafft` raises `MafftError` naming the configured executable and explaining
  that MAFFT can be installed or alignment disabled.
- A nonzero process exit raises sanitized `MafftError` with exit code and at most
  2,000 normalized stderr characters.
- Missing or empty stdout output, malformed FASTA, membership mismatches, reference
  reversal, sequence mutation, unequal lengths, or invalid symbols fail explicitly.
- Invalid configuration, record membership, empty sequences, and explicit-reference
  errors occur before subprocess execution and final artifact mutation.

## Testing

The standard suite remains deterministic and offline:

- configuration tests cover defaults, YAML, unknown fields, types, bounds, blank
  references, and direct construction;
- selection tests cover explicit reference, automatic ambiguity/length/order ties,
  empty records, and references absent from Discovery;
- parser/service tests use realistic `_R_` output and cover original-ID restoration,
  Discovery ordering, reverse complements, strict composition, coordinate gaps,
  empty/singleton/disabled behavior, artifacts, atomic replacement, and failures;
- subprocess tests assert the exact structured command, paths with spaces, lazy
  construction, missing executable, bounded stderr, and missing/empty output;
- pipeline tests prove MAFFT receives only Discovery representatives, Evaluation and
  Discovery sets remain unchanged, disabled local/NCBI runs do not require MAFFT, and
  top-level traceability is present.

A real test lives at `integration_tests/test_mafft.py`, outside default test discovery.
It skips unless `mafft` is present, aligns a deterministic small nucleotide set with
one reverse-complemented member, and verifies orientations, coordinates, and artifacts.

MAFFT is not installed on the current Windows host, so this real-tool test is expected
to skip here.

## Acceptance mapping

- MAFFT with automatic strategy and configurable resources: safe runner using
  `--auto` and configurable threads.
- Explicit reference: validated Discovery Set `reference_id`.
- Automatic reference: deterministic ambiguity, length, and order selection with the
  decision persisted.
- Reverse normalization: `--adjustdirectionaccurately`, strict `_R_` handling, and
  per-sequence orientation records.
- Alignment and genomic coordinates: aligned FASTA plus 1-based coordinate TSV and
  immutable coordinate models.
- Real integration: separately discovered conditional MAFFT test.
