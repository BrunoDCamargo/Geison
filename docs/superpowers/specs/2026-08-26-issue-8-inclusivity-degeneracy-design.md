# Issue #8: Inclusivity and bounded IUPAC degeneracy

## Goal

Evaluate every Primer3 assay from issue #7 against every QC-approved target in
the `EvaluationSet`, make each oligo mismatch auditable in synthesis orientation,
and propose bounded IUPAC expansions only when they improve exact binding-site
coverage. The original assay is always preserved and evaluated beside the
proposed version.

This specification implements the acceptance criteria of GitHub issue #8. The
live issue is authoritative: inclusivity uses the complete `EvaluationSet`, not
the clustered `DiscoverySet`. The README wording introduced with issue #7 is
corrected accordingly. Off-target specificity remains issue #9.

## Scope and non-goals

This issue adds immutable inclusivity configuration and results, a deterministic
local oligo search, complete-assay geometry checks, variation and mismatch
classification, bounded degeneracy proposals, normalized TSV/JSON artifacts,
pipeline integration, documentation, and synthetic tests.

It does not scan off-target databases, score thermodynamic consequences of a
degenerate proposal, replace Primer3 assays, rank final assay risk, or add UI.
Specificity belongs to issue #9; consolidated reports and UI belong to later
issues. A proposed sequence is therefore a transparent computational candidate,
not a laboratory recommendation.

## Architecture and inputs

The new `qpcr_pipeline.inclusivity` module exposes one service:

```python
evaluate_inclusivity(
    records: tuple[SequenceRecord, ...],
    evaluation_set: EvaluationSet,
    primer_design: PrimerDesignResult,
    config: InclusivityConfig,
    output_dir: Path,
) -> InclusivityResult
```

`records` contains every and only QC-approved `EvaluationSet` member, once and in
the same order. `primer_design` must be `COMPLETE` when inclusivity is enabled.
These validations happen before issue-#8 artifacts are mutated. The service has
no external executable or new dependency.

When inclusivity is disabled, the service validates only `InclusivityConfig`; it
does not inspect records, the Evaluation Set, or Primer3 assays. This keeps the
stage independently skippable while still publishing its `SKIPPED` report.

The implementation is divided into focused internal units:

- IUPAC support-set and reverse-complement operations;
- bounded local hit enumeration and deterministic hit ordering;
- compatible assay-triplet selection and coordinate conversion;
- variation aggregation and bounded proposal search;
- immutable result construction and atomic publication.

## Configuration

`PipelineConfig` gains a frozen, slotted model:

```python
@dataclass(frozen=True, slots=True)
class InclusivityConfig:
    enabled: bool = False
    search_flank: int = 250
    max_hits_per_oligo: int = 20
    max_primer_mismatches: int = 2
    max_probe_mismatches: int = 1
    reject_primer_3_prime_mismatch: bool = True
    primer_3_prime_bases: int = 5
    max_primer_degeneracy: int = 16
    max_probe_degeneracy: int = 4
    allow_primer_3_prime_degeneracy: bool = False
    max_amplicon_size_delta: int = 20
```

YAML uses an optional `inclusivity` mapping. Unknown fields are rejected.
Booleans are accepted only for boolean fields. Counts and limits exclude
booleans; `max_hits_per_oligo`, `primer_3_prime_bases`, and both degeneracy limits
must be positive, while flank, mismatch limits, and amplicon delta may be zero.
Direct construction and YAML loading share the same validation. Enabled
inclusivity requires enabled primer design, which already requires conservation
and alignment. Inclusivity defaults to disabled for backward compatibility.

The degeneracy limits are total sequence degeneracy, calculated as the product
of the cardinalities of all IUPAC symbols in the oligo. They are not counts of
changed positions.

## IUPAC comparison semantics

Accepted sequence symbols are uppercase or lowercase
`ACGTRYSWKMBDHVN`; values are normalized to uppercase. Any other symbol fails
with the record ID and position but without embedding the full sequence in the
error. Reverse complements support the same alphabet.

Each symbol maps to its set of possible canonical bases. At an oligo position,
the target is an exact match when the target symbol's support set is a subset of
the oligo symbol's support set. This conservative rule means an ambiguous target
base is covered only when the oligo covers every represented possibility. A
mismatch is any position that fails that test.

Mismatch positions are always 1-based in the oligo's 5-prime-to-3-prime
synthesis orientation. For primers, a mismatch is a 3-prime mismatch when its
position is greater than `oligo_length - primer_3_prime_bases`; a mismatch at any
other position is internal or 5-prime as determined by its exact position. The
probe-mismatch flag is true for any probe hit with one or more mismatch positions.

## Directed local hit search

Each raw Evaluation sequence is considered in two possible orientations:
as supplied and as its reverse complement. In either orientation, the reference
coordinates stored on the issue-#7 assay remain the expected coordinates. This
handles records supplied on either strand without a whole-Evaluation alignment.

For each assay, orientation, and oligo:

1. Create a search window by expanding the oligo's expected 1-based reference
   interval by `search_flank` on both sides and clipping it to the record.
2. Slide a same-length target segment over every possible start in that window.
3. Compare the forward primer and probe directly to the oriented target segment.
   Compare the reverse primer to the reverse complement of its target segment so
   mismatch positions remain in reverse-primer synthesis orientation.
4. Retain at most `max_hits_per_oligo` hits in this order: mismatch count,
   3-prime mismatch count, absolute displacement from the expected start,
   oriented start, then target segment lexicographically.

Search retention is deliberately independent of compatibility thresholds. This
preserves the best near-match for variation analysis even when the original
oligo is incompatible. A record too short for an oligo, or an empty search
window, produces no hit rather than an exception.

For each orientation, all retained forward/probe/reverse combinations are tested
for the strict reference-oriented geometry:

```text
forward_end < probe_start <= probe_end < reverse_start
```

The candidate amplicon spans the forward-primer start through the reverse-primer
end. Its size must be within the issue-#7 assay product size plus or minus
`max_amplicon_size_delta`. Plausible triplets are ordered by total mismatches,
total primer 3-prime mismatches, absolute amplicon-size delta, total coordinate
displacement, orientation (`FORWARD` before `REVERSE_COMPLEMENT`), and coordinates.
The first is the selected triplet.

Oriented coordinates are converted back to 1-based inclusive source-record
coordinates for publication. For a reverse-complement hit of oriented interval
`[start, end]` in a record of length `L`, source coordinates are
`[L - end + 1, L - start + 1]`. Both source coordinates and the chosen record
orientation are reported.

If no plausible triplet exists, the result retains the independently best hit
for each oligo when one exists, marks geometry false, and never marks the complete
assay compatible.

## Compatibility results

An original primer hit is compatible when its mismatch count is at most
`max_primer_mismatches` and, when `reject_primer_3_prime_mismatch` is true, it has
no mismatch in the configured 3-prime suffix. A probe hit is compatible when its
mismatch count is at most `max_probe_mismatches`. A complete assay is compatible
only when a plausible triplet exists and all three selected hits are compatible.

Every assay/Evaluation-sequence result reports:

- chosen orientation, source amplicon coordinates and size when geometry exists;
- for each oligo, whether a hit exists, exact match, mismatch count, ordered
  mismatch positions, 3-prime mismatch, probe mismatch, compatibility, source
  coordinates, and displacement;
- geometry status and original complete-assay compatibility;
- the same compatibility fields recomputed with accepted proposals, using the
  original sequence for any oligo without an accepted proposal.

An incompatible or missing result is data, not an error. Empty assay and empty
Evaluation sets are valid complete results with header-only TSVs. Fractions with
a zero denominator are serialized as `null`, never as zero or NaN.

## Variation aggregation and degeneracy proposals

Variation aggregation uses the selected plausible binding site for each
assay/Evaluation-sequence pair. It does not use an unrelated independent hit when
the assay has no plausible geometry. Reverse-primer sites are converted to oligo
synthesis orientation before aggregation.

For every assay and oligo role, a variation row is emitted at each position where
at least one selected target symbol's support differs from the original symbol's
support. The result records the original IUPAC symbol and support, the union of
observed target support, the minimal IUPAC symbol representing that union,
affected Evaluation IDs and count, affected fraction, and whether the position
lies in a primer 3-prime suffix. Raw Evaluation sequences are not copied into
reports.

A position is eligible for expansion only when its observed union adds support
to the original symbol. A primer position in its final configured 3-prime bases
is ineligible by default. Each eligible position has exactly two choices: retain
the original symbol or replace it with the minimal IUPAC symbol for the union of
the original and observed support.

The proposal search enumerates combinations deterministically and prunes any
combination whose total sequence degeneracy exceeds the role-specific limit. For
each surviving sequence, exact binding-site coverage is recalculated over all
selected sites for that assay and role. Ordering is:

1. exact covered site count descending;
2. total sequence degeneracy ascending;
3. number of changed positions ascending;
4. proposed sequence lexicographically.

A proposal is accepted only when it strictly increases exact covered-site count
over the original. Otherwise the original remains effective. The record uses
status `ACCEPTED`, `UNCHANGED`, or `REJECTED` with one of these stable reasons:

- `NO_VARIATION`;
- `NO_GEOMETRIC_SITES`;
- `NO_IMPROVEMENT`;
- `REJECTED_3_PRIME` when useful changes exist only in a protected suffix;
- `REJECTED_LIMIT` when every improving candidate exceeds the configured limit;
- `ACCEPTED_IMPROVEMENT`.

Each proposal records role, original and proposed sequences, status and reason,
original and proposed degeneracy, changed 1-based positions, and original and
proposed exact counts and fractions. It never mutates `AssayCandidate` or silently
substitutes its sequence. Complete-assay inclusivity is then recomputed with the
accepted per-oligo proposals on the same selected geometric triplet and reported
beside the original result. This is sufficient because an IUPAC expansion cannot
invalidate the selected coordinates; issue #8 does not use a proposal to discover
a different binding locus.

## Public result models

The module exposes frozen, slotted immutable models with tuples rather than
mutable collections:

```python
@dataclass(frozen=True, slots=True)
class OligoMatch:
    assay_id: str
    sequence_id: str
    role: Literal["FORWARD", "PROBE", "REVERSE"]
    orientation: Literal["FORWARD", "REVERSE_COMPLEMENT"]
    hit_rank: int
    source_start: int
    source_end: int
    expected_start: int
    expected_end: int
    displacement: int
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    exact_match: bool
    three_prime_mismatch: bool
    probe_mismatch: bool
    compatible: bool
    selected: bool

@dataclass(frozen=True, slots=True)
class ProposedOligoCompatibility:
    role: Literal["FORWARD", "PROBE", "REVERSE"]
    effective_sequence: str
    exact_match: bool
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    three_prime_mismatch: bool
    probe_mismatch: bool
    compatible: bool

@dataclass(frozen=True, slots=True)
class AssayInclusivity:
    assay_id: str
    sequence_id: str
    orientation: Literal["FORWARD", "REVERSE_COMPLEMENT"] | None
    geometry_found: bool
    source_amplicon_start: int | None
    source_amplicon_end: int | None
    amplicon_size: int | None
    forward_match: OligoMatch | None
    probe_match: OligoMatch | None
    reverse_match: OligoMatch | None
    original_compatible: bool
    proposed_forward: ProposedOligoCompatibility | None
    proposed_probe: ProposedOligoCompatibility | None
    proposed_reverse: ProposedOligoCompatibility | None
    proposed_compatible: bool

@dataclass(frozen=True, slots=True)
class OligoVariation:
    assay_id: str
    role: Literal["FORWARD", "PROBE", "REVERSE"]
    oligo_position: int
    original_symbol: str
    original_support: tuple[str, ...]
    observed_symbol: str
    observed_support: tuple[str, ...]
    affected_sequence_ids: tuple[str, ...]
    affected_sequence_count: int
    affected_fraction: float | None
    primer_3_prime_position: bool

@dataclass(frozen=True, slots=True)
class DegeneracyProposal:
    assay_id: str
    role: Literal["FORWARD", "PROBE", "REVERSE"]
    original_sequence: str
    proposed_sequence: str
    status: Literal["ACCEPTED", "UNCHANGED", "REJECTED"]
    reason: str
    original_degeneracy: int
    proposed_degeneracy: int
    changed_positions: tuple[int, ...]
    binding_site_count: int
    original_exact_count: int
    original_exact_fraction: float | None
    proposed_exact_count: int
    proposed_exact_fraction: float | None

@dataclass(frozen=True, slots=True)
class InclusivityResult:
    status: Literal["SKIPPED", "COMPLETE"]
    evaluation_sequence_ids: tuple[str, ...]
    oligo_matches: tuple[OligoMatch, ...]
    assay_results: tuple[AssayInclusivity, ...]
    variations: tuple[OligoVariation, ...]
    proposals: tuple[DegeneracyProposal, ...]
    oligo_matches_path: Path | None
    assay_inclusivity_path: Path | None
    oligo_variations_path: Path | None
    degeneracy_proposals_path: Path | None
    report_path: Path
```

No public model may store an untyped dictionary as a substitute for required
acceptance data.

## Artifacts and publication

Every run publishes:

```text
inclusivity/
  inclusivity_report.json
```

Enabled runs additionally publish these normalized UTF-8 TSVs, including headers
when they have zero rows:

```text
inclusivity/
  oligo_matches.tsv
  assay_inclusivity.tsv
  oligo_variations.tsv
  degeneracy_proposals.tsv
```

`oligo_matches.tsv` contains every retained local hit and marks the selected
triplet or independent fallback hit. `assay_inclusivity.tsv` has one row per
assay/Evaluation sequence with original and proposed complete compatibility.
`oligo_variations.tsv` has one row per observed oligo position variation.
`degeneracy_proposals.tsv` has one row per assay and oligo role, including
rejected and no-op outcomes.

The schema-version-1 JSON report contains status, effective configuration,
Evaluation IDs, counts and nullable fractions, complete typed result data, and
relative artifact paths. It includes oligo proposal sequences but no raw target
sequences.

All writes use unique temporary siblings and atomic replacement. Before any
issue-#8 data artifact is mutated, the prior `inclusivity_report.json` is removed.
Data TSVs are published first and the JSON report last. Disabled publication
removes only the four exact issue-#8 TSVs, preserves unrelated files and all
upstream artifacts, then publishes `SKIPPED` last. Any failure after invalidation
propagates without leaving an apparently current final report.

## Pipeline integration

After `design_primers`, `run_pipeline` calls `evaluate_inclusivity` with the
approved records, exact `EvaluationSet`, Primer3 result, inclusivity config, and
output directory. No test injection is needed because the evaluator is pure
Python.

`qc_report.json` gains:

```json
"inclusivity": {
  "status": "COMPLETE",
  "evaluation_sequence_count": 42,
  "assay_count": 15,
  "assay_evaluation_count": 630,
  "original_compatible_count": 570,
  "proposed_compatible_count": 602
}
```

Disabled runs use zero counts. Existing sequence-set, clustering, alignment,
conservation, Primer3, and run-summary semantics remain unchanged. Invalid direct
configuration and the enabled-inclusivity/disabled-primer-design mismatch fail
before the output directory is created.

README configuration and artifact documentation is extended, and its issue
boundary is corrected to say that issue #8 evaluates inclusivity against the
complete Evaluation Set while issue #9 evaluates specificity against off-target
datasets.

## Error handling

- enabled inclusivity with a non-`COMPLETE` Primer3 result fails before issue-#8
  artifact mutation;
- duplicate, blank, missing, extra, or out-of-order Evaluation IDs fail clearly;
- invalid IUPAC symbols, inconsistent oligo lengths or coordinates, and invalid
  direct configuration fail explicitly;
- no oligo hit, no geometric triplet, no assay, no Evaluation member, protected
  3-prime variation, and an unhelpful or over-limit proposal are valid outcomes;
- publication failures propagate after invalidating the prior report;
- errors identify the relevant assay, record, role, or position without exposing
  complete biological sequences.

## Testing and acceptance mapping

Development follows test-driven implementation. Unit tests cover:

- configuration defaults, YAML overrides, unknown fields, direct validation,
  dependency validation, and validation before output creation;
- IUPAC support, conservative ambiguity comparison, reverse complement, and total
  degeneracy;
- exact, internal, 5-prime, final-five 3-prime, and probe mismatches with 1-based
  synthesis-orientation positions;
- both record orientations, source-coordinate conversion, clipped local windows,
  deterministic repeated hits, hit caps, plausible geometry, amplicon tolerance,
  and no-hit/no-geometry cases;
- complete original compatibility thresholds and default 3-prime rejection;
- positional variation aggregation over the entire Evaluation Set, including a
  member excluded from the Discovery Set;
- accepted IUPAC improvement, original preservation, primer/probe limits,
  protected 3-prime rejection, limit rejection, no improvement, and deterministic
  tie-breaking;
- original and proposed complete-assay results side by side;
- empty sets, header-only artifacts, skipped stale cleanup, unrelated-file
  preservation, JSON-last atomic publication, replacement failures, normalized
  newlines, and report schema;
- end-to-end pipeline/QC integration and README responsibility wording.

The existing full unit, integration, and pytest suites must remain green. No new
external-binary integration test is required for this pure-Python stage.
