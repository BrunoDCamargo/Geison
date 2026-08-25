"""Safe Primer3 Boulder-IO request and response boundary."""

from __future__ import annotations

import math
import shutil
import subprocess
from typing import Protocol

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.primer_design import (
    AssayCandidate,
    CandidateRegion,
    DesignedOligo,
    PrimerDesignError,
)


class Primer3Runner(Protocol):
    def run(self, input_text: str) -> str:
        """Run Primer3 against complete Boulder-IO input."""


class SubprocessPrimer3Runner:
    def __init__(self, executable: str = "primer3_core") -> None:
        self._executable = executable

    def run(self, input_text: str) -> str:
        executable = shutil.which(self._executable)
        if executable is None:
            raise PrimerDesignError(
                f"Primer3 executable '{self._executable}' was not found on PATH."
            )
        try:
            input_bytes = input_text.encode("utf-8", errors="strict")
            completed = subprocess.run(
                [executable, "--strict_tags", "--io_version=4"],
                input=input_bytes,
                capture_output=True,
                check=False,
                shell=False,
            )
        except Exception as error:
            raise PrimerDesignError(
                f"Primer3 execution failed ({type(error).__name__})."
            ) from error
        try:
            stdout = completed.stdout.decode("utf-8", errors="strict")
            stderr = completed.stderr.decode("utf-8", errors="strict")
        except UnicodeDecodeError as error:
            raise PrimerDesignError("Primer3 output was not valid UTF-8.")
        if completed.returncode != 0:
            stderr_excerpt = _redact_consensus(stderr, input_text)[:2000]
            raise PrimerDesignError(
                f"Primer3 exited with status {completed.returncode}: "
                f"{stderr_excerpt}"
            )
        if not stdout:
            raise PrimerDesignError("Primer3 produced empty stdout.")
        return stdout


def _redact_consensus(text: str, input_text: str) -> str:
    for line in input_text.splitlines():
        if line.startswith("SEQUENCE_TEMPLATE="):
            consensus = line.removeprefix("SEQUENCE_TEMPLATE=")
            if consensus:
                text = text.replace(consensus, "<consensus omitted>")
    return text


def build_primer3_input(
    consensus: str,
    candidates: tuple[CandidateRegion, ...],
    config: PrimerDesignConfig,
) -> str:
    """Serialize deterministic Primer3 v4 Boulder-IO input records."""
    if (
        not isinstance(consensus, str)
        or not consensus
        or any(base not in "ACGT" for base in consensus)
    ):
        raise PrimerDesignError(
            "Primer3 consensus must contain only uppercase canonical DNA bases."
        )
    records: list[str] = []
    for candidate in candidates:
        if any(delimiter in candidate.region_id for delimiter in ("=", "\r", "\n")):
            raise PrimerDesignError(
                "Primer3 candidate region ID contains unsafe Boulder-IO delimiters."
            )
        included_start = candidate.reference_start - 1
        included_length = candidate.reference_end - candidate.reference_start + 1
        lines = [
            f"SEQUENCE_ID={candidate.region_id}",
            f"SEQUENCE_TEMPLATE={consensus}",
            f"SEQUENCE_INCLUDED_REGION={included_start},{included_length}",
            "PRIMER_TASK=generic",
            "PRIMER_FIRST_BASE_INDEX=0",
            "PRIMER_PICK_LEFT_PRIMER=1",
            "PRIMER_PICK_INTERNAL_OLIGO=1",
            "PRIMER_PICK_RIGHT_PRIMER=1",
            f"PRIMER_NUM_RETURN={config.assays_per_region}",
            f"PRIMER_PRODUCT_SIZE_RANGE={config.product_size_min}-{config.product_size_max}",
            "PRIMER_EXPLAIN_FLAG=1",
            f"PRIMER_MIN_SIZE={config.primer.min_size}",
            f"PRIMER_OPT_SIZE={config.primer.opt_size}",
            f"PRIMER_MAX_SIZE={config.primer.max_size}",
            f"PRIMER_MIN_TM={config.primer.min_tm}",
            f"PRIMER_OPT_TM={config.primer.opt_tm}",
            f"PRIMER_MAX_TM={config.primer.max_tm}",
            f"PRIMER_MIN_GC={config.primer.min_gc_percent}",
            f"PRIMER_MAX_GC={config.primer.max_gc_percent}",
            f"PRIMER_INTERNAL_MIN_SIZE={config.probe.min_size}",
            f"PRIMER_INTERNAL_OPT_SIZE={config.probe.opt_size}",
            f"PRIMER_INTERNAL_MAX_SIZE={config.probe.max_size}",
            f"PRIMER_INTERNAL_MIN_TM={config.probe.min_tm}",
            f"PRIMER_INTERNAL_OPT_TM={config.probe.opt_tm}",
            f"PRIMER_INTERNAL_MAX_TM={config.probe.max_tm}",
            f"PRIMER_INTERNAL_MIN_GC={config.probe.min_gc_percent}",
            f"PRIMER_INTERNAL_MAX_GC={config.probe.max_gc_percent}",
            "=",
        ]
        records.append("\n".join(lines))
    return "\n".join(records) + ("\n" if records else "")


def parse_primer3_output(
    text: str,
    candidates: tuple[CandidateRegion, ...],
    consensus: str,
) -> tuple[tuple[AssayCandidate, ...], dict[str, dict[str, str]]]:
    """Parse Primer3 v4 Boulder-IO records into typed assay candidates."""
    records = _split_records(text)
    for record in records:
        if "PRIMER_ERROR" in record:
            error = record["PRIMER_ERROR"].replace(consensus, "<consensus omitted>")
            raise PrimerDesignError(f"Primer3 reported an error: {error[:2000]}")
    records_by_id: dict[str, dict[str, str]] = {}
    for record in records:
        record_id = _required_tag(record, "SEQUENCE_ID")
        if record_id in records_by_id:
            raise PrimerDesignError(
                f"Primer3 output contains duplicate SEQUENCE_ID '{record_id}'."
            )
        records_by_id[record_id] = record
    requested_ids = {candidate.region_id for candidate in candidates}
    for record_id in records_by_id:
        if record_id not in requested_ids:
            raise PrimerDesignError(
                f"Primer3 output contains unknown SEQUENCE_ID '{record_id}'."
            )
    for candidate in candidates:
        if candidate.region_id not in records_by_id:
            raise PrimerDesignError(
                "Primer3 output is missing requested "
                f"SEQUENCE_ID '{candidate.region_id}'."
            )
    assays: list[AssayCandidate] = []
    details: dict[str, dict[str, str]] = {}
    for candidate in candidates:
        record = records_by_id[candidate.region_id]
        details[candidate.region_id] = {
            key: value
            for key, value in record.items()
            if key.endswith(("_ERROR", "_WARNING", "_EXPLAIN"))
        }
        returned_counts = {
            kind: _parse_int(record, f"PRIMER_{kind}_NUM_RETURNED")
            for kind in ("LEFT", "INTERNAL", "RIGHT", "PAIR")
        }
        for kind, count in returned_counts.items():
            if count < 0:
                raise PrimerDesignError(
                    f"Primer3 output PRIMER_{kind}_NUM_RETURNED must be nonnegative."
                )
        pair_count = returned_counts["PAIR"]
        for kind in ("LEFT", "INTERNAL", "RIGHT"):
            if pair_count > returned_counts[kind]:
                raise PrimerDesignError(
                    "Primer3 output PRIMER_PAIR_NUM_RETURNED count "
                    f"{pair_count} exceeds PRIMER_{kind}_NUM_RETURNED count "
                    f"{returned_counts[kind]}."
                )
        for index in range(pair_count):
            forward = _parse_oligo(record, "LEFT", index)
            probe = _parse_oligo(record, "INTERNAL", index)
            reverse = _parse_oligo(record, "RIGHT", index)
            _validate_oligo_bounds(forward, candidate, "LEFT", index)
            _validate_oligo_bounds(probe, candidate, "INTERNAL", index)
            _validate_oligo_bounds(reverse, candidate, "RIGHT", index)
            pair_prefix = f"PRIMER_PAIR_{index}_"
            pair_penalty_text = record.get(f"{pair_prefix}PENALTY")
            product_size = _parse_int(record, f"{pair_prefix}PRODUCT_SIZE")
            coordinate_product_size = (
                reverse.reference_end - forward.reference_start + 1
            )
            if product_size <= 0 or coordinate_product_size <= 0:
                raise PrimerDesignError(
                    f"Primer3 output pair {index} product size must be positive."
                )
            if product_size != coordinate_product_size:
                raise PrimerDesignError(
                    f"Primer3 output pair {index} reports product size "
                    f"{product_size}, but its primers imply coordinate-derived "
                    f"size {coordinate_product_size}."
                )
            assays.append(
                AssayCandidate(
                    assay_id=f"{candidate.region_id}-assay-{index + 1:03d}",
                    region_id=candidate.region_id,
                    primer3_index=index,
                    forward_primer=forward,
                    probe=probe,
                    reverse_primer=reverse,
                    product_size=product_size,
                    pair_penalty=(
                        _parse_float(record, f"{pair_prefix}PENALTY")
                        if pair_penalty_text is not None
                        else None
                    ),
                    metrics=_extra_metrics(
                        record,
                        pair_prefix,
                        {"PENALTY", "PRODUCT_SIZE"},
                    ),
                )
            )
    return tuple(assays), details


def _split_records(text: str) -> list[dict[str, str]]:
    if not isinstance(text, str) or not text:
        raise PrimerDesignError("Primer3 output contains malformed Boulder-IO.")
    records: list[dict[str, str]] = []
    current: dict[str, str] = {}
    for line in text.splitlines():
        if line == "=":
            if not current:
                raise PrimerDesignError(
                    "Primer3 output contains malformed Boulder-IO: empty record."
                )
            records.append(current)
            current = {}
            continue
        if "=" not in line:
            raise PrimerDesignError(
                "Primer3 output contains malformed Boulder-IO: "
                "record line has no delimiter."
            )
        key, value = line.split("=", 1)
        if not key or key in current:
            raise PrimerDesignError(
                "Primer3 output contains malformed Boulder-IO: "
                "empty or duplicate tag."
            )
        current[key] = value
    if current:
        raise PrimerDesignError(
            "Primer3 output contains malformed Boulder-IO: unterminated record."
        )
    if not records:
        raise PrimerDesignError("Primer3 output contains malformed Boulder-IO.")
    return records


def _parse_oligo(
    record: dict[str, str], kind: str, index: int
) -> DesignedOligo:
    location_key = f"PRIMER_{kind}_{index}"
    location_parts = _required_tag(record, location_key).split(",")
    if len(location_parts) != 2:
        raise PrimerDesignError(
            f"Primer3 output tag '{location_key}' must be a position,length pair."
        )
    position_text, length_text = location_parts
    position = _parse_int_text(position_text, location_key)
    length = _parse_int_text(length_text, location_key)
    if length <= 0:
        raise PrimerDesignError(
            f"Primer3 output tag '{location_key}' must have a positive length."
        )
    if kind == "RIGHT":
        reference_start = position - length + 2
        reference_end = position + 1
    else:
        reference_start = position + 1
        reference_end = position + length
    prefix = f"{location_key}_"
    sequence_key = f"{prefix}SEQUENCE"
    sequence = _required_tag(record, sequence_key)
    if len(sequence) != length:
        raise PrimerDesignError(
            f"Primer3 output tag '{sequence_key}' length does not match "
            f"the reported oligo length {length}."
        )
    penalty_text = record.get(f"{prefix}PENALTY")
    return DesignedOligo(
        sequence=sequence,
        reference_start=reference_start,
        reference_end=reference_end,
        length=length,
        tm=_parse_float(record, f"{prefix}TM"),
        gc_percent=_parse_float(record, f"{prefix}GC_PERCENT"),
        penalty=(
            _parse_float(record, f"{prefix}PENALTY")
            if penalty_text is not None
            else None
        ),
        metrics=_extra_metrics(
            record,
            prefix,
            {"GC_PERCENT", "PENALTY", "SEQUENCE", "TM"},
        ),
    )


def _extra_metrics(
    record: dict[str, str], prefix: str, known_suffixes: set[str]
) -> tuple[tuple[str, str], ...]:
    return tuple(
        sorted(
            (key, value)
            for key, value in record.items()
            if key.startswith(prefix) and key.removeprefix(prefix) not in known_suffixes
        )
    )


def _validate_oligo_bounds(
    oligo: DesignedOligo,
    candidate: CandidateRegion,
    kind: str,
    index: int,
) -> None:
    if not (
        candidate.reference_start
        <= oligo.reference_start
        <= oligo.reference_end
        <= candidate.reference_end
    ):
        raise PrimerDesignError(
            f"Primer3 output tag 'PRIMER_{kind}_{index}' has coordinates "
            f"{oligo.reference_start}..{oligo.reference_end} outside candidate "
            f"'{candidate.region_id}'."
        )


def _required_tag(record: dict[str, str], key: str) -> str:
    try:
        return record[key]
    except KeyError as error:
        raise PrimerDesignError(
            f"Primer3 output is missing required tag '{key}'."
        ) from error


def _parse_int(record: dict[str, str], key: str) -> int:
    return _parse_int_text(_required_tag(record, key), key)


def _parse_int_text(value: str, key: str) -> int:
    try:
        return int(value)
    except ValueError as error:
        raise PrimerDesignError(
            f"Primer3 output has an invalid numeric value for '{key}'."
        ) from error


def _parse_float(record: dict[str, str], key: str) -> float:
    try:
        value = float(_required_tag(record, key))
    except ValueError as error:
        raise PrimerDesignError(
            f"Primer3 output has an invalid numeric value for '{key}'."
        ) from error
    if not math.isfinite(value):
        raise PrimerDesignError(
            f"Primer3 output tag '{key}' must have a finite numeric value."
        )
    return value
