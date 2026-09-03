"""Deterministic, traceable alignment of a Discovery Set."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import uuid
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Literal, Protocol

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import AlignmentConfig, validate_alignment_config
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet


class MafftError(RuntimeError):
    """Raised when MAFFT cannot produce a valid, traceable alignment."""


class MafftRunner(Protocol):
    """Boundary used by the deterministic alignment service."""

    def run(
        self,
        input_path: Path,
        output_path: Path,
        config: AlignmentConfig,
    ) -> None: ...


class SubprocessMafftRunner:
    """Execute MAFFT using the fixed, non-shell alignment command."""

    def __init__(self, executable: str = "mafft"):
        self.executable = executable

    def run(
        self,
        input_path: Path,
        output_path: Path,
        config: AlignmentConfig,
    ) -> None:
        executable = shutil.which(self.executable)
        if executable is None:
            raise MafftError(
                f"{self.executable!r} was not found on PATH; install MAFFT or disable alignment."
            )
        completed = subprocess.run(
            [
                executable,
                "--auto",
                "--nuc",
                "--inputorder",
                "--adjustdirection",
                "--thread",
                str(config.threads),
                "--threadit",
                "0",
                "--quiet",
                str(input_path),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            stderr = " ".join(completed.stderr.split())[:2_000]
            raise MafftError(f"MAFFT exited with status {completed.returncode}: {stderr}")
        if not completed.stdout.strip():
            raise MafftError("MAFFT produced empty alignment output.")
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            handle.write(completed.stdout)


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


_VALID_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_AUTOMATIC_SELECTION_RULE = "lowest_ambiguity_fraction_then_longest_then_discovery_order"
_COORDINATE_HEADER = "alignment_position\treference_position\treference_base\n"


def align_discovery(
    records: tuple[LocalSequenceRecord, ...],
    discovery_set: DiscoverySet,
    config: AlignmentConfig,
    output_dir: Path,
    *,
    runner: MafftRunner | None = None,
) -> AlignmentResult:
    """Align unchanged Discovery records and publish auditable artifacts.

    Inputs are validated before any output is altered.  A subprocess runner is
    constructed only when enabled alignment needs it and none was injected.
    """
    validate_alignment_config(config)
    output_dir = Path(output_dir)
    records_by_id = _validated_records(records, discovery_set)
    report_path = output_dir / "alignment" / "alignment_report.json"
    fasta_path = output_dir / "alignment" / "discovery_alignment.fasta"
    coordinate_path = output_dir / "alignment" / "coordinate_map.tsv"

    if (
        config.reference_id is not None
        and config.reference_id not in discovery_set.sequence_ids
    ):
        raise MafftError(
            f"Configured alignment reference {config.reference_id!r} is not in the Discovery Set."
        )

    if not config.enabled:
        _publish_skipped(report_path, fasta_path, coordinate_path, discovery_set, config)
        return AlignmentResult(
            status="SKIPPED",
            discovery_set=discovery_set,
            reference_id=None,
            reference_mode=None,
            sequences=(),
            coordinates=(),
            alignment_fasta_path=None,
            coordinate_map_path=None,
            report_path=report_path,
        )

    mafft_executed = False
    if not discovery_set.sequence_ids:
        sequences: tuple[AlignedSequence, ...] = ()
        coordinates: tuple[AlignmentCoordinate, ...] = ()
        reference_id = None
        reference_mode: Literal["explicit", "automatic"] | None = None
    else:
        reference_id, reference_mode = _select_reference(
            records_by_id, discovery_set, config.reference_id
        )
        if len(discovery_set.sequence_ids) == 1:
            only_record = records_by_id[reference_id]
            sequences = (AlignedSequence(reference_id, only_record.sequence, "forward"),)
        else:
            if runner is None:
                runner = SubprocessMafftRunner()
            sequences = _run_and_parse(records_by_id, discovery_set, reference_id, config, runner)
            mafft_executed = True
        reference_sequence = next(
            sequence.aligned_sequence
            for sequence in sequences
            if sequence.sequence_id == reference_id
        )
        coordinates = _coordinates(reference_sequence, records_by_id[reference_id].sequence)

    report = _report(
        status="COMPLETE",
        discovery_set=discovery_set,
        config=config,
        reference_id=reference_id,
        reference_mode=reference_mode,
        sequences=sequences,
        coordinates=coordinates,
        mafft_executed=mafft_executed,
        artifacts={
            "alignment_fasta": "alignment/discovery_alignment.fasta",
            "coordinate_map": "alignment/coordinate_map.tsv",
        },
    )
    _publish_complete(
        fasta_path,
        coordinate_path,
        report_path,
        sequences,
        coordinates,
        report,
    )
    return AlignmentResult(
        status="COMPLETE",
        discovery_set=discovery_set,
        reference_id=reference_id,
        reference_mode=reference_mode,
        sequences=sequences,
        coordinates=coordinates,
        alignment_fasta_path=fasta_path,
        coordinate_map_path=coordinate_path,
        report_path=report_path,
    )


def _validated_records(
    records: tuple[LocalSequenceRecord, ...], discovery_set: DiscoverySet
) -> dict[str, LocalSequenceRecord]:
    if not isinstance(discovery_set, DiscoverySet):
        raise MafftError("Discovery Set must be a DiscoverySet.")
    discovery_ids = discovery_set.sequence_ids
    if len(discovery_ids) != len(set(discovery_ids)):
        raise MafftError("Discovery Set contains duplicate sequence IDs.")
    if any(not isinstance(sequence_id, str) or not sequence_id for sequence_id in discovery_ids):
        raise MafftError("Discovery Set sequence IDs must be non-empty strings.")

    record_ids: list[str] = []
    for record in records:
        if not isinstance(record, LocalSequenceRecord):
            raise MafftError("Discovery records must be LocalSequenceRecord instances.")
        if not isinstance(record.sequence_id, str) or not record.sequence_id:
            raise MafftError("Discovery record IDs must be non-empty strings.")
        if not isinstance(record.sequence, str) or not record.sequence:
            raise MafftError(f"Discovery record {record.sequence_id!r} has an empty sequence.")
        if any(base not in _VALID_IUPAC_DNA for base in record.sequence):
            raise MafftError(
                f"Discovery record {record.sequence_id!r} contains invalid IUPAC DNA symbols."
            )
        record_ids.append(record.sequence_id)
    if len(record_ids) != len(set(record_ids)):
        raise MafftError("Discovery records contain duplicate sequence IDs.")
    if set(record_ids) != set(discovery_ids):
        raise MafftError("Discovery records must contain exactly the Discovery Set sequence IDs.")
    return {record.sequence_id: record for record in records}


def _select_reference(
    records_by_id: dict[str, LocalSequenceRecord],
    discovery_set: DiscoverySet,
    explicit_reference_id: str | None,
) -> tuple[str, Literal["explicit", "automatic"]]:
    if explicit_reference_id is not None:
        if explicit_reference_id not in records_by_id:
            raise MafftError(
                f"Configured alignment reference {explicit_reference_id!r} is not in the Discovery Set."
            )
        return explicit_reference_id, "explicit"

    def selection_key(item: tuple[int, LocalSequenceRecord]) -> tuple[float, int, int]:
        position, record = item
        ambiguous = sum(base not in "ACGT" for base in record.sequence)
        return (ambiguous / len(record.sequence), -len(record.sequence), position)

    position, selected = min(
        (
            (position, records_by_id[sequence_id])
            for position, sequence_id in enumerate(discovery_set.sequence_ids)
        ),
        key=selection_key,
    )
    del position
    return selected.sequence_id, "automatic"


def _run_and_parse(
    records_by_id: dict[str, LocalSequenceRecord],
    discovery_set: DiscoverySet,
    reference_id: str,
    config: AlignmentConfig,
    runner: MafftRunner,
) -> tuple[AlignedSequence, ...]:
    original_ids = discovery_set.sequence_ids
    positions = {sequence_id: position for position, sequence_id in enumerate(original_ids)}
    internal_to_original = {
        f"geison-{position:08d}": sequence_id
        for sequence_id, position in positions.items()
    }
    ordered_original_ids = (reference_id,) + tuple(
        sequence_id for sequence_id in original_ids if sequence_id != reference_id
    )
    with tempfile.TemporaryDirectory(prefix="geison-alignment-") as temporary_directory:
        temporary_path = Path(temporary_directory)
        input_path = temporary_path / "input.fasta"
        output_path = temporary_path / "aligned.fasta"
        _write_fasta(
            input_path,
            tuple(
                (f"geison-{positions[sequence_id]:08d}", records_by_id[sequence_id].sequence)
                for sequence_id in ordered_original_ids
            ),
        )
        try:
            runner.run(input_path, output_path, config)
        except MafftError:
            raise
        except Exception as error:
            raise MafftError(f"MAFFT runner failed: {error}") from error
        if not output_path.is_file() or output_path.stat().st_size == 0:
            raise MafftError("MAFFT runner did not produce a non-empty aligned FASTA output.")
        parsed = _parse_aligned_fasta(output_path, internal_to_original, records_by_id, reference_id)

    return tuple(parsed[sequence_id] for sequence_id in original_ids)


def _parse_aligned_fasta(
    output_path: Path,
    internal_to_original: dict[str, str],
    records_by_id: dict[str, LocalSequenceRecord],
    reference_id: str,
) -> dict[str, AlignedSequence]:
    try:
        output_text = output_path.read_text(encoding="utf-8")
    except Exception as error:
        raise MafftError(f"MAFFT output is not valid FASTA: {error}") from error
    output_text = _normalize_mafft_fasta_case(output_text)
    _validate_raw_aligned_fasta(output_text)
    try:
        output_records = list(SeqIO.parse(StringIO(output_text), "fasta"))
    except Exception as error:
        raise MafftError(f"MAFFT output is not valid FASTA: {error}") from error
    if not output_records:
        raise MafftError("MAFFT output FASTA contains no records.")

    parsed: dict[str, AlignedSequence] = {}
    alignment_length: int | None = None
    for output_record in output_records:
        reversed_sequence, internal_id = _split_mafft_id(output_record.id)
        if internal_id not in internal_to_original:
            raise MafftError(f"MAFFT output contains unknown internal ID {output_record.id!r}.")
        original_id = internal_to_original[internal_id]
        if original_id in parsed:
            raise MafftError("MAFFT output must contain each expected internal ID exactly once.")
        if reversed_sequence and original_id == reference_id:
            raise MafftError("MAFFT reversed the selected reference sequence.")
        sequence = str(output_record.seq)
        if not sequence:
            raise MafftError(f"MAFFT output sequence {output_record.id!r} is empty.")
        if any(base != "-" and base not in _VALID_IUPAC_DNA for base in sequence):
            raise MafftError(f"MAFFT output sequence {output_record.id!r} has invalid IUPAC DNA symbols.")
        if not sequence.replace("-", ""):
            raise MafftError(f"MAFFT output sequence {output_record.id!r} is all gaps.")
        if alignment_length is None:
            alignment_length = len(sequence)
        elif len(sequence) != alignment_length:
            raise MafftError("MAFFT output sequences have unequal alignment lengths.")

        original_sequence = records_by_id[original_id].sequence
        ungapped = sequence.replace("-", "")
        expected_sequence = str(Seq(original_sequence).reverse_complement()) if reversed_sequence else original_sequence
        if ungapped != expected_sequence:
            direction = "reverse-complemented" if reversed_sequence else "forward"
            raise MafftError(
                f"MAFFT {direction} output for {original_id!r} does not match its original sequence."
            )
        parsed[original_id] = AlignedSequence(
            sequence_id=original_id,
            aligned_sequence=sequence,
            orientation="reverse_complemented" if reversed_sequence else "forward",
        )

    if set(parsed) != set(records_by_id):
        raise MafftError("MAFFT output does not contain exactly the expected internal IDs.")
    return parsed


def _normalize_mafft_fasta_case(output_text: str) -> str:
    return "".join(
        line if line.startswith(">") else line.upper()
        for line in output_text.splitlines(keepends=True)
    )


def _validate_raw_aligned_fasta(output_text: str) -> None:
    seen_header = False
    for line_number, line in enumerate(output_text.splitlines(), start=1):
        if line == "":
            continue
        if line.startswith(">"):
            header = line[1:]
            if not header or header[0].isspace():
                raise MafftError(
                    f"MAFFT output FASTA header on line {line_number} has no valid identifier."
                )
            seen_header = True
            continue
        if not seen_header:
            raise MafftError(
                f"MAFFT output FASTA has non-header text before the first header on line {line_number}."
            )
        if any(base != "-" and base not in _VALID_IUPAC_DNA for base in line):
            raise MafftError(
                f"MAFFT output FASTA sequence line {line_number} contains whitespace or invalid symbols."
            )


def _split_mafft_id(output_id: str) -> tuple[bool, str]:
    if output_id.startswith("_R_"):
        internal_id = output_id[3:]
        if internal_id.startswith("_R_"):
            raise MafftError("MAFFT output ID contains multiple _R_ direction prefixes.")
        return True, internal_id
    return False, output_id


def _coordinates(aligned_reference: str, original_reference: str) -> tuple[AlignmentCoordinate, ...]:
    reference_position = 0
    coordinates: list[AlignmentCoordinate] = []
    for alignment_position, base in enumerate(aligned_reference, start=1):
        if base == "-":
            coordinates.append(AlignmentCoordinate(alignment_position, None, None))
        else:
            reference_position += 1
            coordinates.append(AlignmentCoordinate(alignment_position, reference_position, base))
    if reference_position != len(original_reference):
        raise MafftError("Aligned reference length does not match the original reference sequence.")
    return tuple(coordinates)


def _report(
    *,
    status: Literal["SKIPPED", "COMPLETE"],
    discovery_set: DiscoverySet,
    config: AlignmentConfig,
    reference_id: str | None,
    reference_mode: Literal["explicit", "automatic"] | None,
    sequences: tuple[AlignedSequence, ...],
    coordinates: tuple[AlignmentCoordinate, ...],
    mafft_executed: bool,
    artifacts: dict[str, str | None],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "status": status,
        "enabled": config.enabled,
        "tool": {
            "name": "mafft",
            "executed": mafft_executed,
            "parameters": {
                "strategy": "auto",
                "nucleotide_mode": True,
                "input_order": True,
                "adjust_direction": True,
                "threads": config.threads,
                "iterative_refinement_threads": 0,
                "quiet": True,
            },
        },
        "discovery_set_ids": list(discovery_set.sequence_ids),
        "reference": {
            "id": reference_id,
            "mode": reference_mode,
            "automatic_selection_rule": (
                _AUTOMATIC_SELECTION_RULE if reference_mode == "automatic" else None
            ),
        },
        "orientations": [
            {"sequence_id": sequence.sequence_id, "orientation": sequence.orientation}
            for sequence in sequences
        ],
        "counts": {
            "discovery": len(discovery_set.sequence_ids),
            "alignment_length": len(coordinates),
            "reverse_complemented": sum(
                sequence.orientation == "reverse_complemented" for sequence in sequences
            ),
        },
        "artifacts": artifacts,
    }


def _publish_complete(
    fasta_path: Path,
    coordinate_path: Path,
    report_path: Path,
    sequences: tuple[AlignedSequence, ...],
    coordinates: tuple[AlignmentCoordinate, ...],
    report: dict[str, object],
) -> None:
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.unlink(missing_ok=True)
    _atomic_write_text(fasta_path, _fasta_text(sequences))
    _atomic_write_text(coordinate_path, _coordinate_text(coordinates))
    _atomic_write_text(report_path, json.dumps(report, indent=2, sort_keys=True) + "\n")


def _publish_skipped(
    report_path: Path,
    fasta_path: Path,
    coordinate_path: Path,
    discovery_set: DiscoverySet,
    config: AlignmentConfig,
) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.unlink(missing_ok=True)
    fasta_path.unlink(missing_ok=True)
    coordinate_path.unlink(missing_ok=True)
    _atomic_write_text(
        report_path,
        json.dumps(
            _report(
                status="SKIPPED",
                discovery_set=discovery_set,
                config=config,
                reference_id=None,
                reference_mode=None,
                sequences=(),
                coordinates=(),
                mafft_executed=False,
                artifacts={"alignment_fasta": None, "coordinate_map": None},
            ),
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )


def _write_fasta(path: Path, records: tuple[tuple[str, str], ...]) -> None:
    SeqIO.write(
        [SeqRecord(Seq(sequence), id=sequence_id, description="") for sequence_id, sequence in records],
        path,
        "fasta",
    )


def _fasta_text(sequences: tuple[AlignedSequence, ...]) -> str:
    handle = StringIO()
    SeqIO.write(
        [
            SeqRecord(Seq(sequence.aligned_sequence), id=sequence.sequence_id, description="")
            for sequence in sequences
        ],
        handle,
        "fasta",
    )
    return handle.getvalue()


def _coordinate_text(coordinates: tuple[AlignmentCoordinate, ...]) -> str:
    return _COORDINATE_HEADER + "".join(
        f"{coordinate.alignment_position}\t"
        f"{'' if coordinate.reference_position is None else coordinate.reference_position}\t"
        f"{'' if coordinate.reference_base is None else coordinate.reference_base}\n"
        for coordinate in coordinates
    )


def _atomic_write_text(destination: Path, content: str) -> None:
    temporary_path = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.tmp"
    try:
        with temporary_path.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        temporary_path.replace(destination)
    finally:
        temporary_path.unlink(missing_ok=True)