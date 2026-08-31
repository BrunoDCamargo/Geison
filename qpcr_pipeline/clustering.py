"""Deterministic CD-HIT clustering artifacts and strict cluster parsing."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
import uuid
from collections import Counter
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Literal, Protocol

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import ClusteringConfig, validate_clustering_config
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet, EvaluationSet


class CdHitError(RuntimeError):
    """Raised when CD-HIT output cannot produce a traceable clustering result."""


class CdHitRunner(Protocol):
    def run(
        self,
        input_path: Path,
        output_path: Path,
        config: ClusteringConfig,
    ) -> None: ...


class SubprocessCdHitRunner:
    """Run CD-HIT-EST through a structured subprocess boundary."""

    def __init__(self, executable: str = "cd-hit-est") -> None:
        self.executable = executable

    def run(
        self,
        input_path: Path,
        output_path: Path,
        config: ClusteringConfig,
    ) -> None:
        executable = shutil.which(self.executable)
        if executable is None:
            raise CdHitError(
                f"{self.executable!r} was not found on PATH; install CD-HIT or disable clustering."
            )

        word_length = derive_word_length(config.identity)
        args = [
            executable,
            "-i", str(input_path),
            "-o", str(output_path),
            "-c", str(config.identity),
            "-n", str(word_length),
            "-l", str(word_length - 1),
            "-d", "0",
            "-g", "1",
            "-r", "0",
            "-T", str(config.threads),
            "-M", str(config.memory_mb),
        ]
        completed = subprocess.run(
            args,
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            stderr = " ".join(completed.stderr.split())[:2_000]
            detail = f": {stderr}" if stderr else ""
            raise CdHitError(
                f"cd-hit-est failed with exit code {completed.returncode}{detail}"
            )

        raw_cluster_path = Path(str(output_path) + ".clstr")
        if not output_path.is_file() or not raw_cluster_path.is_file():
            raise CdHitError(
                "cd-hit-est did not produce representative FASTA and .clstr output."
            )


@dataclass(frozen=True, slots=True)
class ClusterMember:
    sequence_id: str
    representative: bool
    identity: float | None
    strand: Literal["+", "-"] | None


@dataclass(frozen=True, slots=True)
class SequenceCluster:
    cluster_id: str
    representative_id: str
    members: tuple[ClusterMember, ...]


@dataclass(frozen=True, slots=True)
class ClusteringResult:
    discovery_set: DiscoverySet
    clusters: tuple[SequenceCluster, ...]
    discovery_fasta_path: Path
    report_path: Path
    raw_cluster_path: Path | None


@dataclass(frozen=True, slots=True)
class _ParsedMember:
    internal_id: str
    representative: bool
    identity: float | None
    strand: Literal["+", "-"] | None


@dataclass(frozen=True, slots=True)
class _ParsedCluster:
    raw_cluster_id: str
    members: tuple[_ParsedMember, ...]


_HEADER_PATTERN = re.compile(r"^>Cluster (?P<cluster_id>[0-9]+)$")
_MEMBER_PATTERN = re.compile(
    r"^(?P<ordinal>[0-9]+)[ \t]+(?P<length>[0-9]+)nt, "
    r">(?P<internal_id>[A-Za-z0-9][A-Za-z0-9._-]*)\.\.\. "
    r"(?:(?P<representative>\*)|at (?P<strand>[+-])/(?P<identity>[0-9]+(?:\.[0-9]+)?)%)$"
)


def derive_word_length(identity: float) -> int:
    """Return the CD-HIT-EST word length compatible with ``identity``."""
    if not 0.80 <= identity <= 1.0:
        raise ValueError("CD-HIT identity must be between 0.80 and 1.0.")
    if identity >= 0.95:
        return 10
    if identity >= 0.90:
        return 8
    if identity >= 0.88:
        return 7
    if identity >= 0.85:
        return 6
    return 5


def cluster_sequences(
    records: tuple[LocalSequenceRecord, ...],
    evaluation_set: EvaluationSet,
    config: ClusteringConfig,
    output_dir: Path,
    *,
    runner: CdHitRunner | None = None,
) -> ClusteringResult:
    """Create a Discovery Set and its traceable clustering artifacts."""
    validate_clustering_config(config)
    output_dir = Path(output_dir)
    evaluation_ids = evaluation_set.sequence_ids
    records_by_id = _validated_records(records, evaluation_ids)

    if not config.enabled:
        discovery_ids = evaluation_ids
        clusters: tuple[SequenceCluster, ...] = ()
        raw_cluster_text: str | None = None
    elif not evaluation_ids:
        discovery_ids = evaluation_ids
        clusters = ()
        raw_cluster_text = ""
    else:
        _validate_minimum_lengths(
            records_by_id,
            evaluation_ids,
            derive_word_length(config.identity),
        )
        if runner is None:
            runner = SubprocessCdHitRunner()
        discovery_ids, clusters, raw_cluster_text = _run_and_parse(
            records_by_id, evaluation_ids, config, runner
        )

    discovery_set = DiscoverySet(sequence_ids=discovery_ids)
    discovery_fasta_path = output_dir / "discovery_set.fasta"
    report_path = output_dir / "clustering_report.json"
    raw_cluster_path = (
        output_dir / "clustering" / "cd-hit-est.clstr"
        if raw_cluster_text is not None
        else None
    )
    report = _report(
        evaluation_ids,
        discovery_set,
        clusters,
        config,
        raw_cluster_path is not None,
    )

    _publish_artifacts(
        discovery_fasta_path,
        report_path,
        raw_cluster_path,
        tuple(records_by_id[sequence_id] for sequence_id in discovery_ids),
        raw_cluster_text,
        report,
    )
    return ClusteringResult(
        discovery_set=discovery_set,
        clusters=clusters,
        discovery_fasta_path=discovery_fasta_path,
        report_path=report_path,
        raw_cluster_path=raw_cluster_path,
    )


def _validated_records(
    records: tuple[LocalSequenceRecord, ...], evaluation_ids: tuple[str, ...]
) -> dict[str, LocalSequenceRecord]:
    if len(evaluation_ids) != len(set(evaluation_ids)):
        raise CdHitError("Duplicate approved sequence IDs cannot be clustered.")
    record_ids = tuple(record.sequence_id for record in records)
    if len(record_ids) != len(set(record_ids)):
        raise CdHitError("Duplicate approved sequence IDs cannot be clustered.")
    if set(record_ids) != set(evaluation_ids):
        raise CdHitError(
            "Approved records must contain exactly the Evaluation Set sequence IDs."
        )
    return {record.sequence_id: record for record in records}


def _validate_minimum_lengths(
    records_by_id: dict[str, LocalSequenceRecord],
    evaluation_ids: tuple[str, ...],
    minimum_length: int,
) -> None:
    for sequence_id in evaluation_ids:
        sequence_length = len(records_by_id[sequence_id].sequence)
        if sequence_length < minimum_length:
            raise CdHitError(
                f"Approved sequence {sequence_id!r} has length {sequence_length} nt; "
                f"CD-HIT-EST requires a minimum length of {minimum_length} nt "
                "for the configured identity."
            )


def _run_and_parse(
    records_by_id: dict[str, LocalSequenceRecord],
    evaluation_ids: tuple[str, ...],
    config: ClusteringConfig,
    runner: CdHitRunner,
) -> tuple[tuple[str, ...], tuple[SequenceCluster, ...], str]:
    internal_ids = tuple(f"geison-{position:08d}" for position in range(len(evaluation_ids)))
    internal_to_original = dict(zip(internal_ids, evaluation_ids, strict=True))

    with tempfile.TemporaryDirectory(prefix="geison-clustering-") as temporary_directory:
        temporary_path = Path(temporary_directory)
        input_path = temporary_path / "input.fasta"
        representatives_path = temporary_path / "representatives.fasta"
        _write_fasta(
            input_path,
            tuple(
                (internal_id, records_by_id[original_id].sequence)
                for internal_id, original_id in internal_to_original.items()
            ),
        )
        runner.run(input_path, representatives_path, config)
        raw_cluster_temp_path = Path(str(representatives_path) + ".clstr")
        if not representatives_path.is_file() or not raw_cluster_temp_path.is_file():
            raise CdHitError("CD-HIT runner did not produce representative FASTA and .clstr output.")
        raw_cluster_text = raw_cluster_temp_path.read_text(encoding="utf-8")
        parsed_clusters = _parse_clusters(raw_cluster_text, internal_ids)
        representative_ids = _representative_ids(representatives_path, internal_ids)

    cluster_representative_ids = {
        member.internal_id
        for cluster in parsed_clusters
        for member in cluster.members
        if member.representative
    }
    if representative_ids != cluster_representative_ids:
        raise CdHitError(
            "Representative FASTA IDs do not match representatives in CD-HIT clusters."
        )

    evaluation_positions = {sequence_id: position for position, sequence_id in enumerate(evaluation_ids)}
    sorted_clusters = sorted(
        parsed_clusters,
        key=lambda cluster: min(
            evaluation_positions[internal_to_original[member.internal_id]]
            for member in cluster.members
        ),
    )
    clusters = tuple(
        SequenceCluster(
            cluster_id=f"cluster-{position:05d}",
            representative_id=internal_to_original[
                next(member.internal_id for member in cluster.members if member.representative)
            ],
            members=tuple(
                ClusterMember(
                    sequence_id=internal_to_original[member.internal_id],
                    representative=member.representative,
                    identity=member.identity,
                    strand=member.strand,
                )
                for member in cluster.members
            ),
        )
        for position, cluster in enumerate(sorted_clusters)
    )
    discovery_ids = tuple(
        sequence_id
        for sequence_id in evaluation_ids
        if sequence_id in {
            cluster.representative_id for cluster in clusters
        }
    )
    return discovery_ids, clusters, raw_cluster_text


def _parse_clusters(raw_text: str, expected_internal_ids: tuple[str, ...]) -> tuple[_ParsedCluster, ...]:
    clusters: list[_ParsedCluster] = []
    raw_cluster_ids: set[str] = set()
    current_cluster_id: str | None = None
    current_members: list[_ParsedMember] = []

    for line_number, line in enumerate(raw_text.splitlines(), start=1):
        header = _HEADER_PATTERN.fullmatch(line)
        if header is not None:
            if current_cluster_id is not None:
                clusters.append(_finish_cluster(current_cluster_id, current_members))
            current_cluster_id = header.group("cluster_id")
            if current_cluster_id in raw_cluster_ids:
                raise CdHitError(f"Duplicate CD-HIT cluster header at line {line_number}.")
            raw_cluster_ids.add(current_cluster_id)
            current_members = []
            continue

        member = _MEMBER_PATTERN.fullmatch(line)
        if member is None or current_cluster_id is None:
            raise CdHitError(f"Malformed CD-HIT cluster line {line_number}.")
        representative = member.group("representative") is not None
        current_members.append(
            _ParsedMember(
                internal_id=member.group("internal_id"),
                representative=representative,
                identity=None if representative else float(member.group("identity")),
                strand=None if representative else member.group("strand"),
            )
        )

    if current_cluster_id is None:
        raise CdHitError("CD-HIT cluster output is empty or missing a cluster header.")
    clusters.append(_finish_cluster(current_cluster_id, current_members))

    parsed_ids = [member.internal_id for cluster in clusters for member in cluster.members]
    expected_ids = set(expected_internal_ids)
    occurrences = Counter(parsed_ids)
    if set(parsed_ids) != expected_ids:
        raise CdHitError("CD-HIT cluster members do not match the expected internal IDs.")
    if any(occurrences[internal_id] != 1 for internal_id in expected_ids):
        raise CdHitError("CD-HIT cluster members must occur exactly once.")
    return tuple(clusters)


def _finish_cluster(
    raw_cluster_id: str, members: list[_ParsedMember]
) -> _ParsedCluster:
    representatives = [member for member in members if member.representative]
    if len(representatives) != 1:
        raise CdHitError(
            f"CD-HIT cluster {raw_cluster_id} must contain exactly one representative."
        )
    return _ParsedCluster(raw_cluster_id=raw_cluster_id, members=tuple(members))


def _representative_ids(path: Path, expected_internal_ids: tuple[str, ...]) -> set[str]:
    with path.open(encoding="utf-8") as handle:
        representative_ids = [record.id for record in SeqIO.parse(handle, "fasta")]
    if len(representative_ids) != len(set(representative_ids)):
        raise CdHitError("Representative FASTA contains duplicate internal IDs.")
    unknown_ids = set(representative_ids) - set(expected_internal_ids)
    if unknown_ids:
        raise CdHitError("Representative FASTA contains unknown internal IDs.")
    return set(representative_ids)


def _report(
    evaluation_ids: tuple[str, ...],
    discovery_set: DiscoverySet,
    clusters: tuple[SequenceCluster, ...],
    config: ClusteringConfig,
    has_raw_cluster: bool,
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "status": "COMPLETE",
        "clustering_enabled": config.enabled,
        "tool": "cd-hit-est",
        "parameters": {
            "identity": config.identity,
            "word_length": derive_word_length(config.identity),
            "threads": config.threads,
            "memory_mb": config.memory_mb,
        },
        "evaluation_set": {"sequence_ids": list(evaluation_ids)},
        "discovery_set": {"sequence_ids": list(discovery_set.sequence_ids)},
        "counts": {
            "evaluation": len(evaluation_ids),
            "discovery": len(discovery_set.sequence_ids),
        },
        "clusters": [
            {
                "cluster_id": cluster.cluster_id,
                "representative_id": cluster.representative_id,
                "members": [
                    {
                        "sequence_id": member.sequence_id,
                        "representative": member.representative,
                        "identity": member.identity,
                        "strand": member.strand,
                    }
                    for member in cluster.members
                ],
            }
            for cluster in clusters
        ],
        "artifacts": {
            "discovery_fasta": "discovery_set.fasta",
            "raw_cluster": "clustering/cd-hit-est.clstr" if has_raw_cluster else None,
        },
    }


def _publish_artifacts(
    discovery_fasta_path: Path,
    report_path: Path,
    raw_cluster_path: Path | None,
    discovery_records: tuple[LocalSequenceRecord, ...],
    raw_cluster_text: str | None,
    report: dict[str, object],
) -> None:
    discovery_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    if raw_cluster_path is not None:
        raw_cluster_path.parent.mkdir(parents=True, exist_ok=True)

    fasta_text = _fasta_text(
        tuple((record.sequence_id, record.sequence) for record in discovery_records)
    )
    _atomic_write_text(discovery_fasta_path, fasta_text)
    if raw_cluster_path is not None and raw_cluster_text is not None:
        _atomic_write_text(raw_cluster_path, raw_cluster_text)
    else:
        (report_path.parent / "clustering" / "cd-hit-est.clstr").unlink(
            missing_ok=True
        )
    _atomic_write_text(
        report_path, json.dumps(report, indent=2, sort_keys=True) + "\n"
    )


def _write_fasta(path: Path, records: tuple[tuple[str, str], ...]) -> None:
    SeqIO.write(
        [SeqRecord(Seq(sequence), id=sequence_id, description="") for sequence_id, sequence in records],
        path,
        "fasta",
    )


def _fasta_text(records: tuple[tuple[str, str], ...]) -> str:
    handle = StringIO()
    SeqIO.write(
        [SeqRecord(Seq(sequence), id=sequence_id, description="") for sequence_id, sequence in records],
        handle,
        "fasta",
    )
    return handle.getvalue()


def _atomic_write_text(destination: Path, content: str) -> None:
    temporary_path = destination.parent / (
        f".{destination.name}.{uuid.uuid4().hex}.tmp"
    )
    try:
        with temporary_path.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        temporary_path.replace(destination)
    finally:
        temporary_path.unlink(missing_ok=True)
