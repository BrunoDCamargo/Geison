"""Strict typed JSON codecs for checkpointed pipeline stage state."""

from __future__ import annotations

from collections.abc import Mapping as MappingABC
from dataclasses import fields, is_dataclass
from enum import Enum
import math
from pathlib import Path
import types
from typing import Any, Generic, Literal, Mapping, TypeVar, Union, get_args, get_origin, get_type_hints

from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

from qpcr_pipeline.alignment import AlignmentResult
from qpcr_pipeline.clustering import ClusteringResult
from qpcr_pipeline.conservation import ConservationResult
from qpcr_pipeline.contrastive_conservation import ContrastiveConservationResult
from qpcr_pipeline.inclusivity import InclusivityResult
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.panel_manifest import PanelResult
from qpcr_pipeline.primer_design import PrimerDesignResult
from qpcr_pipeline.qc import QCResult
from qpcr_pipeline.ranking import RankingResult
from qpcr_pipeline.specificity import SpecificityResult


T = TypeVar("T")


def _relative_path(path: Path, outdir: Path) -> str:
    resolved_outdir = Path(outdir).resolve()
    resolved = Path(path).resolve()
    try:
        relative = resolved.relative_to(resolved_outdir)
    except ValueError as error:
        raise ValueError("checkpoint state path escapes the output directory") from error
    return relative.as_posix()


def _resolved_path(value: object, outdir: Path) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError("checkpoint state path must be a non-empty relative string")
    path = Path(value)
    if path.is_absolute():
        raise ValueError("checkpoint state path must be relative")
    resolved_outdir = Path(outdir).resolve()
    resolved = (Path(outdir) / path).resolve()
    try:
        resolved.relative_to(resolved_outdir)
    except ValueError as error:
        raise ValueError("checkpoint state path escapes the output directory") from error
    return resolved


def _feature_payload(feature: SeqFeature) -> dict[str, object]:
    if not isinstance(feature, SeqFeature):
        raise ValueError("record feature must be a SeqFeature")
    parts: list[dict[str, object]] = []
    if feature.location is not None:
        for part in feature.location.parts:
            parts.append(
                {
                    "start": int(part.start),
                    "end": int(part.end),
                    "strand": part.strand,
                    "ref": part.ref,
                }
            )
    qualifiers: dict[str, list[str]] = {}
    raw_qualifiers = feature.qualifiers if isinstance(feature.qualifiers, MappingABC) else {}
    for key in ("gene", "locus_tag", "product"):
        raw = raw_qualifiers.get(key, [])
        if isinstance(raw, str):
            values = [raw]
        elif isinstance(raw, (list, tuple)) and all(isinstance(item, str) for item in raw):
            values = list(raw)
        else:
            values = []
        qualifiers[key] = values
    return {
        "type": str(feature.type),
        "parts": parts,
        "qualifiers": qualifiers,
    }


def _decode_feature(payload: object) -> SeqFeature:
    if not isinstance(payload, dict) or set(payload) != {"type", "parts", "qualifiers"}:
        raise ValueError("record feature fields are invalid")
    feature_type = payload["type"]
    raw_parts = payload["parts"]
    raw_qualifiers = payload["qualifiers"]
    if not isinstance(feature_type, str) or not isinstance(raw_parts, list):
        raise ValueError("record feature structure is invalid")
    if not isinstance(raw_qualifiers, dict) or set(raw_qualifiers) != {"gene", "locus_tag", "product"}:
        raise ValueError("record feature qualifiers are invalid")
    qualifiers: dict[str, list[str]] = {}
    for key, raw in raw_qualifiers.items():
        if not isinstance(raw, list) or any(not isinstance(item, str) for item in raw):
            raise ValueError("record feature qualifier values are invalid")
        qualifiers[key] = list(raw)

    parts: list[SimpleLocation] = []
    for raw in raw_parts:
        if not isinstance(raw, dict) or set(raw) != {"start", "end", "strand", "ref"}:
            raise ValueError("record feature location fields are invalid")
        start = raw["start"]
        end = raw["end"]
        strand = raw["strand"]
        ref = raw["ref"]
        if type(start) is not int or type(end) is not int or start < 0 or end < start:
            raise ValueError("record feature location coordinates are invalid")
        if strand not in {-1, 0, 1, None}:
            raise ValueError("record feature strand is invalid")
        if ref is not None and not isinstance(ref, str):
            raise ValueError("record feature reference is invalid")
        parts.append(SimpleLocation(start, end, strand=strand, ref=ref))

    if not parts:
        location = None
    elif len(parts) == 1:
        location = parts[0]
    else:
        location = CompoundLocation(parts, operator="join")
    return SeqFeature(location=location, type=feature_type, qualifiers=qualifiers)


def _record_payload(record: LocalSequenceRecord) -> dict[str, object]:
    if not isinstance(record, LocalSequenceRecord):
        raise ValueError("input checkpoint must contain LocalSequenceRecord values")
    metadata = record.metadata if isinstance(record.metadata, MappingABC) else {}
    raw_features = metadata.get("features", ())
    if not isinstance(raw_features, (tuple, list)):
        raise ValueError("record features metadata must be a sequence")
    return {
        "sequence_id": record.sequence_id,
        "sequence": record.sequence,
        "metadata": {
            "features": [_feature_payload(feature) for feature in raw_features],
        },
    }


def _decode_record(payload: object) -> LocalSequenceRecord:
    if not isinstance(payload, dict) or set(payload) != {"sequence_id", "sequence", "metadata"}:
        raise ValueError("input record fields are invalid")
    sequence_id = payload["sequence_id"]
    sequence = payload["sequence"]
    metadata = payload["metadata"]
    if not isinstance(sequence_id, str) or not sequence_id:
        raise ValueError("input record sequence_id is invalid")
    if not isinstance(sequence, str):
        raise ValueError("input record sequence is invalid")
    if not isinstance(metadata, dict) or set(metadata) != {"features"}:
        raise ValueError("input record metadata fields are invalid")
    raw_features = metadata["features"]
    if not isinstance(raw_features, list):
        raise ValueError("input record features are invalid")
    return LocalSequenceRecord(
        sequence_id=sequence_id,
        sequence=sequence,
        metadata={"features": tuple(_decode_feature(feature) for feature in raw_features)},
    )


def _encode(value: object, outdir: Path) -> object:
    if isinstance(value, LocalSequenceRecord):
        return _record_payload(value)
    if isinstance(value, Path):
        return _relative_path(value, outdir)
    if isinstance(value, Enum):
        return value.value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("checkpoint state cannot contain non-finite floats")
        return value
    if is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: _encode(getattr(value, field.name), outdir)
            for field in fields(value)
        }
    if isinstance(value, tuple):
        return [_encode(item, outdir) for item in value]
    if isinstance(value, list):
        return [_encode(item, outdir) for item in value]
    if isinstance(value, MappingABC):
        if any(not isinstance(key, str) for key in value):
            raise ValueError("checkpoint state mapping keys must be strings")
        return {key: _encode(item, outdir) for key, item in value.items()}
    raise ValueError(f"unsupported checkpoint state value: {type(value).__name__}")


def _decode_object(payload: object) -> object:
    if payload is None or isinstance(payload, (str, bool, int, float)):
        return payload
    if isinstance(payload, list):
        return [_decode_object(item) for item in payload]
    if isinstance(payload, dict):
        if any(not isinstance(key, str) for key in payload):
            raise ValueError("checkpoint state mapping keys must be strings")
        return {key: _decode_object(item) for key, item in payload.items()}
    raise ValueError("checkpoint state contains a non-JSON value")


def _decode(payload: object, annotation: object, outdir: Path) -> object:
    if annotation is Any or annotation is object:
        return _decode_object(payload)
    if annotation is type(None):
        if payload is not None:
            raise ValueError("checkpoint state expected null")
        return None
    if annotation is LocalSequenceRecord:
        return _decode_record(payload)
    if annotation is Path:
        return _resolved_path(payload, outdir)

    origin = get_origin(annotation)
    args = get_args(annotation)

    if origin is Literal:
        if payload not in args:
            raise ValueError(f"checkpoint state literal value is invalid: {payload!r}")
        return payload

    if origin in (Union, types.UnionType):
        if payload is None and type(None) in args:
            return None
        errors: list[Exception] = []
        for candidate in args:
            if candidate is type(None):
                continue
            try:
                return _decode(payload, candidate, outdir)
            except (ValueError, TypeError) as error:
                errors.append(error)
        raise ValueError("checkpoint state does not match any allowed union type") from (errors[-1] if errors else None)

    if origin is tuple:
        if not isinstance(payload, list):
            raise ValueError("checkpoint state expected a list for tuple field")
        if len(args) == 2 and args[1] is Ellipsis:
            return tuple(_decode(item, args[0], outdir) for item in payload)
        if len(payload) != len(args):
            raise ValueError("checkpoint state tuple length is invalid")
        return tuple(
            _decode(item, item_type, outdir)
            for item, item_type in zip(payload, args, strict=True)
        )

    if origin is list:
        if not isinstance(payload, list):
            raise ValueError("checkpoint state expected a list")
        item_type = args[0] if args else object
        return [_decode(item, item_type, outdir) for item in payload]

    if origin in (dict, Mapping, MappingABC):
        if not isinstance(payload, dict):
            raise ValueError("checkpoint state expected a mapping")
        key_type, value_type = args if len(args) == 2 else (str, object)
        if key_type is not str:
            raise ValueError("checkpoint state only supports string mapping keys")
        return {
            key: _decode(item, value_type, outdir)
            for key, item in payload.items()
            if isinstance(key, str)
        }

    if isinstance(annotation, type) and issubclass(annotation, Enum):
        try:
            return annotation(payload)
        except (ValueError, TypeError) as error:
            raise ValueError(f"checkpoint enum value for {annotation.__name__} is invalid") from error

    if isinstance(annotation, type) and is_dataclass(annotation):
        if not isinstance(payload, dict):
            raise ValueError(f"checkpoint state for {annotation.__name__} must be a mapping")
        expected = {field.name for field in fields(annotation)}
        if set(payload) != expected:
            raise ValueError(
                f"checkpoint state fields for {annotation.__name__} are invalid"
            )
        hints = get_type_hints(annotation)
        values = {
            field.name: _decode(payload[field.name], hints[field.name], outdir)
            for field in fields(annotation)
        }
        return annotation(**values)

    if annotation is bool:
        if type(payload) is not bool:
            raise ValueError("checkpoint state expected bool")
        return payload
    if annotation is int:
        if type(payload) is not int:
            raise ValueError("checkpoint state expected int")
        return payload
    if annotation is float:
        if type(payload) not in (int, float) or isinstance(payload, bool):
            raise ValueError("checkpoint state expected float")
        value = float(payload)
        if not math.isfinite(value):
            raise ValueError("checkpoint state float must be finite")
        return value
    if annotation is str:
        if not isinstance(payload, str):
            raise ValueError("checkpoint state expected string")
        return payload

    raise ValueError(f"unsupported checkpoint state annotation: {annotation!r}")


class _StructuralCodec(Generic[T]):
    def __init__(self, root_type: object):
        self.root_type = root_type

    def encode(self, value: T, outdir: Path) -> object:
        return _encode(value, Path(outdir))

    def decode(self, payload: object, outdir: Path) -> T:
        return _decode(payload, self.root_type, Path(outdir))  # type: ignore[return-value]


INPUT_CODEC = _StructuralCodec[tuple[LocalSequenceRecord, ...]](
    tuple[LocalSequenceRecord, ...]
)
PANEL_CODEC = _StructuralCodec[PanelResult](PanelResult)
QC_CODEC = _StructuralCodec[QCResult](QCResult)
CLUSTERING_CODEC = _StructuralCodec[ClusteringResult](ClusteringResult)
ALIGNMENT_CODEC = _StructuralCodec[AlignmentResult](AlignmentResult)
CONSERVATION_CODEC = _StructuralCodec[ConservationResult](ConservationResult)
CONTRASTIVE_CONSERVATION_CODEC = _StructuralCodec[ContrastiveConservationResult](
    ContrastiveConservationResult
)
PRIMER_DESIGN_CODEC = _StructuralCodec[PrimerDesignResult](PrimerDesignResult)
INCLUSIVITY_CODEC = _StructuralCodec[InclusivityResult](InclusivityResult)
SPECIFICITY_CODEC = _StructuralCodec[SpecificityResult](SpecificityResult)
RANKING_CODEC = _StructuralCodec[RankingResult](RankingResult)
