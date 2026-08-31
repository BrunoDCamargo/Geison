"""Trustworthy stage checkpoints with deterministic fingerprints and integrity checks."""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Generic, Literal, Mapping, Protocol, TypeVar


CHECKPOINT_SCHEMA_VERSION = 1
CheckpointStatus = Literal["RUNNING", "COMPLETE", "FAILED"]
T = TypeVar("T")


class CheckpointInvalidity(str, Enum):
    MISSING_MANIFEST = "MISSING_MANIFEST"
    INVALID_MANIFEST = "INVALID_MANIFEST"
    NON_COMPLETE_STATUS = "NON_COMPLETE_STATUS"
    FINGERPRINT_MISMATCH = "FINGERPRINT_MISMATCH"
    MISSING_STATE = "MISSING_STATE"
    STATE_HASH_MISMATCH = "STATE_HASH_MISMATCH"
    INVALID_STATE = "INVALID_STATE"
    MISSING_OUTPUT = "MISSING_OUTPUT"
    OUTPUT_HASH_MISMATCH = "OUTPUT_HASH_MISMATCH"
    RESULT_FINGERPRINT_MISMATCH = "RESULT_FINGERPRINT_MISMATCH"


class StateCodec(Protocol[T]):
    def encode(self, value: T, outdir: Path) -> object: ...

    def decode(self, payload: object, outdir: Path) -> T: ...


@dataclass(frozen=True, slots=True)
class CheckpointRequest:
    stage: str
    dependencies: Mapping[str, str]
    inputs: Mapping[str, object]
    parameters: Mapping[str, object]
    software: Mapping[str, object]
    tools: Mapping[str, object]

    def fingerprint_payload(self) -> dict[str, object]:
        return {
            "schema_version": CHECKPOINT_SCHEMA_VERSION,
            "stage": self.stage,
            "dependencies": dict(self.dependencies),
            "inputs": dict(self.inputs),
            "parameters": dict(self.parameters),
            "software": dict(self.software),
            "tools": dict(self.tools),
        }


@dataclass(frozen=True, slots=True)
class OutputIdentity:
    path: str
    sha256: str


@dataclass(frozen=True, slots=True)
class StateIdentity:
    path: str
    sha256: str


@dataclass(frozen=True, slots=True)
class CheckpointManifest:
    schema_version: int
    stage: str
    status: CheckpointStatus
    fingerprint: str
    result_fingerprint: str | None
    dependencies: Mapping[str, str]
    inputs: Mapping[str, object]
    parameters: Mapping[str, object]
    software: Mapping[str, object]
    tools: Mapping[str, object]
    outputs: tuple[OutputIdentity, ...]
    state: StateIdentity | None
    error: Mapping[str, str] | None = None


@dataclass(frozen=True, slots=True)
class CheckpointLoad(Generic[T]):
    state: T
    manifest: CheckpointManifest


@dataclass(frozen=True, slots=True)
class CheckpointValidation(Generic[T]):
    valid: bool
    invalidity: CheckpointInvalidity | None
    detail: str | None
    loaded: CheckpointLoad[T] | None


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def canonical_sha256(value: object) -> str:
    return "sha256:" + hashlib.sha256(_canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(Path(path).read_bytes()).hexdigest()


def causal_fingerprint(request: CheckpointRequest) -> str:
    return canonical_sha256(request.fingerprint_payload())


def result_fingerprint(
    causal: str,
    state_sha256: str,
    outputs: tuple[OutputIdentity, ...],
) -> str:
    ordered_outputs = sorted(outputs, key=lambda item: item.path)
    return canonical_sha256(
        {
            "fingerprint": causal,
            "state_sha256": state_sha256,
            "outputs": [
                {"path": output.path, "sha256": output.sha256}
                for output in ordered_outputs
            ],
        }
    )


def _atomic_write_text(destination: Path, text: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        os.close(descriptor)
        temporary.write_text(text, encoding="utf-8")
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _json_text(value: object) -> str:
    return json.dumps(
        value,
        indent=2,
        sort_keys=True,
        ensure_ascii=False,
        allow_nan=False,
    ) + "\n"


def _manifest_payload(manifest: CheckpointManifest) -> dict[str, object]:
    return {
        "schema_version": manifest.schema_version,
        "stage": manifest.stage,
        "status": manifest.status,
        "fingerprint": manifest.fingerprint,
        "result_fingerprint": manifest.result_fingerprint,
        "dependencies": dict(manifest.dependencies),
        "inputs": dict(manifest.inputs),
        "parameters": dict(manifest.parameters),
        "software": dict(manifest.software),
        "tools": dict(manifest.tools),
        "outputs": [
            {"path": output.path, "sha256": output.sha256}
            for output in manifest.outputs
        ],
        "state": (
            {"path": manifest.state.path, "sha256": manifest.state.sha256}
            if manifest.state is not None
            else None
        ),
        "error": dict(manifest.error) if manifest.error is not None else None,
    }


def _parse_manifest(value: object) -> CheckpointManifest:
    if not isinstance(value, dict):
        raise ValueError("checkpoint manifest must be a mapping")
    required = {
        "schema_version",
        "stage",
        "status",
        "fingerprint",
        "result_fingerprint",
        "dependencies",
        "inputs",
        "parameters",
        "software",
        "tools",
        "outputs",
        "state",
        "error",
    }
    if set(value) != required:
        raise ValueError("checkpoint manifest fields are invalid")

    schema_version = value["schema_version"]
    stage = value["stage"]
    status = value["status"]
    fingerprint = value["fingerprint"]
    result = value["result_fingerprint"]
    mapping_keys = ("dependencies", "inputs", "parameters", "software", "tools")

    if schema_version != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError("checkpoint schema version is unsupported")
    if not isinstance(stage, str) or not stage:
        raise ValueError("checkpoint stage is invalid")
    if status not in {"RUNNING", "COMPLETE", "FAILED"}:
        raise ValueError("checkpoint status is invalid")
    if not isinstance(fingerprint, str) or not fingerprint.startswith("sha256:"):
        raise ValueError("checkpoint fingerprint is invalid")
    if result is not None and (
        not isinstance(result, str) or not result.startswith("sha256:")
    ):
        raise ValueError("checkpoint result fingerprint is invalid")
    if any(not isinstance(value[key], dict) for key in mapping_keys):
        raise ValueError("checkpoint request metadata is invalid")

    raw_outputs = value["outputs"]
    if not isinstance(raw_outputs, list):
        raise ValueError("checkpoint outputs are invalid")
    outputs: list[OutputIdentity] = []
    for item in raw_outputs:
        if (
            not isinstance(item, dict)
            or set(item) != {"path", "sha256"}
            or not isinstance(item["path"], str)
            or not item["path"]
            or not isinstance(item["sha256"], str)
            or not item["sha256"].startswith("sha256:")
        ):
            raise ValueError("checkpoint output identity is invalid")
        outputs.append(OutputIdentity(item["path"], item["sha256"]))

    raw_state = value["state"]
    state: StateIdentity | None
    if raw_state is None:
        state = None
    elif (
        isinstance(raw_state, dict)
        and set(raw_state) == {"path", "sha256"}
        and isinstance(raw_state["path"], str)
        and raw_state["path"]
        and isinstance(raw_state["sha256"], str)
        and raw_state["sha256"].startswith("sha256:")
    ):
        state = StateIdentity(raw_state["path"], raw_state["sha256"])
    else:
        raise ValueError("checkpoint state identity is invalid")

    raw_error = value["error"]
    if raw_error is not None and (
        not isinstance(raw_error, dict)
        or any(
            not isinstance(key, str) or not isinstance(item, str)
            for key, item in raw_error.items()
        )
    ):
        raise ValueError("checkpoint error metadata is invalid")

    return CheckpointManifest(
        schema_version=schema_version,
        stage=stage,
        status=status,
        fingerprint=fingerprint,
        result_fingerprint=result,
        dependencies=dict(value["dependencies"]),
        inputs=dict(value["inputs"]),
        parameters=dict(value["parameters"]),
        software=dict(value["software"]),
        tools=dict(value["tools"]),
        outputs=tuple(outputs),
        state=state,
        error=dict(raw_error) if raw_error is not None else None,
    )


class CheckpointManager:
    def __init__(self, outdir: str | Path):
        self.outdir = Path(outdir)

    def _stage_dir(self, stage: str) -> Path:
        if (
            not isinstance(stage, str)
            or not stage
            or "/" in stage
            or "\\" in stage
            or stage in {".", ".."}
        ):
            raise ValueError("checkpoint stage must be a simple non-empty name")
        return self.outdir / ".checkpoints" / stage

    def _manifest_path(self, stage: str) -> Path:
        return self._stage_dir(stage) / "manifest.json"

    def _state_path(self, stage: str) -> Path:
        return self._stage_dir(stage) / "state.json"

    def _base_manifest(
        self, request: CheckpointRequest, status: CheckpointStatus
    ) -> CheckpointManifest:
        return CheckpointManifest(
            schema_version=CHECKPOINT_SCHEMA_VERSION,
            stage=request.stage,
            status=status,
            fingerprint=causal_fingerprint(request),
            result_fingerprint=None,
            dependencies=dict(request.dependencies),
            inputs=dict(request.inputs),
            parameters=dict(request.parameters),
            software=dict(request.software),
            tools=dict(request.tools),
            outputs=(),
            state=None,
            error=None,
        )

    def _write_manifest(self, manifest: CheckpointManifest) -> None:
        _atomic_write_text(
            self._manifest_path(manifest.stage),
            _json_text(_manifest_payload(manifest)),
        )

    def begin(self, request: CheckpointRequest) -> CheckpointManifest:
        manifest = self._base_manifest(request, "RUNNING")
        self._write_manifest(manifest)
        return manifest

    def _relative_output(self, output: Path) -> str:
        resolved_outdir = self.outdir.resolve()
        resolved_output = output.resolve()
        try:
            relative = resolved_output.relative_to(resolved_outdir)
        except ValueError as error:
            raise ValueError(
                "Checkpoint outputs must be inside the pipeline output directory."
            ) from error
        if relative.parts and relative.parts[0] == ".checkpoints":
            raise ValueError(
                "Scientific checkpoint outputs must not point inside .checkpoints."
            )
        return relative.as_posix()

    def complete(
        self,
        request: CheckpointRequest,
        state: T,
        codec: StateCodec[T],
        outputs: tuple[Path, ...],
    ) -> CheckpointManifest:
        stage_dir = self._stage_dir(request.stage)
        stage_dir.mkdir(parents=True, exist_ok=True)
        state_payload = codec.encode(state, self.outdir)
        state_path = self._state_path(request.stage)
        _atomic_write_text(state_path, _json_text(state_payload))
        state_identity = StateIdentity(
            path=state_path.relative_to(self.outdir).as_posix(),
            sha256=file_sha256(state_path),
        )

        identities: list[OutputIdentity] = []
        for raw_output in outputs:
            output = Path(raw_output)
            relative = self._relative_output(output)
            if not output.is_file():
                raise ValueError(
                    f"Checkpoint output does not exist as a file: {relative}"
                )
            identities.append(OutputIdentity(relative, file_sha256(output)))
        ordered = tuple(sorted(identities, key=lambda item: item.path))
        causal = causal_fingerprint(request)
        manifest = CheckpointManifest(
            schema_version=CHECKPOINT_SCHEMA_VERSION,
            stage=request.stage,
            status="COMPLETE",
            fingerprint=causal,
            result_fingerprint=result_fingerprint(
                causal, state_identity.sha256, ordered
            ),
            dependencies=dict(request.dependencies),
            inputs=dict(request.inputs),
            parameters=dict(request.parameters),
            software=dict(request.software),
            tools=dict(request.tools),
            outputs=ordered,
            state=state_identity,
            error=None,
        )
        self._write_manifest(manifest)
        return manifest

    def fail(
        self, request: CheckpointRequest, error: BaseException
    ) -> CheckpointManifest:
        failed = replace(
            self._base_manifest(request, "FAILED"),
            error={"type": type(error).__name__, "message": str(error)},
        )
        try:
            self._write_manifest(failed)
        except Exception:
            pass
        return failed

    def validate(
        self, request: CheckpointRequest, codec: StateCodec[T]
    ) -> CheckpointValidation[T]:
        manifest_path = self._manifest_path(request.stage)
        if not manifest_path.is_file():
            return CheckpointValidation(
                False,
                CheckpointInvalidity.MISSING_MANIFEST,
                "checkpoint manifest is missing",
                None,
            )
        try:
            manifest = _parse_manifest(
                json.loads(manifest_path.read_text(encoding="utf-8"))
            )
        except (
            OSError,
            UnicodeError,
            json.JSONDecodeError,
            ValueError,
            TypeError,
        ) as error:
            return CheckpointValidation(
                False, CheckpointInvalidity.INVALID_MANIFEST, str(error), None
            )
        if manifest.stage != request.stage:
            return CheckpointValidation(
                False,
                CheckpointInvalidity.INVALID_MANIFEST,
                "checkpoint stage does not match request",
                None,
            )
        if manifest.status != "COMPLETE":
            return CheckpointValidation(
                False,
                CheckpointInvalidity.NON_COMPLETE_STATUS,
                f"checkpoint status is {manifest.status}",
                None,
            )

        expected_causal = causal_fingerprint(request)
        if manifest.fingerprint != expected_causal:
            return CheckpointValidation(
                False,
                CheckpointInvalidity.FINGERPRINT_MISMATCH,
                "causal fingerprint changed",
                None,
            )
        if manifest.state is None:
            return CheckpointValidation(
                False,
                CheckpointInvalidity.MISSING_STATE,
                "checkpoint state identity is missing",
                None,
            )

        try:
            state_path = self._validated_internal_path(
                manifest.state.path, allow_checkpoint=True
            )
        except ValueError as error:
            return CheckpointValidation(
                False, CheckpointInvalidity.INVALID_MANIFEST, str(error), None
            )
        if not state_path.is_file():
            return CheckpointValidation(
                False,
                CheckpointInvalidity.MISSING_STATE,
                "checkpoint state file is missing",
                None,
            )
        try:
            actual_state_hash = file_sha256(state_path)
        except OSError as error:
            return CheckpointValidation(
                False, CheckpointInvalidity.MISSING_STATE, str(error), None
            )
        if actual_state_hash != manifest.state.sha256:
            return CheckpointValidation(
                False,
                CheckpointInvalidity.STATE_HASH_MISMATCH,
                "checkpoint state hash changed",
                None,
            )

        for output in manifest.outputs:
            try:
                path = self._validated_internal_path(
                    output.path, allow_checkpoint=False
                )
            except ValueError as error:
                return CheckpointValidation(
                    False,
                    CheckpointInvalidity.INVALID_MANIFEST,
                    str(error),
                    None,
                )
            if not path.is_file():
                return CheckpointValidation(
                    False,
                    CheckpointInvalidity.MISSING_OUTPUT,
                    f"checkpoint output is missing: {output.path}",
                    None,
                )
            try:
                actual_hash = file_sha256(path)
            except OSError as error:
                return CheckpointValidation(
                    False, CheckpointInvalidity.MISSING_OUTPUT, str(error), None
                )
            if actual_hash != output.sha256:
                return CheckpointValidation(
                    False,
                    CheckpointInvalidity.OUTPUT_HASH_MISMATCH,
                    f"checkpoint output hash changed: {output.path}",
                    None,
                )

        expected_result = result_fingerprint(
            expected_causal, actual_state_hash, manifest.outputs
        )
        if manifest.result_fingerprint != expected_result:
            return CheckpointValidation(
                False,
                CheckpointInvalidity.RESULT_FINGERPRINT_MISMATCH,
                "result fingerprint changed",
                None,
            )

        try:
            payload = json.loads(state_path.read_text(encoding="utf-8"))
            state = codec.decode(payload, self.outdir)
        except (
            OSError,
            UnicodeError,
            json.JSONDecodeError,
            ValueError,
            TypeError,
            KeyError,
        ) as error:
            return CheckpointValidation(
                False, CheckpointInvalidity.INVALID_STATE, str(error), None
            )
        return CheckpointValidation(
            True, None, None, CheckpointLoad(state, manifest)
        )

    def _validated_internal_path(
        self, relative: str, *, allow_checkpoint: bool
    ) -> Path:
        candidate = Path(relative)
        if candidate.is_absolute():
            raise ValueError(
                "checkpoint paths must be relative to the output directory"
            )
        resolved_outdir = self.outdir.resolve()
        resolved = (self.outdir / candidate).resolve()
        try:
            rel = resolved.relative_to(resolved_outdir)
        except ValueError as error:
            raise ValueError("checkpoint path escapes the output directory") from error
        if (
            not allow_checkpoint
            and rel.parts
            and rel.parts[0] == ".checkpoints"
        ):
            raise ValueError("scientific output path points inside .checkpoints")
        return resolved
