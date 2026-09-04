"""Deterministic allowlisted packaging of published Geison run evidence."""

from __future__ import annotations

from pathlib import Path
import zipfile


_ROOT_FILES = (
    "report.html",
    "report_error.json",
    "run_manifest.json",
    "run_summary.json",
    "qc_report.json",
    "ncbi_dataset_manifest.json",
)

_SCIENTIFIC_DIRECTORIES = (
    "panel",
    "conservation",
    "contrastive_conservation",
    "primer_design",
    "inclusivity",
    "specificity",
    "ranking",
)

_FIXED_ZIP_TIMESTAMP = (1980, 1, 1, 0, 0, 0)


def _regular_file(path: Path) -> bool:
    return path.is_file() and not path.is_symlink()


def _archive_name(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def _published_files(output_dir: Path) -> tuple[tuple[str, Path], ...]:
    entries: dict[str, Path] = {}

    for relative in _ROOT_FILES:
        path = output_dir / relative
        if _regular_file(path):
            entries[relative] = path

    for directory_name in _SCIENTIFIC_DIRECTORIES:
        directory = output_dir / directory_name
        if not directory.is_dir() or directory.is_symlink():
            continue
        for path in directory.rglob("*"):
            if not _regular_file(path):
                continue
            archive_name = _archive_name(path, output_dir)
            entries[archive_name] = path

    checkpoint_root = output_dir / ".checkpoints"
    if checkpoint_root.is_dir() and not checkpoint_root.is_symlink():
        for path in checkpoint_root.rglob("manifest.json"):
            if not _regular_file(path):
                continue
            archive_name = _archive_name(path, output_dir)
            entries[archive_name] = path

    return tuple(sorted(entries.items(), key=lambda item: item[0]))


def _validated_extra_files(
    extra_files: tuple[Path, ...],
) -> tuple[tuple[str, Path], ...]:
    entries: dict[str, Path] = {}
    for raw_path in extra_files:
        path = Path(raw_path)
        if not _regular_file(path):
            raise ValueError(f"Evidence bundle extra file is not a regular file: {path}")
        archive_name = f"inputs/{path.name}"
        if archive_name in entries:
            raise ValueError(
                f"Evidence bundle extra files contain duplicate name: {path.name}"
            )
        entries[archive_name] = path
    return tuple(sorted(entries.items(), key=lambda item: item[0]))


def _write_entry(
    archive: zipfile.ZipFile,
    archive_name: str,
    path: Path,
) -> None:
    info = zipfile.ZipInfo(archive_name, date_time=_FIXED_ZIP_TIMESTAMP)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = (0o100644 & 0xFFFF) << 16
    archive.writestr(info, path.read_bytes())


def create_evidence_bundle(
    output_dir: Path,
    destination: Path,
    *,
    extra_files: tuple[Path, ...] = (),
) -> Path:
    """Package allowlisted run artifacts without recalculating scientific results."""
    root = Path(output_dir)
    if not root.is_dir():
        raise ValueError(f"Evidence bundle output directory does not exist: {root}")

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    entries: dict[str, Path] = dict(_published_files(root))
    for archive_name, path in _validated_extra_files(tuple(extra_files)):
        if archive_name in entries:
            raise ValueError(f"Evidence bundle archive name collision: {archive_name}")
        entries[archive_name] = path

    ordered = tuple(sorted(entries.items(), key=lambda item: item[0]))
    temporary = destination.with_name(f".{destination.name}.tmp")
    temporary.unlink(missing_ok=True)
    try:
        with zipfile.ZipFile(
            temporary,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=9,
        ) as archive:
            for archive_name, path in ordered:
                _write_entry(archive, archive_name, path)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)

    return destination
