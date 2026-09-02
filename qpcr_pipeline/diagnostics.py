from __future__ import annotations

import platform
import subprocess
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from typing import Mapping, Protocol

from qpcr_pipeline.config import PipelineConfig


@dataclass(frozen=True, slots=True)
class CommandResult:
    returncode: int
    stdout: str
    stderr: str


class CommandRunner(Protocol):
    def run(self, argv: tuple[str, ...]) -> CommandResult: ...


class SubprocessCommandRunner:
    def run(self, argv: tuple[str, ...]) -> CommandResult:
        try:
            completed = subprocess.run(
                argv,
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as error:
            return CommandResult(127, "", str(error))
        return CommandResult(completed.returncode, completed.stdout, completed.stderr)


@dataclass(frozen=True, slots=True)
class ComponentReport:
    name: str
    status: str
    required: bool
    installed: bool
    version: str | None
    commit: str | None = None
    dirty: bool | None = None


@dataclass(frozen=True, slots=True)
class EnvironmentReport:
    python: ComponentReport
    geison: ComponentReport
    git: ComponentReport
    tools: Mapping[str, ComponentReport]

    @property
    def missing_required_tools(self) -> tuple[str, ...]:
        return tuple(
            name
            for name, report in self.tools.items()
            if report.required and not report.installed
        )


TOOL_PROBES: dict[str, tuple[str, ...]] = {
    "cd-hit-est": ("cd-hit-est", "-h"),
    "mafft": ("mafft", "--version"),
    "primer3_core": ("primer3_core", "--about"),
    "blast+": ("blastn", "-version"),
}


def _bounded_version(result: CommandResult) -> str | None:
    text = " ".join(
        line.strip()
        for chunk in (result.stdout, result.stderr)
        for line in chunk.splitlines()
        if line.strip()
    )
    return text[:500] or None


class EnvironmentInspector:
    def __init__(self, *, runner: CommandRunner | None = None) -> None:
        self.runner = runner or SubprocessCommandRunner()

    def inspect(self, config: PipelineConfig | None = None) -> EnvironmentReport:
        python_report = ComponentReport(
            name="Python",
            status="USED",
            required=True,
            installed=True,
            version=platform.python_version(),
        )

        try:
            geison_version = version("geison-qpcr")
        except PackageNotFoundError:
            geison_report = ComponentReport(
                name="Geison",
                status="UNAVAILABLE",
                required=True,
                installed=False,
                version=None,
            )
        else:
            geison_report = ComponentReport(
                name="Geison",
                status="USED",
                required=True,
                installed=True,
                version=geison_version,
            )

        git_commit = self.runner.run(("git", "rev-parse", "HEAD"))
        if git_commit.returncode == 0:
            git_status = self.runner.run(("git", "status", "--porcelain"))
            git_report = ComponentReport(
                name="Git",
                status="AVAILABLE",
                required=False,
                installed=True,
                version=None,
                commit=git_commit.stdout.strip() or None,
                dirty=bool(git_status.stdout.strip()) if git_status.returncode == 0 else None,
            )
        else:
            git_report = ComponentReport(
                name="Git",
                status="UNAVAILABLE",
                required=False,
                installed=False,
                version=None,
            )

        required = {
            "cd-hit-est": bool(config and config.clustering.enabled),
            "mafft": bool(config and config.alignment.enabled),
            "primer3_core": bool(config and config.primer_design.enabled),
            "blast+": False,
        }
        tools: dict[str, ComponentReport] = {}
        for name, probe in TOOL_PROBES.items():
            result = self.runner.run(probe)
            installed = result.returncode == 0
            if name == "blast+":
                status = "NOT_USED"
            elif installed:
                status = "AVAILABLE"
            else:
                status = "UNAVAILABLE"
            tools[name] = ComponentReport(
                name=name,
                status=status,
                required=required[name],
                installed=installed,
                version=_bounded_version(result) if installed else None,
            )

        return EnvironmentReport(
            python=python_report,
            geison=geison_report,
            git=git_report,
            tools=tools,
        )
