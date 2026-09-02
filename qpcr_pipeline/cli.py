import argparse
from pathlib import Path

from qpcr_pipeline.config import load_config
from qpcr_pipeline.diagnostics import EnvironmentInspector, doctor_exit_code, render_environment_report
from qpcr_pipeline.dry_run import dry_run_pipeline
from qpcr_pipeline.execution import ExecutionPolicy, STAGE_ORDER
from qpcr_pipeline.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="qpcr-pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("doctor", help="Inspect the Geison execution environment")

    run_parser = subparsers.add_parser(
        "run", help="Run or validate the qPCR pipeline configuration"
    )
    run_parser.add_argument("config", type=Path)
    run_parser.add_argument("--outdir", type=Path)
    run_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate environment and preview stage actions without executing or writing",
    )
    run_parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse valid checkpoints and recompute only invalid stages",
    )
    run_parser.add_argument(
        "--from-step",
        choices=STAGE_ORDER,
        help="Restart strictly from this stage using valid prerequisite checkpoints",
    )
    run_parser.add_argument(
        "--force-step",
        choices=STAGE_ORDER,
        help="With --resume, force this stage and its dependent subgraph",
    )

    return parser


def _render_dry_run(report) -> str:
    lines = [f"Dry run for target: {report.target_name}", "stage\taction\treason"]
    lines.extend(
        f"{decision.stage}\t{decision.action}\t{decision.reason}"
        for decision in report.decisions
    )
    if report.environment.missing_required_tools:
        lines.append(
            "Missing required tools: "
            + ", ".join(report.environment.missing_required_tools)
        )
    return "\n".join(lines)


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "doctor":
        report = EnvironmentInspector().inspect()
        print(render_environment_report(report))
        return doctor_exit_code(report)

    if args.command == "run":
        config = load_config(args.config)
        resume_control_requested = (
            args.resume or args.from_step is not None or args.force_step is not None
        )

        if args.outdir is None and resume_control_requested:
            parser.error("--resume, --from-step, and --force-step require --outdir")

        try:
            execution = ExecutionPolicy(
                resume=args.resume,
                from_step=args.from_step,
                force_step=args.force_step,
            )
        except ValueError as error:
            parser.error(str(error))

        if args.dry_run:
            report = dry_run_pipeline(config, args.outdir, execution=execution)
            print(_render_dry_run(report))
            return 2 if report.environment.missing_required_tools else 0

        if args.outdir is None:
            print(f"Loaded configuration for target: {config.target_name}")
            return 0

        summary = run_pipeline(config, args.outdir, execution=execution)
        print(
            f"{summary.status}: {summary.target_name} "
            f"({summary.sequence_count} sequences)"
        )
        return 0

    raise RuntimeError(f"Unsupported command: {args.command}")
