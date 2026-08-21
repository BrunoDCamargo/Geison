import argparse
from pathlib import Path

from qpcr_pipeline.config import load_config
from qpcr_pipeline.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="qpcr-pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run the qPCR pipeline")
    run_parser.add_argument("config", type=Path)
    run_parser.add_argument("--outdir", type=Path, required=True)

    return parser


def main() -> int:
    args = build_parser().parse_args()

    if args.command == "run":
        config = load_config(args.config)
        summary = run_pipeline(config, args.outdir)
        print(f"{summary.status}: {summary.target_name} ({summary.sequence_count} sequences)")
        return 0

    raise RuntimeError(f"Unsupported command: {args.command}")
