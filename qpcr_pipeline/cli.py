import argparse
from pathlib import Path

from qpcr_pipeline.config import load_config


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="qpcr-pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Load and validate a pipeline configuration")
    run_parser.add_argument("config", type=Path)

    return parser


def main() -> int:
    args = build_parser().parse_args()

    if args.command == "run":
        config = load_config(args.config)
        print(f"Loaded configuration for target: {config.target_name}")
        return 0

    raise RuntimeError(f"Unsupported command: {args.command}")
