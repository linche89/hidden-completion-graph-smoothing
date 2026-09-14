"""Command-line entry point for P3/P4 graph-smoother experiments."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .provenance import experiment_provenance
from .robust import run_robust_benchmark
from .stress import run_stress_benchmark


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m experiments.smoother",
        description="Run frozen-scope graph-smoother robust and stress experiments.",
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--section",
        choices=("all", "robust", "stress"),
        default="all",
        help="select one configured experiment section",
    )
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    config_path = arguments.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    repository_root = Path(__file__).resolve().parents[2]
    payload: dict[str, Any] = {
        "schema_version": 1,
        "experiment_program": "P3/P4 graph-regularized state smoothing",
        "config_path": str(config_path),
        "provenance": experiment_provenance(repository_root),
    }
    if arguments.section in {"all", "robust"}:
        if "robust" not in config:
            raise ValueError("configuration has no robust section")
        payload["robust"] = run_robust_benchmark(dict(config["robust"]))
    if arguments.section in {"all", "stress"}:
        if "stress" not in config:
            raise ValueError("configuration has no stress section")
        payload["stress"] = run_stress_benchmark(
            dict(config["stress"]), base_directory=config_path.parent
        )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(arguments.output)
