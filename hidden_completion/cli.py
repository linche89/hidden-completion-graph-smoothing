"""Command-line interface for the exact hidden-completion solver."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .scenarios import finite_scenarios, partial_partition_scenarios
from .solver import solve_exact_spectral


def _write_json(payload: object, destination: str | None) -> None:
    rendered = json.dumps(payload, indent=2, sort_keys=False) + "\n"
    if destination:
        Path(destination).write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


def _solve(config_path: str, output: str | None, verbose: bool) -> None:
    config: dict[str, Any] = json.loads(Path(config_path).read_text(encoding="utf-8"))
    completion = config.get("completion", {})
    mode = completion.get("mode", config.get("mode", "finite"))
    h = completion.get("hidden_budget", config.get("h"))
    solver_config = config.get("solver", {})
    result = solve_exact_spectral(
        np.asarray(config["C"], dtype=float),
        q=config.get("q"),
        h=h,
        mode=mode,
        local_rule=config.get("local_rule"),
        prune_vertices=config.get("prune_vertices", True),
        solver=solver_config.get("name"),
        solver_options=solver_config.get("options"),
        active_tolerance=config.get("active_tolerance", 2e-6),
        verbose=verbose,
    )
    _write_json(result.to_dict(), output)


def _enumerate(q: int, h: int | None, vertices_only: bool, output: str | None) -> None:
    if h is None:
        scenarios = partial_partition_scenarios(q)
        mode = "unbounded"
    else:
        scenarios = finite_scenarios(q, h, vertices_only=vertices_only)
        mode = "finite"
    _write_json(
        {
            "mode": mode,
            "q": q,
            "h": h,
            "scenario_count": len(scenarios),
            "scenarios": [scenario.to_dict() for scenario in scenarios],
        },
        output,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m hidden_completion",
        description="Exact finite-h and unbounded hidden-completion scenario solver.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    solve_parser = subparsers.add_parser("solve", help="solve a JSON-configured spectral SDP")
    solve_parser.add_argument("config", help="path to the JSON problem description")
    solve_parser.add_argument("--output", help="write JSON output to this path")
    solve_parser.add_argument("--verbose", action="store_true", help="show CVXPY solver output")

    enumerate_parser = subparsers.add_parser("enumerate", help="enumerate exact scenarios")
    enumerate_parser.add_argument("--q", type=int, required=True, help="number of visible ports")
    enumerate_parser.add_argument(
        "--h",
        type=int,
        help="hidden budget; omit to enumerate unbounded partial partitions",
    )
    enumerate_parser.add_argument(
        "--vertices-only",
        action="store_true",
        help="retain only exact terminal-hull vertices in finite mode",
    )
    enumerate_parser.add_argument("--output", help="write JSON output to this path")
    return parser


def main() -> None:
    parser = build_parser()
    arguments = parser.parse_args()
    if arguments.command == "solve":
        _solve(arguments.config, arguments.output, arguments.verbose)
    else:
        _enumerate(arguments.q, arguments.h, arguments.vertices_only, arguments.output)
