"""Factorial physical-completion sweep under exact finite-h certificates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from hidden_completion.local_rules import LocalRuleSpec
from hidden_completion.solver import solve_exact_spectral

from .provenance import experiment_provenance
from .robust import full_error_operator, generate_completion


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


def _positive_integer(value: Any, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _weight_profiles(raw: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    profiles: list[dict[str, Any]] = []
    names: set[str] = set()
    for entry in raw:
        name = str(entry["name"])
        minimum = float(entry["min"])
        maximum = float(entry["max"])
        if not name or name in names:
            raise ValueError("weight profile names must be nonempty and unique")
        if not 0.0 < minimum <= maximum < np.inf:
            raise ValueError("weight profiles require 0 < min <= max < infinity")
        names.add(name)
        profiles.append({"name": name, "min": minimum, "max": maximum})
    if not profiles:
        raise ValueError("at least one weight profile is required")
    return profiles


def _sample_summary(samples: list[dict[str, Any]], radius: float) -> dict[str, Any]:
    errors = np.asarray([row["operator_2_norm_error"] for row in samples])
    cut_edges = np.asarray([row["port_cut_edges"] for row in samples])
    cut_weights = np.asarray([row["port_cut_weight"] for row in samples])
    return {
        "sample_count": len(samples),
        "maximum_operator_2_norm_error": float(np.max(errors)),
        "mean_operator_2_norm_error": float(np.mean(errors)),
        "maximum_error_over_certificate": float(np.max(errors) / radius)
        if radius > 0.0
        else 0.0,
        "minimum_certificate_slack": float(np.min(radius - errors)),
        "realized_port_cut_edges_min": int(np.min(cut_edges)),
        "realized_port_cut_edges_max": int(np.max(cut_edges)),
        "realized_port_cut_weight_min": float(np.min(cut_weights)),
        "realized_port_cut_weight_max": float(np.max(cut_weights)),
    }


def run_parameter_sweep(config: Mapping[str, Any]) -> dict[str, Any]:
    """Run the configured physical samples and audit the exact certificates."""

    q = _positive_integer(config["q"], "q")
    c_matrix = np.asarray(config["C"], dtype=np.float64)
    if c_matrix.ndim != 2 or c_matrix.shape[1] != q or not np.all(np.isfinite(c_matrix)):
        raise ValueError("C must be a finite matrix with q columns")

    hidden_budgets = sorted({int(value) for value in config["hidden_budgets"]})
    if not hidden_budgets or hidden_budgets[0] < 0:
        raise ValueError("hidden_budgets must be nonempty and nonnegative")
    profiles = _weight_profiles(config["weight_profiles"])
    mean_degrees = sorted({float(value) for value in config["mean_degrees"]})
    if not mean_degrees or not all(0.0 < value < np.inf for value in mean_degrees):
        raise ValueError("mean_degrees must be finite and positive")

    replicates = _positive_integer(config.get("replicates", 1), "replicates")
    seed = int(config.get("seed", 20260914))
    family = str(config.get("family", "random_cyclic"))
    base_graph_options = dict(config.get("graph_options", {}))
    local_rule = LocalRuleSpec.from_dict(config.get("local_rule"))
    solver_config = dict(config.get("solver", {}))
    solver = solver_config.get("name")
    solver_options = solver_config.get("options")
    tolerance = float(config.get("certificate_tolerance", 2e-6))
    if tolerance < 0.0 or not np.isfinite(tolerance):
        raise ValueError("certificate_tolerance must be finite and nonnegative")

    total_samples = (
        len(hidden_budgets) * len(profiles) * len(mean_degrees) * replicates
    )
    seeds = iter(np.random.SeedSequence(seed).spawn(total_samples))
    groups: list[dict[str, Any]] = []
    all_samples: list[dict[str, Any]] = []

    for h in hidden_budgets:
        exact = solve_exact_spectral(
            c_matrix,
            q=q,
            h=h,
            mode="finite",
            local_rule=local_rule,
            prune_vertices=True,
            solver=solver,
            solver_options=solver_options,
        )
        group_samples: list[dict[str, Any]] = []
        for profile in profiles:
            for mean_degree in mean_degrees:
                options = dict(base_graph_options)
                options["mean_degree"] = mean_degree
                for replicate in range(replicates):
                    child_seed = int(
                        next(seeds).generate_state(1, dtype=np.uint32)[0]
                    )
                    completion = generate_completion(
                        q,
                        h,
                        family,
                        child_seed,
                        weight_min=profile["min"],
                        weight_max=profile["max"],
                        graph_options=options,
                    )
                    error = full_error_operator(c_matrix, exact.q_matrix, completion)
                    operator_error = float(np.linalg.norm(error, ord=2))
                    metadata = completion.metadata
                    row = {
                        "h": h,
                        "hidden_count": completion.hidden_count,
                        "family": metadata["family"],
                        "seed": child_seed,
                        "replicate": replicate,
                        "weight_profile": profile["name"],
                        "requested_weight_min": profile["min"],
                        "requested_weight_max": profile["max"],
                        "realized_weight_min": metadata.get("weight_min_realized"),
                        "realized_weight_max": metadata.get("weight_max_realized"),
                        "requested_mean_degree": mean_degree,
                        "nodes": completion.nodes,
                        "edges": metadata["edges"],
                        "cycle_rank": metadata["cycle_rank"],
                        "port_cut_edges": metadata["port_cut_edges"],
                        "port_cut_weight": metadata["port_cut_weight"],
                        "operator_2_norm_error": operator_error,
                        "isotropic_expected_mse_per_output": float(
                            np.sum(np.square(error)) / error.shape[0]
                        ),
                        "exact_finite_h_radius": exact.radius,
                        "certificate_slack": exact.radius - operator_error,
                    }
                    group_samples.append(row)
                    all_samples.append(row)
        groups.append(
            {
                "h": h,
                "exact_design": exact.to_dict(),
                "summary": _sample_summary(group_samples, exact.radius),
                "samples": group_samples,
            }
        )

    positive_hidden_samples = [row for row in all_samples if row["h"] > 0]
    violations = [
        row["operator_2_norm_error"] - row["exact_finite_h_radius"]
        for row in all_samples
    ]
    return {
        "schema_version": 1,
        "classification": "theorem_matching_physical_completion_sweep",
        "semantics": {
            "certificate": "exact finite-h full-input spectral operator norm",
            "samples": (
                "finite-conductance physical completions used only as lower-bound "
                "stress points; they do not define the certificate"
            ),
            "port_cut": "realized visible-to-hidden edge count and total weight",
        },
        "problem": {
            "q": q,
            "C": c_matrix.tolist(),
            "local_rule": local_rule.data,
            "unit_grounding": True,
        },
        "sweep_axes": {
            "hidden_budgets": hidden_budgets,
            "weight_profiles": profiles,
            "mean_degrees": mean_degrees,
            "replicates": replicates,
            "family": family,
            "seed": seed,
            "total_physical_completions": len(all_samples),
        },
        "by_hidden_budget": groups,
        "checks": {
            "all_physical_errors_within_exact_certificate": max(violations)
            <= tolerance,
            "maximum_certificate_violation": max(0.0, float(max(violations))),
            "multiple_realized_port_cut_sizes": len(
                {row["port_cut_edges"] for row in positive_hidden_samples}
            )
            >= 2,
            "all_requested_cells_present": len(all_samples) == total_samples,
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep hidden budget, edge-weight scale and graph density under exact "
            "finite-h certificates."
        )
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    config_path = arguments.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    repository_root = Path(__file__).resolve().parents[2]
    payload = run_parameter_sweep(config)
    payload["config_path"] = str(config_path)
    payload["provenance"] = experiment_provenance(repository_root)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(arguments.output)


if __name__ == "__main__":
    main()
