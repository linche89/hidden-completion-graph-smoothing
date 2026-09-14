"""Independent numerical cross-check for flagged exact-scenario SDP solves.

This audit compares directly re-evaluated radii from two conic solvers.  It is
numerical corroboration, not a rigorous dual lower bound or a uniqueness test
for the returned local-rule matrix.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter
from typing import Any, Mapping

import numpy as np
import scs

from experiments.smoother.provenance import experiment_provenance
from experiments.smoother.robust import scenario_certificate
from hidden_completion.local_rules import LocalRuleSpec
from hidden_completion.scenarios import finite_scenarios
from hidden_completion.solver import solve_exact_spectral


Array = np.ndarray


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _deployed_matrix(matrix: Array, tolerance: float) -> Array:
    value = np.asarray(matrix, dtype=np.float64).copy()
    scale = max(1.0, float(np.max(np.abs(value))))
    value[np.abs(value) <= tolerance * scale] = 0.0
    return value


def _local_rule_residuals(rule: LocalRuleSpec, matrix: Array) -> dict[str, float]:
    rows, columns = matrix.shape
    mask = rule.allowed_mask(rows, columns)
    residuals: dict[str, float] = {
        "forced_zero_max_abs": (
            float(np.max(np.abs(matrix[~mask]))) if np.any(~mask) else 0.0
        )
    }

    shared = 0.0
    for group in rule.data.get("shared_coefficients", []):
        anchor = tuple(group[0])
        for entry in group[1:]:
            shared = max(
                shared,
                abs(float(matrix[tuple(entry)] - matrix[anchor])),
            )
    residuals["shared_coefficient_max_abs"] = shared

    fixed = 0.0
    for item in rule.data.get("fixed_entries", []):
        if isinstance(item, dict):
            row = int(item["row"])
            column = int(item["column"])
            target = float(item["value"])
        else:
            row, column, target = int(item[0]), int(item[1]), float(item[2])
        fixed = max(fixed, abs(float(matrix[row, column]) - target))
    residuals["fixed_entry_max_abs"] = fixed

    linear = 0.0
    for equality in rule.data.get("linear_equalities", []):
        lhs = sum(
            float(coefficient) * float(matrix[int(row), int(column)])
            for row, column, coefficient in equality["terms"]
        )
        linear = max(linear, abs(lhs - float(equality["rhs"])))
    residuals["linear_equality_max_abs"] = linear

    row_sum = 0.0
    if "row_sums" in rule.data:
        targets = rule.data["row_sums"]
        if isinstance(targets, (int, float)):
            targets = [float(targets)] * rows
        row_sum = float(
            np.max(np.abs(np.sum(matrix, axis=1) - np.asarray(targets)))
        )
    residuals["row_sum_max_abs"] = row_sum
    residuals["maximum"] = max(residuals.values(), default=0.0)
    return residuals


def run_crosscheck(
    source: Mapping[str, Any],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    selection_status = str(config.get("selection_status", "optimal_inaccurate"))
    solver_config = dict(config["solver"])
    solver = str(solver_config["name"])
    solver_options = dict(solver_config.get("options", {}))
    deployment_tolerance = float(source["deployment_zero_tolerance"])
    selected = [
        row
        for row in source["multi_instance"]["rows"]
        if row["exact_finite"]["status"] == selection_status
    ]
    records: list[dict[str, Any]] = []
    started = perf_counter()

    for index, row in enumerate(selected, start=1):
        q = int(row["q"])
        h = int(row["h"])
        c_matrix = np.asarray(row["C"], dtype=np.float64)
        rule = LocalRuleSpec.from_dict(row["local_rule"])
        original = row["exact_finite"]
        original_raw_radius = float(
            original.get("raw_verified_radius", original["radius"])
        )
        original_deployed_radius = float(original["radius"])
        cell_started = perf_counter()
        result = solve_exact_spectral(
            c_matrix,
            q=q,
            h=h,
            mode="finite",
            local_rule=rule,
            prune_vertices=True,
            solver=solver,
            solver_options=solver_options,
        )
        deployed_q = _deployed_matrix(result.q_matrix, deployment_tolerance)
        deployed_certificate = scenario_certificate(
            c_matrix,
            deployed_q,
            finite_scenarios(q, h),
        )
        deployed_radius = float(deployed_certificate["certified_radius"])
        raw_relative_delta = (
            result.radius - original_raw_radius
        ) / max(original_raw_radius, np.finfo(float).tiny)
        deployed_relative_delta = (
            deployed_radius - original_deployed_radius
        ) / max(original_deployed_radius, np.finfo(float).tiny)
        original_raw_q = np.asarray(original["raw_Q"], dtype=np.float64)
        records.append(
            {
                "problem_id": row["problem_id"],
                "q": q,
                "h": h,
                "family": row["family"],
                "rule_name": row["rule_name"],
                "primary_status": original["status"],
                "secondary_status": result.status,
                "primary_raw_radius": original_raw_radius,
                "secondary_raw_radius": result.radius,
                "raw_radius_relative_delta": raw_relative_delta,
                "primary_deployed_radius": original_deployed_radius,
                "secondary_deployed_radius": deployed_radius,
                "deployed_radius_relative_delta": deployed_relative_delta,
                "secondary_epigraph_minus_recomputed_squared_risk": (
                    result.optimization_squared_radius - result.squared_radius
                ),
                "secondary_certificate_min_eigenvalue": (
                    result.certificate_min_eigenvalue
                ),
                "local_rule_residuals": _local_rule_residuals(
                    rule, result.q_matrix
                ),
                "q_frobenius_delta": float(
                    np.linalg.norm(result.q_matrix - original_raw_q)
                ),
                "solver_seconds": result.timing_seconds["solver"],
                "wall_seconds": perf_counter() - cell_started,
            }
        )
        print(
            f"[{index:02d}/{len(selected)}] {row['problem_id']} h={h}: "
            f"{result.status}; relative radius delta={raw_relative_delta:.3e}",
            flush=True,
        )

    maximum_raw_delta = max(
        (abs(float(row["raw_radius_relative_delta"])) for row in records),
        default=0.0,
    )
    maximum_deployed_delta = max(
        (abs(float(row["deployed_radius_relative_delta"])) for row in records),
        default=0.0,
    )
    maximum_epigraph_gap = max(
        (
            abs(float(row["secondary_epigraph_minus_recomputed_squared_risk"]))
            for row in records
        ),
        default=0.0,
    )
    maximum_local_residual = max(
        (float(row["local_rule_residuals"]["maximum"]) for row in records),
        default=0.0,
    )
    thresholds = dict(config.get("checks", {}))
    checks = {
        "selected_cell_count_matches_configuration": len(selected)
        == int(config["expected_selected_cells"]),
        "all_secondary_statuses_optimal": all(
            row["secondary_status"] == "optimal" for row in records
        ),
        "relative_radius_agreement": maximum_raw_delta
        <= float(thresholds["maximum_relative_radius_delta"]),
        "secondary_epigraph_agrees_with_direct_reevaluation": maximum_epigraph_gap
        <= float(thresholds["maximum_epigraph_squared_risk_gap"]),
        "secondary_local_rule_residuals_below_tolerance": maximum_local_residual
        <= float(thresholds["maximum_local_rule_residual"]),
    }
    if not all(checks.values()):
        raise RuntimeError(f"cross-solver checks failed: {checks}")
    return {
        "schema_version": 1,
        "classification": "independent_numerical_cross_solver_audit",
        "scope": (
            "all primary exact_finite cells with status "
            f"{selection_status!r}"
        ),
        "interpretation": (
            "numerical cross-solver corroboration only; no rigorous dual "
            "lower bound and no claim that Q is unique"
        ),
        "primary_solver": source["multi_instance"]["rows"][0]["exact_finite"][
            "solver"
        ],
        "secondary_solver": solver,
        "secondary_solver_options": solver_options,
        "deployment_zero_tolerance": deployment_tolerance,
        "summary": {
            "selected_cells": len(selected),
            "secondary_optimal_cells": sum(
                row["secondary_status"] == "optimal" for row in records
            ),
            "maximum_absolute_relative_raw_radius_delta": maximum_raw_delta,
            "maximum_absolute_relative_deployed_radius_delta": maximum_deployed_delta,
            "maximum_absolute_epigraph_squared_risk_gap": maximum_epigraph_gap,
            "maximum_local_rule_residual": maximum_local_residual,
            "minimum_secondary_certificate_eigenvalue": min(
                (
                    float(row["secondary_certificate_min_eigenvalue"])
                    for row in records
                ),
                default=0.0,
            ),
            "maximum_q_frobenius_delta": max(
                (float(row["q_frobenius_delta"]) for row in records),
                default=0.0,
            ),
            "wall_seconds": perf_counter() - started,
        },
        "checks": checks,
        "rows": records,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Cross-check flagged exact-scenario SDP solves."
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    repository_root = Path(__file__).resolve().parents[2]
    config_path = arguments.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    source_path = (repository_root / str(config["source"])).resolve()
    source = json.loads(source_path.read_text(encoding="utf-8"))
    payload = run_crosscheck(source, config)
    provenance = experiment_provenance(repository_root)
    provenance["packages"]["scs"] = scs.__version__
    payload["inputs"] = {
        "config_path": str(config_path),
        "config_sha256": _sha256(config_path),
        "source_path": str(source_path),
        "source_sha256": _sha256(source_path),
    }
    payload["provenance"] = provenance
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(arguments.output)


if __name__ == "__main__":
    main()
