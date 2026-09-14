"""Core theorem-matched evidence for the TSP manuscript.

This module deliberately excludes fixed-known-graph and hardware benchmarks.
Every returned baseline rule is re-evaluated on the complete finite scenario
set, and every path witness is a finite-conductance connected degree-two graph.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from statistics import NormalDist
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from hidden_completion.local_rules import LocalRuleSpec
from hidden_completion.metrics import spectral_squared_risk
from hidden_completion.scenarios import Scenario, finite_scenarios
from hidden_completion.solver import ExactSolveResult, solve_exact_spectral

from .provenance import experiment_provenance, record_path
from .resources import direct_rule_resources
from .robust import (
    CompletionSample,
    _solve_sampled_operators,
    full_error_operator,
    generate_completion_samples,
    scenario_certificate,
    solve_psd_contraction_relaxation,
)


Array = NDArray[np.float64]


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
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def path_adjacency(q: int) -> Array:
    if not isinstance(q, int) or isinstance(q, bool) or q < 1:
        raise ValueError("q must be a positive integer")
    adjacency = np.zeros((q, q), dtype=np.float64)
    for index in range(q - 1):
        adjacency[index, index + 1] = 1.0
        adjacency[index + 1, index] = 1.0
    return adjacency


def _r_hop_rule(q: int, radius: int, *, row_sums: float | None = None) -> dict[str, Any]:
    rule: dict[str, Any] = {
        "r_hop": {
            "adjacency": path_adjacency(q).tolist(),
            "output_nodes": list(range(q)),
            "radius": radius,
        }
    }
    if row_sums is not None:
        rule["row_sums"] = row_sums
    return rule


def frozen_benchmark_definitions(q: int) -> list[dict[str, Any]]:
    """Return the fixed problem families used in the breadth audit."""

    if q < 3:
        raise ValueError("the frozen benchmark uses q >= 3")
    definitions: list[dict[str, Any]] = []

    # Direct graph-state smoothing.  Q 1 = 1 is the natural exact-constant
    # constraint because (I+L_G)^-1 1 = 1 for every completion.
    for label, radius in (("zero_hop", 0), ("one_hop", 1), ("full_port", q - 1)):
        definitions.append(
            {
                "problem_id": f"q{q}_constant_preserving_state_{label}",
                "family": "constant_preserving_state",
                "rule_name": label,
                "rounds": radius,
                "C": np.eye(q).tolist(),
                "local_rule": _r_hop_rule(q, radius, row_sums=1.0),
            }
        )

    # Hardware/protocol-constrained diagonal output gains.
    diagonal_c = np.diag(np.linspace(0.5, 1.0, q))
    diagonal_rule = _r_hop_rule(q, 0)
    definitions.append(
        {
            "problem_id": f"q{q}_heterogeneous_diagonal_independent",
            "family": "heterogeneous_diagonal",
            "rule_name": "independent_opcode",
            "rounds": 0,
            "C": diagonal_c.tolist(),
            "local_rule": diagonal_rule,
        }
    )
    shared_rule = dict(diagonal_rule)
    shared_rule["shared_coefficients"] = [[[index, index] for index in range(q)]]
    definitions.append(
        {
            "problem_id": f"q{q}_heterogeneous_diagonal_shared",
            "family": "heterogeneous_diagonal",
            "rule_name": "shared_opcode",
            "rounds": 0,
            "C": diagonal_c.tolist(),
            "local_rule": shared_rule,
        }
    )

    # Negative control: without constant preservation, one-hop state design
    # collapses to the completion-independent center and need not communicate.
    definitions.append(
        {
            "problem_id": f"q{q}_unconstrained_state_control",
            "family": "unconstrained_state_control",
            "rule_name": "one_hop",
            "rounds": 1,
            "C": np.eye(q).tolist(),
            "local_rule": _r_hop_rule(q, 1),
        }
    )
    return definitions


def _deployed_matrix(matrix: Array, tolerance: float) -> Array:
    value = np.asarray(matrix, dtype=np.float64).copy()
    if tolerance < 0.0:
        raise ValueError("deployment tolerance must be nonnegative")
    value[np.abs(value) <= tolerance * max(1.0, float(np.max(np.abs(value))))] = 0.0
    return value


def _quantile_summary(values: Sequence[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    if array.size == 0:
        raise ValueError("cannot summarize an empty sequence")
    return {
        "minimum": float(np.min(array)),
        "q25": float(np.quantile(array, 0.25)),
        "median": float(np.median(array)),
        "q75": float(np.quantile(array, 0.75)),
        "maximum": float(np.max(array)),
    }


def wilson_interval(successes: int, trials: int, confidence: float = 0.95) -> tuple[float, float]:
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("Wilson inputs require 0 <= successes <= trials and trials > 0")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie strictly between zero and one")
    z = NormalDist().inv_cdf(0.5 + confidence / 2.0)
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    center = (proportion + z * z / (2.0 * trials)) / denominator
    half_width = (
        z
        * math.sqrt(proportion * (1.0 - proportion) / trials + z * z / (4.0 * trials**2))
        / denominator
    )
    return max(0.0, center - half_width), min(1.0, center + half_width)


def _deployed_exact_record(
    result: ExactSolveResult,
    c_matrix: Array,
    scenarios: Sequence[Scenario],
    *,
    tolerance: float,
) -> tuple[dict[str, Any], Array, float]:
    deployed = _deployed_matrix(result.q_matrix, tolerance)
    certificate = scenario_certificate(c_matrix, deployed, scenarios)
    outputs, q = deployed.shape
    identity_outputs = np.eye(outputs)
    identity_ports = np.eye(q)
    verified_squared_radius = float(certificate["certified_squared_radius"])
    deployed_minimum_eigenvalue = math.inf
    for scenario in scenarios:
        x_matrix = scenario.matrix
        top_left = (
            verified_squared_radius * identity_outputs
            - c_matrix @ x_matrix @ c_matrix.T
            + c_matrix @ x_matrix @ deployed.T
            + deployed @ x_matrix @ c_matrix.T
        )
        top_left = (top_left + top_left.T) / 2.0
        block = np.block([[top_left, deployed], [deployed.T, identity_ports]])
        deployed_minimum_eigenvalue = min(
            deployed_minimum_eigenvalue,
            float(np.linalg.eigvalsh((block + block.T) / 2.0)[0]),
        )
    record = result.to_dict()
    record["raw_Q"] = record.pop("Q")
    record["Q"] = deployed.tolist()
    record["raw_verified_radius"] = record["radius"]
    record["radius"] = certificate["certified_radius"]
    record["squared_radius"] = certificate["certified_squared_radius"]
    record["deployment_zero_tolerance"] = tolerance
    record["deployment_radius_change"] = (
        float(certificate["certified_radius"]) - result.radius
    )
    record["raw_certificate_min_eigenvalue"] = record["certificate_min_eigenvalue"]
    record["certificate_min_eigenvalue"] = deployed_minimum_eigenvalue
    record["active_worst_case_scenarios"] = certificate[
        "active_worst_case_scenarios"
    ]
    return record, deployed, float(certificate["certified_radius"])


def _sampling_record(
    c_matrix: Array,
    rule: LocalRuleSpec,
    scenarios: Sequence[Scenario],
    samples: Sequence[CompletionSample],
    *,
    exact_radius: float,
    solver: str | None,
    solver_options: Mapping[str, Any] | None,
    absolute_failure_tolerance: float,
    relative_failure_tolerance: float,
) -> dict[str, Any]:
    sampled = _solve_sampled_operators(
        c_matrix,
        samples,
        local_rule=rule,
        solver=solver,
        solver_options=dict(solver_options or {}),
    )
    q_value = np.asarray(sampled["Q"], dtype=np.float64)
    exact_evaluation = scenario_certificate(c_matrix, q_value, scenarios)
    true_radius = float(exact_evaluation["certified_radius"])
    training_radius = float(sampled["training_radius_recomputed"])
    threshold = max(
        absolute_failure_tolerance,
        relative_failure_tolerance * true_radius,
    )
    return {
        "sample_size": len(samples),
        "status": sampled["status"],
        "training_radius": training_radius,
        "exact_finite_radius": true_radius,
        "absolute_underreport": true_radius - training_radius,
        "relative_underreport": (true_radius - training_radius) / true_radius,
        "relative_true_risk_over_exact_design": true_radius / exact_radius - 1.0,
        "underreport_failure": true_radius - training_radius > threshold,
        "failure_threshold": threshold,
        "Q": sampled["Q"],
        "timing_seconds": sampled["timing_seconds"],
        "active_exact_scenarios": exact_evaluation["active_worst_case_scenarios"],
    }


def run_multi_instance_benchmark(
    config: Mapping[str, Any],
    *,
    solver: str | None,
    solver_options: Mapping[str, Any] | None,
    deployment_tolerance: float,
) -> dict[str, Any]:
    q_values = [int(value) for value in config.get("q_values", [3, 4, 5])]
    h_values = [int(value) for value in config.get("h_values", [1, 2, 4])]
    sample_size = int(config.get("sample_size", 32))
    replicates = int(config.get("sampling_replicates", 5))
    seed = int(config.get("seed", 20260914))
    families = list(config.get("families", ["random_cyclic", "grid", "geometric"]))
    weight_min = float(config.get("weight_min", 0.01))
    weight_max = float(config.get("weight_max", 100.0))
    graph_options = dict(config.get("graph_options", {"mean_degree": 4.0}))
    absolute_tolerance = float(config.get("failure_absolute_tolerance", 2e-5))
    relative_tolerance = float(config.get("failure_relative_tolerance", 1e-3))

    if not q_values or min(q_values) < 3 or not h_values or min(h_values) < 0:
        raise ValueError("benchmark q_values/h_values are invalid")
    if sample_size < 1 or replicates < 1:
        raise ValueError("benchmark sample_size and replicates must be positive")

    work_items = [
        (q, h, definition)
        for q in q_values
        for definition in frozen_benchmark_definitions(q)
        for h in h_values
    ]
    child_seeds = iter(np.random.SeedSequence(seed).spawn(len(work_items) * replicates))
    unbounded_cache: dict[str, dict[str, Any]] = {}
    relaxation_cache: dict[str, dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []

    for q, h, definition in work_items:
        c_matrix = np.asarray(definition["C"], dtype=np.float64)
        rule = LocalRuleSpec.from_dict(definition["local_rule"])
        scenarios = finite_scenarios(q, h)
        exact_raw = solve_exact_spectral(
            c_matrix,
            q=q,
            h=h,
            mode="finite",
            local_rule=rule,
            prune_vertices=True,
            solver=solver,
            solver_options=dict(solver_options or {}),
        )
        exact, deployed_q, exact_radius = _deployed_exact_record(
            exact_raw,
            c_matrix,
            scenarios,
            tolerance=deployment_tolerance,
        )

        cache_key = str(definition["problem_id"])
        if cache_key not in unbounded_cache:
            unbounded_cache[cache_key] = solve_exact_spectral(
                c_matrix,
                q=q,
                mode="unbounded",
                local_rule=rule,
                solver=solver,
                solver_options=dict(solver_options or {}),
            ).to_dict()
        if cache_key not in relaxation_cache:
            relaxation_cache[cache_key] = solve_psd_contraction_relaxation(
                c_matrix,
                local_rule=rule,
                solver=solver,
                solver_options=dict(solver_options or {}),
                projection_audit_samples=int(config.get("projection_audit_samples", 32)),
                seed=seed + q,
            )

        unbounded = unbounded_cache[cache_key]
        relaxation = relaxation_cache[cache_key]
        relaxation_q = np.asarray(relaxation["Q"], dtype=np.float64)
        relaxation_true = scenario_certificate(c_matrix, relaxation_q, scenarios)

        sampled_records: list[dict[str, Any]] = []
        for replicate in range(replicates):
            child_seed = int(next(child_seeds).generate_state(1, dtype=np.uint32)[0])
            samples = generate_completion_samples(
                q,
                h,
                sample_size,
                families=families,
                seed=child_seed,
                weight_min=weight_min,
                weight_max=weight_max,
                graph_options=graph_options,
            )
            record = _sampling_record(
                c_matrix,
                rule,
                scenarios,
                samples,
                exact_radius=exact_radius,
                solver=solver,
                solver_options=solver_options,
                absolute_failure_tolerance=absolute_tolerance,
                relative_failure_tolerance=relative_tolerance,
            )
            record["replicate"] = replicate
            record["seed"] = child_seed
            sampled_records.append(record)

        hop = definition["local_rule"]["r_hop"]
        resources = direct_rule_resources(
            deployed_q,
            np.asarray(hop["adjacency"], dtype=np.float64),
            hop["output_nodes"],
            zero_tolerance=deployment_tolerance,
        )
        sampling_overheads = [
            float(item["relative_true_risk_over_exact_design"])
            for item in sampled_records
        ]
        sampling_underreports = [
            float(item["relative_underreport"]) for item in sampled_records
        ]
        relaxation_true_radius = float(relaxation_true["certified_radius"])
        rows.append(
            {
                "q": q,
                "h": h,
                "problem_id": definition["problem_id"],
                "family": definition["family"],
                "rule_name": definition["rule_name"],
                "rounds": definition["rounds"],
                "C": definition["C"],
                "local_rule": definition["local_rule"],
                "exact_finite": exact,
                "exact_unbounded": unbounded,
                "finite_budget_radius_reduction_from_unbounded": (
                    float(unbounded["radius"]) - exact_radius
                )
                / float(unbounded["radius"]),
                "psd_contraction": {
                    "reported_radius": relaxation["certified_radius"],
                    "exact_finite_radius": relaxation_true_radius,
                    "relative_true_risk_over_exact_design": (
                        relaxation_true_radius / exact_radius - 1.0
                    ),
                    "relative_certificate_margin_over_true_risk": (
                        float(relaxation["certified_radius"]) / relaxation_true_radius
                        - 1.0
                    ),
                    "status": relaxation["status"],
                },
                "random_sampling": {
                    "sample_size": sample_size,
                    "replicates": sampled_records,
                    "true_risk_overhead_summary": _quantile_summary(sampling_overheads),
                    "underreport_summary": _quantile_summary(sampling_underreports),
                    "failure_count": sum(
                        bool(item["underreport_failure"]) for item in sampled_records
                    ),
                },
                "online_resources": resources,
            }
        )

    sampling_overheads = [
        float(sample["relative_true_risk_over_exact_design"])
        for row in rows
        for sample in row["random_sampling"]["replicates"]
    ]
    psd_overheads = [
        float(row["psd_contraction"]["relative_true_risk_over_exact_design"])
        for row in rows
    ]
    finite_gains = [
        float(row["finite_budget_radius_reduction_from_unbounded"]) for row in rows
    ]

    constant_rows = [row for row in rows if row["family"] == "constant_preserving_state"]
    paired: dict[tuple[int, int], dict[str, float]] = {}
    for row in constant_rows:
        key = (int(row["q"]), int(row["h"]))
        paired.setdefault(key, {})[str(row["rule_name"])] = float(
            row["exact_finite"]["radius"]
        )
    round_tradeoffs = []
    for (q, h), radii in sorted(paired.items()):
        required = {"zero_hop", "one_hop", "full_port"}
        if set(radii) != required:
            raise RuntimeError(f"incomplete round comparison for q={q}, h={h}")
        round_tradeoffs.append(
            {
                "q": q,
                "h": h,
                "zero_hop_radius": radii["zero_hop"],
                "one_hop_radius": radii["one_hop"],
                "full_port_radius": radii["full_port"],
                "one_hop_reduction_from_zero": 1.0
                - radii["one_hop"] / radii["zero_hop"],
                "full_port_reduction_from_zero": 1.0
                - radii["full_port"] / radii["zero_hop"],
            }
        )

    return {
        "classification": "theorem_matched_fixed_design_grid",
        "grid": {
            "q_values": q_values,
            "h_values": h_values,
            "problem_families_per_q": 6,
            "problem_count": len(rows),
            "random_sampling_replicates_per_problem": replicates,
            "random_sampling_size": sample_size,
            "random_sampling_distribution": {
                "families": families,
                "weight_range": [weight_min, weight_max],
                "nested": False,
            },
        },
        "design_protocol": {
            "material_radius_gain_threshold": 0.05,
            "material_squared_risk_gain_threshold": 0.10,
            "sampling_failure_definition": (
                "R_exact(Q_N)-R_train > max(2e-5, 1e-3 R_exact(Q_N))"
            ),
            "all_configurations_reported_without_result_based_filtering": True,
        },
        "rows": rows,
        "round_tradeoffs": round_tradeoffs,
        "summary": {
            "sampling_true_risk_overhead": _quantile_summary(sampling_overheads),
            "psd_true_risk_overhead": _quantile_summary(psd_overheads),
            "finite_budget_radius_gain": _quantile_summary(finite_gains),
            "one_hop_radius_reduction_from_zero": _quantile_summary(
                [row["one_hop_reduction_from_zero"] for row in round_tradeoffs]
            ),
            "full_port_radius_reduction_from_zero": _quantile_summary(
                [row["full_port_reduction_from_zero"] for row in round_tradeoffs]
            ),
            "solver_status_counts": dict(
                Counter(
                    str(row["exact_finite"]["status"])
                    for row in rows
                )
            ),
        },
    }


def _sampling_curve_cases() -> list[dict[str, Any]]:
    q_state = 4
    q_opcode = 3
    opcode_rule = _r_hop_rule(q_opcode, 0)
    opcode_rule["shared_coefficients"] = [
        [[index, index] for index in range(q_opcode)]
    ]
    return [
        {
            "case_id": "constant_preserving_state",
            "label": "State smoother, one hop",
            "q": q_state,
            "h": 1,
            "C": np.eye(q_state).tolist(),
            "local_rule": _r_hop_rule(q_state, 1, row_sums=1.0),
        },
        {
            "case_id": "shared_opcode",
            "label": "Shared diagonal opcode",
            "q": q_opcode,
            "h": 2,
            "C": np.diag([0.5, 0.75, 1.0]).tolist(),
            "local_rule": opcode_rule,
        },
    ]


def run_sampling_budget_curve(
    config: Mapping[str, Any],
    *,
    solver: str | None,
    solver_options: Mapping[str, Any] | None,
    deployment_tolerance: float,
) -> dict[str, Any]:
    sample_sizes = sorted({int(value) for value in config.get("sample_sizes", [4, 8, 16, 32, 64, 128, 256])})
    replicates = int(config.get("replicates", 30))
    seed = int(config.get("seed", 20260915))
    families = list(config.get("families", ["random_cyclic", "grid", "geometric"]))
    weight_min = float(config.get("weight_min", 0.01))
    weight_max = float(config.get("weight_max", 100.0))
    graph_options = dict(config.get("graph_options", {"mean_degree": 4.0}))
    absolute_tolerance = float(config.get("failure_absolute_tolerance", 2e-5))
    relative_tolerance = float(config.get("failure_relative_tolerance", 1e-3))
    if not sample_sizes or sample_sizes[0] < 1 or replicates < 1:
        raise ValueError("sampling curve sizes/replicates are invalid")

    case_seeds = iter(
        np.random.SeedSequence(seed).spawn(len(_sampling_curve_cases()) * replicates)
    )
    case_records: list[dict[str, Any]] = []
    for case in _sampling_curve_cases():
        q = int(case["q"])
        h = int(case["h"])
        c_matrix = np.asarray(case["C"], dtype=np.float64)
        rule = LocalRuleSpec.from_dict(case["local_rule"])
        scenarios = finite_scenarios(q, h)
        raw_exact = solve_exact_spectral(
            c_matrix,
            q=q,
            h=h,
            mode="finite",
            local_rule=rule,
            solver=solver,
            solver_options=dict(solver_options or {}),
        )
        exact, _, exact_radius = _deployed_exact_record(
            raw_exact,
            c_matrix,
            scenarios,
            tolerance=deployment_tolerance,
        )
        records: list[dict[str, Any]] = []
        for replicate in range(replicates):
            child_seed = int(next(case_seeds).generate_state(1, dtype=np.uint32)[0])
            nested_samples = generate_completion_samples(
                q,
                h,
                sample_sizes[-1],
                families=families,
                seed=child_seed,
                weight_min=weight_min,
                weight_max=weight_max,
                graph_options=graph_options,
            )
            for sample_size in sample_sizes:
                record = _sampling_record(
                    c_matrix,
                    rule,
                    scenarios,
                    nested_samples[:sample_size],
                    exact_radius=exact_radius,
                    solver=solver,
                    solver_options=solver_options,
                    absolute_failure_tolerance=absolute_tolerance,
                    relative_failure_tolerance=relative_tolerance,
                )
                record["replicate"] = replicate
                record["seed"] = child_seed
                records.append(record)

        summaries = []
        for sample_size in sample_sizes:
            selected = [row for row in records if row["sample_size"] == sample_size]
            failures = sum(bool(row["underreport_failure"]) for row in selected)
            low, high = wilson_interval(failures, len(selected))
            summaries.append(
                {
                    "sample_size": sample_size,
                    "replicates": len(selected),
                    "failure_count": failures,
                    "failure_rate": failures / len(selected),
                    "failure_rate_wilson95_low": low,
                    "failure_rate_wilson95_high": high,
                    "relative_underreport": _quantile_summary(
                        [float(row["relative_underreport"]) for row in selected]
                    ),
                    "relative_true_risk_over_exact_design": _quantile_summary(
                        [
                            float(row["relative_true_risk_over_exact_design"])
                            for row in selected
                        ]
                    ),
                    "solver_status_counts": dict(
                        Counter(str(row["status"]) for row in selected)
                    ),
                }
            )
        case_records.append(
            {
                **case,
                "exact_design": exact,
                "records": records,
                "summary_by_sample_size": summaries,
            }
        )

    return {
        "classification": "theorem_matched_nested_random_sampling",
        "protocol": {
            "sample_sizes": sample_sizes,
            "replicates_per_case": replicates,
            "same_distribution_for_all_sample_sizes": True,
            "nested_samples_within_each_seed": True,
            "families": families,
            "weight_range": [weight_min, weight_max],
            "failure_absolute_tolerance": absolute_tolerance,
            "failure_relative_tolerance": relative_tolerance,
        },
        "cases": case_records,
    }


def connected_path_completion(
    q: int,
    h: int,
    partition: Sequence[Sequence[int]],
    allocation: Sequence[int],
    *,
    strong_conductance: float,
    weak_conductance: float,
) -> CompletionSample:
    """Construct the connected path used by the finite-h converse."""

    blocks = [tuple(int(vertex) for vertex in block) for block in partition]
    allocation_tuple = tuple(int(value) for value in allocation)
    if len(blocks) != len(allocation_tuple):
        raise ValueError("partition and allocation lengths differ")
    if sorted(vertex for block in blocks for vertex in block) != list(range(q)):
        raise ValueError("partition must cover every port exactly once")
    if any(value < 0 for value in allocation_tuple) or sum(allocation_tuple) > h:
        raise ValueError("allocation must be nonnegative and use at most h vertices")
    if not (0.0 < strong_conductance < math.inf):
        raise ValueError("strong conductance must be finite and positive")
    if not (0.0 < weak_conductance < math.inf):
        raise ValueError("weak conductance must be finite and positive")

    next_hidden = q
    segments: list[list[int]] = []
    for block, hidden_count in zip(blocks, allocation_tuple, strict=True):
        hidden = list(range(next_hidden, next_hidden + hidden_count))
        next_hidden += hidden_count
        segments.append(list(block) + hidden)
    while next_hidden < q + h:
        segments.append([next_hidden])
        next_hidden += 1

    n = q + h
    adjacency = np.zeros((n, n), dtype=np.float64)
    edge_kinds: list[dict[str, Any]] = []
    for segment in segments:
        for left, right in zip(segment[:-1], segment[1:], strict=True):
            adjacency[left, right] = strong_conductance
            adjacency[right, left] = strong_conductance
            edge_kinds.append({"edge": [left, right], "kind": "strong"})
    for left_segment, right_segment in zip(segments[:-1], segments[1:], strict=True):
        left = left_segment[-1]
        right = right_segment[0]
        adjacency[left, right] = weak_conductance
        adjacency[right, left] = weak_conductance
        edge_kinds.append({"edge": [left, right], "kind": "weak"})

    degrees = np.count_nonzero(adjacency, axis=1)
    laplacian = np.diag(np.sum(adjacency, axis=1)) - adjacency
    connected = n == 1 or np.linalg.matrix_rank(laplacian, tol=1e-9) == n - 1
    if not connected or int(np.max(degrees, initial=0)) > 2:
        raise RuntimeError("constructed witness is not a connected degree-two path")
    smoother = np.linalg.solve(np.eye(n) + laplacian, np.eye(n))
    return CompletionSample(
        q=q,
        hidden_count=h,
        smoother=smoother,
        adjacency=adjacency,
        metadata={
            "classification": "theorem_matching_connected_path_witness",
            "family": "connected_path_witness",
            "seed": None,
            "nodes": n,
            "edges": max(0, n - 1),
            "cycle_rank": 0,
            "contains_cycle": False,
            "port_cut_edges": int(np.count_nonzero(adjacency[:q, q:])),
            "port_cut_weight": float(np.sum(adjacency[:q, q:])),
            "maximum_degree": int(np.max(degrees, initial=0)),
            "connected": connected,
            "segments": segments,
            "edge_kinds": edge_kinds,
            "strong_conductance": strong_conductance,
            "weak_conductance": weak_conductance,
        },
    )


def _select_active_vertices(
    c_matrix: Array,
    q_matrix: Array,
    scenarios: Sequence[Scenario],
    maximum_count: int,
) -> list[Scenario]:
    risks = [spectral_squared_risk(c_matrix, q_matrix, scenario.matrix) for scenario in scenarios]
    maximum = max(risks)
    tolerance = 2e-6 * max(1.0, maximum)
    active = [
        scenario
        for scenario, risk in zip(scenarios, risks, strict=True)
        if scenario.is_vertex and maximum - risk <= tolerance
    ]
    chosen: list[Scenario] = []
    for block_count in sorted({len(scenario.partition) for scenario in active}):
        candidates = [scenario for scenario in active if len(scenario.partition) == block_count]
        candidates.sort(
            key=lambda scenario: (
                -(sum(scenario.allocation or ())),
                scenario.label,
            )
        )
        chosen.append(candidates[0])
        if len(chosen) == maximum_count:
            break
    if len(chosen) < min(maximum_count, len(active)):
        for scenario in sorted(active, key=lambda item: item.label):
            if scenario not in chosen:
                chosen.append(scenario)
            if len(chosen) == maximum_count:
                break
    return chosen


def run_path_tightness(
    config: Mapping[str, Any],
    *,
    solver: str | None,
    solver_options: Mapping[str, Any] | None,
    deployment_tolerance: float,
) -> dict[str, Any]:
    q = int(config.get("q", 4))
    h = int(config.get("h", 1))
    strong_values = [float(value) for value in config.get("strong_values", np.logspace(1, 6, 11))]
    maximum_witnesses = int(config.get("maximum_witnesses", 3))
    weak_exponent = float(config.get("weak_exponent", -1.0))
    final_tolerance = float(config.get("final_relative_tolerance", 1e-4))
    c_matrix = np.eye(q)
    rule_payload = _r_hop_rule(q, 1, row_sums=1.0)
    rule = LocalRuleSpec.from_dict(rule_payload)
    scenarios = finite_scenarios(q, h)
    raw_exact = solve_exact_spectral(
        c_matrix,
        q=q,
        h=h,
        mode="finite",
        local_rule=rule,
        solver=solver,
        solver_options=dict(solver_options or {}),
    )
    exact, q_matrix, exact_radius = _deployed_exact_record(
        raw_exact,
        c_matrix,
        scenarios,
        tolerance=deployment_tolerance,
    )
    selected = _select_active_vertices(
        c_matrix, q_matrix, scenarios, maximum_witnesses
    )
    witnesses: list[dict[str, Any]] = []
    invariant_checks: list[bool] = []
    final_checks: list[bool] = []
    for witness_index, scenario in enumerate(selected):
        canonical_radius = math.sqrt(
            spectral_squared_risk(c_matrix, q_matrix, scenario.matrix)
        )
        rows = []
        for strong in strong_values:
            weak = strong**weak_exponent
            completion = connected_path_completion(
                q,
                h,
                scenario.partition,
                scenario.allocation or (),
                strong_conductance=strong,
                weak_conductance=weak,
            )
            physical_radius = float(
                np.linalg.norm(full_error_operator(c_matrix, q_matrix, completion), ord=2)
            )
            port_block_distance = float(
                np.linalg.norm(
                    completion.smoother[:q, :q] - scenario.matrix,
                    ord=2,
                )
            )
            relative_radius_gap = abs(physical_radius - canonical_radius) / max(
                canonical_radius, np.finfo(float).tiny
            )
            invariant = (
                completion.hidden_count == h
                and bool(completion.metadata["connected"])
                and int(completion.metadata["maximum_degree"]) <= 2
                and np.all(completion.adjacency[np.nonzero(completion.adjacency)] > 0.0)
            )
            invariant_checks.append(invariant)
            rows.append(
                {
                    "strong_conductance": strong,
                    "weak_conductance": weak,
                    "physical_radius": physical_radius,
                    "canonical_radius": canonical_radius,
                    "physical_over_canonical": physical_radius / canonical_radius,
                    "relative_radius_gap": relative_radius_gap,
                    "port_block_operator_distance": port_block_distance,
                    "connected": completion.metadata["connected"],
                    "maximum_degree": completion.metadata["maximum_degree"],
                    "hidden_count": completion.hidden_count,
                }
            )
        final_checks.append(
            rows[-1]["relative_radius_gap"] <= final_tolerance
            and rows[-1]["port_block_operator_distance"] <= final_tolerance
        )
        witnesses.append(
            {
                "witness_index": witness_index,
                "scenario": scenario.to_dict(include_matrix=True),
                "canonical_radius": canonical_radius,
                "rows": rows,
            }
        )

    return {
        "classification": "theorem_matched_connected_degree_two_path_tightness",
        "problem": {
            "q": q,
            "h": h,
            "C": c_matrix.tolist(),
            "local_rule": rule_payload,
            "exact_design": exact,
            "exact_radius": exact_radius,
        },
        "protocol": {
            "strong_values": strong_values,
            "weak_rule": f"epsilon=t^({weak_exponent:g})",
            "maximum_witnesses": maximum_witnesses,
            "final_relative_tolerance": final_tolerance,
        },
        "witnesses": witnesses,
        "checks": {
            "all_graphs_connected_degree_two_with_exact_hidden_budget": all(
                invariant_checks
            ),
            "all_final_witness_gaps_below_tolerance": all(final_checks),
        },
    }


def run_core_evidence(config: Mapping[str, Any]) -> dict[str, Any]:
    solver_config = dict(config.get("solver", {}))
    solver = solver_config.get("name")
    solver_options = solver_config.get("options")
    deployment_tolerance = float(config.get("deployment_zero_tolerance", 2e-5))
    multi = run_multi_instance_benchmark(
        dict(config.get("multi_instance", {})),
        solver=solver,
        solver_options=solver_options,
        deployment_tolerance=deployment_tolerance,
    )
    sampling = run_sampling_budget_curve(
        dict(config.get("sampling_curve", {})),
        solver=solver,
        solver_options=solver_options,
        deployment_tolerance=deployment_tolerance,
    )
    path = run_path_tightness(
        dict(config.get("path_tightness", {})),
        solver=solver,
        solver_options=solver_options,
        deployment_tolerance=deployment_tolerance,
    )

    numerical_tolerance = float(config.get("numerical_check_tolerance", 2e-5))
    rows = multi["rows"]
    exact_dominates = all(
        float(sample["exact_finite_radius"])
        + numerical_tolerance
        >= float(row["exact_finite"]["radius"])
        for row in rows
        for sample in row["random_sampling"]["replicates"]
    ) and all(
        float(row["psd_contraction"]["exact_finite_radius"])
        + numerical_tolerance
        >= float(row["exact_finite"]["radius"])
        for row in rows
    )
    finite_not_worse = all(
        float(row["exact_finite"]["radius"])
        <= float(row["exact_unbounded"]["radius"]) + numerical_tolerance
        for row in rows
    )
    checks = {
        "fixed_multi_instance_count_is_54": multi["grid"]["problem_count"] == 54,
        "all_baselines_re_evaluated_no_better_than_exact_within_tolerance": exact_dominates,
        "finite_optimum_not_above_unbounded_optimum": finite_not_worse,
        **path["checks"],
    }
    if not all(checks.values()):
        raise RuntimeError(f"core-evidence checks failed: {checks}")
    return {
        "schema_version": 1,
        "classification": "tsp_core_theorem_matched_evidence",
        "semantics": {
            "target": "C E_Gamma (I+L_G)^-1 y",
            "local_rule": "Q E_Gamma y, common to every completion",
            "certified_error": "full-input spectral operator norm",
            "communication_graph": "fixed protocol overlay H, independent of physical completion G",
            "baseline_policy": "every returned Q is re-evaluated on every exact finite scenario",
        },
        "deployment_zero_tolerance": deployment_tolerance,
        "multi_instance": multi,
        "sampling_curve": sampling,
        "path_tightness": path,
        "checks": checks,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the theorem-matched evidence suite used by the TSP draft."
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    config_path = arguments.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    repository_root = Path(__file__).resolve().parents[2]
    payload = run_core_evidence(config)
    payload["config_path"] = record_path(config_path, repository_root)
    payload["provenance"] = experiment_provenance(repository_root)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(arguments.output)


if __name__ == "__main__":
    main()
