"""Fixed-known-graph graph-smoother baselines and external scale tests."""

from __future__ import annotations

from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
from numpy.typing import NDArray
import scipy.sparse as sp
import scipy.sparse.csgraph as csgraph

from .algorithms import (
    chebyshev_apply_coefficients,
    chebyshev_coefficients,
    gershgorin_spectral_interval,
    global_smoother,
    induced_local_operator,
    neumann_apply,
)
from .graphs import (
    WeightedGraph,
    load_pglib_topology,
    load_suitesparse_pattern,
    smoother_matrix,
    synthetic_weighted_graph,
)
from .resources import (
    PeakRssSampler,
    ball_nodes,
    csr_storage_bytes,
    local_collection_resources,
    polynomial_iteration_resources,
)


Array = NDArray[np.float64]


def _resolve_path(value: str, base_directory: Path | None) -> Path:
    path = Path(value)
    if path.is_absolute() or base_directory is None:
        return path
    return (base_directory / path).resolve()


def _load_graph(specification: dict[str, Any], base_directory: Path | None) -> WeightedGraph:
    source = str(specification.get("source", "synthetic"))
    if source == "synthetic":
        return synthetic_weighted_graph(
            str(specification["family"]),
            int(specification["n"]),
            int(specification["seed"]),
            weight_min=float(specification.get("weight_min", 0.25)),
            weight_max=float(specification.get("weight_max", 4.0)),
            options=specification,
        )
    if source == "pglib":
        return load_pglib_topology(
            _resolve_path(str(specification["path"]), base_directory),
            normalize_median=bool(specification.get("normalize_median", True)),
        )
    if source == "suitesparse":
        return load_suitesparse_pattern(
            _resolve_path(str(specification["archive"]), base_directory),
            member=specification.get("member"),
            use_value_magnitudes=bool(specification.get("use_value_magnitudes", False)),
        )
    raise ValueError(f"unknown stress-graph source {source!r}")


def _signal_metrics(approximation: Array, reference: Array) -> dict[str, float]:
    error = np.asarray(approximation) - np.asarray(reference)
    denominator = max(float(np.linalg.norm(reference)), 1e-300)
    node_rms = np.sqrt(np.mean(np.square(error), axis=1))
    return {
        "empirical_mse_per_node": float(np.mean(np.square(error))),
        "empirical_relative_l2_error": float(np.linalg.norm(error) / denominator),
        "maximum_node_rmse": float(np.max(node_rms)),
        "median_node_rmse": float(np.median(node_rms)),
    }


def fixed_operator_metrics(approximation: Array, exact: Array) -> dict[str, float]:
    """Exact fixed-instance metrics; these are not hidden-completion certificates."""

    error = np.asarray(approximation, dtype=np.float64) - np.asarray(exact, dtype=np.float64)
    return {
        "fixed_instance_operator_2_norm_error": float(np.linalg.norm(error, ord=2)),
        "isotropic_expected_mse_per_node": float(np.sum(np.square(error)) / error.shape[0]),
        "fixed_instance_frobenius_error": float(np.linalg.norm(error, ord="fro")),
    }


def _cut_statistics(adjacency: sp.csr_matrix, roots: Array, radius: int) -> dict[str, float | int]:
    counts: list[int] = []
    weights: list[float] = []
    for root in roots:
        inside_nodes = ball_nodes(adjacency, int(root), radius)
        inside = np.zeros(adjacency.shape[0], dtype=bool)
        inside[inside_nodes] = True
        cut = adjacency[inside][:, ~inside]
        counts.append(int(cut.nnz))
        weights.append(float(cut.sum()))
    return {
        "sampled_roots": len(roots),
        "mean_cut_edges": float(np.mean(counts)) if counts else 0.0,
        "maximum_cut_edges": max(counts, default=0),
        "mean_cut_weight": float(np.mean(weights)) if weights else 0.0,
        "maximum_cut_weight": max(weights, default=0.0),
    }


def _graph_summary(graph: WeightedGraph, system: sp.csr_matrix) -> dict[str, Any]:
    adjacency = graph.adjacency
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    components, _ = csgraph.connected_components(adjacency, directed=False)
    upper = sp.triu(adjacency, k=1)
    weights = upper.data
    return {
        "nodes": graph.nodes,
        "edges": graph.edges,
        "connected_components": int(components),
        "weighted_degree_min": float(np.min(degree)),
        "weighted_degree_median": float(np.median(degree)),
        "weighted_degree_max": float(np.max(degree)),
        "edge_weight_min": float(np.min(weights)) if weights.size else 0.0,
        "edge_weight_median": float(np.median(weights)) if weights.size else 0.0,
        "edge_weight_max": float(np.max(weights)) if weights.size else 0.0,
        "adjacency_csr_bytes": csr_storage_bytes(adjacency),
        "system_nnz": int(system.nnz),
        "system_csr_bytes": csr_storage_bytes(system),
    }


def _one_case(
    specification: dict[str, Any],
    config: dict[str, Any],
    base_directory: Path | None,
    case_index: int,
) -> dict[str, Any]:
    case_started = perf_counter()
    graph_started = perf_counter()
    with PeakRssSampler() as graph_memory:
        graph = _load_graph(specification, base_directory)
        edge_scale = float(specification.get("edge_scale", config.get("edge_scale", 1.0)))
        system = smoother_matrix(graph, edge_scale=edge_scale)
    graph_seconds = perf_counter() - graph_started
    n = graph.nodes
    signal_count = int(specification.get("signal_samples", config.get("signal_samples", 16)))
    seed = int(specification.get("signal_seed", int(config.get("seed", 20260914)) + 1000 * case_index))
    rng = np.random.default_rng(seed)
    signals = rng.standard_normal((n, signal_count))
    exact_operator_max_n = int(config.get("exact_operator_max_n", 96))
    audit_operator = n <= exact_operator_max_n
    with PeakRssSampler() as global_memory:
        global_result = global_smoother(
            system,
            signals,
            direct_max_n=int(config.get("direct_max_n", 50_000)),
            relative_tolerance=float(config.get("reference_relative_tolerance", 1e-10)),
            max_iterations=config.get("reference_max_iterations"),
        )
    exact_signals = np.asarray(global_result.solution)
    exact_operator = None
    operator_reference_seconds = 0.0
    if audit_operator:
        operator_reference_started = perf_counter()
        exact_operator = np.asarray(
            global_smoother(
                system,
                np.eye(n),
                direct_max_n=int(config.get("direct_max_n", 50_000)),
                relative_tolerance=float(config.get("reference_relative_tolerance", 1e-10)),
                max_iterations=config.get("reference_max_iterations"),
            ).solution
        )
        operator_reference_seconds = perf_counter() - operator_reference_started
    global_record: dict[str, Any] = {
        "access_class": "full_graph_and_full_input",
        "rounds": None,
        "signal_metrics": _signal_metrics(exact_signals, exact_signals),
        "offline": {
            "factorization_seconds": global_result.metadata["factorization_seconds"],
        },
        "online": {
            "solve_seconds": global_result.metadata["solve_seconds"],
            "right_hand_sides_in_timed_solve": signal_count,
        },
        "relative_residual_max": global_result.metadata["relative_residual_max"],
        "memory": global_memory.to_dict(),
    }
    if exact_operator is not None:
        global_record["operator_audit"] = fixed_operator_metrics(exact_operator, exact_operator)
        global_record["operator_audit_reference_seconds"] = operator_reference_seconds

    lower, upper = gershgorin_spectral_interval(graph, edge_scale=edge_scale)
    root_count = min(n, int(config.get("cut_root_samples", 32)))
    roots = np.sort(rng.choice(n, size=root_count, replace=False)).astype(np.int64)
    scalar_bytes = int(config.get("scalar_bytes", 8))
    local_max_n = int(config.get("local_truncation_max_n", 128))
    methods_by_round: list[dict[str, Any]] = []
    for rounds in [int(value) for value in config["rounds"]]:
        row: dict[str, Any] = {
            "rounds": rounds,
            "cut": _cut_statistics(graph.adjacency, roots, rounds),
            "methods": {},
        }

        if bool(config.get("enable_local_truncation", True)) and n <= local_max_n:
            with PeakRssSampler() as local_memory:
                construction_started = perf_counter()
                local_operator = induced_local_operator(graph.adjacency, rounds)
                construction_seconds = perf_counter() - construction_started
                apply_started = perf_counter()
                local_values = np.asarray(local_operator @ signals)
                apply_seconds = perf_counter() - apply_started
            local_record: dict[str, Any] = {
                "access_class": "known_induced_radius_r_subgraph_and_local_inputs",
                "boundary_condition": "unit-grounded smoother on the induced subgraph",
                "signal_metrics": _signal_metrics(local_values, exact_signals),
                "offline": {
                    "all_root_operator_construction_seconds": construction_seconds,
                    "operator_nnz": int(local_operator.nnz),
                    "operator_storage_bytes": csr_storage_bytes(local_operator),
                },
                "online": {
                    "batch_sparse_apply_seconds": apply_seconds,
                    "signal_count": signal_count,
                },
                "resources": local_collection_resources(
                    graph.adjacency,
                    range(n),
                    rounds,
                    right_hand_sides=signal_count,
                    scalar_bytes=scalar_bytes,
                ),
                "memory": local_memory.to_dict(),
            }
            if exact_operator is not None:
                local_record["operator_audit"] = fixed_operator_metrics(
                    local_operator.toarray(), exact_operator
                )
            row["methods"]["induced_local_truncation"] = local_record
        else:
            row["methods"]["induced_local_truncation"] = {
                "skipped": True,
                "reason": (
                    "disabled by configuration"
                    if not bool(config.get("enable_local_truncation", True))
                    else f"n={n} exceeds local_truncation_max_n={local_max_n}"
                ),
            }

        with PeakRssSampler() as neumann_memory:
            neumann_started = perf_counter()
            neumann_values = np.asarray(neumann_apply(system, signals, rounds, upper))
            neumann_seconds = perf_counter() - neumann_started
        neumann_record: dict[str, Any] = {
            "access_class": "degree_r_local_polynomial_with_global_Gershgorin_bound",
            "signal_metrics": _signal_metrics(neumann_values, exact_signals),
            "offline": {"spectral_upper_bound": upper, "coefficient_seconds": 0.0},
            "online": {"batch_apply_seconds": neumann_seconds, "signal_count": signal_count},
            "resources": polynomial_iteration_resources(
                graph.adjacency,
                system,
                rounds,
                right_hand_sides=signal_count,
                scalar_bytes=scalar_bytes,
                recurrence_vectors=2,
            ),
            "memory": neumann_memory.to_dict(),
        }
        if exact_operator is not None:
            audit_started = perf_counter()
            neumann_operator = neumann_apply(system, np.eye(n), rounds, upper)
            neumann_record["operator_audit"] = fixed_operator_metrics(
                neumann_operator, exact_operator
            )
            neumann_record["operator_audit_seconds"] = perf_counter() - audit_started
        row["methods"]["neumann"] = neumann_record

        coefficient_started = perf_counter()
        coefficients = chebyshev_coefficients(rounds, lower, upper)
        coefficient_seconds = perf_counter() - coefficient_started
        with PeakRssSampler() as chebyshev_memory:
            chebyshev_started = perf_counter()
            chebyshev_values = np.asarray(
                chebyshev_apply_coefficients(system, signals, coefficients, lower, upper)
            )
            chebyshev_seconds = perf_counter() - chebyshev_started
        chebyshev_record: dict[str, Any] = {
            "access_class": "degree_r_local_polynomial_with_global_Gershgorin_interval",
            "signal_metrics": _signal_metrics(chebyshev_values, exact_signals),
            "offline": {
                "spectral_interval": [lower, upper],
                "coefficient_seconds": coefficient_seconds,
                "coefficient_count": int(coefficients.size),
            },
            "online": {"batch_apply_seconds": chebyshev_seconds, "signal_count": signal_count},
            "resources": polynomial_iteration_resources(
                graph.adjacency,
                system,
                rounds,
                right_hand_sides=signal_count,
                scalar_bytes=scalar_bytes,
                recurrence_vectors=3,
            ),
            "memory": chebyshev_memory.to_dict(),
        }
        if exact_operator is not None:
            audit_started = perf_counter()
            chebyshev_operator = chebyshev_apply_coefficients(
                system, np.eye(n), coefficients, lower, upper
            )
            chebyshev_record["operator_audit"] = fixed_operator_metrics(
                chebyshev_operator, exact_operator
            )
            chebyshev_record["operator_audit_seconds"] = perf_counter() - audit_started
        row["methods"]["chebyshev"] = chebyshev_record
        methods_by_round.append(row)

    checks = {
        "reference_residual_below_tolerance": global_result.metadata["relative_residual_max"]
        <= max(1e-9, 20 * float(config.get("reference_relative_tolerance", 1e-10))),
        "global_empirical_error_is_zero": global_record["signal_metrics"]["empirical_mse_per_node"]
        == 0.0,
        "all_reported_errors_finite": all(
            np.isfinite(method["signal_metrics"]["empirical_mse_per_node"])
            for row in methods_by_round
            for method in row["methods"].values()
            if not method.get("skipped", False)
        ),
    }
    if not all(checks.values()):
        raise RuntimeError(f"fixed-graph stress checks failed: {checks}")
    return {
        "case_index": case_index,
        "classification": (
            "external_application_stress"
            if graph.metadata.get("classification") == "external_application_stress"
            else "fixed_graph_synthetic_stress"
        ),
        "source": graph.metadata,
        "graph": _graph_summary(graph, system),
        "edge_scale": edge_scale,
        "signal_seed": seed,
        "signal_samples": signal_count,
        "spectral_interval_used": [lower, upper],
        "construction": {"seconds": graph_seconds, "memory": graph_memory.to_dict()},
        "global_exact": global_record,
        "round_sweep": methods_by_round,
        "checks": checks,
        "case_total_seconds": perf_counter() - case_started,
    }


def run_stress_benchmark(
    config: dict[str, Any],
    *,
    base_directory: str | Path | None = None,
) -> dict[str, Any]:
    """Run global, truncation, Neumann, and Chebyshev on fixed graphs."""

    base = Path(base_directory) if base_directory is not None else None
    cases = [
        _one_case(dict(specification), config, base, index)
        for index, specification in enumerate(config["cases"])
    ]
    return {
        "schema_version": 1,
        "classification": "graph_smoother_application_stress",
        "semantics": {
            "objective": "0.5||x-y||^2 + 0.5 sum_edges w_ij(x_i-x_j)^2",
            "target": "(I+L)^-1 y on one fully specified graph",
            "operator_audit": "exact fixed-instance error, not a hidden-completion certificate",
            "empirical_mse": "synthetic graph-signal error, not WLS/Kalman measurement-noise MSE",
        },
        "config": config,
        "cases": cases,
        "checks": {
            "all_cases_passed": all(all(case["checks"].values()) for case in cases),
            "contains_external_data": any(
                case["classification"] == "external_application_stress" for case in cases
            ),
            "contains_synthetic_data": any(
                case["classification"] == "fixed_graph_synthetic_stress" for case in cases
            ),
        },
    }
