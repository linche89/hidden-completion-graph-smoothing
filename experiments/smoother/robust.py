"""Theorem-matching robust-design experiments and deliberately weaker baselines."""

from __future__ import annotations

from dataclasses import dataclass
import math
from time import perf_counter
import tracemalloc
from typing import Any, Iterable

import cvxpy as cp
import numpy as np
from numpy.typing import NDArray
import scipy.sparse as sp

from hidden_completion.local_rules import LocalRuleSpec
from hidden_completion.metrics import spectral_squared_risk
from hidden_completion.scenarios import (
    Scenario,
    finite_scenarios,
    partial_partition_scenarios,
    scenario_count,
    vertex_count,
)
from hidden_completion.solver import ExactSolverError, solve_exact_spectral

from .graphs import WeightedGraph, smoother_matrix, synthetic_weighted_graph
from .resources import PeakRssSampler, direct_rule_resources


Array = NDArray[np.float64]


@dataclass(slots=True)
class CompletionSample:
    """One physical finite-conductance graph completion."""

    q: int
    hidden_count: int
    smoother: Array
    adjacency: Array
    metadata: dict[str, Any]

    @property
    def nodes(self) -> int:
        return int(self.smoother.shape[0])


def _choose_solver(requested: str | None) -> str:
    installed = set(cp.installed_solvers())
    if requested is not None:
        normalized = requested.upper()
        if normalized not in installed:
            raise ExactSolverError(
                f"requested solver {normalized!r} unavailable; installed={sorted(installed)}"
            )
        return normalized
    for candidate in ("CLARABEL", "SCS", "CVXOPT", "MOSEK"):
        if candidate in installed:
            return candidate
    raise ExactSolverError("no SDP-capable CVXPY solver is installed")


def generate_completion(
    q: int,
    hidden_count: int,
    family: str,
    seed: int,
    *,
    weight_min: float,
    weight_max: float,
    graph_options: dict[str, Any] | None = None,
) -> CompletionSample:
    """Generate a unit-grounded physical completion with ports ``0,...,q-1``."""

    if not isinstance(q, int) or q < 1:
        raise ValueError("q must be a positive integer")
    if not isinstance(hidden_count, int) or hidden_count < 0:
        raise ValueError("hidden_count must be a nonnegative integer")
    n = q + hidden_count
    if n == 1:
        graph = WeightedGraph(
            sp.csr_matrix((1, 1), dtype=np.float64),
            {
                "classification": "synthetic",
                "family": "singleton",
                "seed": seed,
                "nodes": 1,
                "edges": 0,
                "weight_min_requested": weight_min,
                "weight_max_requested": weight_max,
            },
        )
    else:
        graph = synthetic_weighted_graph(
            family,
            n,
            seed,
            weight_min=weight_min,
            weight_max=weight_max,
            options=graph_options,
        )
    system = smoother_matrix(graph).toarray()
    smoother = np.linalg.solve(system, np.eye(n))
    adjacency = np.asarray(graph.adjacency.toarray(), dtype=np.float64)
    port_hidden = adjacency[:q, q:]
    port_cut_edges = int(np.count_nonzero(port_hidden))
    port_cut_weight = float(np.sum(port_hidden))
    cycle_rank = graph.edges - graph.nodes + 1
    eigenvalues = np.linalg.eigvalsh((smoother + smoother.T) / 2.0)
    metadata = dict(graph.metadata)
    metadata.update(
        {
            "classification": "theorem_matching",
            "unit_grounding": True,
            "q": q,
            "hidden_count": hidden_count,
            "port_cut_edges": port_cut_edges,
            "port_cut_weight": port_cut_weight,
            "cycle_rank": cycle_rank,
            "contains_cycle": cycle_rank > 0,
            "smoother_eigenvalue_min": float(eigenvalues[0]),
            "smoother_eigenvalue_max": float(eigenvalues[-1]),
        }
    )
    return CompletionSample(q, hidden_count, smoother, adjacency, metadata)


def generate_completion_samples(
    q: int,
    h: int,
    count: int,
    *,
    families: Iterable[str],
    seed: int,
    weight_min: float,
    weight_max: float,
    graph_options: dict[str, Any] | None = None,
) -> list[CompletionSample]:
    if h < 0 or count < 1:
        raise ValueError("h must be nonnegative and count positive")
    family_tuple = tuple(families)
    if not family_tuple:
        raise ValueError("at least one completion family is required")
    samples: list[CompletionSample] = []
    seed_sequence = np.random.SeedSequence(seed).spawn(count)
    for index, child in enumerate(seed_sequence):
        hidden_count = index % (h + 1)
        family = family_tuple[index % len(family_tuple)]
        child_seed = int(child.generate_state(1, dtype=np.uint32)[0])
        samples.append(
            generate_completion(
                q,
                hidden_count,
                family,
                child_seed,
                weight_min=weight_min,
                weight_max=weight_max,
                graph_options=graph_options,
            )
        )
    return samples


def full_error_operator(c_matrix: Array, q_matrix: Array, sample: CompletionSample) -> Array:
    c_value = np.asarray(c_matrix, dtype=np.float64)
    q_value = np.asarray(q_matrix, dtype=np.float64)
    if c_value.shape[1] != sample.q or q_value.shape != c_value.shape:
        raise ValueError("C/Q dimensions do not match the completion ports")
    target = c_value @ sample.smoother[: sample.q, :]
    local = np.zeros_like(target)
    local[:, : sample.q] = q_value
    return target - local


def completion_evaluation(
    c_matrix: Array,
    q_matrix: Array,
    samples: Iterable[CompletionSample],
    *,
    input_draws: int,
    seed: int,
) -> dict[str, Any]:
    sample_tuple = tuple(samples)
    if input_draws < 1 or not sample_tuple:
        raise ValueError("input_draws and the completion sample set must be nonempty")
    rng = np.random.default_rng(seed)
    operator_errors: list[float] = []
    expected_mse: list[float] = []
    empirical_squared_errors: list[float] = []
    records: list[dict[str, Any]] = []
    for sample in sample_tuple:
        error = full_error_operator(c_matrix, q_matrix, sample)
        radius = float(np.linalg.norm(error, ord=2))
        per_output_expected = float(np.sum(np.square(error)) / error.shape[0])
        inputs = rng.standard_normal((sample.nodes, input_draws))
        empirical = float(np.mean(np.square(error @ inputs)))
        operator_errors.append(radius)
        expected_mse.append(per_output_expected)
        empirical_squared_errors.append(empirical)
        records.append(
            {
                "family": sample.metadata["family"],
                "seed": sample.metadata["seed"],
                "hidden_count": sample.hidden_count,
                "nodes": sample.nodes,
                "cycle_rank": sample.metadata["cycle_rank"],
                "port_cut_edges": sample.metadata["port_cut_edges"],
                "port_cut_weight": sample.metadata["port_cut_weight"],
                "operator_2_norm_error": radius,
                "isotropic_expected_mse_per_output": per_output_expected,
                "empirical_mse_per_output": empirical,
            }
        )
    worst_index = int(np.argmax(operator_errors))
    return {
        "completion_count": len(sample_tuple),
        "input_draws_per_completion": input_draws,
        "maximum_operator_2_norm_error": max(operator_errors),
        "mean_operator_2_norm_error": float(np.mean(operator_errors)),
        "mean_isotropic_expected_mse_per_output": float(np.mean(expected_mse)),
        "mean_empirical_mse_per_output": float(np.mean(empirical_squared_errors)),
        "worst_completion": records[worst_index],
        "by_completion": records,
    }


def _solve_sampled_operators(
    c_matrix: Array,
    samples: Iterable[CompletionSample],
    *,
    local_rule: LocalRuleSpec,
    solver: str | None,
    solver_options: dict[str, Any] | None,
    verbose: bool = False,
) -> dict[str, Any]:
    sample_tuple = tuple(samples)
    if not sample_tuple:
        raise ValueError("at least one training completion is required")
    c_value = np.asarray(c_matrix, dtype=np.float64)
    outputs, q = c_value.shape
    if any(sample.q != q for sample in sample_tuple):
        raise ValueError("all training completions must use C's port count")
    chosen_solver = _choose_solver(solver)
    q_variable = cp.Variable((outputs, q), name="Q_sampled")
    radius = cp.Variable(nonneg=True, name="sampled_radius")
    constraints, _ = local_rule.cvxpy_constraints(q_variable)
    identity_outputs = np.eye(outputs)
    started = perf_counter()
    for sample in sample_tuple:
        target = c_value @ sample.smoother[:q, :]
        selector = np.zeros((q, sample.nodes))
        selector[:, :q] = np.eye(q)
        error = target - q_variable @ selector
        constraints.append(
            cp.bmat(
                [
                    [radius * identity_outputs, error],
                    [error.T, radius * np.eye(sample.nodes)],
                ]
            )
            >> 0
        )
    modeled_at = perf_counter()
    problem = cp.Problem(cp.Minimize(radius), constraints)
    try:
        problem.solve(solver=chosen_solver, verbose=verbose, **dict(solver_options or {}))
    except cp.error.SolverError as error:
        raise ExactSolverError(f"sampled SDP failed: {error}") from error
    solved_at = perf_counter()
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE} or q_variable.value is None:
        raise ExactSolverError(f"sampled SDP ended with status {problem.status!r}")
    q_value = np.asarray(q_variable.value, dtype=np.float64)
    risks = [float(np.linalg.norm(full_error_operator(c_value, q_value, sample), ord=2)) for sample in sample_tuple]
    maximum = max(risks)
    threshold = 2e-6 * max(1.0, maximum)
    active = [
        {
            "training_index": index,
            "family": sample.metadata["family"],
            "seed": sample.metadata["seed"],
            "hidden_count": sample.hidden_count,
            "operator_2_norm_error": risk,
        }
        for index, (sample, risk) in enumerate(zip(sample_tuple, risks, strict=True))
        if maximum - risk <= threshold
    ]
    return {
        "Q": q_value.tolist(),
        "training_radius_recomputed": maximum,
        "optimization_radius": max(0.0, float(radius.value)),
        "status": str(problem.status),
        "solver": chosen_solver,
        "active_training_completions": active,
        "timing_seconds": {
            "model_construction": modeled_at - started,
            "solver": solved_at - modeled_at,
            "total": solved_at - started,
        },
    }


def _solve_projection_subset(
    c_matrix: Array,
    projections: list[tuple[str, Array]],
    *,
    local_rule: LocalRuleSpec,
    solver: str | None,
    solver_options: dict[str, Any] | None,
) -> dict[str, Any]:
    c_value = np.asarray(c_matrix, dtype=np.float64)
    outputs, q = c_value.shape
    chosen_solver = _choose_solver(solver)
    q_variable = cp.Variable((outputs, q), name="Q_projection_subset")
    radius = cp.Variable(nonneg=True, name="subset_radius")
    constraints, _ = local_rule.cvxpy_constraints(q_variable)
    started = perf_counter()
    for _, projection in projections:
        error = c_value @ projection - q_variable
        constraints.append(
            cp.bmat(
                [
                    [radius * np.eye(outputs), error],
                    [error.T, radius * np.eye(q)],
                ]
            )
            >> 0
        )
    modeled_at = perf_counter()
    problem = cp.Problem(cp.Minimize(radius), constraints)
    problem.solve(solver=chosen_solver, **dict(solver_options or {}))
    solved_at = perf_counter()
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE} or q_variable.value is None:
        raise ExactSolverError(f"projection-subset SDP ended with status {problem.status!r}")
    q_value = np.asarray(q_variable.value, dtype=np.float64)
    risks = [float(np.linalg.norm(c_value @ projection - q_value, ord=2)) for _, projection in projections]
    maximum = max(risks)
    return {
        "Q": q_value.tolist(),
        "training_radius_recomputed": maximum,
        "optimization_radius": max(0.0, float(radius.value)),
        "status": str(problem.status),
        "solver": chosen_solver,
        "checked_projections": [label for label, _ in projections],
        "timing_seconds": {
            "model_construction": modeled_at - started,
            "solver": solved_at - modeled_at,
            "total": solved_at - started,
        },
    }


def solve_psd_contraction_relaxation(
    c_matrix: Array,
    *,
    local_rule: LocalRuleSpec | dict[str, Any] | None = None,
    solver: str | None = None,
    solver_options: dict[str, Any] | None = None,
    projection_audit_samples: int = 256,
    seed: int = 0,
) -> dict[str, Any]:
    """Solve the exact spectral relaxation over every ``0 <= X <= I``.

    With ``X=(I+Delta)/2``, testing the spectral norm only uses the image of
    one vector.  Symmetric and unstructured contractions have the same image
    ball, so the single-full-block Petersen LMI below is lossless.
    """

    c_value = np.asarray(c_matrix, dtype=np.float64)
    if c_value.ndim != 2 or not np.all(np.isfinite(c_value)):
        raise ValueError("C must be a finite matrix")
    outputs, q = c_value.shape
    rule = local_rule if isinstance(local_rule, LocalRuleSpec) else LocalRuleSpec.from_dict(local_rule)
    chosen_solver = _choose_solver(solver)
    q_variable = cp.Variable((outputs, q), name="Q_contraction")
    radius = cp.Variable(nonneg=True, name="contraction_radius")
    multiplier = cp.Variable(nonneg=True, name="petersen_multiplier")
    constraints, _ = rule.cvxpy_constraints(q_variable)
    gain = 0.5 * c_value
    nominal = 0.5 * c_value - q_variable
    zero_output_port = np.zeros((outputs, q))
    lmi = cp.bmat(
        [
            [
                radius * np.eye(outputs) - multiplier * (gain @ gain.T),
                nominal,
                zero_output_port,
            ],
            [nominal.T, radius * np.eye(q), np.eye(q)],
            [zero_output_port.T, np.eye(q), multiplier * np.eye(q)],
        ]
    )
    constraints.append(lmi >> 0)
    started = perf_counter()
    problem = cp.Problem(cp.Minimize(radius), constraints)
    modeled_at = perf_counter()
    try:
        problem.solve(solver=chosen_solver, **dict(solver_options or {}))
    except cp.error.SolverError as error:
        raise ExactSolverError(f"contraction-relaxation SDP failed: {error}") from error
    solved_at = perf_counter()
    if (
        problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
        or q_variable.value is None
        or radius.value is None
        or multiplier.value is None
    ):
        raise ExactSolverError(f"contraction SDP ended with status {problem.status!r}")
    q_value = np.asarray(q_variable.value, dtype=np.float64)
    radius_value = max(0.0, float(radius.value))
    multiplier_value = max(0.0, float(multiplier.value))
    concrete_lmi = np.block(
        [
            [
                radius_value * np.eye(outputs) - multiplier_value * (gain @ gain.T),
                0.5 * c_value - q_value,
                zero_output_port,
            ],
            [
                (0.5 * c_value - q_value).T,
                radius_value * np.eye(q),
                np.eye(q),
            ],
            [zero_output_port.T, np.eye(q), multiplier_value * np.eye(q)],
        ]
    )
    rng = np.random.default_rng(seed)
    audited: list[tuple[str, Array]] = [("rank_0", np.zeros((q, q))), ("rank_q", np.eye(q))]
    for index in range(projection_audit_samples):
        rank = 1 + index % max(1, q - 1) if q > 1 else 1
        basis, _ = np.linalg.qr(rng.standard_normal((q, rank)))
        audited.append((f"random_rank_{rank}_{index}", basis @ basis.T))
    audit_values = [float(np.linalg.norm(c_value @ projection - q_value, ord=2)) for _, projection in audited]
    worst_index = int(np.argmax(audit_values))
    return {
        "Q": q_value.tolist(),
        "certified_radius": radius_value,
        "status": str(problem.status),
        "solver": chosen_solver,
        "petersen_multiplier": multiplier_value,
        "certificate_min_eigenvalue": float(np.linalg.eigvalsh((concrete_lmi + concrete_lmi.T) / 2.0)[0]),
        "random_projection_audit": {
            "projection_count": len(audited),
            "maximum_observed_radius": max(audit_values),
            "slack_to_certificate": radius_value - max(audit_values),
            "worst_projection_label": audited[worst_index][0],
        },
        "timing_seconds": {
            "model_construction": modeled_at - started,
            "solver": solved_at - modeled_at,
            "postsolve_audit": perf_counter() - solved_at,
            "total": perf_counter() - started,
        },
        "semantics": "lossless spectral robust design over the outer set 0 <= X <= I",
    }


def scenario_certificate(
    c_matrix: Array,
    q_matrix: Array,
    scenarios: Iterable[Scenario],
    *,
    active_tolerance: float = 2e-6,
) -> dict[str, Any]:
    scenario_tuple = tuple(scenarios)
    risks = [spectral_squared_risk(c_matrix, q_matrix, scenario.matrix) for scenario in scenario_tuple]
    maximum = max(risks)
    threshold = active_tolerance * max(1.0, maximum)
    active = [
        {
            "label": scenario.label,
            "radius": math.sqrt(max(0.0, risk)),
            "squared_risk": risk,
        }
        for scenario, risk in zip(scenario_tuple, risks, strict=True)
        if maximum - risk <= threshold
    ]
    return {
        "scenario_count": len(scenario_tuple),
        "certified_radius": math.sqrt(maximum),
        "certified_squared_radius": maximum,
        "active_worst_case_scenarios": active,
    }


def _independent_parameter_count(rule: LocalRuleSpec, rows: int, columns: int) -> int:
    summary = rule.summary(rows, columns)
    count = int(summary["free_entries_before_equalities"])
    count -= sum(max(0, len(group) - 1) for group in rule.data.get("shared_coefficients", []))
    count -= len(rule.data.get("fixed_entries", []))
    count -= len(rule.data.get("linear_equalities", []))
    if "row_sums" in rule.data:
        count -= rows
    return max(0, count)


def _online_resources(q_matrix: Array, rule: LocalRuleSpec, scalar_bytes: int) -> dict[str, Any] | None:
    hop = rule.data.get("r_hop")
    if not isinstance(hop, dict):
        return None
    return direct_rule_resources(
        q_matrix,
        np.asarray(hop["adjacency"], dtype=np.float64),
        hop["output_nodes"],
        scalar_bytes=scalar_bytes,
        independent_coefficients=_independent_parameter_count(
            rule, q_matrix.shape[0], q_matrix.shape[1]
        ),
    )


def _attach_method_evaluation(
    method: dict[str, Any],
    c_matrix: Array,
    h: int,
    training: list[CompletionSample],
    test: list[CompletionSample],
    *,
    input_draws: int,
    seed: int,
    rule: LocalRuleSpec,
    scalar_bytes: int,
) -> dict[str, Any]:
    answer = dict(method)
    q_value = np.asarray(answer["Q"], dtype=np.float64)
    finite = scenario_certificate(c_matrix, q_value, finite_scenarios(c_matrix.shape[1], h))
    unbounded = scenario_certificate(c_matrix, q_value, partial_partition_scenarios(c_matrix.shape[1]))
    answer["exact_finite_h_evaluation"] = finite
    answer["exact_unbounded_evaluation"] = unbounded
    answer["training_completion_evaluation"] = completion_evaluation(
        c_matrix, q_value, training, input_draws=input_draws, seed=seed
    )
    answer["held_out_completion_evaluation"] = completion_evaluation(
        c_matrix, q_value, test, input_draws=input_draws, seed=seed + 1
    )
    training_radius = answer.get("training_radius_recomputed")
    if training_radius is not None:
        answer["finite_certificate_gap_over_training"] = max(
            0.0, finite["certified_radius"] - float(training_radius)
        )
    answer["online_resources"] = _online_resources(q_value, rule, scalar_bytes)
    return answer


def _timed_exact_solve(**kwargs: Any) -> tuple[dict[str, Any], Array]:
    with PeakRssSampler() as memory:
        result = solve_exact_spectral(**kwargs)
    record = result.to_dict()
    record["memory"] = memory.to_dict()
    return record, result.q_matrix


def _budget_sweep(
    c_matrix: Array,
    h_values: list[int],
    unbounded_q: Array,
    unbounded_radius: float,
    *,
    local_rule: LocalRuleSpec,
    solver: str | None,
    solver_options: dict[str, Any] | None,
    solve_optima: bool,
) -> list[dict[str, Any]]:
    q = c_matrix.shape[1]
    coefficient = float(
        np.linalg.norm(c_matrix, ord=2) ** 2
        + 2.0 * np.linalg.norm(c_matrix, ord=2) * np.linalg.norm(unbounded_q, ord=2)
    )
    rows: list[dict[str, Any]] = []
    for h in h_values:
        finite_certificate = scenario_certificate(c_matrix, unbounded_q, finite_scenarios(q, h))
        squared_gap = max(
            0.0,
            unbounded_radius**2 - float(finite_certificate["certified_squared_radius"]),
        )
        row: dict[str, Any] = {
            "h": h,
            "candidate_scenarios": scenario_count(q, h),
            "retained_vertices": vertex_count(q, h),
            "exact_hausdorff_distance": q / (q + h),
            "fixed_unbounded_Q_finite_radius": finite_certificate["certified_radius"],
            "fixed_unbounded_Q_unbounded_radius": unbounded_radius,
            "fixed_Q_squared_risk_gap": squared_gap,
            "fixed_Q_squared_gap_divided_by_hausdorff": squared_gap / (q / (q + h)),
            "theorem_squared_risk_gap_upper_bound": coefficient * q / (q + h),
        }
        if solve_optima:
            exact, _ = _timed_exact_solve(
                c_matrix=c_matrix,
                h=h,
                mode="finite",
                local_rule=local_rule,
                prune_vertices=True,
                solver=solver,
                solver_options=solver_options,
            )
            row["optimized_finite_radius"] = exact["radius"]
            row["optimization_seconds"] = exact["timing_seconds"]["total"]
        rows.append(row)
    return rows


def _scenario_scaling(config: dict[str, Any], solver: str | None, solver_options: dict[str, Any] | None) -> list[dict[str, Any]]:
    q_values = [int(value) for value in config.get("q_values", range(2, 9))]
    h = int(config.get("h", 2))
    solve_max_q = int(config.get("solve_max_q", 0))
    rows: list[dict[str, Any]] = []
    for q in q_values:
        tracemalloc.start()
        started = perf_counter()
        candidates = finite_scenarios(q, h)
        generated_candidates_at = perf_counter()
        vertices = [scenario for scenario in candidates if scenario.is_vertex]
        partial = partial_partition_scenarios(q)
        finished = perf_counter()
        _, python_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        row: dict[str, Any] = {
            "q": q,
            "h": h,
            "finite_candidates": len(candidates),
            "finite_vertices": len(vertices),
            "pruned_scenarios": len(candidates) - len(vertices),
            "unbounded_partial_partitions": len(partial),
            "finite_generation_seconds": generated_candidates_at - started,
            "unbounded_generation_seconds": finished - generated_candidates_at,
            "python_peak_tracemalloc_bytes": int(python_peak),
        }
        if q <= solve_max_q:
            c_value = np.diag(np.linspace(1.0, 2.0, q))
            zeros = [[row_index, column] for row_index in range(q) for column in range(q) if row_index != column]
            shared = [[[index, index] for index in range(q)]] if q > 1 else []
            rule = LocalRuleSpec.from_dict(
                {"zeros": zeros, "shared_coefficients": shared}
            )
            finite_result, _ = _timed_exact_solve(
                c_matrix=c_value,
                h=h,
                mode="finite",
                local_rule=rule,
                solver=solver,
                solver_options=solver_options,
            )
            unbounded_result, _ = _timed_exact_solve(
                c_matrix=c_value,
                mode="unbounded",
                local_rule=rule,
                solver=solver,
                solver_options=solver_options,
            )
            row["finite_sdp"] = {
                "radius": finite_result["radius"],
                "timing_seconds": finite_result["timing_seconds"],
                "memory": finite_result["memory"],
                "active": finite_result["scenario_counts"]["active"],
            }
            row["unbounded_sdp"] = {
                "radius": unbounded_result["radius"],
                "timing_seconds": unbounded_result["timing_seconds"],
                "memory": unbounded_result["memory"],
                "active": unbounded_result["scenario_counts"]["active"],
            }
        else:
            row["sdp_skipped_reason"] = (
                f"q={q} exceeds configured solve_max_q={solve_max_q}; generation is still measured"
            )
        rows.append(row)
        del candidates, vertices, partial
    return rows


def run_robust_benchmark(config: dict[str, Any]) -> dict[str, Any]:
    """Run exact and baseline robust designs under one frozen semantic contract."""

    q = int(config["q"])
    h = int(config["h"])
    c_matrix = np.asarray(config["C"], dtype=np.float64)
    if c_matrix.ndim != 2 or c_matrix.shape[1] != q:
        raise ValueError("C must be a matrix with q columns")
    rule = LocalRuleSpec.from_dict(config.get("local_rule"))
    seed = int(config.get("seed", 20260914))
    scalar_bytes = int(config.get("scalar_bytes", 8))
    solver_config = dict(config.get("solver", {}))
    solver = solver_config.get("name")
    solver_options = solver_config.get("options")
    sampling = dict(config.get("random_sampling", {}))
    families = sampling.get(
        "families", ["random_cyclic", "grid", "geometric", "small_world"]
    )
    train = generate_completion_samples(
        q,
        h,
        int(sampling.get("training_count", 16)),
        families=families,
        seed=seed + 10,
        weight_min=float(sampling.get("training_weight_min", 0.1)),
        weight_max=float(sampling.get("training_weight_max", 10.0)),
        graph_options=sampling.get("graph_options"),
    )
    test = generate_completion_samples(
        q,
        h,
        int(sampling.get("test_count", 32)),
        families=families,
        seed=seed + 20,
        weight_min=float(sampling.get("test_weight_min", 0.01)),
        weight_max=float(sampling.get("test_weight_max", 100.0)),
        graph_options=sampling.get("graph_options"),
    )
    input_draws = int(sampling.get("input_draws", 32))

    finite_record, finite_q = _timed_exact_solve(
        c_matrix=c_matrix,
        q=q,
        h=h,
        mode="finite",
        local_rule=rule,
        prune_vertices=True,
        solver=solver,
        solver_options=solver_options,
    )
    unbounded_record, unbounded_q = _timed_exact_solve(
        c_matrix=c_matrix,
        q=q,
        mode="unbounded",
        local_rule=rule,
        solver=solver,
        solver_options=solver_options,
    )
    with PeakRssSampler() as sampled_memory:
        sampled_record = _solve_sampled_operators(
            c_matrix,
            train,
            local_rule=rule,
            solver=solver,
            solver_options=solver_options,
        )
    sampled_record["memory"] = sampled_memory.to_dict()
    with PeakRssSampler() as relaxation_memory:
        relaxation_record = solve_psd_contraction_relaxation(
            c_matrix,
            local_rule=rule,
            solver=solver,
            solver_options=solver_options,
            projection_audit_samples=int(config.get("projection_audit_samples", 256)),
            seed=seed + 30,
        )
    relaxation_record["memory"] = relaxation_memory.to_dict()
    with PeakRssSampler() as endpoint_memory:
        endpoint_record = _solve_projection_subset(
            c_matrix,
            [("zero", np.zeros((q, q))), ("identity", np.eye(q))],
            local_rule=rule,
            solver=solver,
            solver_options=solver_options,
        )
    endpoint_record["memory"] = endpoint_memory.to_dict()

    methods: dict[str, dict[str, Any]] = {}
    for name, base in (
        ("exact_finite_h", finite_record),
        ("exact_unbounded", unbounded_record),
        ("random_completion_sampling", sampled_record),
        ("psd_contraction_relaxation", relaxation_record),
        ("endpoint_only_no_go", endpoint_record),
    ):
        methods[name] = _attach_method_evaluation(
            base,
            c_matrix,
            h,
            train,
            test,
            input_draws=input_draws,
            seed=seed + 100,
            rule=rule,
            scalar_bytes=scalar_bytes,
        )

    pruning_ablation: dict[str, Any] | None = None
    if bool(config.get("run_pruning_ablation", True)):
        unpruned_record, _ = _timed_exact_solve(
            c_matrix=c_matrix,
            q=q,
            h=h,
            mode="finite",
            local_rule=rule,
            prune_vertices=False,
            solver=solver,
            solver_options=solver_options,
        )
        pruning_ablation = {
            "pruned_radius": finite_record["radius"],
            "unpruned_radius": unpruned_record["radius"],
            "absolute_radius_difference": abs(
                float(finite_record["radius"]) - float(unpruned_record["radius"])
            ),
            "pruned_lmis": finite_record["scenario_counts"]["retained_lmis"],
            "unpruned_lmis": unpruned_record["scenario_counts"]["retained_lmis"],
            "pruned_total_seconds": finite_record["timing_seconds"]["total"],
            "unpruned_total_seconds": unpruned_record["timing_seconds"]["total"],
            "lmi_reduction_fraction": 1.0
            - finite_record["scenario_counts"]["retained_lmis"]
            / unpruned_record["scenario_counts"]["retained_lmis"],
        }

    budget_config = dict(config.get("budget_sweep", {}))
    h_values = [int(value) for value in budget_config.get("h_values", [0, 1, h, 2 * h + 1])]
    budget_c = np.asarray(budget_config.get("C", c_matrix), dtype=np.float64)
    if budget_c.ndim != 2 or budget_c.shape[1] < 1:
        raise ValueError("budget_sweep C must be a nonempty matrix")
    budget_rule = (
        LocalRuleSpec.from_dict(budget_config.get("local_rule"))
        if "C" in budget_config or "local_rule" in budget_config
        else rule
    )
    if budget_c.shape == c_matrix.shape and np.array_equal(budget_c, c_matrix) and budget_rule.data == rule.data:
        budget_unbounded_q = unbounded_q
        budget_unbounded_radius = float(unbounded_record["radius"])
        budget_unbounded_design_seconds = float(unbounded_record["timing_seconds"]["total"])
    else:
        budget_unbounded_record, budget_unbounded_q = _timed_exact_solve(
            c_matrix=budget_c,
            mode="unbounded",
            local_rule=budget_rule,
            solver=solver,
            solver_options=solver_options,
        )
        budget_unbounded_radius = float(budget_unbounded_record["radius"])
        budget_unbounded_design_seconds = float(
            budget_unbounded_record["timing_seconds"]["total"]
        )
    budget = _budget_sweep(
        budget_c,
        sorted(set(h_values)),
        budget_unbounded_q,
        budget_unbounded_radius,
        local_rule=budget_rule,
        solver=solver,
        solver_options=solver_options,
        solve_optima=bool(budget_config.get("solve_optima", True)),
    )

    result = {
        "schema_version": 1,
        "classification": "theorem_matching",
        "semantics": {
            "target": "C E_Gamma (I+L_G)^-1 y",
            "local_rule": "Q E_Gamma y, common to every completion",
            "certified_error": "full-input spectral operator norm",
            "empirical_mse": "graph-signal error for isotropic synthetic y; not WLS/Kalman MSE",
        },
        "problem": {
            "q": q,
            "h": h,
            "C": c_matrix.tolist(),
            "local_rule": rule.data,
            "unit_grounding": True,
            "undirected_nonnegative_edges": True,
        },
        "random_completion_distribution": {
            "families": list(families),
            "training_count": len(train),
            "test_count": len(test),
            "training_weight_range": [
                float(sampling.get("training_weight_min", 0.1)),
                float(sampling.get("training_weight_max", 10.0)),
            ],
            "test_weight_range": [
                float(sampling.get("test_weight_min", 0.01)),
                float(sampling.get("test_weight_max", 100.0)),
            ],
            "training_hidden_count_histogram": {
                str(value): sum(sample.hidden_count == value for sample in train)
                for value in range(h + 1)
            },
            "test_hidden_count_histogram": {
                str(value): sum(sample.hidden_count == value for sample in test)
                for value in range(h + 1)
            },
        },
        "methods": methods,
        "vertex_pruning_ablation": pruning_ablation,
        "budget_sweep": {
            "problem": {
                "q": budget_c.shape[1],
                "C": budget_c.tolist(),
                "local_rule": budget_rule.data,
                "fixed_unbounded_Q": budget_unbounded_q.tolist(),
                "fixed_unbounded_radius": budget_unbounded_radius,
                "unbounded_design_seconds": budget_unbounded_design_seconds,
            },
            "rows": budget,
        },
        "checks": {
            "pruned_unpruned_radius_agree": pruning_ablation is None
            or pruning_ablation["absolute_radius_difference"] <= 2e-5,
            "finite_exact_training_never_exceeds_certificate": methods["exact_finite_h"][
                "training_completion_evaluation"
            ]["maximum_operator_2_norm_error"]
            <= methods["exact_finite_h"]["exact_finite_h_evaluation"]["certified_radius"]
            + 2e-6,
            "finite_exact_test_never_exceeds_certificate": methods["exact_finite_h"][
                "held_out_completion_evaluation"
            ]["maximum_operator_2_norm_error"]
            <= methods["exact_finite_h"]["exact_finite_h_evaluation"]["certified_radius"]
            + 2e-6,
            "contraction_outer_radius_not_below_exact_unbounded_for_same_Q": methods[
                "psd_contraction_relaxation"
            ]["certified_radius"]
            + 2e-5
            >= methods["psd_contraction_relaxation"]["exact_unbounded_evaluation"][
                "certified_radius"
            ],
        },
    }
    scaling_config = config.get("scenario_scaling")
    if isinstance(scaling_config, dict):
        result["scenario_scaling"] = _scenario_scaling(
            scaling_config, solver, solver_options
        )
    if not all(result["checks"].values()):
        raise RuntimeError(f"robust benchmark checks failed: {result['checks']}")
    return result
