"""Lossless spectral-norm SDP for exact hidden-completion scenarios."""

from __future__ import annotations

from dataclasses import dataclass
import math
from time import perf_counter
from typing import Any, Literal

import cvxpy as cp
import numpy as np
from numpy.typing import NDArray

from .local_rules import LocalRuleSpec
from .metrics import spectral_squared_risk
from .scenarios import Scenario, finite_scenarios, partial_partition_scenarios


Array = NDArray[np.float64]


class ExactSolverError(RuntimeError):
    """Raised when the exact conic formulation cannot be solved."""


@dataclass(slots=True)
class ExactSolveResult:
    mode: Literal["finite", "unbounded"]
    q: int
    hidden_budget: int | None
    q_matrix: Array
    radius: float
    squared_radius: float
    optimization_squared_radius: float
    status: str
    solver: str
    total_scenario_count: int
    retained_scenario_count: int
    active_scenarios: list[dict[str, object]]
    local_rule_summary: dict[str, object]
    timing_seconds: dict[str, float]
    certificate_min_eigenvalue: float

    def to_dict(self) -> dict[str, object]:
        return {
            "mode": self.mode,
            "q": self.q,
            "hidden_budget": self.hidden_budget,
            "Q": self.q_matrix.tolist(),
            "radius": self.radius,
            "squared_radius": self.squared_radius,
            "optimization_squared_radius": self.optimization_squared_radius,
            "status": self.status,
            "solver": self.solver,
            "scenario_counts": {
                "total_exact_candidates": self.total_scenario_count,
                "retained_lmis": self.retained_scenario_count,
                "pruned": self.total_scenario_count - self.retained_scenario_count,
                "active": len(self.active_scenarios),
            },
            "active_worst_case_scenarios": self.active_scenarios,
            "local_rule": self.local_rule_summary,
            "timing_seconds": self.timing_seconds,
            "certificate_min_eigenvalue": self.certificate_min_eigenvalue,
        }


def _choose_solver(requested: str | None) -> str:
    installed = set(cp.installed_solvers())
    if requested is not None:
        normalized = requested.upper()
        if normalized not in installed:
            raise ExactSolverError(
                f"requested solver {normalized!r} is unavailable; installed solvers: {sorted(installed)}"
            )
        return normalized
    for candidate in ("CLARABEL", "SCS", "CVXOPT", "MOSEK"):
        if candidate in installed:
            return candidate
    raise ExactSolverError(
        "no SDP-capable solver found; install cvxpy with CLARABEL, SCS, CVXOPT, or MOSEK"
    )


def _validate_c(c_matrix: Array, q: int | None) -> tuple[Array, int, int]:
    matrix = np.asarray(c_matrix, dtype=float)
    if matrix.ndim != 2 or not np.all(np.isfinite(matrix)):
        raise ValueError("C must be a finite two-dimensional matrix")
    outputs, inferred_q = matrix.shape
    if outputs < 1 or inferred_q < 1:
        raise ValueError("C must have at least one row and one column")
    if q is not None:
        if not isinstance(q, int) or isinstance(q, bool) or q < 1:
            raise ValueError("q must be a positive integer when supplied")
        if q != inferred_q:
            raise ValueError(f"q={q} does not match C's {inferred_q} columns")
    return matrix, outputs, inferred_q


def _scenario_record(scenario: Scenario, squared_risk: float, maximum: float) -> dict[str, object]:
    record = scenario.to_dict(include_matrix=True)
    record["squared_risk"] = squared_risk
    record["radius"] = math.sqrt(max(0.0, squared_risk))
    record["squared_risk_slack_from_maximum"] = maximum - squared_risk
    return record


def solve_exact_spectral(
    c_matrix: Array,
    *,
    q: int | None = None,
    h: int | None = None,
    mode: Literal["finite", "unbounded"] = "finite",
    local_rule: LocalRuleSpec | dict[str, Any] | None = None,
    prune_vertices: bool = True,
    solver: str | None = None,
    solver_options: dict[str, Any] | None = None,
    active_tolerance: float = 2e-6,
    verbose: bool = False,
) -> ExactSolveResult:
    """Solve the exact common-local-rule spectral minimax problem.

    The conic formulation is lossless.  Returned floating-point values inherit
    the selected SDP solver's numerical tolerance; ``squared_radius`` is
    recomputed directly from the returned Q over every exact candidate.
    """

    c_matrix, outputs, q = _validate_c(c_matrix, q)
    if mode not in {"finite", "unbounded"}:
        raise ValueError("mode must be 'finite' or 'unbounded'")
    if mode == "finite":
        if not isinstance(h, int) or isinstance(h, bool) or h < 0:
            raise ValueError("finite mode requires a nonnegative integer h")
    elif h is not None:
        raise ValueError("h must be omitted in unbounded mode")
    if active_tolerance < 0:
        raise ValueError("active_tolerance must be nonnegative")

    rule = local_rule if isinstance(local_rule, LocalRuleSpec) else LocalRuleSpec.from_dict(local_rule)
    chosen_solver = _choose_solver(solver)
    solver_options = dict(solver_options or {})

    start = perf_counter()
    if mode == "finite":
        all_scenarios = finite_scenarios(q, h, vertices_only=False)
        retained = (
            [scenario for scenario in all_scenarios if scenario.is_vertex]
            if prune_vertices
            else all_scenarios
        )
    else:
        all_scenarios = partial_partition_scenarios(q)
        retained = all_scenarios
    generated_at = perf_counter()

    q_variable = cp.Variable((outputs, q), name="Q")
    rho = cp.Variable(nonneg=True, name="rho")
    constraints, _ = rule.cvxpy_constraints(q_variable)
    identity_outputs = np.eye(outputs)
    identity_ports = np.eye(q)
    for scenario in retained:
        x_matrix = scenario.matrix
        top_left = (
            rho * identity_outputs
            - c_matrix @ x_matrix @ c_matrix.T
            + c_matrix @ x_matrix @ q_variable.T
            + q_variable @ x_matrix @ c_matrix.T
        )
        top_left = (top_left + top_left.T) / 2.0
        lmi = cp.bmat([[top_left, q_variable], [q_variable.T, identity_ports]])
        constraints.append(lmi >> 0)
    problem = cp.Problem(cp.Minimize(rho), constraints)
    modeled_at = perf_counter()

    try:
        problem.solve(solver=chosen_solver, verbose=verbose, **solver_options)
    except cp.error.SolverError as error:
        raise ExactSolverError(f"{chosen_solver} failed to solve the SDP: {error}") from error
    solved_at = perf_counter()

    accepted_statuses = {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
    if problem.status not in accepted_statuses or q_variable.value is None or rho.value is None:
        raise ExactSolverError(f"SDP ended with status {problem.status!r}")

    q_value = np.asarray(q_variable.value, dtype=float)
    risks = [spectral_squared_risk(c_matrix, q_value, scenario.matrix) for scenario in all_scenarios]
    verified_squared_radius = max(risks)
    optimization_squared_radius = max(0.0, float(rho.value))
    threshold = active_tolerance * max(1.0, verified_squared_radius)
    active = [
        _scenario_record(scenario, risk, verified_squared_radius)
        for scenario, risk in zip(all_scenarios, risks, strict=True)
        if verified_squared_radius - risk <= threshold
    ]

    # Re-evaluate every retained LMI at the verified risk.  The Schur
    # complement implies nonnegativity; the minimum eigenvalue reports only
    # floating-point residual, not an additional relaxation.
    minimum_eigenvalue = math.inf
    for scenario in retained:
        x_matrix = scenario.matrix
        top_left = (
            verified_squared_radius * identity_outputs
            - c_matrix @ x_matrix @ c_matrix.T
            + c_matrix @ x_matrix @ q_value.T
            + q_value @ x_matrix @ c_matrix.T
        )
        block = np.block([[top_left, q_value], [q_value.T, identity_ports]])
        minimum_eigenvalue = min(
            minimum_eigenvalue,
            float(np.linalg.eigvalsh((block + block.T) / 2.0)[0]),
        )
    evaluated_at = perf_counter()

    return ExactSolveResult(
        mode=mode,
        q=q,
        hidden_budget=h,
        q_matrix=q_value,
        radius=math.sqrt(verified_squared_radius),
        squared_radius=verified_squared_radius,
        optimization_squared_radius=optimization_squared_radius,
        status=str(problem.status),
        solver=chosen_solver,
        total_scenario_count=len(all_scenarios),
        retained_scenario_count=len(retained),
        active_scenarios=active,
        local_rule_summary=rule.summary(outputs, q),
        timing_seconds={
            "scenario_generation": generated_at - start,
            "model_construction": modeled_at - generated_at,
            "solver": solved_at - modeled_at,
            "postsolve_verification": evaluated_at - solved_at,
            "total": evaluated_at - start,
        },
        certificate_min_eigenvalue=minimum_eigenvalue,
    )
