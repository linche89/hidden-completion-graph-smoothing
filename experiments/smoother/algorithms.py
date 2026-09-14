"""Fixed-known-graph graph-smoother baselines."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Any

import numpy as np
from numpy.typing import NDArray
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .graphs import WeightedGraph, smoother_matrix
from .resources import ball_nodes


Array = NDArray[np.float64]


def _matrix_rhs(values: Array) -> tuple[Array, bool]:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1:
        return array[:, None], True
    if array.ndim != 2:
        raise ValueError("right-hand side must be a vector or a two-dimensional array")
    return array, False


def _restore(values: Array, was_vector: bool) -> Array:
    return values[:, 0] if was_vector else values


@dataclass(slots=True)
class GlobalSolve:
    solution: Array
    metadata: dict[str, Any]


def global_smoother(
    system: sp.spmatrix,
    signals: Array,
    *,
    direct_max_n: int = 50_000,
    relative_tolerance: float = 1e-10,
    max_iterations: int | None = None,
) -> GlobalSolve:
    """Solve the exact reference system by sparse LU or tightly converged CG."""

    matrix = sp.csr_matrix(system, dtype=np.float64)
    rhs, was_vector = _matrix_rhs(signals)
    if rhs.shape[0] != matrix.shape[0]:
        raise ValueError("signal dimension and system size disagree")
    started = perf_counter()
    iterations: list[int] = []
    if matrix.shape[0] <= direct_max_n:
        factor_started = perf_counter()
        factor = spla.splu(matrix.tocsc())
        factor_seconds = perf_counter() - factor_started
        solve_started = perf_counter()
        answer = np.asarray(factor.solve(rhs), dtype=np.float64)
        solve_seconds = perf_counter() - solve_started
        method = "sparse_lu"
    else:
        factor_seconds = 0.0
        columns: list[Array] = []
        solve_started = perf_counter()
        diagonal = matrix.diagonal()
        if np.any(diagonal <= 0):
            raise ValueError("Jacobi preconditioner requires a positive diagonal")
        preconditioner = spla.LinearOperator(
            matrix.shape, matvec=lambda vector: vector / diagonal, dtype=np.float64
        )
        for column in rhs.T:
            counter = [0]

            def callback(_: Array) -> None:
                counter[0] += 1

            solved, info = spla.cg(
                matrix,
                column,
                rtol=relative_tolerance,
                atol=0.0,
                maxiter=max_iterations,
                M=preconditioner,
                callback=callback,
            )
            if info != 0:
                raise RuntimeError(f"CG reference solve failed with info={info}")
            columns.append(np.asarray(solved))
            iterations.append(counter[0])
        answer = np.column_stack(columns)
        solve_seconds = perf_counter() - solve_started
        method = "jacobi_pcg"
    residuals = np.linalg.norm(matrix @ answer - rhs, axis=0) / np.maximum(
        np.linalg.norm(rhs, axis=0), 1e-300
    )
    return GlobalSolve(
        _restore(answer, was_vector),
        {
            "method": method,
            "factorization_seconds": factor_seconds,
            "solve_seconds": solve_seconds,
            "total_seconds": perf_counter() - started,
            "relative_residual_max": float(np.max(residuals)),
            "relative_residual_mean": float(np.mean(residuals)),
            "iterations_per_rhs": iterations,
            "relative_tolerance": relative_tolerance,
        },
    )


def gershgorin_spectral_interval(graph: WeightedGraph, *, edge_scale: float = 1.0) -> tuple[float, float]:
    degree = np.asarray(graph.adjacency.sum(axis=1)).ravel()
    return 1.0, float(1.0 + 2.0 * edge_scale * np.max(degree, initial=0.0))


def neumann_apply(
    system: sp.spmatrix,
    signals: Array,
    rounds: int,
    spectral_upper: float,
) -> Array:
    """Degree-``rounds`` scaled Neumann polynomial, with exactly that many SpMVs."""

    if rounds < 0:
        raise ValueError("rounds must be nonnegative")
    matrix = sp.csr_matrix(system, dtype=np.float64)
    rhs, was_vector = _matrix_rhs(signals)
    if not np.isfinite(spectral_upper) or spectral_upper <= 0:
        raise ValueError("spectral_upper must be positive and finite")
    term = rhs / spectral_upper
    answer = term.copy()
    for _ in range(rounds):
        term = term - (matrix @ term) / spectral_upper
        answer += term
    return _restore(np.asarray(answer), was_vector)


def chebyshev_coefficients(degree: int, lower: float, upper: float) -> Array:
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    if not (0 < lower <= upper < np.inf):
        raise ValueError("spectral interval must satisfy 0 < lower <= upper")
    if lower == upper:
        coefficients = np.zeros(degree + 1)
        coefficients[0] = 1.0 / lower
        return coefficients
    midpoint = 0.5 * (lower + upper)
    halfwidth = 0.5 * (upper - lower)
    return np.asarray(
        np.polynomial.chebyshev.chebinterpolate(
            lambda point: 1.0 / (midpoint + halfwidth * point), degree
        ),
        dtype=np.float64,
    )


def chebyshev_apply(
    system: sp.spmatrix,
    signals: Array,
    degree: int,
    lower: float,
    upper: float,
) -> Array:
    """Apply the interpolating Chebyshev inverse polynomial using ``degree`` SpMVs."""

    coefficients = chebyshev_coefficients(degree, lower, upper)
    return chebyshev_apply_coefficients(system, signals, coefficients, lower, upper)


def chebyshev_apply_coefficients(
    system: sp.spmatrix,
    signals: Array,
    coefficients: Array,
    lower: float,
    upper: float,
) -> Array:
    """Apply precomputed coefficients, separating offline design from online work."""

    matrix = sp.csr_matrix(system, dtype=np.float64)
    rhs, was_vector = _matrix_rhs(signals)
    coefficients = np.asarray(coefficients, dtype=np.float64)
    if coefficients.ndim != 1 or coefficients.size < 1:
        raise ValueError("coefficients must be a nonempty vector")
    degree = coefficients.size - 1
    if lower == upper:
        return _restore(coefficients[0] * rhs, was_vector)
    midpoint = 0.5 * (lower + upper)
    halfwidth = 0.5 * (upper - lower)
    scaled = (matrix - midpoint * sp.eye(matrix.shape[0], format="csr")) / halfwidth
    previous = rhs.copy()
    answer = coefficients[0] * previous
    if degree >= 1:
        current = scaled @ rhs
        answer += coefficients[1] * current
        for order in range(2, degree + 1):
            following = 2.0 * (scaled @ current) - previous
            answer += coefficients[order] * following
            previous, current = current, following
    return _restore(np.asarray(answer), was_vector)


def induced_local_operator(adjacency: sp.spmatrix, radius: int) -> sp.csr_matrix:
    """Build the ordinary induced-subgraph smoother row at every root."""

    graph = sp.csr_matrix(adjacency, dtype=np.float64)
    n = graph.shape[0]
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for root in range(n):
        nodes = ball_nodes(graph, root, radius)
        local_adjacency = graph[nodes][:, nodes].tocsr()
        degree = np.asarray(local_adjacency.sum(axis=1)).ravel()
        local_system = sp.eye(nodes.size, format="csc") + sp.diags(degree) - local_adjacency
        location = int(np.flatnonzero(nodes == root)[0])
        selector = np.zeros(nodes.size)
        selector[location] = 1.0
        coefficients = np.asarray(spla.spsolve(local_system, selector)).ravel()
        nonzero = np.abs(coefficients) > 1e-14
        rows.extend([root] * int(np.count_nonzero(nonzero)))
        columns.extend(nodes[nonzero].tolist())
        values.extend(coefficients[nonzero].tolist())
    return sp.coo_matrix((values, (rows, columns)), shape=(n, n)).tocsr()


def apply_induced_local(operator: sp.spmatrix, signals: Array) -> Array:
    return np.asarray(sp.csr_matrix(operator) @ np.asarray(signals, dtype=np.float64))


def error_metrics(reference: Array, approximation: Array) -> dict[str, float]:
    exact, _ = _matrix_rhs(reference)
    estimate, _ = _matrix_rhs(approximation)
    if exact.shape != estimate.shape:
        raise ValueError("reference and approximation shapes disagree")
    difference = estimate - exact
    return {
        "empirical_mse_per_node": float(np.mean(np.square(difference))),
        "empirical_rmse_per_node": float(np.sqrt(np.mean(np.square(difference)))),
        "relative_frobenius_signal_error": float(
            np.linalg.norm(difference) / max(np.linalg.norm(exact), 1e-300)
        ),
        "maximum_absolute_node_error": float(np.max(np.abs(difference))),
    }


def fixed_operator_metrics(reference: Array, approximation: Array) -> dict[str, float]:
    exact = np.asarray(reference, dtype=np.float64)
    estimate = np.asarray(approximation, dtype=np.float64)
    difference = estimate - exact
    return {
        "fixed_instance_operator_2_norm": float(np.linalg.norm(difference, ord=2)),
        "isotropic_input_expected_mse_per_node": float(
            np.sum(np.square(difference)) / difference.shape[0]
        ),
        "operator_frobenius_norm": float(np.linalg.norm(difference, ord="fro")),
    }


def exact_smoother_operator(graph: WeightedGraph, *, edge_scale: float = 1.0) -> Array:
    matrix = smoother_matrix(graph, edge_scale=edge_scale)
    factor = spla.splu(matrix.tocsc())
    return np.asarray(factor.solve(np.eye(matrix.shape[0])), dtype=np.float64)
