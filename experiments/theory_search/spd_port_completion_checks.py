"""Checks for the arbitrary-SPD port completion theorem candidate.

The theorem concerns

    T = H^(1/2) (H + L)^(-1) H^(1/2),
    H = diag(Lambda, D_hidden),

where Lambda is an arbitrary positive-definite port block and L is a graph
Laplacian on ports plus an arbitrary finite hidden graph.  These small tests
target the delicate points that disappear when Lambda is diagonal:

* graph components are no longer orthogonal after whitening;
* a second DPP/volume-sampling mixture is needed at the port boundary;
* the lower witness must kill inactive ports while retaining the full hidden
  input columns; and
* the full-input Schatten-p reduction has the sharp threshold p >= 2.

The tests are numerical checks, not proofs.
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.optimize import linprog, minimize_scalar


def symmetric_sqrt(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    if eigenvalues.min() <= 0.0:
        raise ValueError("matrix must be positive definite")
    return (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T


def set_partitions(items: tuple[int, ...]) -> list[tuple[tuple[int, ...], ...]]:
    """Enumerate unlabeled set partitions in a canonical order."""

    if not items:
        return [tuple()]
    first, rest = items[0], items[1:]
    result: list[tuple[tuple[int, ...], ...]] = []
    for partition in set_partitions(rest):
        result.append(((first,),) + partition)
        for block_index in range(len(partition)):
            blocks = [tuple(block) for block in partition]
            blocks[block_index] = tuple(sorted((first,) + blocks[block_index]))
            result.append(tuple(sorted(blocks)))
    return list(dict.fromkeys(result))


def spd_partial_partition_projections(port_base: np.ndarray) -> list[np.ndarray]:
    """Return P = Lambda^(1/2) U (U' Lambda U)^-1 U' Lambda^(1/2)."""

    port_count = port_base.shape[0]
    root = symmetric_sqrt(port_base)
    ports = tuple(range(port_count))
    projections: list[np.ndarray] = []
    for active_count in range(port_count + 1):
        for active in itertools.combinations(ports, active_count):
            for partition in set_partitions(active):
                if not partition:
                    projections.append(np.zeros_like(port_base))
                    continue
                indicator = np.zeros((port_count, len(partition)))
                for column, block in enumerate(partition):
                    indicator[list(block), column] = 1.0
                vectors = root @ indicator
                projection = vectors @ np.linalg.solve(
                    vectors.T @ vectors, vectors.T
                )
                projections.append((projection + projection.T) / 2.0)
    return projections


def laplacian(
    node_count: int, edges: list[tuple[int, int, float]]
) -> np.ndarray:
    matrix = np.zeros((node_count, node_count))
    for left, right, weight in edges:
        difference = np.zeros(node_count)
        difference[left] = 1.0
        difference[right] = -1.0
        matrix += weight * np.outer(difference, difference)
    return matrix


def normalized_smoother(
    port_base: np.ndarray,
    hidden_base: np.ndarray,
    edges: list[tuple[int, int, float]],
) -> np.ndarray:
    port_count = port_base.shape[0]
    hidden_count = hidden_base.size
    base = np.zeros((port_count + hidden_count, port_count + hidden_count))
    base[:port_count, :port_count] = port_base
    if hidden_count:
        base[port_count:, port_count:] = np.diag(hidden_base)
    root = symmetric_sqrt(base)
    graph_laplacian = laplacian(base.shape[0], edges)
    smoother = root @ np.linalg.solve(base + graph_laplacian, root)
    return (smoother + smoother.T) / 2.0


def column_subset_projection_mixture(
    vectors: np.ndarray, hidden_masses: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """DPP mean of projections onto subsets of the columns of vectors."""

    scaled = vectors / np.sqrt(hidden_masses)[None, :]
    subsets = []
    raw_weights = []
    projections = []
    for subset_size in range(vectors.shape[1] + 1):
        for subset in itertools.combinations(range(vectors.shape[1]), subset_size):
            subsets.append(subset)
            if not subset:
                raw_weights.append(1.0)
                projections.append(np.zeros((vectors.shape[0], vectors.shape[0])))
                continue
            indices = np.asarray(subset, dtype=int)
            selected_scaled = scaled[:, indices]
            raw_weights.append(float(np.linalg.det(selected_scaled.T @ selected_scaled)))
            selected = vectors[:, indices]
            projections.append(
                selected @ np.linalg.solve(selected.T @ selected, selected.T)
            )
    del subsets
    probabilities = np.maximum(np.asarray(raw_weights), 0.0)
    probabilities /= probabilities.sum()
    mean = sum(
        probability * projection
        for probability, projection in zip(probabilities, projections)
    )
    return mean, probabilities


def schatten_norm(matrix: np.ndarray, exponent: float) -> float:
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    if np.isinf(exponent):
        return float(singular_values.max(initial=0.0))
    return float(np.sum(singular_values**exponent) ** (1.0 / exponent))


def full_error(
    smoother: np.ndarray,
    port_count: int,
    output: np.ndarray,
    local_rule: np.ndarray,
) -> np.ndarray:
    hidden_count = smoother.shape[0] - port_count
    padded_rule = np.hstack(
        (local_rule, np.zeros((local_rule.shape[0], hidden_count)))
    )
    return output @ smoother[:port_count, :] - padded_rule


def in_convex_hull(matrix: np.ndarray, vertices: list[np.ndarray]) -> bool:
    coordinates = np.column_stack([vertex.reshape(-1) for vertex in vertices])
    constraints = np.vstack((coordinates, np.ones((1, len(vertices)))))
    target = np.concatenate((matrix.reshape(-1), [1.0]))
    result = linprog(
        np.zeros(len(vertices)),
        A_eq=constraints,
        b_eq=target,
        bounds=(0.0, None),
        method="highs",
        options={"dual_feasibility_tolerance": 1.0e-9},
    )
    return bool(result.success and np.linalg.norm(constraints @ result.x - target) < 1.0e-7)


def check_scenario_count_and_projection_geometry() -> None:
    port_base = np.array(
        [[2.0, 0.55, -0.20], [0.55, 1.7, 0.35], [-0.20, 0.35, 1.4]]
    )
    scenarios = spd_partial_partition_projections(port_base)
    assert len(scenarios) == 15  # Bell number B_4
    for projection in scenarios:
        np.testing.assert_allclose(projection, projection.T, atol=2.0e-13)
        np.testing.assert_allclose(projection @ projection, projection, atol=8.0e-13)
    assert min(
        np.linalg.norm(left - right)
        for index, left in enumerate(scenarios)
        for right in scenarios[index + 1 :]
    ) > 1.0e-4
    assert any(np.linalg.norm(item) < 1.0e-13 for item in scenarios)
    assert any(np.linalg.norm(item - np.eye(3)) < 1.0e-12 for item in scenarios)


def check_nonorthogonal_component_dpp_mixture() -> None:
    port_base = np.array(
        [
            [2.3, 0.45, -0.25, 0.15],
            [0.45, 1.8, 0.30, -0.10],
            [-0.25, 0.30, 1.6, 0.35],
            [0.15, -0.10, 0.35, 1.5],
        ]
    )
    root = symmetric_sqrt(port_base)
    indicator = np.zeros((4, 3))
    indicator[[0, 2], 0] = 1.0
    indicator[1, 1] = 1.0
    indicator[3, 2] = 1.0
    vectors = root @ indicator
    gram = vectors.T @ vectors
    assert np.max(np.abs(gram - np.diag(np.diag(gram)))) > 0.1

    hidden_masses = np.array([0.7, 2.0, 1.1])
    compression = vectors @ np.linalg.solve(
        gram + np.diag(hidden_masses), vectors.T
    )
    mixture, probabilities = column_subset_projection_mixture(
        vectors, hidden_masses
    )
    np.testing.assert_allclose(probabilities.sum(), 1.0, atol=2.0e-14)
    np.testing.assert_allclose(mixture, compression, rtol=2.0e-11, atol=2.0e-12)

    # Components with no hidden vertices correspond to zero h_j.  The fixed
    # finite projection hull is closed, so h_j -> 0 is the correct limiting
    # interpretation of the same DPP formula.
    semidefinite_masses = np.array([0.0, 1.3, 0.0])
    target = vectors @ np.linalg.solve(
        gram + np.diag(semidefinite_masses), vectors.T
    )
    errors = []
    for epsilon in (1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8):
        regularized = np.maximum(semidefinite_masses, epsilon)
        mixture, _ = column_subset_projection_mixture(vectors, regularized)
        errors.append(np.linalg.norm(mixture - target, ord=2))
    assert all(right < left for left, right in zip(errors, errors[1:]))
    assert errors[-1] < 2.0e-8


def check_random_spd_completion_bounds() -> None:
    rng = np.random.default_rng(20260914)
    port_count = 3
    raw = rng.normal(size=(port_count, port_count))
    port_base = raw @ raw.T + 0.8 * np.eye(port_count)
    hidden_base = np.array([0.55, 1.8])
    scenarios = spd_partial_partition_projections(port_base)
    output = rng.normal(size=(2, port_count))
    local_rule = rng.normal(size=(2, port_count))

    for _ in range(40):
        edges = []
        for left, right in itertools.combinations(range(5), 2):
            if rng.random() < 0.55:
                edges.append((left, right, float(np.exp(rng.normal()))))
        smoother = normalized_smoother(port_base, hidden_base, edges)
        compression = smoother[:port_count, :port_count]
        assert in_convex_hull(compression, scenarios)

        error = full_error(smoother, port_count, output, local_rule)
        for exponent in (2.0, 3.0, 4.0, np.inf):
            actual = schatten_norm(error, exponent)
            scenario_bound = max(
                schatten_norm(output @ projection - local_rule, exponent)
                for projection in scenarios
            )
            assert actual <= scenario_bound + 2.0e-8


def check_degree_two_lower_witness() -> None:
    port_base = np.array(
        [[2.0, 0.65, -0.15], [0.65, 1.5, 0.25], [-0.15, 0.25, 1.3]]
    )
    root = symmetric_sqrt(port_base)
    active_indicator = np.array([[1.0], [0.0], [1.0]])
    vector = root @ active_indicator
    target = vector @ np.linalg.solve(vector.T @ vector, vector.T)

    errors = []
    for scale in (10.0, 100.0, 1000.0, 10000.0):
        # Active block {0,2} is fused by one edge.  Inactive port 1 is fused
        # to one increasingly grounded hidden vertex.  This is a path forest
        # of maximum degree one (hence certainly at most two).
        hidden_base = np.array([scale])
        conductance = scale**3
        edges = [(0, 2, conductance), (1, 3, conductance)]
        smoother = normalized_smoother(port_base, hidden_base, edges)
        padded_target = np.hstack((target, np.zeros((3, 1))))
        errors.append(np.linalg.norm(smoother[:3, :] - padded_target, ord=2))
    assert all(right < left for left, right in zip(errors, errors[1:]))
    assert errors[-1] < 0.012


def check_schatten_threshold_is_sharp() -> None:
    # q=1, one hidden node, unit bases, and an infinitely conductive edge.
    # The limiting forest projection is 11'/2.  Take orthonormal endpoint
    # residuals Q=e1 and C-Q=e2, i.e. C=(1,1)' and Q=(1,0)'.
    full_projection = np.ones((2, 2)) / 2.0
    output = np.array([[1.0], [1.0]])
    local_rule = np.array([[1.0], [0.0]])
    limit_error = output @ full_projection[:1, :] - np.hstack(
        (local_rule, np.zeros((2, 1)))
    )
    np.testing.assert_allclose(
        np.linalg.svd(limit_error, compute_uv=False),
        np.array([1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0)]),
        atol=2.0e-14,
    )
    endpoint_value = 1.0
    for exponent in (1.0, 1.2, 1.5, 1.9):
        actual = schatten_norm(limit_error, exponent)
        expected = 2.0 ** (1.0 / exponent - 0.5)
        np.testing.assert_allclose(actual, expected, atol=2.0e-14)
        assert actual > endpoint_value + 1.0e-3
    np.testing.assert_allclose(schatten_norm(limit_error, 2.0), 1.0)
    assert schatten_norm(limit_error, 3.0) < 1.0
    assert schatten_norm(limit_error, np.inf) < 1.0


def check_structured_opcode_needs_intermediate_partition() -> None:
    # With a shared opcode Q=alpha I, the connected two-port projection is
    # genuinely active.  Omitting it gives alpha=1 and radius 1; including
    # all five scenarios gives the exact larger values below.
    port_base = np.eye(2)
    scenarios = spd_partial_partition_projections(port_base)
    output = np.diag([1.0, 2.0])

    objective = lambda alpha: max(
        np.linalg.norm(output @ projection - alpha * np.eye(2), ord=2)
        for projection in scenarios
    )
    optimum = minimize_scalar(
        objective, bounds=(-0.5, 2.5), method="bounded", options={"xatol": 1.0e-13}
    )
    exact_alpha = 4.0 / 3.0 - 2.0 * np.sqrt(10.0) / 15.0
    exact_radius = (10.0 + 2.0 * np.sqrt(10.0)) / 15.0
    np.testing.assert_allclose(optimum.x, exact_alpha, atol=3.0e-8)
    np.testing.assert_allclose(optimum.fun, exact_radius, atol=3.0e-8)
    assert exact_radius > 1.08

    connected = np.ones((2, 2)) / 2.0
    without_connected = [
        projection
        for projection in scenarios
        if np.linalg.norm(projection - connected) > 1.0e-12
    ]
    endpoint_objective = lambda alpha: max(
        np.linalg.norm(output @ projection - alpha * np.eye(2), ord=2)
        for projection in without_connected
    )
    endpoint_optimum = minimize_scalar(
        endpoint_objective,
        bounds=(-0.5, 2.5),
        method="bounded",
        options={"xatol": 1.0e-13},
    )
    np.testing.assert_allclose(endpoint_optimum.x, 1.0, atol=3.0e-8)
    np.testing.assert_allclose(endpoint_optimum.fun, 1.0, atol=3.0e-8)


def main() -> None:
    checks = [
        check_scenario_count_and_projection_geometry,
        check_nonorthogonal_component_dpp_mixture,
        check_random_spd_completion_bounds,
        check_degree_two_lower_witness,
        check_schatten_threshold_is_sharp,
        check_structured_opcode_needs_intermediate_partition,
    ]
    for check in checks:
        check()
        print(f"PASS {check.__name__}")


if __name__ == "__main__":
    main()
