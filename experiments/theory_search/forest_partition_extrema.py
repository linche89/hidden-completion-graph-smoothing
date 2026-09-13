"""Numerical checks for the grounded-network partition extremum theorem.

The script is deliberately small and dependency-light.  It checks:

1. the matrix-forest convex-mixture identity on a weighted triangle;
2. convergence of high-conductance path components to block-averaging
   projections;
3. an optimized common-rule counterexample in which an intermediate cut of a
   three-port path is active and the two naive endpoint completions are not
   sufficient; and
4. a strict gap between grounded-network completions and the enclosing
   self-adjoint contraction/Petersen uncertainty set;
5. the general L-ensemble mean-projection identity for arbitrary factor
   columns, of which the graph-forest formula is a special case;
6. an explicit two-sparse raw-WLS realization of a prescribed Schur block;
7. the obstruction to a Euclidean-input hidden-node extension with
   heterogeneous grounding; and
8. the strict loss of the flat/projection converse when edge weights have a
   finite upper bound.
"""

from __future__ import annotations

import itertools
from fractions import Fraction

import numpy as np


def laplacian(n: int, edges: list[tuple[int, int, float]]) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=float)
    for left, right, weight in edges:
        matrix[left, left] += weight
        matrix[right, right] += weight
        matrix[left, right] -= weight
        matrix[right, left] -= weight
    return matrix


def components(n: int, edges: list[tuple[int, int, float]]) -> list[list[int]]:
    parent = list(range(n))

    def find(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left: int, right: int) -> bool:
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            return False
        parent[root_right] = root_left
        return True

    for left, right, _ in edges:
        if not union(left, right):
            raise ValueError("edge set is not a forest")

    groups: dict[int, list[int]] = {}
    for node in range(n):
        groups.setdefault(find(node), []).append(node)
    return list(groups.values())


def averaging_projection(n: int, blocks: list[list[int]]) -> np.ndarray:
    projection = np.zeros((n, n), dtype=float)
    for block in blocks:
        projection[np.ix_(block, block)] = 1.0 / len(block)
    return projection


def check_matrix_forest_identity() -> None:
    edges = [(0, 1, 0.7), (1, 2, 1.2), (0, 2, 0.4)]
    target = np.linalg.inv(np.eye(3) + laplacian(3, edges))

    numerator = np.zeros((3, 3), dtype=float)
    denominator = 0.0
    for mask in range(1 << len(edges)):
        chosen = [edges[index] for index in range(len(edges)) if mask >> index & 1]
        try:
            blocks = components(3, chosen)
        except ValueError:
            continue
        edge_weight = np.prod([edge[2] for edge in chosen], dtype=float)
        rooted_multiplicity = np.prod([len(block) for block in blocks], dtype=float)
        total_weight = edge_weight * rooted_multiplicity
        numerator += total_weight * averaging_projection(3, blocks)
        denominator += total_weight

    mixture = numerator / denominator
    np.testing.assert_allclose(mixture, target, rtol=2e-15, atol=2e-15)


def check_high_conductance_limit() -> None:
    # Components {0, 2} and {1, 3, 4}; the latter is a path.
    expected = averaging_projection(5, [[0, 2], [1, 3, 4]])
    edges = [(0, 2, 1.0e8), (1, 3, 1.0e8), (3, 4, 1.0e8)]
    actual = np.linalg.inv(np.eye(5) + laplacian(5, edges))
    assert np.linalg.norm(actual - expected, ord=2) < 1.0e-7


def connected_path_projections() -> list[np.ndarray]:
    # Host path 0--2--1.  Its connected partitions are obtained by deleting
    # any subset of the two edges.
    return [
        averaging_projection(3, [[0], [1], [2]]),
        averaging_projection(3, [[0, 2], [1]]),
        averaging_projection(3, [[1, 2], [0]]),
        averaging_projection(3, [[0, 1, 2]]),
    ]


def squared_row_loss(
    projection: np.ndarray,
    c_row: np.ndarray,
    f_matrix: np.ndarray,
    q_row: np.ndarray,
) -> float:
    visible = c_row @ projection @ f_matrix - q_row
    # projection is symmetric and idempotent, hence
    # ||c P||^2 = c P c^T.
    invisible_squared = c_row @ projection @ c_row
    return float(visible @ visible + invisible_squared)


def check_optimized_intermediate_cut_counterexample() -> None:
    c_row = np.array([-1.0, -1.0, 2.0])
    f_matrix = np.array(
        [
            [-1.0, 3.0, 3.0],
            [-1.0, 2.0, -1.0],
            [2.0, 3.0, -3.0],
        ]
    )
    # C is the first column of F transposed, as required by
    # C = R_i F^T with R_i selecting the first local coordinate.
    np.testing.assert_array_equal(c_row, f_matrix[:, 0])

    q_star = np.array(
        [
            float(Fraction(317, 84)),
            float(Fraction(20, 21)),
            -float(Fraction(157, 42)),
        ]
    )
    projections = connected_path_projections()
    losses = np.array(
        [squared_row_loss(item, c_row, f_matrix, q_star) for item in projections]
    )
    expected = np.array(
        [
            float(Fraction(9785, 336)),
            float(Fraction(9785, 336)),
            float(Fraction(3485, 336)),
            float(Fraction(9785, 336)),
        ]
    )
    np.testing.assert_allclose(losses, expected, rtol=2e-15, atol=2e-15)

    # Positive KKT weights for the three active scenarios prove optimality of
    # q_star for the full four-scenario minimax problem.
    active_centers = np.array(
        [[6.0, 1.0, -8.0], [1.5, 1.0, 1.0], [0.0, 0.0, 0.0]]
    )
    multipliers = np.array([197.0 / 378.0, 163.0 / 378.0, 1.0 / 21.0])
    np.testing.assert_allclose(
        multipliers @ (q_star[None, :] - active_centers),
        np.zeros(3),
        rtol=0.0,
        atol=2e-15,
    )
    np.testing.assert_allclose(multipliers.sum(), 1.0, rtol=0.0, atol=2e-15)

    # If only the all-disconnected and all-fused endpoints are checked, the
    # exact two-scenario center and value are smaller.
    endpoint_q = (107.0 / 202.0) * np.array([6.0, 1.0, -8.0])
    endpoint_value = float(Fraction(11449, 404))
    endpoint_losses = [
        squared_row_loss(projections[index], c_row, f_matrix, endpoint_q)
        for index in (0, 3)
    ]
    np.testing.assert_allclose(
        endpoint_losses, [endpoint_value, endpoint_value], rtol=2e-15, atol=2e-15
    )
    assert float(Fraction(9785, 336)) > endpoint_value


def check_petersen_relaxation_gap() -> None:
    # Network extreme set for two ports with hidden-node dilution allowed.
    zero = np.zeros((2, 2))
    e_11 = np.diag([1.0, 0.0])
    e_22 = np.diag([0.0, 1.0])
    identity = np.eye(2)
    fused = np.full((2, 2), 0.5)
    network_extremes = [zero, e_11, e_22, identity, fused]

    # This version also obeys the WLS compatibility C=R_i F^T: take
    # F=diag(1,sqrt(5)) and R_i=[1,0].
    c_row = np.array([1.0, 0.0])
    f_column = np.diag([1.0, np.sqrt(5.0)])
    q_row = np.zeros(2)
    network_squared = max(
        squared_row_loss(item, c_row, f_column, q_row)
        for item in network_extremes
    )
    np.testing.assert_allclose(network_squared, 2.0)

    # A rank-one projection in the enclosing interval 0 <= T <= I attains the
    # larger Petersen/full-block value 9/4.
    full_block_projection = np.array(
        [[3.0, np.sqrt(3.0)], [np.sqrt(3.0), 1.0]]
    ) / 4.0
    np.testing.assert_allclose(
        full_block_projection @ full_block_projection,
        full_block_projection,
        rtol=2e-15,
        atol=2e-15,
    )
    full_block_squared = squared_row_loss(
        full_block_projection, c_row, f_column, q_row
    )
    np.testing.assert_allclose(full_block_squared, 9.0 / 4.0)
    assert full_block_squared > network_squared


def kernel_projection(columns: np.ndarray) -> np.ndarray:
    """Orthogonal projector onto ker(columns.T)."""
    dimension = columns.shape[0]
    if columns.shape[1] == 0:
        return np.eye(dimension)
    gram = columns.T @ columns
    return np.eye(dimension) - columns @ np.linalg.inv(gram) @ columns.T


def check_general_dpp_projection_identity() -> None:
    # General (not incidence) two-local factors are covered by the same
    # identity.  Dependent subsets have zero L-ensemble mass and can be
    # skipped safely.
    factor_columns = np.array(
        [
            [1.0, -0.2, 0.0, 0.7, -0.4],
            [0.3, 1.1, -0.8, 0.0, 0.6],
            [0.0, 0.5, 1.2, -0.9, 0.2],
        ]
    )
    weights = np.array([0.4, 1.3, 0.8, 2.1, 0.6])
    target = np.linalg.inv(
        np.eye(3) + factor_columns @ np.diag(weights) @ factor_columns.T
    )

    numerator = np.zeros((3, 3))
    denominator = 0.0
    for size in range(4):
        for subset in itertools.combinations(range(5), size):
            selected = factor_columns[:, subset]
            if size == 0:
                gram_determinant = 1.0
            else:
                gram_determinant = float(np.linalg.det(selected.T @ selected))
                if gram_determinant < 1.0e-12:
                    continue
            mass = gram_determinant * float(np.prod(weights[list(subset)]))
            numerator += mass * kernel_projection(selected)
            denominator += mass

    np.testing.assert_allclose(
        denominator,
        np.linalg.det(
            np.eye(3) + factor_columns @ np.diag(weights) @ factor_columns.T
        ),
        rtol=3e-14,
        atol=3e-14,
    )
    np.testing.assert_allclose(numerator / denominator, target, rtol=3e-14, atol=3e-14)


def check_raw_pairwise_wls_realization() -> None:
    # Prescribe F, choose A=aI and E=aF^T, then realize every entry of E by
    # a two-coordinate measurement row.  A sufficiently small a leaves a
    # diagonally-dominant exterior residual, which is itself a sum of unary
    # and two-coordinate squares.
    f_matrix = np.array(
        [
            [-1.0, 3.0, 3.0],
            [-1.0, 2.0, -1.0],
            [2.0, 3.0, -3.0],
        ]
    )
    local_dimension, port_dimension = f_matrix.shape[1], f_matrix.shape[0]
    mu = 1.0
    a_scale = 1.0e-3
    a_block = a_scale * np.eye(local_dimension)
    e_block = a_scale * f_matrix.T

    factor_rows: list[np.ndarray] = []
    local_diagonal_used = np.zeros(local_dimension)
    exterior_diagonal_used = np.zeros(port_dimension)
    for local in range(local_dimension):
        support = np.flatnonzero(f_matrix[:, local])
        share = a_scale / (2.0 * max(1, len(support)))
        for port in support:
            left = np.sqrt(share)
            cross = e_block[local, port]
            right = cross / left
            row = np.zeros(local_dimension + port_dimension)
            row[local], row[local_dimension + port] = left, right
            factor_rows.append(row)
            local_diagonal_used[local] += left * left
            exterior_diagonal_used[port] += right * right
        residual = a_scale - local_diagonal_used[local]
        if residual > 0.0:
            row = np.zeros(local_dimension + port_dimension)
            row[local] = np.sqrt(residual)
            factor_rows.append(row)

    # The fixed exterior target before adding completion factors is
    # mu*I + E^T A^{-1} E.  Subtract what the cross-factor rows already used
    # and decompose the remaining SDD matrix into two-coordinate squares.
    target_exterior = mu * np.eye(port_dimension) + e_block.T @ np.linalg.solve(
        a_block, e_block
    )
    residual_exterior = target_exterior - np.diag(exterior_diagonal_used)
    assert np.all(
        np.diag(residual_exterior)
        >= np.sum(np.abs(residual_exterior), axis=1)
        - np.abs(np.diag(residual_exterior))
        - 1.0e-14
    )
    used = np.zeros(port_dimension)
    for left_port in range(port_dimension):
        for right_port in range(left_port + 1, port_dimension):
            entry = residual_exterior[left_port, right_port]
            if abs(entry) < 1.0e-15:
                continue
            magnitude = abs(entry)
            row = np.zeros(local_dimension + port_dimension)
            row[local_dimension + left_port] = np.sqrt(magnitude)
            row[local_dimension + right_port] = np.sign(entry) * np.sqrt(magnitude)
            factor_rows.append(row)
            used[left_port] += magnitude
            used[right_port] += magnitude
    for port in range(port_dimension):
        residual = residual_exterior[port, port] - used[port]
        assert residual >= -1.0e-14
        if residual > 0.0:
            row = np.zeros(local_dimension + port_dimension)
            row[local_dimension + port] = np.sqrt(residual)
            factor_rows.append(row)

    measurement_matrix = np.vstack(factor_rows)
    information = measurement_matrix.T @ measurement_matrix
    np.testing.assert_allclose(information[:local_dimension, :local_dimension], a_block)
    np.testing.assert_allclose(information[:local_dimension, local_dimension:], e_block)
    np.testing.assert_allclose(
        information[local_dimension:, local_dimension:], target_exterior
    )
    schur = information[local_dimension:, local_dimension:] - (
        information[local_dimension:, :local_dimension]
        @ np.linalg.solve(a_block, information[:local_dimension, local_dimension:])
    )
    np.testing.assert_allclose(schur, mu * np.eye(port_dimension), atol=2.0e-14)
    np.testing.assert_allclose(
        information[local_dimension:, :local_dimension] @ np.linalg.inv(a_block),
        f_matrix,
    )


def check_raw_measurement_metric_orthogonalization() -> None:
    # The nontrivial score-space loss must not be confused with the raw
    # whitened-measurement metric.  For b=H^T xi and J=H^T H, the latter
    # exactly orthogonalizes into Q A Q^T + C Sigma^{-1} C^T.
    rng = np.random.default_rng(17092026)
    local_dimension, port_dimension, output_dimension = 4, 3, 2
    a_seed = rng.normal(size=(local_dimension, local_dimension))
    a_block = a_seed @ a_seed.T + 0.8 * np.eye(local_dimension)
    e_block = rng.normal(size=(local_dimension, port_dimension))
    sigma_seed = rng.normal(size=(port_dimension, port_dimension))
    sigma = sigma_seed @ sigma_seed.T + 0.6 * np.eye(port_dimension)
    d_block = sigma + e_block.T @ np.linalg.solve(a_block, e_block)
    information = np.block([[a_block, e_block], [e_block.T, d_block]])

    f_matrix = e_block.T @ np.linalg.inv(a_block)
    c_matrix = rng.normal(size=(output_dimension, port_dimension))
    q_matrix = rng.normal(size=(output_dimension, local_dimension))
    transfer = c_matrix @ np.linalg.inv(sigma)
    score_error = np.hstack(
        (transfer @ f_matrix - q_matrix, -transfer)
    )
    expected_bias = np.hstack(
        (-q_matrix @ a_block, -q_matrix @ e_block - c_matrix)
    )
    np.testing.assert_allclose(
        score_error @ information,
        expected_bias,
        rtol=2.0e-13,
        atol=2.0e-13,
    )
    raw_covariance = score_error @ information @ score_error.T
    expected = (
        q_matrix @ a_block @ q_matrix.T
        + c_matrix @ np.linalg.solve(sigma, c_matrix.T)
    )
    np.testing.assert_allclose(raw_covariance, expected, rtol=2.0e-13, atol=2.0e-13)

    whitening_design = np.linalg.cholesky(information).T
    raw_operator = score_error @ whitening_design.T
    np.testing.assert_allclose(
        np.linalg.norm(raw_operator, ord=2) ** 2,
        np.linalg.eigvalsh(expected)[-1],
        rtol=3.0e-13,
        atol=3.0e-13,
    )
    q_penalty = q_matrix @ a_block @ q_matrix.T
    assert np.linalg.eigvalsh(q_penalty)[0] >= -1.0e-12


def check_heterogeneous_hidden_euclidean_obstruction() -> None:
    # WLS-compatible one-port example: C=F=R_i=1 and Q=1/4.  One hidden node
    # with lower grounding creates a raw-Euclidean score loss larger than
    # both the inactive and no-hidden endpoint scenarios.
    conductance = 1.0e9
    precision = np.diag([3.0, 1.0]) + laplacian(
        2, [(0, 1, conductance)]
    )
    inverse_row = np.linalg.inv(precision)[0, :]
    q_value = 0.25
    hidden_squared = (inverse_row[0] - q_value) ** 2 + np.dot(
        inverse_row, inverse_row
    )
    inactive_squared = q_value**2
    no_hidden_squared = (1.0 / 3.0 - q_value) ** 2 + (1.0 / 3.0) ** 2
    np.testing.assert_allclose(hidden_squared, 1.0 / 8.0, rtol=3.0e-8)
    np.testing.assert_allclose(no_hidden_squared, 17.0 / 144.0)
    assert hidden_squared > max(inactive_squared, no_hidden_squared)


def check_bounded_weight_projection_gap() -> None:
    # q=1, one factor of weight w in [0,1], C=F=1 and Q=2.  The attainable
    # normalized inverse is T in [1/2,1].  Its exact worst squared row loss is
    # 5/2, whereas adding the unattainable infinite-weight projection T=0
    # gives the strictly conservative value 4.
    def loss(normalized_inverse: float) -> float:
        return (normalized_inverse - 2.0) ** 2 + normalized_inverse**2

    bounded_exact = max(loss(1.0), loss(0.5))
    projection_outer_bound = max(loss(1.0), loss(0.0))
    np.testing.assert_allclose(bounded_exact, 2.5)
    np.testing.assert_allclose(projection_outer_bound, 4.0)
    assert projection_outer_bound > bounded_exact


def main() -> None:
    checks = [
        check_matrix_forest_identity,
        check_high_conductance_limit,
        check_optimized_intermediate_cut_counterexample,
        check_petersen_relaxation_gap,
        check_general_dpp_projection_identity,
        check_raw_pairwise_wls_realization,
        check_raw_measurement_metric_orthogonalization,
        check_heterogeneous_hidden_euclidean_obstruction,
        check_bounded_weight_projection_gap,
    ]
    for check in checks:
        check()
        print(f"PASS {check.__name__}")


if __name__ == "__main__":
    main()
