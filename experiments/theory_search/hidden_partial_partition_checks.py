"""Numerical stress tests for the hidden-node partial-partition reduction.

These checks do not prove the theorem.  They target the two steps most likely
to conceal an algebraic mistake:

1. a finite graph with hidden inputs must be bounded by the finite family of
   port partial-partition scenarios; and
2. an inactive port can be diluted by a long, highly conducting hidden path
   without silently discarding the hidden-input columns of the error operator.
"""

from __future__ import annotations

import itertools

import numpy as np


def laplacian(node_count: int, edges: list[tuple[int, int, float]]) -> np.ndarray:
    matrix = np.zeros((node_count, node_count), dtype=float)
    for left, right, weight in edges:
        matrix[left, left] += weight
        matrix[right, right] += weight
        matrix[left, right] -= weight
        matrix[right, left] -= weight
    return matrix


def set_partitions(items: tuple[int, ...]) -> list[tuple[tuple[int, ...], ...]]:
    """Enumerate set partitions in a canonical order (small sets only)."""

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


def partial_partition_projections(port_count: int) -> list[np.ndarray]:
    projections: list[np.ndarray] = []
    ports = tuple(range(port_count))
    for active_count in range(port_count + 1):
        for active in itertools.combinations(ports, active_count):
            for partition in set_partitions(active):
                projection = np.zeros((port_count, port_count), dtype=float)
                for block in partition:
                    projection[np.ix_(block, block)] = 1.0 / len(block)
                projections.append(projection)
    return projections


def weighted_partial_partition_projections(
    port_weights: np.ndarray,
) -> list[np.ndarray]:
    """Partial-partition projections in the grounding-normalized geometry."""

    port_count = len(port_weights)
    projections: list[np.ndarray] = []
    ports = tuple(range(port_count))
    for active_count in range(port_count + 1):
        for active in itertools.combinations(ports, active_count):
            for partition in set_partitions(active):
                projection = np.zeros((port_count, port_count), dtype=float)
                for block in partition:
                    indices = np.asarray(block, dtype=int)
                    direction = np.sqrt(port_weights[indices])
                    projection[np.ix_(indices, indices)] = (
                        np.outer(direction, direction) / np.dot(direction, direction)
                    )
                projections.append(projection)
    return projections


def full_error_norm(
    inverse: np.ndarray,
    port_count: int,
    c_matrix: np.ndarray,
    f_matrix: np.ndarray,
    q_matrix: np.ndarray,
) -> float:
    hidden_count = inverse.shape[0] - port_count
    c_full = np.hstack((c_matrix, np.zeros((c_matrix.shape[0], hidden_count))))
    f_full = np.vstack((f_matrix, np.zeros((hidden_count, f_matrix.shape[1]))))
    error = np.hstack((c_full @ inverse @ f_full - q_matrix, -c_full @ inverse))
    return float(np.linalg.norm(error, ord=2))


def partial_scenario_norm(
    projection: np.ndarray,
    mu: float,
    c_matrix: np.ndarray,
    f_matrix: np.ndarray,
    q_matrix: np.ndarray,
) -> float:
    error = np.hstack(
        (
            (c_matrix @ projection @ f_matrix) / mu - q_matrix,
            -(c_matrix @ projection) / mu,
        )
    )
    return float(np.linalg.norm(error, ord=2))


def normalized_full_error_norm(
    normalized_inverse: np.ndarray,
    port_count: int,
    c_matrix: np.ndarray,
    f_matrix: np.ndarray,
    q_matrix: np.ndarray,
) -> float:
    """Error norm when the exterior variables use grounding whitening."""

    hidden_count = normalized_inverse.shape[0] - port_count
    c_full = np.hstack((c_matrix, np.zeros((c_matrix.shape[0], hidden_count))))
    f_full = np.vstack((f_matrix, np.zeros((hidden_count, f_matrix.shape[1]))))
    error = np.hstack(
        (
            c_full @ normalized_inverse @ f_full - q_matrix,
            -c_full @ normalized_inverse,
        )
    )
    return float(np.linalg.norm(error, ord=2))


def normalized_partial_scenario_norm(
    projection: np.ndarray,
    c_matrix: np.ndarray,
    f_matrix: np.ndarray,
    q_matrix: np.ndarray,
) -> float:
    error = np.hstack(
        (
            c_matrix @ projection @ f_matrix - q_matrix,
            -c_matrix @ projection,
        )
    )
    return float(np.linalg.norm(error, ord=2))


def check_random_hidden_graph_upper_bound() -> None:
    rng = np.random.default_rng(9132026)
    port_count, hidden_count = 3, 3
    node_count = port_count + hidden_count
    mu = 0.7
    c_matrix = rng.normal(size=(2, port_count))
    f_matrix = rng.normal(size=(port_count, 2))
    q_matrix = rng.normal(size=(2, 2))
    scenarios = partial_partition_projections(port_count)
    assert len(scenarios) == 15  # B_4
    scenario_maximum = max(
        partial_scenario_norm(item, mu, c_matrix, f_matrix, q_matrix)
        for item in scenarios
    )

    complete_edges = list(itertools.combinations(range(node_count), 2))
    for _ in range(100):
        selected = rng.random(len(complete_edges)) < 0.35
        weights = np.exp(rng.uniform(-3.0, 3.0, size=len(complete_edges)))
        edges = [
            (left, right, float(weight))
            for (left, right), keep, weight in zip(
                complete_edges, selected, weights, strict=True
            )
            if keep
        ]
        inverse = np.linalg.inv(mu * np.eye(node_count) + laplacian(node_count, edges))
        actual = full_error_norm(
            inverse, port_count, c_matrix, f_matrix, q_matrix
        )
        assert actual <= scenario_maximum + 2.0e-12


def diluted_path_inverse(
    port_count: int, inactive_port: int, hidden_count: int, mu: float
) -> np.ndarray:
    """Fuse ports 0 and 2; dilute port 1 down a hidden path."""

    assert port_count == 3 and inactive_port == 1
    node_count = port_count + hidden_count
    conductance = mu * hidden_count**3
    edges: list[tuple[int, int, float]] = [(0, 2, conductance)]
    previous = inactive_port
    for hidden in range(port_count, node_count):
        edges.append((previous, hidden, conductance))
        previous = hidden
    return np.linalg.inv(mu * np.eye(node_count) + laplacian(node_count, edges))


def check_degree_two_dilution_converse() -> None:
    rng = np.random.default_rng(14092026)
    port_count = 3
    mu = 1.3
    c_matrix = rng.normal(size=(2, port_count))
    f_matrix = rng.normal(size=(port_count, 2))
    q_matrix = rng.normal(size=(2, 2))

    target_projection = np.zeros((port_count, port_count), dtype=float)
    target_projection[np.ix_([0, 2], [0, 2])] = 0.5
    target = partial_scenario_norm(
        target_projection, mu, c_matrix, f_matrix, q_matrix
    )

    errors = []
    for hidden_count in (8, 16, 32, 64, 128):
        inverse = diluted_path_inverse(
            port_count, inactive_port=1, hidden_count=hidden_count, mu=mu
        )
        errors.append(
            abs(
                full_error_norm(
                    inverse, port_count, c_matrix, f_matrix, q_matrix
                )
                - target
            )
        )

    # The hidden-input columns remain in full_error_norm.  Their contribution
    # decays only as the fused component grows, so monotonicity is not assumed;
    # the final error must nevertheless be small and much lower than initially.
    assert errors[-1] < 0.12 * errors[0]
    assert errors[-1] < 1.0e-2


def check_mixed_state_raw_score_upper_bound() -> None:
    """Check the bias + raw-noise + score-disturbance corollary."""

    rng = np.random.default_rng(18092026)
    local_count, port_count, hidden_count, output_count = 2, 3, 3, 2
    node_count = port_count + hidden_count
    mu, score_radius = 0.8, 0.45

    a_seed = rng.normal(size=(local_count, local_count))
    a_block = a_seed @ a_seed.T + 0.7 * np.eye(local_count)
    e_block = rng.normal(size=(local_count, port_count))
    f_matrix = e_block.T @ np.linalg.inv(a_block)
    output_selector = rng.normal(size=(output_count, local_count))
    c_matrix = output_selector @ np.linalg.solve(a_block, e_block)
    q_matrix = rng.normal(size=(output_count, local_count))
    pi_seed = rng.normal(size=(local_count + port_count,) * 2)
    state_shape = pi_seed @ pi_seed.T + 0.2 * np.eye(local_count + port_count)
    bias = np.hstack(
        (-q_matrix @ a_block, -q_matrix @ e_block - c_matrix)
    )
    fixed_covariance = (
        bias @ state_shape @ bias.T + q_matrix @ a_block @ q_matrix.T
    )

    scenario_maximum = -np.inf
    for projection in partial_partition_projections(port_count):
        delta = np.hstack(
            (
                c_matrix @ projection @ f_matrix / mu - q_matrix,
                -c_matrix @ projection / mu,
            )
        )
        covariance = (
            fixed_covariance
            + c_matrix @ projection @ c_matrix.T / mu
            + score_radius**2 * delta @ delta.T
        )
        scenario_maximum = max(
            scenario_maximum, float(np.linalg.eigvalsh(covariance)[-1])
        )

    complete_edges = list(itertools.combinations(range(node_count), 2))
    for _ in range(100):
        selected = rng.random(len(complete_edges)) < 0.4
        weights = np.exp(rng.uniform(-3.0, 3.0, size=len(complete_edges)))
        edges = [
            (left, right, float(weight))
            for (left, right), keep, weight in zip(
                complete_edges, selected, weights, strict=True
            )
            if keep
        ]
        inverse = np.linalg.inv(
            mu * np.eye(node_count) + laplacian(node_count, edges)
        )
        c_full = np.hstack((c_matrix, np.zeros((output_count, hidden_count))))
        f_full = np.vstack((f_matrix, np.zeros((hidden_count, local_count))))
        delta = np.hstack(
            (c_full @ inverse @ f_full - q_matrix, -c_full @ inverse)
        )
        covariance = (
            fixed_covariance
            + c_full @ inverse @ c_full.T
            + score_radius**2 * delta @ delta.T
        )
        actual = float(np.linalg.eigvalsh(covariance)[-1])
        assert actual <= scenario_maximum + 5.0e-11


def check_heterogeneous_grounding_upper_bound() -> None:
    """Random test of the weighted theorem after diagonal whitening."""

    rng = np.random.default_rng(15092026)
    port_count, hidden_count = 3, 4
    node_count = port_count + hidden_count
    grounding = np.exp(rng.uniform(-1.5, 1.5, size=node_count))
    square_root = np.diag(np.sqrt(grounding))
    c_matrix = rng.normal(size=(2, port_count))
    f_matrix = rng.normal(size=(port_count, 3))
    q_matrix = rng.normal(size=(2, 3))

    scenarios = weighted_partial_partition_projections(grounding[:port_count])
    assert len(scenarios) == 15
    scenario_maximum = max(
        normalized_partial_scenario_norm(
            projection, c_matrix, f_matrix, q_matrix
        )
        for projection in scenarios
    )

    complete_edges = list(itertools.combinations(range(node_count), 2))
    for _ in range(100):
        selected = rng.random(len(complete_edges)) < 0.4
        weights = np.exp(rng.uniform(-4.0, 4.0, size=len(complete_edges)))
        edges = [
            (left, right, float(weight))
            for (left, right), keep, weight in zip(
                complete_edges, selected, weights, strict=True
            )
            if keep
        ]
        precision = np.diag(grounding) + laplacian(node_count, edges)
        normalized_inverse = square_root @ np.linalg.inv(precision) @ square_root
        actual = normalized_full_error_norm(
            normalized_inverse, port_count, c_matrix, f_matrix, q_matrix
        )
        assert actual <= scenario_maximum + 3.0e-12


def check_heterogeneous_grounding_degree_two_converse() -> None:
    """A weighted active block plus one port diluted down a degree-two path."""

    rng = np.random.default_rng(16092026)
    port_count = 3
    port_grounding = np.array([2.0, 0.4, 3.5])
    c_matrix = rng.normal(size=(2, port_count))
    f_matrix = rng.normal(size=(port_count, 2))
    q_matrix = rng.normal(size=(2, 2))

    target_projection = np.zeros((port_count, port_count))
    active = np.array([0, 2])
    direction = np.sqrt(port_grounding[active])
    target_projection[np.ix_(active, active)] = (
        np.outer(direction, direction) / np.dot(direction, direction)
    )
    target = normalized_partial_scenario_norm(
        target_projection, c_matrix, f_matrix, q_matrix
    )

    errors = []
    for hidden_count in (8, 16, 32, 64, 128):
        node_count = port_count + hidden_count
        grounding = np.concatenate((port_grounding, np.ones(hidden_count)))
        conductance = hidden_count**3
        edges: list[tuple[int, int, float]] = [(0, 2, conductance)]
        previous = 1
        for hidden in range(port_count, node_count):
            edges.append((previous, hidden, conductance))
            previous = hidden
        precision = np.diag(grounding) + laplacian(node_count, edges)
        root = np.diag(np.sqrt(grounding))
        normalized_inverse = root @ np.linalg.inv(precision) @ root
        actual = normalized_full_error_norm(
            normalized_inverse, port_count, c_matrix, f_matrix, q_matrix
        )
        errors.append(abs(actual - target))

    assert errors[-1] < 0.15 * errors[0]
    assert errors[-1] < 1.5e-2


def main() -> None:
    checks = [
        check_random_hidden_graph_upper_bound,
        check_degree_two_dilution_converse,
        check_mixed_state_raw_score_upper_bound,
        check_heterogeneous_grounding_upper_bound,
        check_heterogeneous_grounding_degree_two_converse,
    ]
    for check in checks:
        check()
        print(f"PASS {check.__name__}")


if __name__ == "__main__":
    main()
