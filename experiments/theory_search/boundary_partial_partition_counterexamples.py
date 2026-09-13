"""Counterexamples delimiting the hidden-node partition theorem.

The exact finite reduction is valid for the port compression of the
grounding-normalized inverse and for the spectral/Frobenius Schur losses
proved in the accompanying audit.  These checks show why it must not be
advertised for arbitrary convex losses on the full hidden-input matrix, for
Schatten-p losses with p < 2, or for an unwhitened Euclidean RHS under
heterogeneous grounding.
"""

from __future__ import annotations

import numpy as np


def partial_projections_two_ports() -> list[np.ndarray]:
    return [
        np.zeros((2, 2)),
        np.diag([1.0, 0.0]),
        np.diag([0.0, 1.0]),
        np.eye(2),
        np.ones((2, 2)) / 2.0,
    ]


def nuclear_schur_loss(
    projection: np.ndarray,
    c_matrix: np.ndarray,
    f_matrix: np.ndarray,
    q_matrix: np.ndarray,
) -> float:
    visible = c_matrix @ projection @ f_matrix - q_matrix
    gram = visible @ visible.T + c_matrix @ projection @ c_matrix.T
    return float(np.sqrt(np.maximum(np.linalg.eigvalsh(gram), 0.0)).sum())


def check_arbitrary_full_matrix_convex_loss_fails() -> None:
    # One port and one hidden node fused by an infinitely conductive edge.
    # The full normalized inverse tends to the averaging projection below.
    full_projection = np.ones((2, 2)) / 2.0
    full_invisible_row = -full_projection[0]

    # Padding the two one-port partial scenarios into the same two-column
    # space produces no hidden-column response at all.
    partial_rows = [np.array([0.0, 0.0]), np.array([-1.0, 0.0])]
    convex_loss = lambda row: abs(float(row[1]))
    assert convex_loss(full_invisible_row) == 0.5
    assert max(convex_loss(row) for row in partial_rows) == 0.0


def check_schatten_one_reduction_fails() -> None:
    # Two disjoint port--hidden pairs, each fused strongly, give the port
    # compression X = I/2.  The matrices below are rational decimals and
    # also satisfy C = R F^T for a suitable R (F is nonsingular).
    x_matrix = np.eye(2) / 2.0
    c_matrix = np.array([[-1.4, 0.0], [-1.2, -0.1]])
    f_matrix = np.array([[-1.0, 0.0], [-0.36, 0.1]])
    q_matrix = np.array([[6.3, 0.0], [0.5, 0.0]])
    r_matrix = c_matrix @ np.linalg.inv(f_matrix.T)
    np.testing.assert_allclose(c_matrix, r_matrix @ f_matrix.T, atol=1.0e-14)

    hidden_limit = nuclear_schur_loss(
        x_matrix, c_matrix, f_matrix, q_matrix
    )
    finite_scenarios = max(
        nuclear_schur_loss(item, c_matrix, f_matrix, q_matrix)
        for item in partial_projections_two_ports()
    )
    assert hidden_limit > finite_scenarios + 0.07


def check_direct_smoother_schatten_one_reduction_fails() -> None:
    # For direct approximation of the Tikhonov map, the full error is
    # C_full T - [Q, 0].  At a forest projection its Gram is
    # C X C^T - C X Q^T - Q X C^T + Q Q^T.  Schatten-1 can increase when
    # fractional hidden blocks mix singular directions, even though the
    # spectral and Frobenius reductions remain exact.
    x_matrix = np.eye(2) / 2.0
    c_matrix = np.array([[-2.0, 4.0], [-1.0, -6.0]])
    q_matrix = np.array([[-1.0, -1.0], [-1.0, -5.0]])
    gram = (
        c_matrix @ x_matrix @ c_matrix.T
        - c_matrix @ x_matrix @ q_matrix.T
        - q_matrix @ x_matrix @ c_matrix.T
        + q_matrix @ q_matrix.T
    )
    hidden_limit = float(
        np.sqrt(np.maximum(np.linalg.eigvalsh(gram), 0.0)).sum()
    )
    finite_scenarios = max(
        float(np.linalg.svd(c_matrix @ item - q_matrix, compute_uv=False).sum())
        for item in partial_projections_two_ports()
    )
    assert hidden_limit > finite_scenarios + 1.0


def check_direct_smoother_schatten_threshold_is_two() -> None:
    """One port disproves every Schatten-p reduction below p=2."""

    full_projection = np.ones((2, 2)) / 2.0
    c_matrix = np.array([[1.0], [1.0]])
    q_matrix = np.array([[1.0], [0.0]])
    selector = np.array([[1.0, 0.0]])
    full_error = c_matrix @ selector @ full_projection - q_matrix @ selector
    np.testing.assert_allclose(
        np.linalg.svd(full_error, compute_uv=False),
        np.ones(2) / np.sqrt(2.0),
        atol=1.0e-14,
    )

    partial_errors = [-q_matrix, c_matrix - q_matrix]
    for exponent in (1.0, 1.25, 1.5, 1.75):
        full_value = float(
            np.linalg.norm(
                np.linalg.svd(full_error, compute_uv=False), ord=exponent
            )
        )
        partial_value = max(
            float(
                np.linalg.norm(
                    np.linalg.svd(error, compute_uv=False), ord=exponent
                )
            )
            for error in partial_errors
        )
        assert full_value > partial_value + 1.0e-10

    np.testing.assert_allclose(
        np.linalg.norm(full_error, ord="fro"), 1.0, atol=1.0e-14
    )


def check_shared_opcode_connected_scenario_is_active() -> None:
    """A connected two-port block changes the exact shared-gain optimum."""

    c_matrix = np.diag([1.0, 2.0])
    identity = np.eye(2)
    alpha = 4.0 / 3.0 - 2.0 * np.sqrt(10.0) / 15.0
    radius = (10.0 + 2.0 * np.sqrt(10.0)) / 15.0
    values = [
        float(np.linalg.norm(c_matrix @ item - alpha * identity, ord=2))
        for item in partial_projections_two_ports()
    ]
    np.testing.assert_allclose(max(values), radius, atol=1.0e-14)
    np.testing.assert_allclose(values[2:], [radius, radius, radius], atol=1.0e-14)

    # Without the connected averaging projection, the four remaining
    # scenarios have the misleading optimum alpha=1 and radius=1.
    reduced_at_one = max(
        float(np.linalg.norm(c_matrix @ item - identity, ord=2))
        for item in partial_projections_two_ports()[:-1]
    )
    np.testing.assert_allclose(reduced_at_one, 1.0, atol=1.0e-14)
    assert radius > reduced_at_one + 0.08


def check_unwhitened_heterogeneous_rhs_can_be_unbounded() -> None:
    # Fuse one port (grounding 1) to h hidden nodes, each with grounding
    # 1/h^2.  In the infinite-conductance limit, the unnormalized inverse
    # row is constant with value 1 / (1 + 1/h), so its Euclidean norm grows
    # like sqrt(h).  Grounding whitening removes this pathology.
    values = []
    for hidden_count in (10, 100, 1000):
        hidden_grounding = 1.0 / hidden_count**2
        total_grounding = 1.0 + hidden_count * hidden_grounding
        values.append(np.sqrt(hidden_count + 1.0) / total_grounding)
    assert values[1] > 2.5 * values[0]
    assert values[2] > 2.5 * values[1]


def main() -> None:
    checks = [
        check_arbitrary_full_matrix_convex_loss_fails,
        check_schatten_one_reduction_fails,
        check_direct_smoother_schatten_one_reduction_fails,
        check_direct_smoother_schatten_threshold_is_two,
        check_shared_opcode_connected_scenario_is_active,
        check_unwhitened_heterogeneous_rhs_can_be_unbounded,
    ]
    for check in checks:
        check()
        print(f"PASS {check.__name__}")


if __name__ == "__main__":
    main()
