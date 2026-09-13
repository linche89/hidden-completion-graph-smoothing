"""Reproduce the 2x2 counterexample for spectral-window robust boundaries.

This script only evaluates a two-dimensional analytic construction.  It does not
write files and is not used as a proof; the exact calculation is in the report.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize_scalar


ALPHA = 0.5
BETA = 1.0
F_SCALE = 4.0


def loss(q: float, theta: float) -> float:
    """Norm at S=alpha*I+(beta-alpha)*u*u.T, u=(cos,sin)."""
    u = np.array([math.cos(theta), math.sin(theta)], dtype=np.float64)
    s = ALPHA * np.eye(2) + (BETA - ALPHA) * np.outer(u, u)
    c = np.array([[1.0, 0.0]])
    f = np.array([[0.0], [F_SCALE]])
    block = np.concatenate((c @ s @ f - q, -(c @ s)), axis=1)
    return float(np.linalg.norm(block, ord=2))


def worst_projection(q: float) -> float:
    result = minimize_scalar(
        lambda theta: -loss(q, theta),
        bounds=(0.0, math.pi),
        method="bounded",
        options={"xatol": 1e-14},
    )
    return -float(result.fun)


def main() -> None:
    robust = minimize_scalar(
        worst_projection,
        bounds=(-3.0, 3.0),
        method="bounded",
        options={"xatol": 1e-13},
    )
    exact = 5.0 * math.sqrt(17.0) / 16.0
    y_star = 19.0 / 32.0  # cos(theta)^2
    isotropic_only = 1.0
    # Exact continuous-window SDP certificate from the report, at q=0,
    # tau=R^2=425/256 and lambda=17.
    h = (ALPHA + BETA) / 2.0
    c = np.array([[1.0, 0.0]])
    f = np.array([[0.0], [F_SCALE]])
    d_matrix = np.block(
        [
            [-np.eye(2), h * c.T],
            [h * c, -ALPHA * BETA * (c @ c.T)],
        ]
    )
    e_matrix = np.diag([0.0, 0.0, 1.0])
    l_matrix = np.block(
        [[f.T, np.zeros((1, 1))], [np.eye(2), np.zeros((2, 1))]]
    )
    tau = 425.0 / 256.0
    certificate = np.block(
        [
            [-tau * e_matrix + 17.0 * d_matrix, l_matrix.T],
            [l_matrix, -np.eye(3)],
        ]
    )
    certificate_max_eigenvalue = float(np.linalg.eigvalsh(certificate)[-1])
    print(
        {
            "q_numeric": float(robust.x),
            "projection_worst_numeric": float(robust.fun),
            "projection_worst_exact": exact,
            "cos2_theta_exact": y_star,
            "mI_MI_only_minimax": isotropic_only,
            "strict_gap": exact - isotropic_only,
            "exact_sdp_certificate_max_eigenvalue": certificate_max_eigenvalue,
        }
    )
    if (
        abs(robust.x) > 1e-6
        or abs(robust.fun - exact) > 1e-9
        or certificate_max_eigenvalue > 1e-10
    ):
        raise RuntimeError("numeric search did not reproduce the analytic result")


if __name__ == "__main__":
    main()
