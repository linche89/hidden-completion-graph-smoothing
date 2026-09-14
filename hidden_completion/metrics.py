"""Canonical full-input errors associated with terminal scenarios."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray


Array = NDArray[np.float64]


def _as_matrix(value: Array, name: str) -> Array:
    matrix = np.asarray(value, dtype=float)
    if matrix.ndim != 2 or not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} must be a finite two-dimensional array")
    return matrix


def positive_semidefinite_sqrt(matrix: Array, *, tolerance: float = 1e-9) -> Array:
    matrix = _as_matrix(matrix, "matrix")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("matrix must be square")
    symmetric = (matrix + matrix.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    scale = max(1.0, float(np.linalg.norm(symmetric, ord=2)))
    if float(eigenvalues[0]) < -tolerance * scale:
        raise ValueError(f"matrix is not positive semidefinite: lambda_min={eigenvalues[0]:.3e}")
    return (eigenvectors * np.sqrt(np.maximum(eigenvalues, 0.0))) @ eigenvectors.T


def scenario_error_gram(c_matrix: Array, q_matrix: Array, x_matrix: Array) -> Array:
    """Return M_Q(X), the exact left Gram matrix of a scenario error."""

    c_matrix = _as_matrix(c_matrix, "C")
    q_matrix = _as_matrix(q_matrix, "Q")
    x_matrix = _as_matrix(x_matrix, "X")
    if c_matrix.shape != q_matrix.shape:
        raise ValueError("C and Q must have the same shape")
    if x_matrix.shape != (c_matrix.shape[1], c_matrix.shape[1]):
        raise ValueError("X must be square with dimension equal to the number of ports")
    gram = (
        c_matrix @ x_matrix @ c_matrix.T
        - c_matrix @ x_matrix @ q_matrix.T
        - q_matrix @ x_matrix @ c_matrix.T
        + q_matrix @ q_matrix.T
    )
    return (gram + gram.T) / 2.0


def canonical_error(c_matrix: Array, q_matrix: Array, x_matrix: Array) -> Array:
    """Return ``[CX-Q, C(X-X^2)^(1/2)]``.

    It has the singular values of the full hidden-input limiting error, while
    using a fixed ``2q`` column dimension independent of the completion.
    """

    c_matrix = _as_matrix(c_matrix, "C")
    q_matrix = _as_matrix(q_matrix, "Q")
    x_matrix = _as_matrix(x_matrix, "X")
    if c_matrix.shape != q_matrix.shape:
        raise ValueError("C and Q must have the same shape")
    if x_matrix.shape != (c_matrix.shape[1], c_matrix.shape[1]):
        raise ValueError("X must be square with dimension equal to the number of ports")
    hidden_factor = positive_semidefinite_sqrt(x_matrix - x_matrix @ x_matrix)
    return np.concatenate((c_matrix @ x_matrix - q_matrix, c_matrix @ hidden_factor), axis=1)


def full_input_schatten_norm(
    c_matrix: Array,
    q_matrix: Array,
    x_matrix: Array,
    exponent: float = math.inf,
) -> float:
    """Evaluate the exact finite-scenario full-input Schatten norm."""

    if math.isnan(exponent) or exponent < 1 or (math.isinf(exponent) and exponent < 0):
        raise ValueError("the Schatten exponent must be at least one")
    singular_values = np.linalg.svd(canonical_error(c_matrix, q_matrix, x_matrix), compute_uv=False)
    if math.isinf(exponent):
        return float(singular_values[0]) if singular_values.size else 0.0
    return float(np.sum(singular_values**exponent) ** (1.0 / exponent))


def spectral_squared_risk(c_matrix: Array, q_matrix: Array, x_matrix: Array) -> float:
    """Evaluate lambda_max(M_Q(X)), clipped only at roundoff level."""

    gram = scenario_error_gram(c_matrix, q_matrix, x_matrix)
    return max(0.0, float(np.linalg.eigvalsh(gram)[-1]))
