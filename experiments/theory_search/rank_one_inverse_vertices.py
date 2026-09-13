"""Checks for the rank-one inverse-information corner theorem.

For a positive definite base matrix H and scalar-factor vectors a_e, the map

    T(w) = (H + sum_e w_e a_e a_e.T)^(-1)

lies in the convex hull of its box-corner values.  The key fact is that, with
all other weights fixed, a Sherman--Morrison rank-one update makes T(w_e)
trace a *line segment* between its two endpoint inverses.

The checks below verify a constructive recursive decomposition, the related
variable-size volume-sampling identity for unbounded weights, and a rank-two
counterexample showing why independently weighted scalar factors matter.
"""

from __future__ import annotations

import itertools

import numpy as np


def inverse_information(
    base: np.ndarray, factors: np.ndarray, weights: np.ndarray
) -> np.ndarray:
    information = base.copy()
    for index, weight in enumerate(weights):
        vector = factors[:, index]
        information += weight * np.outer(vector, vector)
    return np.linalg.inv(information)


def recursive_corner_decomposition(
    base: np.ndarray,
    factors: np.ndarray,
    weights: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    coordinate: int = 0,
) -> dict[tuple[float, ...], float]:
    """Construct convex coefficients for the box-corner inverse matrices."""

    if coordinate == weights.size:
        return {tuple(float(item) for item in weights): 1.0}

    low, high = float(lower[coordinate]), float(upper[coordinate])
    current = float(weights[coordinate])
    if high == low:
        return recursive_corner_decomposition(
            base, factors, weights, lower, upper, coordinate + 1
        )

    low_weights = weights.copy()
    low_weights[coordinate] = low
    low_inverse = inverse_information(base, factors, low_weights)
    vector = factors[:, coordinate]
    alpha = float(vector @ low_inverse @ vector)
    displacement = current - low
    width = high - low
    current_coefficient = displacement / (1.0 + displacement * alpha)
    endpoint_coefficient = width / (1.0 + width * alpha)
    theta = current_coefficient / endpoint_coefficient

    high_weights = weights.copy()
    high_weights[coordinate] = high
    low_parts = recursive_corner_decomposition(
        base, factors, low_weights, lower, upper, coordinate + 1
    )
    high_parts = recursive_corner_decomposition(
        base, factors, high_weights, lower, upper, coordinate + 1
    )
    result: dict[tuple[float, ...], float] = {}
    for corner, coefficient in low_parts.items():
        result[corner] = result.get(corner, 0.0) + (1.0 - theta) * coefficient
    for corner, coefficient in high_parts.items():
        result[corner] = result.get(corner, 0.0) + theta * coefficient
    return result


def check_finite_box_corner_decomposition() -> None:
    rng = np.random.default_rng(16092026)
    dimension, factor_count = 4, 5
    matrix = rng.normal(size=(dimension, dimension))
    base = matrix @ matrix.T + 0.8 * np.eye(dimension)
    factors = rng.normal(size=(dimension, factor_count))
    lower = rng.uniform(0.0, 0.4, size=factor_count)
    upper = lower + rng.uniform(0.2, 3.0, size=factor_count)

    for _ in range(20):
        weights = rng.uniform(lower, upper)
        decomposition = recursive_corner_decomposition(
            base, factors, weights, lower, upper
        )
        assert len(decomposition) == 2**factor_count
        coefficients = np.array(list(decomposition.values()))
        assert coefficients.min() >= -2.0e-14
        np.testing.assert_allclose(coefficients.sum(), 1.0, atol=2.0e-14)

        reconstructed = np.zeros_like(base)
        for corner, coefficient in decomposition.items():
            reconstructed += coefficient * inverse_information(
                base, factors, np.asarray(corner)
            )
        target = inverse_information(base, factors, weights)
        np.testing.assert_allclose(reconstructed, target, rtol=2e-13, atol=2e-13)


def nullspace_projection(columns: np.ndarray) -> np.ndarray:
    dimension = columns.shape[0]
    if columns.shape[1] == 0:
        return np.eye(dimension)
    return np.eye(dimension) - columns @ np.linalg.pinv(columns)


def check_volume_sampling_projection_identity() -> None:
    rng = np.random.default_rng(17092026)
    dimension, factor_count = 4, 7
    factors = rng.normal(size=(dimension, factor_count))
    weights = np.exp(rng.uniform(-2.0, 2.0, size=factor_count))
    scaled = factors * np.sqrt(weights)[None, :]

    numerator = np.zeros((dimension, dimension), dtype=float)
    normalizer = 0.0
    for subset_size in range(dimension + 1):
        for subset in itertools.combinations(range(factor_count), subset_size):
            columns = scaled[:, subset]
            volume_squared = (
                float(np.linalg.det(columns.T @ columns)) if subset else 1.0
            )
            # Roundoff can make a theoretically zero Gram determinant tiny and
            # negative for dependent subsets.
            if volume_squared <= 1.0e-12:
                continue
            numerator += volume_squared * nullspace_projection(columns)
            normalizer += volume_squared

    target = np.linalg.inv(np.eye(dimension) + scaled @ scaled.T)
    determinant_normalizer = float(
        np.linalg.det(np.eye(dimension) + scaled @ scaled.T)
    )
    np.testing.assert_allclose(normalizer, determinant_normalizer, rtol=2e-12)
    np.testing.assert_allclose(
        numerator / normalizer, target, rtol=2e-12, atol=2e-12
    )


def check_unbounded_corner_limits() -> None:
    rng = np.random.default_rng(18092026)
    dimension, factor_count = 5, 3
    factors = rng.normal(size=(dimension, factor_count))
    expected = nullspace_projection(factors)
    for scale, tolerance in ((1.0e4, 2.0e-3), (1.0e7, 3.0e-6)):
        actual = np.linalg.inv(
            np.eye(dimension) + scale * factors @ factors.T
        )
        assert np.linalg.norm(actual - expected, ord=2) < tolerance


def check_shared_rank_two_endpoint_failure() -> None:
    # A single shared scalar multiplies a rank-two update.  The inverse no
    # longer travels on a line segment, and even a linear (hence convex) loss
    # can be strictly maximized in the interior.
    def inverse_at(weight: float) -> np.ndarray:
        return np.diag([1.0 / (1.0 + weight), 1.0 / (1.0 + 4.0 * weight)])

    def convex_loss(matrix: np.ndarray) -> float:
        return float(matrix[0, 0] - matrix[1, 1])

    lower, upper, interior = 0.0, 10.0, 0.5
    endpoint_maximum = max(convex_loss(inverse_at(lower)), convex_loss(inverse_at(upper)))
    interior_value = convex_loss(inverse_at(interior))
    np.testing.assert_allclose(interior_value, 1.0 / 3.0)
    assert interior_value > endpoint_maximum + 0.25


def check_unbounded_full_block_relaxation_gap() -> None:
    """The direction-free contraction relaxation can be arbitrarily conservative."""

    scale = 40.0
    c_row = np.array([[1.0, 0.0]])
    f_matrix = np.diag([1.0, scale])
    q_row = np.array([[1.0, 0.0]])

    def row_loss(matrix: np.ndarray, correction: np.ndarray) -> float:
        error = np.hstack(
            (c_row @ matrix @ f_matrix - correction, -c_row @ matrix)
        )
        return float(np.linalg.norm(error, ord=2))

    # One uncertain scalar factor a=e_1 generates only the segment
    # diag(t,1), 0<=t<=1.  Convexity puts its maximum at t=0 or t=1.
    actual_extremes = [np.eye(2), np.diag([0.0, 1.0])]
    actual_value = max(row_loss(item, q_row) for item in actual_extremes)
    np.testing.assert_allclose(actual_value, 1.0)

    # The direction-free outer set 0<=T<=I contains two rotated projections.
    # Their visible second coordinates are +/-scale/2.  For every common
    # correction Q, at least one is at distance >=scale/2 in that coordinate.
    plus = 0.5 * np.array([[1.0, 1.0], [1.0, 1.0]])
    minus = 0.5 * np.array([[1.0, -1.0], [-1.0, 1.0]])
    for second_coordinate in np.linspace(-scale, scale, 101):
        arbitrary_q = np.array([[0.37, second_coordinate]])
        relaxed_pair_value = max(
            row_loss(plus, arbitrary_q), row_loss(minus, arbitrary_q)
        )
        assert relaxed_pair_value >= scale / 2.0 - 1.0e-12

    assert scale / (2.0 * actual_value) == 20.0


def main() -> None:
    checks = [
        check_finite_box_corner_decomposition,
        check_volume_sampling_projection_identity,
        check_unbounded_corner_limits,
        check_shared_rank_two_endpoint_failure,
        check_unbounded_full_block_relaxation_gap,
    ]
    for check in checks:
        check()
        print(f"PASS {check.__name__}")


if __name__ == "__main__":
    main()
