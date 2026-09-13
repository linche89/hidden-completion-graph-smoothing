"""Numerical regression tests for the phase-3 mathematical statements.

These tests do not prove the theorems.  They protect the paper and experiment
code against sign, block-orientation, norm, and normalization mistakes.
"""

from __future__ import annotations

import unittest

import numpy as np
import scipy.linalg as la
import scipy.optimize as opt
import scipy.special as special


def spectral_norm(matrix: np.ndarray) -> float:
    return float(la.svdvals(matrix)[0]) if matrix.size else 0.0


def weighted_path_information(
    n: int, *, edge_weight: float = 1.0, root_ground: float = 1.0
) -> np.ndarray:
    """Return L_W + diag(kappa) for a path grounded only at node zero."""

    if n < 2:
        raise ValueError("n must be at least two")
    laplacian = np.zeros((n, n), dtype=float)
    for node in range(n - 1):
        laplacian[node, node] += edge_weight
        laplacian[node + 1, node + 1] += edge_weight
        laplacian[node, node + 1] -= edge_weight
        laplacian[node + 1, node] -= edge_weight
    laplacian[0, 0] += root_ground
    return laplacian


class TestSchurIdentities(unittest.TestCase):
    def test_exact_factorization_and_norm_sandwich(self) -> None:
        rng = np.random.default_rng(20260913)
        for inside_dim, outside_dim in ((3, 2), (6, 4), (8, 5)):
            a_seed = rng.normal(size=(inside_dim, inside_dim))
            a = a_seed.T @ a_seed + 0.75 * np.eye(inside_dim)
            coupling = rng.normal(size=(inside_dim, outside_dim))
            t_seed = rng.normal(size=(outside_dim, outside_dim))
            schur = t_seed.T @ t_seed + 0.5 * np.eye(outside_dim)
            outside = schur + coupling.T @ la.solve(
                a, coupling, assume_a="pos"
            )
            information = np.block(
                [[a, coupling], [coupling.T, outside]]
            )

            selector = np.zeros((1, inside_dim + outside_dim))
            selector[0, 0] = 1.0
            full_row = la.solve(
                information, selector.T, assume_a="pos"
            ).T
            local_inside_row = la.solve(
                a, selector[:, :inside_dim].T, assume_a="pos"
            ).T
            local_row = np.hstack(
                [local_inside_row, np.zeros((1, outside_dim))]
            )

            y = full_row[:, inside_dim:]
            f = coupling.T @ la.solve(a, np.eye(inside_dim), assume_a="pos")
            factorized_error = y @ np.hstack([-f, np.eye(outside_dim)])
            actual_error = full_row - local_row

            np.testing.assert_allclose(
                actual_error, factorized_error, rtol=2e-11, atol=2e-11
            )

            oracle_tail = spectral_norm(y)
            dirichlet_error = spectral_norm(actual_error)
            gamma = np.sqrt(1.0 + spectral_norm(f) ** 2)
            self.assertLessEqual(oracle_tail, dirichlet_error + 1e-11)
            self.assertLessEqual(
                dirichlet_error, gamma * oracle_tail + 1e-11
            )
            self.assertAlmostEqual(
                spectral_norm(np.hstack([-f, np.eye(outside_dim)])),
                gamma,
                places=11,
            )

    def test_sharp_spectral_window_constant(self) -> None:
        lower, upper = 0.2, 8.2
        cosine = np.sqrt(upper / (lower + upper))
        sine = np.sqrt(lower / (lower + upper))
        rotation = np.asarray([[cosine, -sine], [sine, cosine]])
        information = rotation @ np.diag([lower, upper]) @ rotation.T
        inside = float(information[0, 0])
        coupling = float(information[1, 0])
        observed = abs(coupling / inside)
        sharp = (upper - lower) / (2.0 * np.sqrt(lower * upper))
        self.assertAlmostEqual(observed, sharp, places=12)
        self.assertAlmostEqual(
            np.sqrt(1.0 + observed**2),
            (upper + lower) / (2.0 * np.sqrt(lower * upper)),
            places=12,
        )

    def test_scalar_schur_window_robust_minimax(self) -> None:
        transfer = 0.7
        alpha, beta = 0.2, 1.5

        def check(local_factor: float) -> None:
            threshold = (alpha + beta) / (beta - alpha)
            if local_factor**2 <= threshold:
                optimum = abs(transfer) * beta
                correction = transfer * beta * local_factor
            else:
                imbalance = (
                    local_factor**2 * (beta - alpha) - (alpha + beta)
                )
                optimum = abs(transfer) * np.sqrt(
                    beta**2 + imbalance**2 / (4.0 * local_factor**2)
                )
                correction = transfer * (
                    (local_factor**2 + 1.0)
                    * (alpha + beta)
                    / (2.0 * local_factor)
                )

            def worst_squared(q: float) -> float:
                return max(
                    (transfer * local_factor * s - q) ** 2
                    + (transfer * s) ** 2
                    for s in (alpha, beta)
                )

            numerical = opt.minimize_scalar(
                worst_squared,
                bracket=(-10.0, 10.0),
                method="brent",
                options={"xtol": 1e-13},
            )
            self.assertTrue(numerical.success)
            self.assertAlmostEqual(np.sqrt(numerical.fun), optimum, places=9)
            self.assertAlmostEqual(np.sqrt(worst_squared(correction)), optimum, places=10)
            dirichlet = np.sqrt(worst_squared(0.0))
            self.assertGreater(dirichlet, optimum)

        check(0.6)
        check(2.0)


class TestRobustBoundaryGeometry(unittest.TestCase):
    def test_general_matrix_interval_image_is_an_ellipsoid(self) -> None:
        """Construct S in [S_minus,S_plus] for an arbitrary feasible image Sv."""

        rng = np.random.default_rng(31051984)
        dimension = 5
        delta_seed = rng.normal(size=(dimension, dimension))
        delta = delta_seed.T @ delta_seed + 0.5 * np.eye(dimension)
        eigvals, eigvecs = la.eigh(delta)
        delta_sqrt = (eigvecs * np.sqrt(eigvals)) @ eigvecs.T
        delta_inv_sqrt = (eigvecs * (1.0 / np.sqrt(eigvals))) @ eigvecs.T

        midpoint = 3.0 * delta + np.eye(dimension)
        lower = midpoint - delta
        upper = midpoint + delta
        vector = rng.normal(size=dimension)
        source = delta_sqrt @ vector

        target_direction = rng.normal(size=dimension)
        target_direction /= la.norm(target_direction)
        contraction_ratio = 0.63
        target = contraction_ratio * la.norm(source) * target_direction

        source_direction = source / la.norm(source)
        if np.allclose(source_direction, target_direction):
            reflection = np.eye(dimension)
        elif np.allclose(source_direction, -target_direction):
            reflection = -np.eye(dimension)
        else:
            normal = source_direction - target_direction
            reflection = np.eye(dimension) - 2.0 * np.outer(
                normal, normal
            ) / float(normal @ normal)
        contraction = contraction_ratio * reflection
        matrix = midpoint + delta_sqrt @ contraction @ delta_sqrt
        image = midpoint @ vector + delta_sqrt @ target

        np.testing.assert_allclose(matrix, matrix.T, rtol=0.0, atol=2e-12)
        np.testing.assert_allclose(matrix @ vector, image, rtol=2e-11, atol=2e-11)
        self.assertGreaterEqual(float(la.eigvalsh(matrix - lower)[0]), -2e-11)
        self.assertGreaterEqual(float(la.eigvalsh(upper - matrix)[0]), -2e-11)

        displacement = image - midpoint @ vector
        lhs = float(displacement @ la.solve(delta, displacement, assume_a="pos"))
        rhs = float(vector @ delta @ vector)
        self.assertLessEqual(lhs, rhs + 2e-11)
        np.testing.assert_allclose(
            delta_inv_sqrt @ displacement,
            target,
            rtol=2e-11,
            atol=2e-11,
        )

    def test_anisotropic_projection_beats_isotropic_endpoint_check(self) -> None:
        alpha, beta = 0.5, 1.0
        coupling_row = np.array([[1.0, 0.0]])
        local_transfer = np.array([[0.0], [4.0]])
        isotropic_worst = max(
            spectral_norm(
                np.hstack(
                    [
                        coupling_row @ (scale * np.eye(2)) @ local_transfer,
                        -coupling_row @ (scale * np.eye(2)),
                    ]
                )
            )
            for scale in (alpha, beta)
        )

        cosine_squared = 19.0 / 32.0
        direction = np.array(
            [np.sqrt(cosine_squared), np.sqrt(1.0 - cosine_squared)]
        )
        extreme = alpha * np.eye(2) + (beta - alpha) * np.outer(
            direction, direction
        )
        anisotropic_loss = spectral_norm(
            np.hstack(
                [
                    coupling_row @ extreme @ local_transfer,
                    -coupling_row @ extreme,
                ]
            )
        )

        self.assertAlmostEqual(isotropic_worst, 1.0, places=12)
        self.assertAlmostEqual(
            anisotropic_loss, 5.0 * np.sqrt(17.0) / 16.0, places=12
        )
        self.assertGreater(anisotropic_loss, isotropic_worst)


class TestMajoritySpectralTail(unittest.TestCase):
    def test_arbitrary_noncommuting_stable_wls_rule_obeys_tail_converse(
        self,
    ) -> None:
        """Check the low-spectrum converse without assuming L commutes with J."""

        rng = np.random.default_rng(16092026)
        dimension = 11
        eigenvalues = np.geomspace(0.015, 3.0, dimension)
        eigenvectors, _ = la.qr(rng.normal(size=(dimension, dimension)))
        information = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T
        inverse = eigenvectors @ np.diag(1.0 / eigenvalues) @ eigenvectors.T
        square_root = (
            eigenvectors @ np.diag(np.sqrt(eigenvalues)) @ eigenvectors.T
        )
        inverse_square_root = (
            eigenvectors @ np.diag(1.0 / np.sqrt(eigenvalues)) @ eigenvectors.T
        )

        local_rule = rng.normal(size=(dimension, dimension))
        # This is almost surely noncommuting, so the regression guards the
        # genuinely non-spectral step in the converse.
        self.assertGreater(
            spectral_norm(local_rule @ information - information @ local_rule),
            1e-3,
        )

        cutoff = 0.19
        projector = (
            eigenvectors[:, eigenvalues < cutoff]
            @ eigenvectors[:, eigenvalues < cutoff].T
        )
        root = 6
        low_spectrum_tail = la.norm(
            inverse_square_root[root, :] @ projector
        )
        wls_error = la.norm(
            (inverse[root, :] - local_rule[root, :]) @ square_root
        )
        row_gain = la.norm(local_rule[root, :])

        self.assertLessEqual(
            low_spectrum_tail,
            wls_error + row_gain * np.sqrt(cutoff) + 2e-12,
        )

    def test_power_weight_christoffel_formula(self) -> None:
        """The exact endpoint formula behind the candidate round-rate law."""

        exponent = 0.75
        degree = 6
        indices = np.arange(degree + 1)
        moment = 1.0 / (
            exponent + 1.0 + indices[:, None] + indices[None, :]
        )
        inverse_moment = la.inv(moment)
        kernel_at_zero = float(inverse_moment[0, 0])
        exact_kernel = np.exp(
            -np.log(exponent + 1.0)
            + 2.0
            * (
                special.gammaln(degree + exponent + 2.0)
                - special.gammaln(exponent + 1.0)
                - special.gammaln(degree + 1.0)
            )
        )
        self.assertAlmostEqual(
            kernel_at_zero / exact_kernel, 1.0, places=8
        )

        # The constrained quadratic minimum c^T M c, c[0] = 1, equals 1/K.
        optimal_coefficients = (
            inverse_moment[:, 0] / inverse_moment[0, 0]
        )
        optimal_value = float(
            optimal_coefficients @ moment @ optimal_coefficients
        )
        np.testing.assert_allclose(
            optimal_value * kernel_at_zero, 1.0, rtol=5e-9, atol=5e-12
        )

    def test_neumann_row_error_equals_local_spectral_integral(self) -> None:
        rng = np.random.default_rng(1172026)
        dimension = 9
        eigenvalues = np.linspace(0.07, 2.0, dimension)
        orthogonal, _ = la.qr(rng.normal(size=(dimension, dimension)))
        information = orthogonal @ np.diag(eigenvalues) @ orthogonal.T
        spectral_upper = 2.0
        steps = 13

        approximation = np.zeros_like(information)
        power = np.eye(dimension)
        iteration = np.eye(dimension) - information / spectral_upper
        for _ in range(steps):
            approximation += power / spectral_upper
            power = power @ iteration
        error = la.solve(information, np.eye(dimension), assume_a="pos") - approximation

        root = 4
        row_error_squared = float(error[root, :] @ error[root, :])
        local_weights = orthogonal[root, :] ** 2
        integral = float(
            np.sum(
                local_weights
                * (1.0 - eigenvalues / spectral_upper) ** (2 * steps)
                / eigenvalues**2
            )
        )
        self.assertAlmostEqual(row_error_squared, integral, places=11)

        cutoff = 0.4
        inverse_square_tail = float(
            np.sum(local_weights[eigenvalues < cutoff] / eigenvalues[eigenvalues < cutoff] ** 2)
        )
        high_spectrum_bound = (
            (1.0 - cutoff / spectral_upper) ** (2 * steps) / cutoff**2
        )
        self.assertLessEqual(integral, inverse_square_tail + high_spectrum_bound)

    def test_whitened_wls_changes_inverse_square_to_inverse_first_power(self) -> None:
        rng = np.random.default_rng(20112026)
        measurements = 14
        states = 7
        factor = rng.normal(size=(measurements, states))
        information = factor.T @ factor + 0.4 * np.eye(states)
        # Treat the ridge rows as virtual whitened measurements, so G G^T = J.
        whitened_map = np.hstack(
            [factor.T, np.sqrt(0.4) * np.eye(states)]
        )
        np.testing.assert_allclose(
            whitened_map @ whitened_map.T,
            information,
            rtol=2e-12,
            atol=2e-12,
        )

        eigenvalues, eigenvectors = la.eigh(information)
        spectral_upper = float(eigenvalues[-1])
        steps = 8
        polynomial = np.zeros_like(information)
        power = np.eye(states)
        iteration = np.eye(states) - information / spectral_upper
        for _ in range(steps):
            polynomial += power / spectral_upper
            power = power @ iteration

        error_operator = (
            la.solve(information, np.eye(states), assume_a="pos") - polynomial
        ) @ whitened_map
        root = 3
        actual = float(error_operator[root, :] @ error_operator[root, :])
        local_weights = eigenvectors[root, :] ** 2
        spectral_formula = float(
            np.sum(
                local_weights
                * (1.0 - eigenvalues / spectral_upper) ** (2 * steps)
                / eigenvalues
            )
        )
        self.assertAlmostEqual(actual, spectral_formula, places=11)

    def test_localized_and_delocalized_small_modes_differ_for_most_nodes(self) -> None:
        dimension = 32
        tiny = dimension ** -2
        cutoff = 2.0 * tiny

        localized = np.ones(dimension)
        localized[0] = tiny
        localized_tail = np.where(localized < cutoff, localized ** -2, 0.0)
        self.assertEqual(int(np.count_nonzero(localized_tail)), 1)

        laplacian = np.zeros((dimension, dimension))
        for node in range(dimension):
            neighbor = (node + 1) % dimension
            laplacian[node, node] += 1.0
            laplacian[neighbor, neighbor] += 1.0
            laplacian[node, neighbor] -= 1.0
            laplacian[neighbor, node] -= 1.0
        delocalized = laplacian + tiny * np.eye(dimension)
        eigenvalues, eigenvectors = la.eigh(delocalized)
        low = eigenvalues < cutoff
        local_inverse_square_mass = (
            eigenvectors[:, low] ** 2 @ (eigenvalues[low] ** -2)
        )
        np.testing.assert_allclose(
            local_inverse_square_mass,
            np.full(dimension, dimension**3),
            rtol=2e-8,
            atol=2e-8,
        )


class TestGroundedLaplacianIdentities(unittest.TestCase):
    def test_oracle_dirichlet_neumann_ordering(self) -> None:
        information = weighted_path_information(
            9, edge_weight=1.3, root_ground=0.8
        )
        strength = np.diag(information)
        weighted_adjacency = np.diag(strength) - information
        transition = weighted_adjacency / strength[:, None]
        fundamental = la.solve(
            np.eye(information.shape[0]) - transition,
            np.eye(information.shape[0]),
        )

        root = 2
        radius = 2
        inside = np.arange(0, 5)
        outside = np.arange(5, information.shape[0])
        local_fundamental = la.solve(
            np.eye(inside.size) - transition[np.ix_(inside, inside)],
            np.eye(inside.size),
        )

        oracle = float(fundamental[root, outside].sum())
        local_row_sum = float(local_fundamental[root, :].sum())
        dirichlet = float(fundamental[root, :].sum() - local_row_sum)
        neumann_tail = float(
            (np.linalg.matrix_power(transition, radius + 1) @ fundamental)[
                root, :
            ].sum()
        )

        exit_probability = float(
            local_fundamental[root, :]
            @ transition[np.ix_(inside, outside)]
            @ np.ones(outside.size)
        )
        lifetime = fundamental @ np.ones(information.shape[0])
        exit_states = outside[
            np.asarray(transition[np.ix_(inside, outside)].sum(axis=0) > 0)
        ]
        lifetime_envelope = float(lifetime[exit_states].max())

        self.assertGreater(oracle, 0.0)
        self.assertLessEqual(oracle, dirichlet + 1e-11)
        self.assertLessEqual(dirichlet, neumann_tail + 1e-11)
        self.assertLessEqual(exit_probability, oracle + 1e-11)
        self.assertLessEqual(
            dirichlet, lifetime_envelope * exit_probability + 1e-10
        )

    def test_single_ground_tentacle_has_unbounded_row_tail(self) -> None:
        radius = 2
        observed = []
        for n in (8, 16, 32, 64):
            information = weighted_path_information(
                n, edge_weight=1.0, root_ground=2.0
            )
            root_row = la.solve(
                information, np.eye(n)[:, 0], assume_a="pos"
            )
            np.testing.assert_allclose(
                root_row, np.full(n, 0.5), rtol=2e-11, atol=2e-11
            )
            outside_tail = float(la.norm(root_row[radius + 1 :], 2))
            expected = 0.5 * np.sqrt(n - radius - 1)
            self.assertAlmostEqual(outside_tail, expected, places=10)
            observed.append(outside_tail)

        self.assertTrue(all(a < b for a, b in zip(observed, observed[1:])))


if __name__ == "__main__":
    unittest.main(verbosity=2)
