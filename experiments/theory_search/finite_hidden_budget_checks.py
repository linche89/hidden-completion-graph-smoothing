"""Deterministic checks for the finite-hidden-budget completion theorem.

The script is deliberately small-scale.  It checks the identities that are most
likely to be implemented incorrectly:

* the edge-DPP/forest reconstruction of ``(I + L)^{-1}``;
* the finite scenario formula for a forest compression;
* the full-input Gram term ``X - X**2`` (including nuclear norm);
* high-conductance, connected path witnesses;
* the exact vertex pruning rule for the *terminal hull*;
* the exact ``q/(q+h)`` Hausdorff distance for every ``q``; and
* the analytic ``q=2, h=1`` robust-smoother improvement example.

These are numerical regression checks, not substitutes for the proofs in
``research/agent_reports/finite_hidden_budget_theorem.md``.
"""

from __future__ import annotations

import itertools
import math
import unittest

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import linprog


Array = NDArray[np.float64]


def set_partitions(items: tuple[int, ...]) -> list[tuple[tuple[int, ...], ...]]:
    """Return each set partition once, with a deterministic block order."""

    if not items:
        return [tuple()]
    first, rest = items[0], items[1:]
    answer: list[tuple[tuple[int, ...], ...]] = []
    for partition in set_partitions(rest):
        answer.append(((first,),) + partition)
        for j in range(len(partition)):
            blocks = list(partition)
            blocks[j] = (first,) + blocks[j]
            answer.append(tuple(blocks))
    return answer


def allocations_at_most(h: int, blocks: int) -> list[tuple[int, ...]]:
    return [
        values
        for values in itertools.product(range(h + 1), repeat=blocks)
        if sum(values) <= h
    ]


def scenario_matrix(q: int, partition: tuple[tuple[int, ...], ...], k: tuple[int, ...]) -> Array:
    x = np.zeros((q, q))
    for block, hidden_count in zip(partition, k, strict=True):
        idx = np.asarray(block)
        x[np.ix_(idx, idx)] = 1.0 / (len(block) + hidden_count)
    return x


def scenarios(q: int, h: int) -> list[tuple[tuple[tuple[int, ...], ...], tuple[int, ...], Array]]:
    answer = []
    for partition in set_partitions(tuple(range(q))):
        for k in allocations_at_most(h, len(partition)):
            answer.append((partition, k, scenario_matrix(q, partition, k)))
    return answer


def laplacian(n: int, edges: list[tuple[int, int, float]]) -> Array:
    result = np.zeros((n, n))
    for u, v, weight in edges:
        d = np.zeros(n)
        d[u], d[v] = 1.0, -1.0
        result += weight * np.outer(d, d)
    return result


def schatten(matrix: Array, p: float) -> float:
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    if math.isinf(p):
        return float(singular_values[0])
    return float(np.sum(singular_values**p) ** (1.0 / p))


def positive_sqrt(matrix: Array) -> Array:
    eigenvalues, eigenvectors = np.linalg.eigh((matrix + matrix.T) / 2.0)
    if eigenvalues.min() < -2e-9:
        raise AssertionError(f"matrix is not PSD: lambda_min={eigenvalues.min()}")
    return (eigenvectors * np.sqrt(np.maximum(eigenvalues, 0.0))) @ eigenvectors.T


def canonical_error(c: Array, q_map: Array, x: Array) -> Array:
    # A full forest projection with port compression X has precisely this Gram.
    residual = c @ x - q_map
    hidden_energy = c @ positive_sqrt(x - x @ x)
    return np.concatenate((residual, hidden_energy), axis=1)


def union_find_components(n: int, selected_edges: list[tuple[int, int, float]]) -> list[list[int]]:
    parent = list(range(n))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a

    for u, v, _ in selected_edges:
        union(u, v)
    groups: dict[int, list[int]] = {}
    for vertex in range(n):
        groups.setdefault(find(vertex), []).append(vertex)
    return list(groups.values())


class FiniteHiddenBudgetChecks(unittest.TestCase):
    def test_forest_dpp_reconstructs_resolvent_and_scenarios(self) -> None:
        q, h, n = 3, 2, 5
        edges = [
            (0, 1, 0.7),
            (1, 2, 1.2),
            (2, 3, 0.4),
            (3, 4, 1.1),
            (0, 4, 0.9),
            (1, 3, 0.6),
        ]
        incidence_columns = []
        for u, v, weight in edges:
            column = np.zeros(n)
            column[u], column[v] = math.sqrt(weight), -math.sqrt(weight)
            incidence_columns.append(column)
        a = np.column_stack(incidence_columns)
        target = np.linalg.inv(np.eye(n) + a @ a.T)
        normalizer = np.linalg.det(np.eye(len(edges)) + a.T @ a)

        reconstructed = np.zeros((n, n))
        probability_sum = 0.0
        for mask in range(1 << len(edges)):
            indices = [j for j in range(len(edges)) if mask & (1 << j)]
            if indices:
                a_s = a[:, indices]
                gram = a_s.T @ a_s
                determinant = float(np.linalg.det(gram))
                if determinant <= 1e-11:  # a cycle: zero DPP mass
                    continue
                projection = np.eye(n) - a_s @ np.linalg.solve(gram, a_s.T)
            else:
                determinant = 1.0
                projection = np.eye(n)

            probability = determinant / normalizer
            probability_sum += probability
            reconstructed += probability * projection

            selected = [edges[j] for j in indices]
            components = union_find_components(n, selected)
            terminal_components = [component for component in components if any(v < q for v in component)]
            partition = tuple(tuple(v for v in component if v < q) for component in terminal_components)
            k = tuple(sum(v >= q for v in component) for component in terminal_components)
            x = scenario_matrix(q, partition, k)
            np.testing.assert_allclose(projection[:q, :q], x, rtol=2e-10, atol=2e-10)
            self.assertLessEqual(sum(k), h)

        self.assertAlmostEqual(probability_sum, 1.0, places=10)
        np.testing.assert_allclose(reconstructed, target, rtol=3e-10, atol=3e-10)

    def test_full_input_upper_bound_for_all_schatten_norms(self) -> None:
        q, h, n = 3, 2, 5
        edges = [
            (0, 1, 0.7),
            (1, 2, 1.2),
            (2, 3, 0.4),
            (3, 4, 1.1),
            (0, 4, 0.9),
            (1, 3, 0.6),
        ]
        t = np.linalg.inv(np.eye(n) + laplacian(n, edges))
        select = np.zeros((q, n))
        select[:, :q] = np.eye(q)
        c = np.array([[1.1, -0.3, 0.8], [0.2, 0.9, -0.5]])
        q_map = np.array([[0.1, 0.4, -0.2], [-0.3, 0.2, 0.5]])
        actual_error = c @ select @ t - q_map @ select

        for p in (1.0, 2.0, 4.0, math.inf):
            scenario_bound = max(
                schatten(canonical_error(c, q_map, x), p) for _, _, x in scenarios(q, h)
            )
            self.assertLessEqual(schatten(actual_error, p), scenario_bound + 2e-10)

    def test_connected_path_witness_and_hidden_gram(self) -> None:
        q, h, n = 3, 2, 5
        partition = ((0, 2), (1,))
        k = (1, 0)
        target_x = scenario_matrix(q, partition, k)

        # Strong segments [0,3,2], [1], [4] are concatenated by weak edges.
        strong, weak = 2e7, 2e-8
        edges = [(0, 3, strong), (3, 2, strong), (2, 1, weak), (1, 4, weak)]
        degrees = [0] * n
        for u, v, _ in edges:
            degrees[u] += 1
            degrees[v] += 1
        self.assertLessEqual(max(degrees), 2)

        t = np.linalg.inv(np.eye(n) + laplacian(n, edges))
        np.testing.assert_allclose(t[:q, :q], target_x, rtol=0.0, atol=2e-7)

        c = np.array([[1.0, -0.2, 0.7], [0.3, 0.8, -0.4]])
        q_map = np.array([[0.1, 0.2, -0.3], [0.5, -0.2, 0.1]])
        select = np.zeros((q, n))
        select[:, :q] = np.eye(q)
        actual_error = c @ select @ t - q_map @ select
        gram = c @ target_x @ c.T - c @ target_x @ q_map.T - q_map @ target_x @ c.T + q_map @ q_map.T
        np.testing.assert_allclose(actual_error @ actual_error.T, gram, rtol=0.0, atol=5e-7)
        for p in (1.0, 2.0, 4.0, math.inf):
            self.assertAlmostEqual(
                schatten(actual_error, p), schatten(canonical_error(c, q_map, target_x), p), places=6
            )

    def test_candidate_count_and_exact_terminal_vertices(self) -> None:
        q, h = 3, 2
        candidates = scenarios(q, h)
        # S(3,1) C(3,1) + S(3,2) C(4,2) + S(3,3) C(5,3) = 31.
        self.assertEqual(len(candidates), 31)

        upper = np.triu_indices(q)
        points = np.column_stack([x[upper] for _, _, x in candidates])
        numerical_vertices: set[int] = set()
        for target in range(len(candidates)):
            others = [j for j in range(len(candidates)) if j != target]
            equality_matrix = np.vstack((points[:, others], np.ones(len(others))))
            equality_rhs = np.concatenate((points[:, target], [1.0]))
            result = linprog(
                np.zeros(len(others)),
                A_eq=equality_matrix,
                b_eq=equality_rhs,
                bounds=(0.0, None),
                method="highs",
            )
            if not result.success:
                numerical_vertices.add(target)

        predicted_vertices = {
            j
            for j, (_, k, _) in enumerate(candidates)
            if all(value == 0 for value in k) or sum(k) == h
        }
        self.assertEqual(numerical_vertices, predicted_vertices)
        self.assertEqual(len(numerical_vertices), 21)

    def test_hausdorff_distance_is_exact_for_every_q(self) -> None:
        q, h = 4, 3
        bound = q / (q + h)
        for active_tuple in itertools.chain.from_iterable(
            itertools.combinations(range(q), size) for size in range(q + 1)
        ):
            active = set(active_tuple)
            for active_partition in set_partitions(tuple(sorted(active))):
                partial = np.zeros((q, q))
                for block in active_partition:
                    idx = np.asarray(block)
                    partial[np.ix_(idx, idx)] = 1.0 / len(block)
                inactive = tuple(i for i in range(q) if i not in active)
                if inactive:
                    full_partition = active_partition + (inactive,)
                    allocation = (0,) * len(active_partition) + (h,)
                else:
                    full_partition = active_partition
                    allocation = (0,) * len(active_partition)
                approximation = scenario_matrix(q, full_partition, allocation)
                self.assertLessEqual(np.linalg.norm(approximation - partial, 2), bound + 2e-12)

        # The missing endpoint 0 is exactly q/(q+h) away from K_h.  The common
        # Rayleigh direction u=1/sqrt(q) proves the lower bound for every
        # candidate, hence also for every convex combination.
        u = np.ones(q) / math.sqrt(q)
        for _, _, x in scenarios(q, h):
            self.assertGreaterEqual(float(u @ x @ u), bound - 2e-12)
        one_block = scenario_matrix(q, (tuple(range(q)),), (h,))
        self.assertAlmostEqual(np.linalg.norm(one_block, 2), bound, places=12)

    def test_two_port_one_hidden_improvement(self) -> None:
        q, h = 2, 1
        p_average = np.ones((2, 2)) / 2.0
        c = p_average
        q_map = (2.0 / 3.0) * p_average
        values = {
            (partition, k): schatten(canonical_error(c, q_map, x), math.inf)
            for partition, k, x in scenarios(q, h)
        }
        finite_radius = max(values.values())
        expected = math.sqrt(2.0) / 3.0
        self.assertAlmostEqual(finite_radius, expected, places=12)
        self.assertAlmostEqual(values[(((0, 1),), (1,))], expected, places=12)
        self.assertAlmostEqual(values[(((0,), (1,)), (0, 0))], 1.0 / 3.0, places=12)
        self.assertAlmostEqual(values[(((0, 1),), (0,))], 1.0 / 3.0, places=12)
        self.assertAlmostEqual(values[(((0,), (1,)), (1, 0))], math.sqrt(7.0) / 6.0, places=12)
        self.assertAlmostEqual(values[(((0,), (1,)), (0, 1))], math.sqrt(7.0) / 6.0, places=12)

        # The X=(2/3)P scenario contributes an unavoidable hidden-input term
        # C(X-X^2)C^T=(2/9)P, proving global optimality of this Q.
        x_hard = (2.0 / 3.0) * p_average
        intrinsic = c @ (x_hard - x_hard @ x_hard) @ c.T
        self.assertAlmostEqual(math.sqrt(np.linalg.eigvalsh(intrinsic).max()), expected, places=12)
        self.assertLess(finite_radius, 0.5)
        self.assertAlmostEqual((0.5 - finite_radius) / 0.5, 1.0 - 2.0 * math.sqrt(2.0) / 3.0, places=12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
