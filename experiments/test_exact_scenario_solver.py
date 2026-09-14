"""Regression tests for the exact hidden-completion scenario solver."""

from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from hidden_completion import (
    canonical_error,
    finite_scenarios,
    partial_partition_scenarios,
    r_hop_allowed_mask,
    scenario_count,
    scenario_error_gram,
    solve_exact_spectral,
    vertex_count,
)
from hidden_completion.scenarios import bell_number, set_partitions


ROOT = Path(__file__).resolve().parents[1]


class ScenarioEnumerationTests(unittest.TestCase):
    def test_set_partitions_are_unique_and_have_bell_counts(self) -> None:
        for q, expected in enumerate((1, 1, 2, 5, 15, 52)):
            partitions = list(set_partitions(tuple(range(q))))
            self.assertEqual(len(partitions), expected)
            self.assertEqual(len(set(partitions)), expected)
            self.assertEqual(bell_number(q), expected)

    def test_finite_counts_pruning_and_no_duplicate_matrices(self) -> None:
        scenarios = finite_scenarios(3, 2)
        vertices = finite_scenarios(3, 2, vertices_only=True)
        self.assertEqual(len(scenarios), 31)
        self.assertEqual(len(vertices), 21)
        self.assertEqual(scenario_count(3, 2), 31)
        self.assertEqual(vertex_count(3, 2), 21)
        encoded = {tuple(np.round(item.matrix, 14).flat) for item in scenarios}
        self.assertEqual(len(encoded), len(scenarios))
        self.assertTrue(all(not any(item.allocation or ()) or sum(item.allocation or ()) == 2 for item in vertices))

    def test_small_q_h_grid_matches_closed_counts(self) -> None:
        for q in range(1, 5):
            for h in range(4):
                with self.subTest(q=q, h=h):
                    candidates = finite_scenarios(q, h)
                    vertices = finite_scenarios(q, h, vertices_only=True)
                    self.assertEqual(len(candidates), scenario_count(q, h))
                    self.assertEqual(len(vertices), vertex_count(q, h))
                    self.assertEqual(
                        len({tuple(np.round(item.matrix, 14).flat) for item in candidates}),
                        len(candidates),
                    )
                    self.assertTrue(all(item.is_vertex for item in vertices))

    def test_unbounded_count_is_bell_q_plus_one(self) -> None:
        for q in range(1, 6):
            scenarios = partial_partition_scenarios(q)
            self.assertEqual(len(scenarios), bell_number(q + 1))
            self.assertEqual(len(scenarios), scenario_count(q, None))
            self.assertTrue(all(np.allclose(item.matrix @ item.matrix, item.matrix) for item in scenarios))

    def test_canonical_error_has_the_proved_gram(self) -> None:
        rng = np.random.default_rng(20260914)
        c_matrix = rng.normal(size=(2, 3))
        q_matrix = rng.normal(size=(2, 3))
        x_matrix = finite_scenarios(3, 2)[17].matrix
        error = canonical_error(c_matrix, q_matrix, x_matrix)
        np.testing.assert_allclose(
            error @ error.T,
            scenario_error_gram(c_matrix, q_matrix, x_matrix),
            rtol=2e-12,
            atol=2e-12,
        )


class LocalRuleTests(unittest.TestCase):
    def test_r_hop_mask(self) -> None:
        adjacency = np.array(
            [
                [0, 1, 0, 0],
                [1, 0, 1, 0],
                [0, 1, 0, 1],
                [0, 0, 1, 0],
            ]
        )
        expected = np.array([[True, True, False, False], [False, False, True, True]])
        np.testing.assert_array_equal(r_hop_allowed_mask(adjacency, [0, 3], 1), expected)

    def test_zero_r_hop_and_shared_constraints_are_enforced(self) -> None:
        c_matrix = np.diag([1.0, 1.5, 2.0])
        local_rule = {
            "r_hop": {
                "adjacency": [[0, 1, 0], [1, 0, 1], [0, 1, 0]],
                "output_nodes": [0, 1, 2],
                "radius": 1,
            },
            "zeros": [[1, 0]],
            "shared_coefficients": [[[0, 0], [2, 2]]],
        }
        result = solve_exact_spectral(c_matrix, h=0, local_rule=local_rule)
        q_matrix = result.q_matrix
        self.assertAlmostEqual(q_matrix[0, 2], 0.0, places=7)
        self.assertAlmostEqual(q_matrix[1, 0], 0.0, places=7)
        self.assertAlmostEqual(q_matrix[2, 0], 0.0, places=7)
        self.assertAlmostEqual(q_matrix[0, 0], q_matrix[2, 2], places=7)

    def test_unknown_constraint_field_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "unsupported local_rule"):
            solve_exact_spectral(np.eye(2), h=0, local_rule={"zero_patern": []})


class ExactSdpTests(unittest.TestCase):
    def test_finite_two_port_analytic_solution(self) -> None:
        average = np.ones((2, 2)) / 2.0
        result = solve_exact_spectral(average, h=1)
        np.testing.assert_allclose(result.q_matrix, (2.0 / 3.0) * average, atol=8e-7)
        self.assertAlmostEqual(result.radius, math.sqrt(2.0) / 3.0, places=7)
        self.assertEqual(result.total_scenario_count, 5)
        self.assertEqual(result.retained_scenario_count, 5)
        self.assertEqual(
            [item["label"] for item in result.active_scenarios],
            ["finite:{0,1};k=(1)"],
        )

    def test_pruned_and_unpruned_finite_sdps_agree(self) -> None:
        c_matrix = np.array([[1.0, -0.2, 0.4], [0.3, 0.8, -0.5]])
        pruned = solve_exact_spectral(c_matrix, h=2, prune_vertices=True)
        unpruned = solve_exact_spectral(c_matrix, h=2, prune_vertices=False)
        self.assertEqual(pruned.retained_scenario_count, 21)
        self.assertEqual(unpruned.retained_scenario_count, 31)
        self.assertAlmostEqual(pruned.radius, unpruned.radius, places=6)

    def test_unbounded_shared_opcode_analytic_solution(self) -> None:
        c_matrix = np.diag([1.0, 2.0])
        result = solve_exact_spectral(
            c_matrix,
            mode="unbounded",
            local_rule={
                "zeros": [[0, 1], [1, 0]],
                "shared_coefficients": [[[0, 0], [1, 1]]],
            },
        )
        alpha = 4.0 / 3.0 - 2.0 * math.sqrt(10.0) / 15.0
        radius = (10.0 + 2.0 * math.sqrt(10.0)) / 15.0
        np.testing.assert_allclose(result.q_matrix, alpha * np.eye(2), atol=8e-7)
        self.assertAlmostEqual(result.radius, radius, places=7)
        self.assertEqual(result.total_scenario_count, 5)
        labels = {item["label"] for item in result.active_scenarios}
        self.assertEqual(
            labels,
            {
                "partial:{1};inactive=(0)",
                "partial:{0}|{1};inactive=()",
                "partial:{0,1};inactive=()",
            },
        )

    def test_cli_emits_machine_readable_active_scenarios(self) -> None:
        config = ROOT / "experiments" / "configs" / "exact_finite_q2_h1.json"
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "result.json"
            completed = subprocess.run(
                [sys.executable, "-m", "hidden_completion", "solve", str(config), "--output", str(output)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.stdout, "")
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["status"], "optimal")
            self.assertEqual(payload["scenario_counts"]["active"], 1)
            self.assertIn("active_worst_case_scenarios", payload)


if __name__ == "__main__":
    unittest.main(verbosity=2)
