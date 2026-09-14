"""Deterministic regression tests for the P3/P4 graph-smoother pipeline."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from experiments.smoother.algorithms import (
    chebyshev_apply,
    exact_smoother_operator,
    induced_local_operator,
    neumann_apply,
)
from experiments.smoother.graphs import (
    load_pglib_topology,
    load_suitesparse_pattern,
    smoother_matrix,
    synthetic_weighted_graph,
)
from experiments.smoother.paper_data import _export_budget
from experiments.smoother.parameter_sweep import run_parameter_sweep
from experiments.smoother.resources import direct_rule_resources
from experiments.smoother.robust import (
    full_error_operator,
    generate_completion,
    run_robust_benchmark,
    scenario_certificate,
    solve_psd_contraction_relaxation,
)
from experiments.smoother.stress import run_stress_benchmark
from hidden_completion.scenarios import finite_scenarios


ROOT = Path(__file__).resolve().parents[1]


class GraphAndResourceTests(unittest.TestCase):
    def test_required_synthetic_families_are_connected_positive_smoothers(self) -> None:
        for index, family in enumerate(
            ("random_cyclic", "grid", "geometric", "small_world")
        ):
            with self.subTest(family=family):
                graph = synthetic_weighted_graph(
                    family,
                    25,
                    100 + index,
                    weight_min=0.2,
                    weight_max=3.0,
                    options={"geometric_radius": 0.27, "small_world_k": 4},
                )
                np.testing.assert_allclose(graph.adjacency.toarray(), graph.adjacency.toarray().T)
                self.assertGreaterEqual(float(np.min(graph.adjacency.data)), 0.0)
                eigenvalues = np.linalg.eigvalsh(smoother_matrix(graph).toarray())
                self.assertGreaterEqual(eigenvalues[0], 1.0 - 1e-10)

    def test_direct_rule_accounting_matches_path_distances(self) -> None:
        adjacency = np.array(
            [[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]],
            dtype=float,
        )
        q_matrix = np.array([[2.0, 0.0, 3.0, 0.0], [0.0, 1.0, 0.0, 4.0]])
        resources = direct_rule_resources(q_matrix, adjacency, [0, 3])
        self.assertEqual(resources["rounds"], 2)
        self.assertEqual(resources["scalar_hops"], 4)
        self.assertEqual(resources["online_multiplies"], 4)
        self.assertEqual(resources["online_additions"], 2)

    def test_external_loaders_record_transformations(self) -> None:
        pglib = load_pglib_topology(
            ROOT / "datasets/raw/pglib_opf_v23.07/pglib_opf_case14_ieee.m"
        )
        self.assertEqual(pglib.nodes, 14)
        self.assertEqual(pglib.metadata["classification"], "external_application_stress")
        suitesparse = load_suitesparse_pattern(
            ROOT / "datasets/raw/suitesparse_ccby4/bcsstk18.tar.gz"
        )
        self.assertEqual(suitesparse.adjacency.shape[0], suitesparse.adjacency.shape[1])
        self.assertGreater(suitesparse.edges, 0)
        self.assertIn("sparsity", suitesparse.metadata["transformation"])


class FixedGraphAlgorithmTests(unittest.TestCase):
    def test_full_radius_local_truncation_is_exact(self) -> None:
        graph = synthetic_weighted_graph("path", 8, 7, weight_min=1.0, weight_max=1.0)
        exact = exact_smoother_operator(graph)
        local = induced_local_operator(graph.adjacency, 7).toarray()
        np.testing.assert_allclose(local, exact, rtol=2e-12, atol=2e-12)

    def test_polynomials_converge_on_a_path(self) -> None:
        graph = synthetic_weighted_graph("path", 12, 8, weight_min=1.0, weight_max=1.0)
        system = smoother_matrix(graph)
        exact = exact_smoother_operator(graph)
        upper = 1.0 + 2.0 * np.max(np.asarray(graph.adjacency.sum(axis=1)).ravel())
        identity = np.eye(graph.nodes)
        neumann_0 = neumann_apply(system, identity, 0, upper)
        neumann_20 = neumann_apply(system, identity, 20, upper)
        chebyshev_0 = chebyshev_apply(system, identity, 0, 1.0, upper)
        chebyshev_12 = chebyshev_apply(system, identity, 12, 1.0, upper)
        self.assertLess(np.linalg.norm(neumann_20 - exact, ord=2), np.linalg.norm(neumann_0 - exact, ord=2))
        self.assertLess(np.linalg.norm(chebyshev_12 - exact, ord=2), np.linalg.norm(chebyshev_0 - exact, ord=2))


class RobustDesignTests(unittest.TestCase):
    def test_physical_completion_is_below_exact_finite_certificate(self) -> None:
        average = np.ones((2, 2)) / 2.0
        q_matrix = (2.0 / 3.0) * average
        sample = generate_completion(
            2,
            1,
            "random_cyclic",
            17,
            weight_min=0.01,
            weight_max=100.0,
        )
        physical = np.linalg.norm(full_error_operator(average, q_matrix, sample), ord=2)
        certificate = scenario_certificate(average, q_matrix, finite_scenarios(2, 1))
        self.assertLessEqual(physical, certificate["certified_radius"] + 1e-10)

    def test_contraction_relaxation_recovers_free_identity_center(self) -> None:
        result = solve_psd_contraction_relaxation(
            np.eye(2), projection_audit_samples=32, seed=9
        )
        np.testing.assert_allclose(result["Q"], 0.5 * np.eye(2), atol=2e-6)
        self.assertAlmostEqual(result["certified_radius"], 0.5, places=6)
        self.assertGreaterEqual(result["random_projection_audit"]["slack_to_certificate"], -2e-6)

    def test_robust_benchmark_runs_all_certificate_baselines(self) -> None:
        average = np.ones((2, 2)) / 2.0
        result = run_robust_benchmark(
            {
                "q": 2,
                "h": 1,
                "C": average.tolist(),
                "seed": 31,
                "solver": {"name": "CLARABEL"},
                "random_sampling": {
                    "families": ["random_cyclic", "grid"],
                    "training_count": 4,
                    "test_count": 4,
                    "input_draws": 4,
                    "training_weight_min": 0.1,
                    "training_weight_max": 10.0,
                    "test_weight_min": 0.01,
                    "test_weight_max": 100.0
                },
                "projection_audit_samples": 16,
                "run_pruning_ablation": True,
                "budget_sweep": {"h_values": [0, 1], "solve_optima": False}
            }
        )
        self.assertTrue(all(result["checks"].values()))
        self.assertEqual(
            set(result["methods"]),
            {
                "exact_finite_h",
                "exact_unbounded",
                "random_completion_sampling",
                "psd_contraction_relaxation",
                "endpoint_only_no_go",
            },
        )
        self.assertEqual(len(result["budget_sweep"]["rows"]), 2)
        self.assertAlmostEqual(
            result["methods"]["exact_finite_h"]["exact_finite_h_evaluation"][
                "certified_radius"
            ],
            math.sqrt(2.0) / 3.0,
            places=6,
        )

    def test_parameter_sweep_keeps_physical_errors_below_certificates(self) -> None:
        result = run_parameter_sweep(
            {
                "q": 2,
                "C": [[0.5, 0.5]],
                "hidden_budgets": [0, 1],
                "weight_profiles": [
                    {"name": "test", "min": 0.1, "max": 10.0}
                ],
                "mean_degrees": [2.0],
                "replicates": 1,
                "family": "random_cyclic",
                "seed": 73,
                "solver": {"name": "CLARABEL"},
            }
        )
        self.assertTrue(
            result["checks"]["all_physical_errors_within_exact_certificate"]
        )
        self.assertTrue(result["checks"]["all_requested_cells_present"])
        self.assertEqual(result["sweep_axes"]["total_physical_completions"], 2)
        self.assertEqual(
            [group["h"] for group in result["by_hidden_budget"]], [0, 1]
        )


class StressBenchmarkTests(unittest.TestCase):
    def test_stress_runner_separates_fixed_graph_semantics(self) -> None:
        result = run_stress_benchmark(
            {
                "seed": 41,
                "rounds": [0, 7],
                "signal_samples": 3,
                "exact_operator_max_n": 16,
                "local_truncation_max_n": 16,
                "cut_root_samples": 4,
                "enable_local_truncation": True,
                "cases": [
                    {
                        "source": "synthetic",
                        "family": "path",
                        "n": 8,
                        "seed": 42,
                        "weight_min": 1.0,
                        "weight_max": 1.0,
                    }
                ],
            }
        )
        self.assertTrue(result["checks"]["all_cases_passed"])
        case = result["cases"][0]
        full_radius = case["round_sweep"][-1]["methods"]["induced_local_truncation"]
        self.assertLess(
            full_radius["operator_audit"]["fixed_instance_operator_2_norm_error"], 1e-10
        )
        self.assertIn("not WLS/Kalman", result["semantics"]["empirical_mse"])


class PaperDataTests(unittest.TestCase):
    def test_budget_export_preserves_the_exact_rate_data(self) -> None:
        robust = {
            "budget_sweep": {
                "rows": [
                    {
                        "h": 2,
                        "exact_hausdorff_distance": 0.5,
                        "fixed_Q_squared_risk_gap": 0.5 - 2e-13,
                        "candidate_scenarios": 9,
                        "retained_vertices": 6,
                        "optimized_finite_radius": math.sqrt(0.5),
                    }
                ]
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            _export_budget(robust, output)
            with (output / "budget_pareto.csv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["h"], "2")
        self.assertEqual(rows[0]["candidates"], "9")
        self.assertAlmostEqual(float(rows[0]["hausdorff"]), 0.5)
        self.assertAlmostEqual(float(rows[0]["risk_gap"]), 0.5, places=12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
