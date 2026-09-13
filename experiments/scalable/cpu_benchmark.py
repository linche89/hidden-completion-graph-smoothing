"""CPU correctness and scale benchmark for archived PGLib DC-WLS cases."""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import time
from pathlib import Path
from typing import Any

import numpy as np
import psutil
import scipy
import scipy.sparse.linalg as spla

from .pglib_dc import load_dc_wls
from .wls_core import sampled_scalar_tails


DEFAULT_CASES = [
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case118_ieee.m",
    "pglib_opf_case300_ieee.m",
    "pglib_opf_case1354_pegase.m",
    "pglib_opf_case2869_pegase.m",
    "pglib_opf_case10000_goc.m",
]


def dkw_halfwidth(sample_count: int, alpha: float) -> float:
    return math.sqrt(math.log(2.0 / alpha) / (2.0 * sample_count))


def benchmark_case(path: Path, seed: int, root_samples: int,
                   radii: list[int]) -> dict[str, Any]:
    process = psutil.Process(os.getpid())
    rss_points = {"start": process.memory_info().rss}
    instance, timings = load_dc_wls(
        path,
        seed=seed,
        flow_fraction=1.0,
        injection_fraction=0.25,
        angle_fraction=0.05,
        flow_holder="from",
        truth_mode="synthetic_dc",
    )
    rss_points["after_build"] = process.memory_info().rss

    factor_start = time.perf_counter()
    factor = spla.splu(instance.j.tocsc())
    factor_seconds = time.perf_counter() - factor_start
    rss_points["after_factor"] = process.memory_info().rss

    solve_start = time.perf_counter()
    estimate = factor.solve(instance.rhs)
    solve_seconds = time.perf_counter() - solve_start
    residual = float(
        np.linalg.norm(instance.j @ estimate - instance.rhs)
        / max(np.linalg.norm(instance.rhs), 1e-300)
    )

    rng = np.random.default_rng(seed + 991)
    eligible = np.asarray(list(instance.node_to_state), dtype=int)
    sampled_roots = rng.choice(eligible, size=root_samples, replace=True)
    tails_start = time.perf_counter()
    tails: dict[int, dict[int, dict[str, float]]] = {}
    for radius in radii:
        tails[radius] = sampled_scalar_tails(instance, factor, sampled_roots, radius)
    tail_seconds = time.perf_counter() - tails_start
    rss_points["after_tails"] = process.memory_info().rss

    monotone = True
    for root in set(int(x) for x in sampled_roots):
        previous = math.inf
        for radius in sorted(radii):
            value = tails[radius][root]["whitened_z_tail_l2"]
            monotone = monotone and value <= previous + 1e-10 * max(1.0, previous)
            previous = value

    tail_summary = []
    thresholds = (1e-4, 1e-3, 1e-2)
    for radius in radii:
        whitened = np.asarray(
            [tails[radius][int(root)]["whitened_z_tail_l2"] for root in sampled_roots]
        )
        physical_raw = np.asarray(
            [tails[radius][int(root)]["physical_raw_z_tail_l2"] for root in sampled_roots]
        )
        rhs_tail = np.asarray(
            [tails[radius][int(root)]["state_rhs_tail_l2"] for root in sampled_roots]
        )
        tail_summary.append(
            {
                "radius": radius,
                "whitened_z_tail_mean": float(np.mean(whitened)),
                "whitened_z_tail_median": float(np.median(whitened)),
                "whitened_z_tail_max": float(np.max(whitened)),
                "physical_raw_z_tail_mean": float(np.mean(physical_raw)),
                "physical_raw_z_tail_median": float(np.median(physical_raw)),
                "physical_raw_z_tail_max": float(np.max(physical_raw)),
                "state_rhs_tail_mean": float(np.mean(rhs_tail)),
                "whitened_z_tail_exceedance": {
                    str(threshold): float(np.mean(whitened > threshold))
                    for threshold in thresholds
                },
            }
        )

    diameter_check = None
    if instance.graph.number_of_nodes() <= 118:
        diameter = max(
            max(lengths.values())
            for component in (
                instance.graph.subgraph(c).copy()
                for c in __import__("networkx").connected_components(instance.graph)
            )
            for _, lengths in __import__("networkx").all_pairs_shortest_path_length(component)
        )
        full_tails = sampled_scalar_tails(instance, factor, sampled_roots, int(diameter))
        diameter_check = {
            "diameter": int(diameter),
            "max_whitened_z_tail": float(
                max(v["whitened_z_tail_l2"] for v in full_tails.values())
            ),
            "max_physical_raw_z_tail": float(
                max(v["physical_raw_z_tail_l2"] for v in full_tails.values())
            ),
            "max_state_rhs_tail": float(
                max(v["state_rhs_tail_l2"] for v in full_tails.values())
            ),
        }

    return {
        "case": path.name,
        "source_sha256": instance.metadata["source_sha256"],
        "nodes": instance.graph.number_of_nodes(),
        "state_dimension": instance.j.shape[0],
        "graph_edges": instance.graph.number_of_edges(),
        "measurements": instance.h.shape[0],
        "H_nnz": int(instance.h.nnz),
        "J_nnz": int(instance.j.nnz),
        "reference_bus_ids": instance.metadata["reference_bus_ids"],
        "measurement_counts": instance.metadata["counts_by_kind"],
        "flow_holder_semantics": instance.metadata["flow_holder_semantics"],
        "timing_seconds": {
            **timings,
            "factorization": factor_seconds,
            "central_solve": solve_seconds,
            "sampled_tail_solves_and_masks": tail_seconds,
        },
        "rss_bytes_at_checkpoints": rss_points,
        "rss_checkpoint_max": max(rss_points.values()),
        "central_relative_residual": residual,
        "sampled_roots_with_replacement": sampled_roots.tolist(),
        "root_sample_count": root_samples,
        "dkw_95_halfwidth_per_fixed_radius": dkw_halfwidth(root_samples, 0.05),
        "whitened_z_tail_monotone_over_sorted_radii": bool(monotone),
        "diameter_zero_tail_check": diameter_check,
        "tail_summary": tail_summary,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", nargs="*", default=DEFAULT_CASES)
    parser.add_argument("--seed", type=int, default=20260913)
    parser.add_argument("--root-samples", type=int, default=16)
    parser.add_argument("--radii", nargs="*", type=int, default=[0, 1, 2, 4])
    args = parser.parse_args()

    started = time.perf_counter()
    results: list[dict[str, Any]] = []
    for offset, case in enumerate(args.cases):
        path = args.data_dir / case
        case_started = time.perf_counter()
        try:
            result = benchmark_case(path, args.seed + offset, args.root_samples, args.radii)
            result["case_total_seconds"] = time.perf_counter() - case_started
            result["status"] = "ok"
        except Exception as exc:  # preserve failures in the machine-readable audit
            result = {
                "case": case,
                "status": "error",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "case_total_seconds": time.perf_counter() - case_started,
            }
        results.append(result)

    payload = {
        "schema_version": 1,
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "logical_cpus": os.cpu_count(),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
        },
        "benchmark_semantics": {
            "matrix_assembly": "direct COO/CSR; no dense H/J",
            "central_solver": "SciPy SuperLU factorization and triangular solve",
            "tail_solver": "sampled transpose sparse solves; no complete inverse",
            "tail_input": (
                "both physical raw z and dimensionless whitened z=R^{-1/2}(z-offset); "
                "exceedance curves use whitened z"
            ),
            "dkw_scope": "95% per fixed graph and fixed radius, not simultaneous over radii",
            "mpi": False,
        },
        "total_seconds": time.perf_counter() - started,
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
