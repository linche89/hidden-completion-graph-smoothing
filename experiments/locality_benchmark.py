"""Small reproducible benchmark for spatial locality in block WLS.

The code is intentionally dependency-light (NumPy, SciPy, NetworkX, psutil), never
forms a dense inverse, and keeps the Schur complement as a small-instance oracle.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import psutil
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def make_graph(kind: str, n: int, seed: int) -> nx.Graph:
    if kind == "path":
        graph = nx.path_graph(n)
    elif kind == "cycle":
        graph = nx.cycle_graph(n)
    elif kind == "grid":
        side = max(2, math.ceil(math.sqrt(n)))
        graph = nx.convert_node_labels_to_integers(nx.grid_2d_graph(side, side))
    elif kind == "torus":
        side = max(3, math.ceil(math.sqrt(n)))
        graph = nx.convert_node_labels_to_integers(
            nx.grid_2d_graph(side, side, periodic=True)
        )
    elif kind == "random_regular":
        degree = min(4, n - 1)
        if (n * degree) % 2:
            n += 1
        graph = nx.random_regular_graph(degree, n, seed=seed)
    elif kind == "barbell":
        clique = max(2, n // 4)
        bridge = max(0, n - 2 * clique)
        graph = nx.barbell_graph(clique, bridge)
    elif kind == "cover":
        # A random lift of a fixed four-node base graph.  It creates genuine cycles
        # whose local balls can agree while their global covers differ.
        base_edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]
        sheets = max(1, math.ceil(n / 4))
        rng = np.random.default_rng(seed)
        graph = nx.Graph()
        graph.add_nodes_from(range(4 * sheets))
        for u, v in base_edges:
            perm = rng.permutation(sheets)
            for sheet in range(sheets):
                graph.add_edge(u * sheets + sheet, v * sheets + int(perm[sheet]))
    else:
        raise ValueError(f"unknown graph family: {kind}")

    if not nx.is_connected(graph):
        components = [list(c) for c in nx.connected_components(graph)]
        for left, right in zip(components, components[1:]):
            graph.add_edge(left[0], right[0])
    return nx.convert_node_labels_to_integers(graph)


def bad_region(graph: nx.Graph, fraction: float) -> set[int]:
    target = max(0, min(graph.number_of_nodes(), round(fraction * graph.number_of_nodes())))
    if target == 0:
        return set()
    return set(list(nx.bfs_tree(graph, source=0).nodes())[:target])


def random_orthogonal(dim: int, rng: np.random.Generator) -> np.ndarray:
    q, r = np.linalg.qr(rng.standard_normal((dim, dim)))
    signs = np.sign(np.diag(r))
    signs[signs == 0] = 1.0
    return q * signs


def synthesize_block_wls(config: dict[str, Any]) -> dict[str, Any]:
    graph = make_graph(config["graph"], int(config["n"]), int(config["seed"]))
    dim = int(config["block_dim"])
    rng = np.random.default_rng(int(config["seed"]))
    bad = bad_region(graph, float(config["bad_fraction"]))
    n = graph.number_of_nodes()

    # Keep factor support and initial data holders separate.  They coincide for
    # node measurements, but need not coincide for edge measurements.
    row_blocks: list[tuple[np.ndarray, list[int], list[int]]] = []
    anchor_weight = float(config["anchor_weight"])
    edge_weight = float(config["edge_weight"])
    anchor_fraction = float(config.get("anchor_fraction", 1.0))
    if not 0.0 < anchor_fraction <= 1.0:
        raise ValueError("anchor_fraction must lie in (0, 1]")
    anchor_count = max(1, min(n, round(anchor_fraction * n)))
    anchor_rng = np.random.default_rng(int(config["seed"]) + 101)
    anchors = set(int(v) for v in anchor_rng.choice(n, size=anchor_count, replace=False))

    for node in sorted(anchors):
        scale = float(config["bad_anchor_scale"]) if node in bad else 1.0
        block = np.zeros((dim, n * dim))
        block[:, node * dim : (node + 1) * dim] = math.sqrt(anchor_weight * scale) * np.eye(dim)
        row_blocks.append((block, [node], [node]))

    for u, v in sorted(graph.edges()):
        scale = float(config["bad_edge_scale"]) if (u in bad or v in bad) else 1.0
        weight = math.sqrt(edge_weight * scale)
        rotation = random_orthogonal(dim, rng)
        block = np.zeros((dim, n * dim))
        block[:, u * dim : (u + 1) * dim] = weight * np.eye(dim)
        block[:, v * dim : (v + 1) * dim] = -weight * rotation
        holder_mode = str(config.get("edge_holder_mode", "both"))
        if holder_mode == "both":
            holders = [u, v]
        elif holder_mode == "lower":
            holders = [min(u, v)]
        elif holder_mode == "upper":
            holders = [max(u, v)]
        else:
            raise ValueError("edge_holder_mode must be one of: both, lower, upper")
        row_blocks.append((block, [u, v], holders))

    h_dense = np.vstack([block for block, _, _ in row_blocks])
    h = sp.csr_matrix(h_dense)
    measurement_variable_support = [
        support for block, support, _ in row_blocks for _ in range(block.shape[0])
    ]
    measurement_holders = [
        holders for block, _, holders in row_blocks for _ in range(block.shape[0])
    ]
    x_true = rng.standard_normal(n * dim)
    noise = float(config["noise_std"]) * rng.standard_normal(h.shape[0])
    z = h @ x_true + noise
    j = (h.T @ h).tocsr()
    rhs = np.asarray(h.T @ z).ravel()
    return {
        "graph": graph,
        "bad_nodes": bad,
        "anchor_nodes": anchors,
        "H": h,
        "J": j,
        "rhs": rhs,
        "z": z,
        "x_true": x_true,
        "measurement_variable_support": measurement_variable_support,
        "measurement_holders": measurement_holders,
    }


def nodes_in_ball(graph: nx.Graph, root: int, radius: int) -> list[int]:
    return sorted(nx.single_source_shortest_path_length(graph, root, cutoff=radius))


def coordinate_index(nodes: list[int], dim: int) -> np.ndarray:
    return np.asarray([node * dim + k for node in nodes for k in range(dim)], dtype=int)


def root_index(root: int, dim: int) -> np.ndarray:
    return np.arange(root * dim, (root + 1) * dim, dtype=int)


def local_dirichlet(j: sp.csr_matrix, rhs: np.ndarray, graph: nx.Graph, root: int,
                    radius: int, dim: int) -> np.ndarray:
    ball = nodes_in_ball(graph, root, radius)
    idx = coordinate_index(ball, dim)
    local = spla.spsolve(j[idx][:, idx].tocsc(), rhs[idx])
    position = ball.index(root)
    return np.asarray(local[position * dim : (position + 1) * dim])


def schur_oracle(j: sp.csr_matrix, rhs: np.ndarray, graph: nx.Graph, root: int,
                 radius: int, dim: int) -> np.ndarray:
    """Exact root value after eliminating the exterior; verification only."""
    ball = nodes_in_ball(graph, root, radius)
    inside = coordinate_index(ball, dim)
    outside = np.setdiff1d(np.arange(j.shape[0]), inside, assume_unique=True)
    if outside.size == 0:
        reduced = spla.spsolve(j[inside][:, inside].tocsc(), rhs[inside])
    else:
        j_bb = j[inside][:, inside].tocsc()
        j_bo = j[inside][:, outside].tocsc()
        j_ob = j[outside][:, inside].tocsc()
        factor = spla.splu(j[outside][:, outside].tocsc())
        eliminated_columns = factor.solve(j_ob.toarray())
        eliminated_rhs = factor.solve(rhs[outside])
        schur = j_bb.toarray() - j_bo @ eliminated_columns
        reduced_rhs = rhs[inside] - j_bo @ eliminated_rhs
        reduced = np.linalg.solve(schur, reduced_rhs)
    position = ball.index(root)
    return np.asarray(reduced[position * dim : (position + 1) * dim])


def spectral_bounds(j: sp.csr_matrix) -> tuple[float, float]:
    smallest = float(spla.eigsh(j, k=1, which="SA", return_eigenvectors=False, tol=1e-5)[0])
    largest = float(spla.eigsh(j, k=1, which="LA", return_eigenvectors=False, tol=1e-5)[0])
    if smallest <= 0:
        raise RuntimeError(f"information matrix is not positive definite: lambda_min={smallest}")
    return smallest, largest


def chebyshev_inverse_apply(j: sp.csr_matrix, rhs: np.ndarray, degree: int,
                            lambda_min: float, lambda_max: float) -> np.ndarray:
    """Apply a fixed degree-k Chebyshev approximation to 1/x."""
    if degree == 0:
        return rhs * (2.0 / (lambda_min + lambda_max))
    midpoint = 0.5 * (lambda_min + lambda_max)
    halfwidth = 0.5 * (lambda_max - lambda_min)
    if halfwidth == 0:
        return rhs / midpoint
    coefficients = np.polynomial.chebyshev.chebinterpolate(
        lambda t: 1.0 / (midpoint + halfwidth * t), degree
    )
    scaled = (j - midpoint * sp.eye(j.shape[0], format="csr")) * (1.0 / halfwidth)
    t_prev = rhs.copy()
    result = coefficients[0] * t_prev
    t_curr = scaled @ rhs
    result = result + coefficients[1] * t_curr
    for k in range(2, degree + 1):
        t_next = 2.0 * (scaled @ t_curr) - t_prev
        result = result + coefficients[k] * t_next
        t_prev, t_curr = t_curr, t_next
    return np.asarray(result)


def sampled_operator_tails(
    factor_t: spla.SuperLU,
    h: sp.csr_matrix,
    measurement_holders: list[list[int]],
    graph: nx.Graph,
    roots: np.ndarray,
    radius: int,
    dim: int,
) -> dict[int, dict[str, float]]:
    """Compute exact sampled block-row tails using transpose solves, never J^{-1}."""
    cache: dict[int, dict[str, float]] = {}
    measurement_holder_sets = [set(s) for s in measurement_holders]
    for root in sorted(set(int(r) for r in roots)):
        selectors = np.zeros((h.shape[1], dim))
        selectors[root_index(root, dim), np.arange(dim)] = 1.0
        solved = factor_t.solve(selectors)
        inverse_block_row = solved.T
        ball = set(nodes_in_ball(graph, root, radius))

        state_outside = [node for node in graph.nodes if node not in ball]
        state_outside_idx = coordinate_index(state_outside, dim)
        state_tail = inverse_block_row[:, state_outside_idx]

        # T_i = E_i J^{-1} H^T, so H @ solved is T_i^T.
        wls_block_row = np.asarray(h @ solved).T
        measurement_outside = np.asarray(
            [not bool(holders & ball) for holders in measurement_holder_sets], dtype=bool
        )
        wls_tail = wls_block_row[:, measurement_outside]
        cache[root] = {
            "state_rhs_tail_2to2": float(np.linalg.norm(state_tail, ord=2)) if state_tail.size else 0.0,
            "raw_wls_tail_2to2": float(np.linalg.norm(wls_tail, ord=2)) if wls_tail.size else 0.0,
        }
    return cache


def dkw_interval(sample_fraction: float, sample_count: int, alpha: float) -> list[float]:
    halfwidth = math.sqrt(math.log(2.0 / alpha) / (2.0 * sample_count))
    return [max(0.0, sample_fraction - halfwidth), min(1.0, sample_fraction + halfwidth)]


def benchmark(config: dict[str, Any]) -> dict[str, Any]:
    process = psutil.Process(os.getpid())
    rss_start = process.memory_info().rss
    start = time.perf_counter()
    instance = synthesize_block_wls(config)
    graph: nx.Graph = instance["graph"]
    h: sp.csr_matrix = instance["H"]
    j: sp.csr_matrix = instance["J"]
    rhs: np.ndarray = instance["rhs"]
    dim = int(config["block_dim"])

    factor_start = time.perf_counter()
    factor = spla.splu(j.tocsc())
    central = factor.solve(rhs)
    factor_seconds = time.perf_counter() - factor_start
    relative_residual = float(np.linalg.norm(j @ central - rhs) / max(np.linalg.norm(rhs), 1e-30))
    lambda_min, lambda_max = spectral_bounds(j)

    rng = np.random.default_rng(int(config["seed"]) + 1)
    sample_count = int(config["root_samples"])
    # Sampling with replacement makes the standard DKW band directly applicable.
    roots = rng.integers(0, graph.number_of_nodes(), size=sample_count)
    factor_t = spla.splu(j.T.tocsc())
    radii_results: list[dict[str, Any]] = []
    previous_tails: dict[int, float] = {}
    tail_monotonic = True

    for radius in [int(r) for r in config["radii"]]:
        radius_start = time.perf_counter()
        tails_by_root = sampled_operator_tails(
            factor_t,
            h,
            instance["measurement_holders"],
            graph,
            roots,
            radius,
            dim,
        )
        state_tails = np.asarray([tails_by_root[int(root)]["state_rhs_tail_2to2"] for root in roots])
        wls_tails = np.asarray([tails_by_root[int(root)]["raw_wls_tail_2to2"] for root in roots])

        dirichlet_errors = []
        schur_errors = []
        for root in sorted(set(int(r) for r in roots)):
            exact = central[root_index(root, dim)]
            local = local_dirichlet(j, rhs, graph, root, radius, dim)
            dirichlet_errors.append(float(np.linalg.norm(local - exact)))
            if j.shape[0] <= int(config["schur_oracle_max_state_dim"]):
                oracle = schur_oracle(j, rhs, graph, root, radius, dim)
                schur_errors.append(float(np.linalg.norm(oracle - exact)))

            current = tails_by_root[root]["raw_wls_tail_2to2"]
            if root in previous_tails and current > previous_tails[root] + 1e-10:
                tail_monotonic = False
            previous_tails[root] = current

        poly = chebyshev_inverse_apply(j, rhs, radius, lambda_min, lambda_max)
        poly_rel_error = float(np.linalg.norm(poly - central) / max(np.linalg.norm(central), 1e-30))

        raw_z_operator_tail_exceedance = []
        for threshold in [float(x) for x in config["failure_thresholds"]]:
            fraction = float(np.mean(wls_tails > threshold))
            raw_z_operator_tail_exceedance.append(
                {
                    "epsilon": threshold,
                    "sample_exceedance_fraction": fraction,
                    "dkw_interval": dkw_interval(
                        fraction, sample_count, float(config["confidence_alpha"])
                    ),
                }
            )

        directed_edges = 2 * graph.number_of_edges()
        message_scalars = radius * directed_edges * dim
        radii_results.append(
            {
                "radius_rounds": radius,
                "message_dimension_per_directed_edge_round": dim,
                "polynomial_bit_hop_proxy": int(message_scalars * int(config["bits_per_scalar"])),
                "sampled_state_rhs_tail": {
                    "mean": float(np.mean(state_tails)),
                    "median": float(np.median(state_tails)),
                    "max": float(np.max(state_tails)),
                },
                "sampled_raw_wls_tail": {
                    "mean": float(np.mean(wls_tails)),
                    "median": float(np.median(wls_tails)),
                    "max": float(np.max(wls_tails)),
                },
                "dirichlet_root_error": {
                    "mean": float(np.mean(dirichlet_errors)),
                    "max": float(np.max(dirichlet_errors)),
                },
                "schur_oracle_root_error_max": (
                    float(np.max(schur_errors)) if schur_errors else None
                ),
                "chebyshev_global_relative_error": poly_rel_error,
                "raw_z_operator_tail_exceedance": raw_z_operator_tail_exceedance,
                "seconds": time.perf_counter() - radius_start,
            }
        )

    rss_end = process.memory_info().rss
    result = {
        "schema_version": 2,
        "environment": {
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "networkx": nx.__version__,
            "logical_cpus": os.cpu_count(),
        },
        "config": config,
        "instance": {
            "nodes": graph.number_of_nodes(),
            "edges": graph.number_of_edges(),
            "state_dimension": j.shape[0],
            "measurements": h.shape[0],
            "bad_nodes": sorted(instance["bad_nodes"]),
            "anchor_nodes": sorted(instance["anchor_nodes"]),
            "anchor_count": len(instance["anchor_nodes"]),
            "nodes_without_local_anchor": graph.number_of_nodes() - len(instance["anchor_nodes"]),
            "local_anchor_rank_histogram": {
                "0": graph.number_of_nodes() - len(instance["anchor_nodes"]),
                str(dim): len(instance["anchor_nodes"]),
            },
            "information_matrix_positive_definite": lambda_min > 0.0,
            "edge_measurement_holder_mode": str(config.get("edge_holder_mode", "both")),
            "J_nnz": int(j.nnz),
            "H_nnz": int(h.nnz),
            "lambda_min": lambda_min,
            "lambda_max": lambda_max,
            "condition_number_estimate": lambda_max / lambda_min,
        },
        "central": {
            "factor_and_solve_seconds": factor_seconds,
            "relative_residual": relative_residual,
        },
        "sampling": {
            "root_draws_with_replacement": roots.tolist(),
            "sample_count": sample_count,
            "confidence_alpha": float(config["confidence_alpha"]),
            "dkw_halfwidth": math.sqrt(
                math.log(2.0 / float(config["confidence_alpha"])) / (2.0 * sample_count)
            ),
        },
        "checks": {
            "raw_wls_tail_nonincreasing_on_sampled_unique_roots": tail_monotonic,
            "central_residual_below_1e-9": relative_residual < 1e-9,
            "schur_oracle_below_1e-8": all(
                row["schur_oracle_root_error_max"] is None
                or row["schur_oracle_root_error_max"] < 1e-8
                for row in radii_results
            ),
        },
        "radii": radii_results,
        "resources": {
            "rss_start_bytes": rss_start,
            "rss_end_bytes": rss_end,
            "wall_seconds": time.perf_counter() - start,
        },
    }
    if not all(result["checks"].values()):
        raise RuntimeError(f"smoke correctness check failed: {result['checks']}")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = benchmark(config)
    encoded = json.dumps(_jsonable(result), indent=2, sort_keys=True)
    if args.output is None:
        print(encoded)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded + "\n", encoding="utf-8")
        print(args.output)


if __name__ == "__main__":
    main()
