"""Explicit communication, arithmetic, storage, and sampled-RSS accounting."""

from __future__ import annotations

from collections import deque
from contextlib import AbstractContextManager
from dataclasses import dataclass, field
import os
import threading
import time
from typing import Any, Iterable

import numpy as np
from numpy.typing import NDArray
import psutil
import scipy.sparse as sp


Array = NDArray[np.float64]


def csr_storage_bytes(matrix: sp.spmatrix) -> int:
    csr = sp.csr_matrix(matrix)
    return int(csr.data.nbytes + csr.indices.nbytes + csr.indptr.nbytes)


def bfs_distances(
    adjacency: sp.spmatrix,
    source: int,
    radius: int | None = None,
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Return deterministic distances and parents in a sorted CSR graph."""

    graph = sp.csr_matrix(adjacency)
    graph.sort_indices()
    n = graph.shape[0]
    if not 0 <= source < n:
        raise ValueError("source lies outside the graph")
    if radius is not None and (not isinstance(radius, int) or radius < 0):
        raise ValueError("radius must be a nonnegative integer or None")
    distance = np.full(n, -1, dtype=np.int64)
    parent = np.full(n, -1, dtype=np.int64)
    distance[source] = 0
    parent[source] = source
    queue: deque[int] = deque([source])
    while queue:
        vertex = queue.popleft()
        if radius is not None and distance[vertex] >= radius:
            continue
        begin, end = graph.indptr[vertex], graph.indptr[vertex + 1]
        for raw_neighbor in graph.indices[begin:end]:
            neighbor = int(raw_neighbor)
            if distance[neighbor] < 0:
                distance[neighbor] = distance[vertex] + 1
                parent[neighbor] = vertex
                queue.append(neighbor)
    return distance, parent


def ball_nodes(adjacency: sp.spmatrix, root: int, radius: int) -> NDArray[np.int64]:
    distance, _ = bfs_distances(adjacency, root, radius)
    return np.flatnonzero(distance >= 0).astype(np.int64)


def _sparse_dot_flops(matrix: sp.spmatrix) -> int:
    csr = sp.csr_matrix(matrix)
    rows_with_entries = int(np.count_nonzero(np.diff(csr.indptr)))
    return int(csr.nnz + max(0, csr.nnz - rows_with_entries))


def direct_rule_resources(
    q_matrix: Array,
    adjacency: sp.spmatrix | Array,
    output_nodes: Iterable[int],
    *,
    scalar_bytes: int = 8,
    index_bytes: int = 4,
    zero_tolerance: float = 1e-10,
    independent_coefficients: int | None = None,
) -> dict[str, Any]:
    """Account for the documented shortest-path-unicast implementation of Q."""

    matrix = np.asarray(q_matrix, dtype=float)
    graph = sp.csr_matrix(adjacency)
    output_tuple = tuple(int(node) for node in output_nodes)
    if matrix.ndim != 2 or matrix.shape[1] != graph.shape[0]:
        raise ValueError("Q columns and communication graph size must agree")
    if len(output_tuple) != matrix.shape[0]:
        raise ValueError("one output node is required per row of Q")
    if scalar_bytes < 1 or index_bytes < 1:
        raise ValueError("byte widths must be positive")

    mask = np.abs(matrix) > zero_tolerance
    scalar_hops = 0
    maximum_distance = 0
    edge_load: dict[tuple[int, int], int] = {}
    unreachable: list[list[int]] = []
    for row, output_node in enumerate(output_tuple):
        distance, parent = bfs_distances(graph, output_node)
        for column in np.flatnonzero(mask[row]):
            column_int = int(column)
            hop_count = int(distance[column_int])
            if hop_count < 0:
                unreachable.append([row, column_int])
                continue
            scalar_hops += hop_count
            maximum_distance = max(maximum_distance, hop_count)
            vertex = column_int
            while vertex != output_node:
                next_vertex = int(parent[vertex])
                edge = (vertex, next_vertex)
                edge_load[edge] = edge_load.get(edge, 0) + 1
                vertex = next_vertex
    if unreachable:
        raise ValueError(f"nonzero Q entries are not communicable: {unreachable}")

    nonzeros = int(np.count_nonzero(mask))
    row_counts = np.count_nonzero(mask, axis=1)
    additions = int(np.sum(np.maximum(row_counts - 1, 0)))
    parameter_count = nonzeros if independent_coefficients is None else int(independent_coefficients)
    return {
        "protocol": "lexicographic_shortest_path_unicast",
        "rounds": maximum_distance,
        "nonzero_dependencies": nonzeros,
        "message_count_scalar_packets": scalar_hops,
        "scalar_hops": scalar_hops,
        "byte_hop": scalar_hops * scalar_bytes,
        "maximum_directed_edge_scalar_load": max(edge_load.values(), default=0),
        "online_multiplies": nonzeros,
        "online_additions": additions,
        "online_flops": nonzeros + additions,
        "coefficient_scalars_stored": nonzeros,
        "independent_coefficient_parameters": parameter_count,
        "sparse_coefficient_storage_bytes": (
            nonzeros * (scalar_bytes + index_bytes) + (matrix.shape[0] + 1) * index_bytes
        ),
        "working_vector_bytes": (matrix.shape[0] + matrix.shape[1]) * scalar_bytes,
        "scalar_bytes": scalar_bytes,
        "claim": "cost of one explicit schedule; not a communication lower bound",
    }


def polynomial_iteration_resources(
    adjacency: sp.spmatrix,
    system: sp.spmatrix,
    rounds: int,
    *,
    right_hand_sides: int = 1,
    scalar_bytes: int = 8,
    recurrence_vectors: int = 3,
) -> dict[str, Any]:
    if rounds < 0 or right_hand_sides < 1:
        raise ValueError("rounds must be nonnegative and right_hand_sides positive")
    graph = sp.csr_matrix(adjacency)
    matrix = sp.csr_matrix(system)
    undirected_edges = int(sp.triu(graph, k=1).nnz)
    directed_messages = 2 * undirected_edges * rounds * right_hand_sides
    spmv_flops = rounds * right_hand_sides * _sparse_dot_flops(matrix)
    recurrence_flops = rounds * right_hand_sides * 4 * matrix.shape[0]
    return {
        "protocol": "neighbor_sparse_polynomial",
        "rounds": rounds,
        "spmv_count": rounds * right_hand_sides,
        "directed_message_count": directed_messages,
        "scalar_hops": directed_messages,
        "byte_hop": directed_messages * scalar_bytes,
        "estimated_spmv_flops": spmv_flops,
        "estimated_recurrence_flops": recurrence_flops,
        "estimated_online_flops": spmv_flops + recurrence_flops,
        "matrix_storage_bytes": csr_storage_bytes(matrix),
        "working_vector_bytes": recurrence_vectors
        * matrix.shape[0]
        * right_hand_sides
        * scalar_bytes,
        "scalar_bytes": scalar_bytes,
    }


def local_collection_resources(
    adjacency: sp.spmatrix,
    roots: Iterable[int],
    radius: int,
    *,
    right_hand_sides: int = 1,
    scalar_bytes: int = 8,
    index_bytes: int = 4,
) -> dict[str, Any]:
    """Account for independently gathering each induced ball at its root."""

    graph = sp.csr_matrix(adjacency)
    root_tuple = tuple(int(root) for root in roots)
    per_root: list[dict[str, int]] = []
    upper = sp.triu(graph, k=1).tocoo()
    edge_list = list(zip(upper.row.tolist(), upper.col.tolist(), strict=True))
    for root in root_tuple:
        distance, _ = bfs_distances(graph, root, radius)
        inside = distance >= 0
        nodes = np.flatnonzero(inside)
        input_scalar_hops = int(np.sum(distance[nodes])) * right_hand_sides
        edge_record_hops = 0
        induced_edges = 0
        for u, v in edge_list:
            if inside[u] and inside[v]:
                induced_edges += 1
                edge_record_hops += int(min(distance[u], distance[v]))
        payload_byte_hop = input_scalar_hops * scalar_bytes + edge_record_hops * (
            scalar_bytes + 2 * index_bytes
        )
        ball_size = int(nodes.size)
        factor_flops = int(round((2.0 / 3.0) * ball_size**3))
        solve_flops = int(2 * ball_size**2 * right_hand_sides)
        per_root.append(
            {
                "root": root,
                "ball_nodes": ball_size,
                "ball_edges": induced_edges,
                "input_scalar_hops": input_scalar_hops,
                "edge_record_hops": edge_record_hops,
                "message_item_hops": input_scalar_hops + edge_record_hops,
                "byte_hop": payload_byte_hop,
                "estimated_dense_factor_flops": factor_flops,
                "estimated_dense_solve_flops": solve_flops,
                "estimated_dense_matrix_bytes": ball_size * ball_size * scalar_bytes,
            }
        )

    def aggregate(key: str, operation: str) -> float | int:
        values = [row[key] for row in per_root]
        if not values:
            return 0
        if operation == "sum":
            return int(sum(values))
        if operation == "max":
            return int(max(values))
        return float(np.mean(values))

    return {
        "protocol": "independent_root_induced_ball_gather",
        "rounds": radius,
        "roots": len(root_tuple),
        "total_message_item_hops": aggregate("message_item_hops", "sum"),
        "total_byte_hop": aggregate("byte_hop", "sum"),
        "mean_byte_hop_per_root": aggregate("byte_hop", "mean"),
        "maximum_byte_hop_per_root": aggregate("byte_hop", "max"),
        "mean_ball_nodes": aggregate("ball_nodes", "mean"),
        "maximum_ball_nodes": aggregate("ball_nodes", "max"),
        "total_estimated_factor_flops": aggregate("estimated_dense_factor_flops", "sum"),
        "total_estimated_solve_flops": aggregate("estimated_dense_solve_flops", "sum"),
        "maximum_dense_matrix_bytes_per_root": aggregate(
            "estimated_dense_matrix_bytes", "max"
        ),
        "claim": "independent collection/factorization reference implementation",
    }


@dataclass(slots=True)
class PeakRssSampler(AbstractContextManager["PeakRssSampler"]):
    """Sample process RSS in a background thread during a code region."""

    interval_seconds: float = 0.005
    baseline_bytes: int = field(init=False, default=0)
    peak_bytes: int = field(init=False, default=0)
    _stop: threading.Event = field(init=False, default_factory=threading.Event)
    _thread: threading.Thread | None = field(init=False, default=None)

    def __enter__(self) -> "PeakRssSampler":
        if self.interval_seconds <= 0:
            raise ValueError("sampling interval must be positive")
        process = psutil.Process(os.getpid())
        self.baseline_bytes = int(process.memory_info().rss)
        self.peak_bytes = self.baseline_bytes
        self._stop.clear()

        def sample() -> None:
            while not self._stop.is_set():
                self.peak_bytes = max(self.peak_bytes, int(process.memory_info().rss))
                self._stop.wait(self.interval_seconds)

        self._thread = threading.Thread(target=sample, name="rss-sampler", daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=max(0.1, 5 * self.interval_seconds))
        self.peak_bytes = max(
            self.peak_bytes, int(psutil.Process(os.getpid()).memory_info().rss)
        )
        return None

    def to_dict(self) -> dict[str, int | float]:
        return {
            "sampling_interval_seconds": self.interval_seconds,
            "baseline_rss_bytes": self.baseline_bytes,
            "peak_sampled_rss_bytes": self.peak_bytes,
            "peak_increment_bytes": max(0, self.peak_bytes - self.baseline_bytes),
        }
