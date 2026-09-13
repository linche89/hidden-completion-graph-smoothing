"""Sparse-only WLS data structures, synthetic generator, and locality metrics.

Large matrices are assembled directly as COO/CSR.  The module never constructs a
dense H, dense J, or a complete inverse.  Small dense arrays are limited to local
state blocks and sampled inverse rows obtained from sparse solves.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import networkx as nx
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass
class WLSInstance:
    graph: nx.Graph
    h: sp.csr_matrix
    precision: np.ndarray
    measurement: np.ndarray
    corrected_measurement: np.ndarray
    known_offset: np.ndarray
    j: sp.csr_matrix
    rhs: np.ndarray
    truth: np.ndarray
    state_nodes: np.ndarray
    node_to_state: dict[int, int]
    reference_nodes: tuple[int, ...]
    measurement_kind: tuple[str, ...]
    measurement_variable_support: tuple[tuple[int, ...], ...]
    measurement_physical_support: tuple[tuple[int, ...], ...]
    measurement_holders: tuple[tuple[int, ...], ...]
    metadata: dict[str, Any]


def information_system(h: sp.csr_matrix, precision: np.ndarray,
                       corrected_measurement: np.ndarray,
                       prior_precision: float = 0.0) -> tuple[sp.csr_matrix, np.ndarray]:
    weighted_h = h.multiply(precision[:, None])
    j = (h.T @ weighted_h).tocsr()
    if prior_precision:
        j = j + float(prior_precision) * sp.eye(j.shape[0], format="csr")
    rhs = np.asarray(h.T @ (precision * corrected_measurement)).ravel()
    return j, rhs


def _append_row(rows: list[int], cols: list[int], data: list[float], row: int,
                entries: Iterable[tuple[int, float]]) -> None:
    for col, value in entries:
        if value != 0.0:
            rows.append(row)
            cols.append(int(col))
            data.append(float(value))


def sparse_synthetic_block_wls(
    graph: nx.Graph,
    block_dim: int = 2,
    anchor_fraction: float = 0.1,
    anchor_rank: int | None = None,
    anchor_weight: float = 0.2,
    edge_weight: float = 1.0,
    noise_std: float = 0.02,
    seed: int = 0,
) -> WLSInstance:
    """Build a block WLS model directly in COO.

    Edge factors are ``x_u - Q_e x_v`` with a small orthogonal block ``Q_e``.
    Only a configurable fraction of nodes receives a local anchor.  Consequently
    this generator can test local rank deficiency rather than silently imposing a
    full-rank measurement at every node.
    """
    graph = nx.convert_node_labels_to_integers(graph)
    n = graph.number_of_nodes()
    d = int(block_dim)
    rank = d if anchor_rank is None else int(anchor_rank)
    if not (0 <= rank <= d):
        raise ValueError("anchor_rank must lie between zero and block_dim")
    rng = np.random.default_rng(seed)
    anchor_count = min(n, max(1 if rank else 0, round(anchor_fraction * n)))
    anchors = set(int(x) for x in rng.choice(n, size=anchor_count, replace=False))
    if rank == d and n and 0 not in anchors:
        anchors.discard(next(iter(anchors)))
        anchors.add(0)

    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    kinds: list[str] = []
    variable_support: list[tuple[int, ...]] = []
    physical_support: list[tuple[int, ...]] = []
    holders: list[tuple[int, ...]] = []
    row = 0

    for node in sorted(anchors):
        local_q, _ = np.linalg.qr(rng.standard_normal((d, d)))
        block = np.sqrt(anchor_weight) * local_q[:rank, :]
        for local_row in range(rank):
            _append_row(rows, cols, values, row,
                        ((node * d + k, block[local_row, k]) for k in range(d)))
            kinds.append("anchor")
            variable_support.append((node,))
            physical_support.append((node,))
            holders.append((node,))
            row += 1

    for u, v in sorted(graph.edges()):
        q, r = np.linalg.qr(rng.standard_normal((d, d)))
        signs = np.sign(np.diag(r))
        signs[signs == 0] = 1.0
        q = q * signs
        scale = np.sqrt(edge_weight)
        for local_row in range(d):
            entries = [(u * d + local_row, scale)]
            entries.extend((v * d + k, -scale * q[local_row, k]) for k in range(d))
            _append_row(rows, cols, values, row, entries)
            kinds.append("edge")
            variable_support.append((u, v))
            physical_support.append((u, v))
            holders.append((u, v))
            row += 1

    h = sp.coo_matrix((values, (rows, cols)), shape=(row, n * d)).tocsr()
    truth = rng.standard_normal(n * d)
    offset = np.zeros(row)
    measurement = np.asarray(h @ truth).ravel() + noise_std * rng.standard_normal(row)
    precision = np.full(row, 1.0 / (noise_std * noise_std))
    j, rhs = information_system(h, precision, measurement)
    return WLSInstance(
        graph=graph,
        h=h,
        precision=precision,
        measurement=measurement,
        corrected_measurement=measurement,
        known_offset=offset,
        j=j,
        rhs=rhs,
        truth=truth,
        state_nodes=np.repeat(np.arange(n), d),
        node_to_state={node: node * d for node in range(n)},
        reference_nodes=(),
        measurement_kind=tuple(kinds),
        measurement_variable_support=tuple(variable_support),
        measurement_physical_support=tuple(physical_support),
        measurement_holders=tuple(holders),
        metadata={
            "kind": "synthetic_block",
            "block_dim": d,
            "anchors": sorted(anchors),
            "anchor_rank": rank,
            "local_full_rank_fraction": len(anchors) / n if rank == d and n else 0.0,
        },
    )


def graph_ball(graph: nx.Graph, root: int, radius: int) -> set[int]:
    return set(nx.single_source_shortest_path_length(graph, root, cutoff=radius))


def sampled_scalar_tails(instance: WLSInstance, factor: spla.SuperLU,
                         root_nodes: Iterable[int], radius: int) -> dict[int, dict[str, float]]:
    """Exact scalar target tails from sparse transpose solves.

    This routine is for one scalar state per non-reference node (the DC model).
    It materializes one sampled row at a time, never the complete inverse.
    """
    result: dict[int, dict[str, float]] = {}
    state_nodes = np.asarray(instance.state_nodes)
    for root in sorted(set(int(x) for x in root_nodes)):
        if root not in instance.node_to_state:
            continue
        state_index = instance.node_to_state[root]
        selector = np.zeros(instance.j.shape[0])
        selector[state_index] = 1.0
        inverse_row = np.asarray(factor.solve(selector, trans="T")).ravel()
        ball = graph_ball(instance.graph, root, radius)

        state_mask = np.asarray([node not in ball for node in state_nodes], dtype=bool)
        state_tail = float(np.linalg.norm(inverse_row[state_mask]))

        h_inverse = np.asarray(instance.h @ inverse_row).ravel()
        # Physical raw-data map: T_i^T = R^{-1} H J^{-T} E_i^T.
        raw_row = instance.precision * h_inverse
        # Dimensionless whitened-data map for w=R^{-1/2}(z-offset):
        # T_i^T = R^{-1/2} H J^{-T} E_i^T.
        whitened_row = np.sqrt(instance.precision) * h_inverse
        measurement_mask = np.asarray(
            [not bool(set(owner) & ball) for owner in instance.measurement_holders],
            dtype=bool,
        )
        raw_tail = float(np.linalg.norm(raw_row[measurement_mask]))
        whitened_tail = float(np.linalg.norm(whitened_row[measurement_mask]))
        result[root] = {
            "state_rhs_tail_l2": state_tail,
            "physical_raw_z_tail_l2": raw_tail,
            "whitened_z_tail_l2": whitened_tail,
        }
    return result
