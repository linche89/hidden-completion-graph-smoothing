"""Weighted graph generators and external-topology loaders."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
import math
from pathlib import Path
import tarfile
from typing import Any

import networkx as nx
import numpy as np
from numpy.typing import NDArray
import scipy.io
import scipy.sparse as sp
import scipy.sparse.csgraph as csgraph

from experiments.scalable.pglib_dc import parse_matpower


Array = NDArray[np.float64]


@dataclass(slots=True)
class WeightedGraph:
    """A nonnegative symmetric adjacency and provenance metadata."""

    adjacency: sp.csr_matrix
    metadata: dict[str, Any]

    @property
    def nodes(self) -> int:
        return int(self.adjacency.shape[0])

    @property
    def edges(self) -> int:
        return int(sp.triu(self.adjacency, k=1).nnz)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validated_adjacency(adjacency: sp.spmatrix | Array) -> sp.csr_matrix:
    matrix = sp.csr_matrix(adjacency, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] < 1:
        raise ValueError("adjacency must be a nonempty square matrix")
    matrix.sum_duplicates()
    matrix.setdiag(0.0)
    matrix.eliminate_zeros()
    if matrix.data.size and (not np.all(np.isfinite(matrix.data)) or np.min(matrix.data) < 0):
        raise ValueError("adjacency weights must be finite and nonnegative")
    difference = matrix - matrix.T
    if difference.nnz and np.max(np.abs(difference.data)) > 1e-11:
        raise ValueError("adjacency must be symmetric")
    matrix = ((matrix + matrix.T) * 0.5).tocsr()
    matrix.sort_indices()
    return matrix


def _rectangular_grid(n: int) -> nx.Graph:
    rows = max(1, int(math.floor(math.sqrt(n))))
    columns = int(math.ceil(n / rows))
    full = nx.grid_2d_graph(rows, columns)
    ordered = sorted(full.nodes())[:n]
    return nx.convert_node_labels_to_integers(full.subgraph(ordered).copy())


def _connect_components_geometrically(graph: nx.Graph, positions: Array) -> None:
    while not nx.is_connected(graph):
        components = [sorted(component) for component in nx.connected_components(graph)]
        left = components[0]
        best: tuple[float, int, int] | None = None
        for right in components[1:]:
            for u in left:
                differences = positions[np.asarray(right)] - positions[u]
                index = int(np.argmin(np.sum(differences * differences, axis=1)))
                v = right[index]
                candidate = (float(np.sum((positions[u] - positions[v]) ** 2)), u, v)
                if best is None or candidate < best:
                    best = candidate
        assert best is not None
        graph.add_edge(best[1], best[2])


def _base_graph(kind: str, n: int, seed: int, options: dict[str, Any]) -> tuple[nx.Graph, dict[str, Any]]:
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError("n must be an integer at least two")
    rng = np.random.default_rng(seed)
    metadata: dict[str, Any] = {"family": kind, "seed": seed}
    if kind == "path":
        graph = nx.path_graph(n)
    elif kind == "cycle":
        graph = nx.cycle_graph(n)
    elif kind == "grid":
        graph = _rectangular_grid(n)
    elif kind == "random_cyclic":
        mean_degree = float(options.get("mean_degree", 4.0))
        target_edges = min(n * (n - 1) // 2, max(n, int(round(mean_degree * n / 2))))
        graph = nx.cycle_graph(n)
        nonedges = list(nx.non_edges(graph))
        rng.shuffle(nonedges)
        graph.add_edges_from(nonedges[: max(0, target_edges - graph.number_of_edges())])
        metadata["target_mean_degree"] = mean_degree
    elif kind == "geometric":
        radius = float(options.get("geometric_radius", 1.7 / math.sqrt(n)))
        positions = rng.random((n, 2))
        position_map = {node: positions[node] for node in range(n)}
        graph = nx.random_geometric_graph(n, radius, pos=position_map)
        _connect_components_geometrically(graph, positions)
        if graph.number_of_edges() == n - 1 and n >= 3:
            candidates = sorted(
                (
                    float(np.sum((positions[u] - positions[v]) ** 2)),
                    u,
                    v,
                )
                for u, v in nx.non_edges(graph)
            )
            if candidates:
                _, u, v = candidates[0]
                graph.add_edge(u, v)
        metadata["geometric_radius"] = radius
        metadata["positions"] = positions.tolist()
    elif kind == "small_world":
        requested = int(options.get("small_world_k", 4))
        k = min(requested, n - 1)
        if k % 2:
            k -= 1
        k = max(2, k) if n >= 3 else 1
        probability = float(options.get("rewire_probability", 0.2))
        if not 0.0 <= probability <= 1.0:
            raise ValueError("rewire_probability must lie in [0,1]")
        graph = nx.connected_watts_strogatz_graph(n, k, probability, tries=200, seed=seed)
        metadata.update({"small_world_k": k, "rewire_probability": probability})
    else:
        raise ValueError(f"unknown synthetic graph family {kind!r}")
    if not nx.is_connected(graph):
        raise RuntimeError(f"generator returned a disconnected {kind} graph")
    return nx.convert_node_labels_to_integers(graph), metadata


def synthetic_weighted_graph(
    kind: str,
    n: int,
    seed: int,
    *,
    weight_min: float = 0.25,
    weight_max: float = 4.0,
    options: dict[str, Any] | None = None,
) -> WeightedGraph:
    """Build a connected synthetic graph with deterministic log-uniform weights."""

    if not (0.0 < weight_min <= weight_max < math.inf):
        raise ValueError("weights must satisfy 0 < weight_min <= weight_max < infinity")
    graph, metadata = _base_graph(kind, n, seed, dict(options or {}))
    rng = np.random.default_rng(seed + 104729)
    edges = sorted((min(u, v), max(u, v)) for u, v in graph.edges())
    logs = rng.uniform(math.log(weight_min), math.log(weight_max), size=len(edges))
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    weights: list[float] = []
    for (u, v), log_weight in zip(edges, logs, strict=True):
        weight = float(math.exp(log_weight))
        rows.extend((u, v))
        columns.extend((v, u))
        values.extend((weight, weight))
        weights.append(weight)
    adjacency = _validated_adjacency(
        sp.coo_matrix((values, (rows, columns)), shape=(n, n))
    )
    metadata.update(
        {
            "classification": "synthetic",
            "nodes": n,
            "edges": len(edges),
            "weight_min_requested": weight_min,
            "weight_max_requested": weight_max,
            "weight_min_realized": min(weights),
            "weight_max_realized": max(weights),
            "weight_median_realized": float(np.median(weights)),
        }
    )
    return WeightedGraph(adjacency, metadata)


def load_pglib_topology(path: str | Path, *, normalize_median: bool = True) -> WeightedGraph:
    """Load active PGLib branches as a positive weighted topology.

    This is a graph-signal transformation, not an OPF or WLS benchmark.
    Parallel branch weights are added.
    """

    source = Path(path)
    parsed = parse_matpower(source)
    bus = np.asarray(parsed["bus"])
    branch = np.asarray(parsed["branch"])
    bus_ids = [int(value) for value in bus[:, 0]]
    lookup = {bus_id: index for index, bus_id in enumerate(bus_ids)}
    if len(lookup) != len(bus_ids):
        raise ValueError("duplicate PGLib bus identifiers")
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    skipped = 0
    for record in branch:
        if int(record[10]) <= 0 or abs(float(record[3])) <= 1e-12:
            skipped += 1
            continue
        u, v = lookup[int(record[0])], lookup[int(record[1])]
        tap = float(record[8]) if float(record[8]) != 0.0 else 1.0
        weight = abs(float(parsed["base_mva"]) / (float(record[3]) * tap))
        rows.extend((u, v))
        columns.extend((v, u))
        values.extend((weight, weight))
    adjacency = _validated_adjacency(
        sp.coo_matrix((values, (rows, columns)), shape=(len(bus_ids), len(bus_ids)))
    )
    positive = adjacency.data[adjacency.data > 0]
    scale = float(np.median(positive)) if normalize_median and positive.size else 1.0
    if scale <= 0:
        scale = 1.0
    adjacency = (adjacency / scale).tocsr()
    components, _ = csgraph.connected_components(adjacency, directed=False)
    return WeightedGraph(
        adjacency,
        {
            "classification": "external_application_stress",
            "source_kind": "PGLib topology",
            "source_path": str(source.resolve()),
            "source_sha256": parsed["source_sha256"],
            "transformation": (
                "active branches mapped to abs(baseMVA/(x*tap)), parallel weights summed, "
                "then divided by their median; graph smoother only"
            ),
            "nodes": adjacency.shape[0],
            "edges": int(sp.triu(adjacency, k=1).nnz),
            "connected_components": int(components),
            "skipped_inactive_or_zero_reactance_branches": skipped,
            "normalization_scale": scale,
        },
    )


def load_suitesparse_pattern(
    archive: str | Path,
    *,
    member: str | None = None,
    use_value_magnitudes: bool = False,
) -> WeightedGraph:
    """Load a square Matrix Market archive as a nonnegative sparsity graph."""

    source = Path(archive)
    with tarfile.open(source, "r:gz") as bundle:
        candidates = sorted(name for name in bundle.getnames() if name.lower().endswith(".mtx"))
        if member is None:
            non_coordinate = [name for name in candidates if "coord" not in name.lower()]
            if len(non_coordinate) != 1:
                raise ValueError(f"archive needs an explicit member; candidates={candidates}")
            member = non_coordinate[0]
        if member not in candidates:
            raise ValueError(f"Matrix Market member {member!r} not found")
        extracted = bundle.extractfile(member)
        if extracted is None:
            raise ValueError(f"cannot extract {member!r}")
        raw = extracted.read()
    matrix = scipy.io.mmread(io.BytesIO(raw))
    matrix = sp.coo_matrix(matrix, dtype=np.float64)
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("SuiteSparse stress matrix must be square")
    matrix.setdiag(0.0)
    matrix.eliminate_zeros()
    if use_value_magnitudes:
        matrix.data = np.abs(matrix.data)
    else:
        matrix.data = np.ones_like(matrix.data)
    adjacency = matrix.maximum(matrix.T).tocsr()
    adjacency = _validated_adjacency(adjacency)
    positive = adjacency.data[adjacency.data > 0]
    scale = float(np.median(positive)) if use_value_magnitudes and positive.size else 1.0
    if scale <= 0:
        scale = 1.0
    adjacency = (adjacency / scale).tocsr()
    components, _ = csgraph.connected_components(adjacency, directed=False)
    return WeightedGraph(
        adjacency,
        {
            "classification": "external_application_stress",
            "source_kind": "SuiteSparse sparsity graph",
            "source_path": str(source.resolve()),
            "source_sha256": _sha256(source),
            "archive_member": member,
            "transformation": (
                "symmetric off-diagonal sparsity pattern with unit positive weights; "
                "original matrix is not treated as a graph smoother"
                if not use_value_magnitudes
                else "absolute symmetric off-diagonal values normalized by their median"
            ),
            "nodes": adjacency.shape[0],
            "edges": int(sp.triu(adjacency, k=1).nnz),
            "connected_components": int(components),
            "uses_original_value_magnitudes": bool(use_value_magnitudes),
            "normalization_scale": scale,
        },
    )


def smoother_matrix(graph: WeightedGraph, *, edge_scale: float = 1.0) -> sp.csr_matrix:
    """Return ``I + edge_scale * L`` for the graph smoother objective."""

    if not np.isfinite(edge_scale) or edge_scale < 0:
        raise ValueError("edge_scale must be finite and nonnegative")
    adjacency = graph.adjacency
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    laplacian = sp.diags(degree) - adjacency
    return (sp.eye(graph.nodes, format="csr") + edge_scale * laplacian).tocsr()


def dense_adjacency(graph: WeightedGraph) -> Array:
    return np.asarray(graph.adjacency.toarray(), dtype=np.float64)
