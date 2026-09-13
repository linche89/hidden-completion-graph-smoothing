"""Minimal MATPOWER/PGLib parser and sparse DC-WLS constructor.

The parser intentionally supports the numeric ``mpc.baseMVA``, ``mpc.bus`` and
``mpc.branch`` blocks used by the archived PGLib cases.  It is not a MATLAB
interpreter and fails loudly if those blocks are absent or ragged.
"""

from __future__ import annotations

import hashlib
import math
import re
import time
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .wls_core import WLSInstance, _append_row, information_system


def _without_comments(text: str) -> str:
    return "\n".join(line.split("%", 1)[0] for line in text.splitlines())


def _matrix(text: str, name: str) -> np.ndarray:
    match = re.search(rf"mpc\.{re.escape(name)}\s*=\s*\[(.*?)\]\s*;", text, re.S)
    if not match:
        raise ValueError(f"MATPOWER block mpc.{name} not found")
    parsed: list[list[float]] = []
    for raw_row in match.group(1).split(";"):
        tokens = raw_row.split()
        if tokens:
            parsed.append([float(token) for token in tokens])
    if not parsed or len({len(row) for row in parsed}) != 1:
        raise ValueError(f"MATPOWER block mpc.{name} is empty or ragged")
    return np.asarray(parsed, dtype=np.float64)


def parse_matpower(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    text = _without_comments(raw.decode("utf-8", errors="strict"))
    scalar = re.search(r"mpc\.baseMVA\s*=\s*([^;]+);", text)
    if not scalar:
        raise ValueError("mpc.baseMVA not found")
    return {
        "base_mva": float(scalar.group(1).strip()),
        "bus": _matrix(text, "bus"),
        "branch": _matrix(text, "branch"),
        "source_path": str(path.resolve()),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
    }


def _choose_references(bus: np.ndarray, graph: nx.Graph) -> tuple[int, ...]:
    slack = {idx for idx, row in enumerate(bus) if int(row[1]) == 3}
    refs: list[int] = []
    for component in nx.connected_components(graph):
        choices = sorted(slack & set(component))
        refs.append(choices[0] if choices else min(component))
    return tuple(sorted(refs))


def _append_bus_coeff(entries: list[tuple[int, float]], node_to_state: dict[int, int],
                      node: int, value: float) -> None:
    if node in node_to_state and value != 0.0:
        entries.append((node_to_state[node], value))


def build_dc_wls(
    parsed: dict[str, Any],
    seed: int = 0,
    flow_fraction: float = 1.0,
    injection_fraction: float = 0.25,
    angle_fraction: float = 0.05,
    flow_sigma_mw: float = 1.0,
    injection_sigma_mw: float = 1.5,
    angle_sigma_rad: float = 0.002,
    flow_holder: str = "from",
    truth_mode: str = "synthetic_dc",
) -> WLSInstance:
    bus = parsed["bus"]
    branch = parsed["branch"]
    base_mva = float(parsed["base_mva"])
    if bus.shape[1] < 13 or branch.shape[1] < 13:
        raise ValueError("expected MATPOWER v2 bus (13+) and branch (13+) columns")
    if flow_holder not in {"from", "both"}:
        raise ValueError("flow_holder must be 'from' or 'both'")
    for name, fraction in (
        ("flow_fraction", flow_fraction),
        ("injection_fraction", injection_fraction),
        ("angle_fraction", angle_fraction),
    ):
        if not 0.0 <= fraction <= 1.0:
            raise ValueError(f"{name} must lie in [0, 1]")
    for name, sigma in (
        ("flow_sigma_mw", flow_sigma_mw),
        ("injection_sigma_mw", injection_sigma_mw),
        ("angle_sigma_rad", angle_sigma_rad),
    ):
        if sigma <= 0.0:
            raise ValueError(f"{name} must be positive")

    bus_ids = [int(x) for x in bus[:, 0]]
    if len(set(bus_ids)) != len(bus_ids):
        raise ValueError("duplicate MATPOWER bus identifiers")
    bus_id_to_node = {bus_id: node for node, bus_id in enumerate(bus_ids)}
    graph = nx.Graph()
    graph.add_nodes_from(range(len(bus_ids)))
    active: list[dict[str, float | int]] = []
    skipped_zero_x = 0
    for row in branch:
        if int(row[10]) <= 0:
            continue
        f_id, t_id = int(row[0]), int(row[1])
        if f_id not in bus_id_to_node or t_id not in bus_id_to_node:
            raise ValueError(f"branch refers to unknown bus {f_id}->{t_id}")
        reactance = float(row[3])
        if abs(reactance) <= 1e-12:
            skipped_zero_x += 1
            continue
        f, t = bus_id_to_node[f_id], bus_id_to_node[t_id]
        tap = float(row[8]) if float(row[8]) != 0.0 else 1.0
        shift = math.radians(float(row[9]))
        coefficient = base_mva / (reactance * tap)
        active.append({"f": f, "t": t, "coefficient": coefficient, "shift": shift})
        graph.add_edge(f, t)

    references = _choose_references(bus, graph)
    state_nodes = np.asarray([node for node in graph.nodes if node not in references], dtype=int)
    node_to_state = {node: idx for idx, node in enumerate(state_nodes)}
    rng = np.random.default_rng(seed)

    # Branch Laplacian for a smooth, reproducible non-flat angle field.
    lap_rows: list[int] = []
    lap_cols: list[int] = []
    lap_values: list[float] = []
    for edge in active:
        f, t = int(edge["f"]), int(edge["t"])
        coeff = abs(float(edge["coefficient"]))
        # Accumulate Laplacian entries without forming a dense incidence matrix.
        for ni, si in ((f, 1.0), (t, -1.0)):
            if ni not in node_to_state:
                continue
            ii = node_to_state[ni]
            for nj, sj in ((f, 1.0), (t, -1.0)):
                if nj in node_to_state:
                    lap_rows.append(ii)
                    lap_cols.append(node_to_state[nj])
                    lap_values.append(coeff * si * sj)
    state_dim = len(state_nodes)
    if state_dim == 0:
        raise ValueError("DC model has no free angle state after reference elimination")
    laplacian = sp.coo_matrix(
        (lap_values, (lap_rows, lap_cols)), shape=(state_dim, state_dim)
    ).tocsr()
    if truth_mode == "case_va" and np.max(np.abs(bus[:, 8])) > 1e-10:
        full_angles = np.radians(bus[:, 8].astype(float))
        for component in nx.connected_components(graph):
            ref = next(x for x in references if x in component)
            full_angles[list(component)] -= full_angles[ref]
        truth = full_angles[state_nodes]
    elif truth_mode == "synthetic_dc":
        forcing = rng.standard_normal(state_dim)
        regularized = laplacian + 1e-8 * sp.eye(state_dim, format="csr")
        truth = np.asarray(spla.spsolve(regularized.tocsc(), forcing)).ravel()
        scale = np.max(np.abs(truth))
        if scale:
            truth *= 0.10 / scale
    else:
        raise ValueError("truth_mode must be 'synthetic_dc' or a non-flat 'case_va'")

    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    means: list[float] = []
    offsets: list[float] = []
    sigmas: list[float] = []
    kinds: list[str] = []
    variable_support: list[tuple[int, ...]] = []
    physical_support: list[tuple[int, ...]] = []
    holders: list[tuple[int, ...]] = []
    row_index = 0

    selected_flow = rng.random(len(active)) < flow_fraction
    for selected, edge in zip(selected_flow, active):
        if not selected:
            continue
        f, t = int(edge["f"]), int(edge["t"])
        coeff, shift = float(edge["coefficient"]), float(edge["shift"])
        entries: list[tuple[int, float]] = []
        _append_bus_coeff(entries, node_to_state, f, coeff)
        _append_bus_coeff(entries, node_to_state, t, -coeff)
        _append_row(rows, cols, values, row_index, entries)
        means.append(sum(value * truth[col] for col, value in entries))
        offsets.append(-coeff * shift)
        sigmas.append(flow_sigma_mw)
        kinds.append("branch_flow_from")
        variable_support.append(tuple(sorted(node for node in (f, t) if node in node_to_state)))
        physical_support.append(tuple(sorted((f, t))))
        holders.append((f,) if flow_holder == "from" else tuple(sorted((f, t))))
        row_index += 1

    incident: dict[int, list[tuple[int, int, float, float]]] = {node: [] for node in graph.nodes}
    for edge in active:
        f, t = int(edge["f"]), int(edge["t"])
        incident[f].append((f, t, float(edge["coefficient"]), -float(edge["coefficient"]) * float(edge["shift"])))
        incident[t].append((t, f, float(edge["coefficient"]), float(edge["coefficient"]) * float(edge["shift"])))

    selected_injection = rng.random(graph.number_of_nodes()) < injection_fraction
    for node in graph.nodes:
        if not selected_injection[node]:
            continue
        coefficient_by_node: dict[int, float] = {}
        known_offset = 0.0
        physical = {node}
        for center, neighbor, coeff, edge_offset in incident[node]:
            coefficient_by_node[center] = coefficient_by_node.get(center, 0.0) + coeff
            coefficient_by_node[neighbor] = coefficient_by_node.get(neighbor, 0.0) - coeff
            known_offset += edge_offset
            physical.add(neighbor)
        entries = []
        for variable_node, value in coefficient_by_node.items():
            _append_bus_coeff(entries, node_to_state, variable_node, value)
        _append_row(rows, cols, values, row_index, entries)
        means.append(sum(value * truth[col] for col, value in entries))
        offsets.append(known_offset)
        sigmas.append(injection_sigma_mw)
        kinds.append("bus_injection")
        variable_support.append(tuple(sorted(x for x in physical if x in node_to_state)))
        physical_support.append(tuple(sorted(physical)))
        holders.append((node,))
        row_index += 1

    selected_angle = rng.random(graph.number_of_nodes()) < angle_fraction
    for node in state_nodes:
        if not selected_angle[node]:
            continue
        state_index = node_to_state[int(node)]
        _append_row(rows, cols, values, row_index, [(state_index, 1.0)])
        means.append(float(truth[state_index]))
        offsets.append(0.0)
        sigmas.append(angle_sigma_rad)
        kinds.append("bus_angle")
        variable_support.append((int(node),))
        physical_support.append((int(node),))
        holders.append((int(node),))
        row_index += 1

    h = sp.coo_matrix((values, (rows, cols)), shape=(row_index, state_dim)).tocsr()
    if row_index == 0:
        raise ValueError("measurement sampling produced no measurements")
    sigma = np.asarray(sigmas)
    precision = 1.0 / np.square(sigma)
    known_offset_array = np.asarray(offsets)
    corrected_mean = np.asarray(means)
    measurement = corrected_mean + known_offset_array + sigma * rng.standard_normal(row_index)
    corrected = measurement - known_offset_array
    j, rhs = information_system(h, precision, corrected)
    return WLSInstance(
        graph=graph,
        h=h,
        precision=precision,
        measurement=measurement,
        corrected_measurement=corrected,
        known_offset=known_offset_array,
        j=j,
        rhs=rhs,
        truth=truth,
        state_nodes=state_nodes,
        node_to_state=node_to_state,
        reference_nodes=references,
        measurement_kind=tuple(kinds),
        measurement_variable_support=tuple(variable_support),
        measurement_physical_support=tuple(physical_support),
        measurement_holders=tuple(holders),
        metadata={
            "kind": "pglib_dc_wls",
            "source_path": parsed["source_path"],
            "source_sha256": parsed["source_sha256"],
            "base_mva": base_mva,
            "bus_ids": bus_ids,
            "reference_bus_ids": [bus_ids[node] for node in references],
            "active_modeled_branches": len(active),
            "skipped_zero_reactance_branches": skipped_zero_x,
            "flow_holder_semantics": flow_holder,
            "truth_mode": truth_mode,
            "counts_by_kind": {
                kind: kinds.count(kind) for kind in sorted(set(kinds))
            },
        },
    )


def load_dc_wls(path: Path, **kwargs: Any) -> tuple[WLSInstance, dict[str, float]]:
    start = time.perf_counter()
    parsed = parse_matpower(path)
    parse_seconds = time.perf_counter() - start
    start = time.perf_counter()
    instance = build_dc_wls(parsed, **kwargs)
    build_seconds = time.perf_counter() - start
    return instance, {"parse_seconds": parse_seconds, "build_seconds": build_seconds}
