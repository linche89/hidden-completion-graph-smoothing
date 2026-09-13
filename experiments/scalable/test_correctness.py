"""Fast correctness tests; run with ``python -m experiments.scalable.test_correctness``."""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np
import scipy.sparse.linalg as spla

from .pglib_dc import load_dc_wls, parse_matpower
from .wls_core import sampled_scalar_tails, sparse_synthetic_block_wls


ROOT = Path(__file__).resolve().parents[2]
CASE14 = ROOT / "datasets" / "raw" / "pglib_opf_v23.07" / "pglib_opf_case14_ieee.m"


def main() -> None:
    parsed = parse_matpower(CASE14)
    assert parsed["bus"].shape[0] == 14
    assert parsed["branch"].shape[0] == 20

    instance, _ = load_dc_wls(CASE14, seed=11, flow_holder="from")
    assert instance.j.shape == (13, 13)
    assert instance.h.shape[1] == 13
    assert len(instance.measurement_holders) == instance.h.shape[0]
    assert len(instance.measurement_variable_support) == instance.h.shape[0]
    flow_index = instance.measurement_kind.index("branch_flow_from")
    assert len(instance.measurement_holders[flow_index]) == 1
    assert len(instance.measurement_physical_support[flow_index]) == 2

    factor = spla.splu(instance.j.tocsc())
    estimate = factor.solve(instance.rhs)
    residual = np.linalg.norm(instance.j @ estimate - instance.rhs) / np.linalg.norm(instance.rhs)
    assert residual < 1e-10
    roots = [node for node in instance.node_to_state][:3]
    small = sampled_scalar_tails(instance, factor, roots, 0)
    diameter = nx.diameter(instance.graph)
    full = sampled_scalar_tails(instance, factor, roots, diameter)
    assert all(full[root]["whitened_z_tail_l2"] < 1e-12 for root in roots)
    assert all(full[root]["physical_raw_z_tail_l2"] < 1e-12 for root in roots)
    assert all(full[root]["state_rhs_tail_l2"] < 1e-12 for root in roots)
    assert all(
        small[root]["whitened_z_tail_l2"] >= full[root]["whitened_z_tail_l2"]
        for root in roots
    )

    synthetic = sparse_synthetic_block_wls(
        nx.cycle_graph(40), block_dim=2, anchor_fraction=0.05, anchor_rank=2, seed=5
    )
    assert synthetic.h.shape[1] == 80
    assert synthetic.metadata["local_full_rank_fraction"] < 1.0
    synthetic_factor = spla.splu(synthetic.j.tocsc())
    synthetic_estimate = synthetic_factor.solve(synthetic.rhs)
    synthetic_residual = (
        np.linalg.norm(synthetic.j @ synthetic_estimate - synthetic.rhs)
        / np.linalg.norm(synthetic.rhs)
    )
    assert synthetic_residual < 1e-10
    print("all scalable correctness tests passed")


if __name__ == "__main__":
    main()
