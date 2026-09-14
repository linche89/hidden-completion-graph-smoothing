# Locality benchmark skeleton

For the exact checked environment and the one-command paper verification
workflow, see [`REPRODUCIBILITY.md`](../REPRODUCIBILITY.md).

This directory is isolated from the literature survey. It contains a small CPU-only
benchmark that checks mathematical locality quantities for synthetic block WLS
instances. It is not yet the production HPC implementation.

The benchmark deliberately does **not** form `J^{-1}`. For sampled target nodes it
solves `J^T y=e_i` and measures the target block row outside an `r`-hop ball. It also
compares:

- the centralized sparse WLS solve;
- an `r`-ball Dirichlet/truncated solve;
- a Schur-complement oracle (small instances only);
- a degree-`r` Chebyshev inverse polynomial.

It records residuals, state error, exact sampled operator-row tails, sampled raw-data
operator-tail exceedance fractions with a Dvoretzky–Kiefer–Wolfowitz confidence band,
rounds, message width,
bit-hop proxies, wall time, sparse nonzeros, and process RSS.

The smoke configuration gives full-rank node measurements to only 20% of nodes; the
remaining nodes have no local anchor. Edge measurements are initially held by both
endpoints. Factor variable support and measurement holders are represented separately.
The current generator intentionally builds dense temporary factor rows before converting
to CSR, so it is a mathematical interface test only, not the large-scale generator.

Run the seconds-scale smoke case from the repository root:

```powershell
python experiments/locality_benchmark.py --config experiments/configs/smoke.json
```

A stricter smoke case gives a full-rank node measurement to only one of 72 nodes:

```powershell
python experiments/locality_benchmark.py --config experiments/configs/single_anchor_smoke.json
```

The theorem-level numerical regression tests check the Schur factorization and
norm sandwich, the sharp spectral-window constant, the scalar robust Schur-window
minimax formula, the grounded-Laplacian oracle/Dirichlet/Neumann ordering, and the
single-ground path counterexample:

```powershell
python -m unittest experiments.test_theory_identities -v
```

The grounded-network completion checks verify the matrix-forest/DPP projection
identities, high-conductance limits, raw pairwise-WLS realization, optimized
intermediate-cut and bounded-weight counterexamples, the strict Petersen gap, and
the raw-measurement metric orthogonalization audit plus uniform/heterogeneous
hidden-node partial-partition bounds with degree-two converses:

```powershell
python experiments/theory_search/forest_partition_extrema.py
python experiments/theory_search/hidden_partial_partition_checks.py
```

Results are JSON and may be redirected by the caller. Large sweeps should use a
separate output directory and the cluster workflow described in
`research/agent_reports/hpc_experiments_phase2.md`.

The exact hidden-completion solver enumerates set partitions and hidden-budget
allocations, applies the proved vertex pruning rule, solves the lossless
spectral SDP, and reports active worst-case topology scenarios:

```powershell
python -m pip install -r hidden_completion/requirements.txt
python -m hidden_completion solve experiments/configs/exact_finite_q2_h1.json
python -m hidden_completion solve experiments/configs/exact_unbounded_shared_opcode.json
python -m unittest experiments.test_exact_scenario_solver -v
```

See `hidden_completion/README.md` for the JSON schema, including zero patterns,
`r`-hop masks, and shared-coefficient constraints.

## P3/P4 graph-smoother experiments

The primary-application runner keeps theorem-matched robust design separate
from fixed-known-graph stress tests:

```powershell
python -m experiments.smoother `
  --config experiments/configs/p3_p4_smoke.json `
  --output experiments/results/p3_p4_smoke.json `
  --section all

python -m experiments.smoother `
  --config experiments/configs/p4_external_scale.json `
  --output experiments/results/p4_external_scale.json `
  --section stress
```

The first command runs exact finite-`h` and unbounded certificates, random
completion sampling, the PSD-contraction outer relaxation, the endpoint NO-GO
ablation, vertex-pruning and scenario-count sweeps, and four 36-node synthetic
graph families.  The second command treats PGLib topologies and SuiteSparse
sparsity patterns only as external graph-smoother stress transformations; it
does not relabel them as raw WLS or Kalman models.

Export the checked JSON records to backend-neutral CSV files for TikZ with:

```powershell
python -m experiments.smoother.paper_data `
  --p3-p4 experiments/results/p3_p4_smoke.json `
  --external experiments/results/p4_external_scale.json `
  --gpu experiments/results/smoother_gpu_scale.json `
        experiments/results/smoother_gpu_scale_recheck1.json `
        experiments/results/smoother_gpu_scale_recheck2.json
```

The lightweight factorial audit varies hidden budget, edge-weight scale and
graph density (recording the realized visible/hidden cut for every sample):

```powershell
python -m experiments.smoother.parameter_sweep `
  --config experiments/configs/p4_parameter_sweep.json `
  --output experiments/results/p4_parameter_sweep.json
```

Physical samples are lower-bound stress points only.  Each row is checked
against the exact finite-`h` certificate produced independently of sampling.
