# Result files

- `pglib_cpu_final.json`: final six-case CPU correctness/scale run. OpenBLAS and
  OpenMP were fixed to one thread. Tails use 32 uniform root draws with replacement.
- `pglib_gpu_case10000_final.json`: final same-matrix CPU/CuPy comparison using
  Jacobi-preconditioned CG at relative tolerance `1e-8` and 1,000 repeated SpMVs.
- `pglib_cpu_recheck.json` and `pglib_gpu_case10000_recheck.json`: independent
  reruns by the coordinating process.  Structural counts, checksums, convergence
  status, iteration counts, and residuals reproduce the final runs; wall times are
  retained separately rather than expected to match exactly.
- `grid_gpu_scale_v4.json`: final six-case grid throughput run produced by the
  implementation task. `grid_gpu_scale_v5_recheck.json` is the coordinating
  process's independent rerun after the exact degree/radius wording correction;
  it reproduces all structural and numerical invariants while retaining timing
  variation instead of selecting only the fastest run.
- `pglib_gpu_case10000.json`: retained negative result from unpreconditioned CG;
  both devices reached 20,000 iterations without meeting tolerance.
- `grid_gpu_scale_v4.json`: final direct-sparse 2-D grid CPU/GPU scale run through
  about one million states, including separate SpMV, batched SpMM, and 40-round
  local Richardson timings. `grid_gpu_scale_v3.json` is the preceding independent
  run with a colder CUDA/JIT path; earlier `grid_gpu_scale*.json` files are retained
  development runs.
- `pglib_cpu_small.json`, `pglib_cpu_large.json`, and
  `pglib_gpu_case10000_jacobi.json`: development runs retained for auditability;
  use the two `final` files for reported numbers.

CPU RSS is the maximum of explicit process checkpoints, not an operating-system
peak-memory measurement. DKW intervals are conditional on one fixed graph and valid
per fixed radius, not simultaneously over all radii.

## Graph-smoother application records

- `p3_p4_smoke.json`: canonical theorem-matched robust-design experiment plus
  four fixed-graph synthetic stress cases.  It records the exact finite-`h`
  and unbounded scenarios, sampled and relaxed baselines, pruning ablation,
  active worst cases, resource models and the sharp budget sweep.
- `p4_external_scale.json`: six PGLib topology transformations and three
  SuiteSparse sparsity-graph transformations, reaching 129,164 nodes.  These
  are explicitly external graph-smoother stress tests, not raw WLS instances.
- `p4_parameter_sweep.json`: 135 theorem-matched physical completions spanning
  five hidden budgets, three edge-weight ranges, three requested graph
  densities, and three repetitions per cell.  Every row records its realized
  port-cut size/weight and is checked against an independently solved exact
  finite-`h` certificate.  Its provenance points to clean commit `6d1e183`.
- `smoother_gpu_scale.json`, `smoother_gpu_scale_recheck1.json` and
  `smoother_gpu_scale_recheck2.json`: three independent CPU/RTX 5080 timing
  repetitions for a fixed 40-round local Richardson smoother on grids from
  10,000 to 1,000,000 states.  Report medians and retain the full timing spread.
