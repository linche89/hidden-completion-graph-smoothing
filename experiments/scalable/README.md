# Scalable sparse WLS implementation

Exact package pins for the checked CPU and RTX 5080 runs are listed in
[`REPRODUCIBILITY.md`](../../REPRODUCIBILITY.md).  The shorter requirement file
below remains the portable GPU environment specification.

All large matrices are assembled directly as COO/CSR. No complete inverse or dense
`H/J` is constructed.

CPU correctness and scale benchmark:

```powershell
python -m experiments.scalable.test_correctness
python -m experiments.scalable.cpu_benchmark `
  --data-dir datasets/raw/pglib_opf_v23.07 `
  --output experiments/results/pglib_cpu.json
```

Optional project-isolated GPU environment:

```powershell
python -m venv experiments/.venv-gpu
experiments/.venv-gpu/Scripts/python -m pip install -r experiments/scalable/requirements-gpu.txt
experiments/.venv-gpu/Scripts/python -m experiments.scalable.gpu_sparse_benchmark `
  --case datasets/raw/pglib_opf_v23.07/pglib_opf_case2869_pegase.m `
  --output experiments/results/pglib_gpu_case2869.json
```

The GPU benchmark separately records CPU preprocessing, host/device transfers,
warmup, SpMV/CG algorithm time, and device-to-host transfer. It does not use or claim
MPI.

Large synthetic 2-D grid CPU/GPU scale benchmark (about `1e5` and `1e6` states,
stable-mass and single-ground families, no NetworkX):

```powershell
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
experiments/.venv-gpu/Scripts/python -m experiments.scalable.grid_gpu_scale `
  --config experiments/scalable/grid_gpu_scale_config.json `
  --output experiments/results/grid_gpu_scale_rerun.json
```

This benchmark uses the same float64 CSR on both devices and C-contiguous dense
batches. Construction, CUDA context initialization, matrix/vector transfers, separate
SpMV/SpMM/local-iteration warmups, kernels, D2H, and end-to-end times are recorded.
The `single_ground` case is the grid Laplacian with one reference row/column removed;
the fixed-step Richardson recurrence is a local Neumann polynomial, not a claimed
optimal estimator.

The checked-in `grid_gpu_scale_v4.json` is the final warm-cache rerun;
`grid_gpu_scale_v3.json` retains a substantially colder CUDA/cuSPARSE startup. Use a
new output name for additional reruns so the recorded results remain auditable.
