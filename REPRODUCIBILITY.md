# Reproducibility guide

The repository separates theorem-matched hidden-completion experiments from
fixed-known-graph stress tests.  The former produce exact finite-scenario or
unbounded certificates.  The latter measure graph-smoother signal error and
resource use; they are not general WLS or Kalman MSE experiments.

## Checked environments

The committed CPU results were produced with Python 3.13.5 on Windows 11
(build 26100). Install the pinned numerical stack in the current project
environment:

```powershell
python -m pip install -r experiments/requirements-repro-cpu.txt
```

The paper figures were checked with MiKTeX 25.12 (XeTeX 4.16), Latexmk 4.87,
and Poppler 24.04.0.  The manuscript itself uses pdfLaTeX; each standalone TikZ
figure uses XeLaTeX and is included as vector PDF.

The committed accelerator crossover was measured on an NVIDIA GeForce RTX
5080 with Python 3.13.5, CuPy 14.2.0, CUDA runtime 12.9 (locally installed
toolkit 12.8), and NVIDIA driver API 13.02. Its optional dependencies are:

```powershell
python -m pip install -r experiments/scalable/requirements-repro-gpu.txt
```

GPU wall-clock values are hardware-specific.  The recorded comparison fixes
OpenMP and OpenBLAS to one CPU thread, uses the same float64 CSR operator on
both devices, and reports both resident-kernel and transfer-inclusive timings.

## One-command TSP verification

From the repository root, the checked TSP build runs the regression suite,
regenerates the paper CSV and macro files from committed JSON, builds the
standalone TikZ figures and PNG previews, and compiles both submission PDFs:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/build_tsp.ps1
```

This writes `output/pdf/tsp_manuscript.pdf` and
`output/pdf/tsp_supplement.pdf`. Use `-SkipTests` for a formatting-only build
or `-SkipData` to leave committed CSV/table artifacts untouched.

## Re-run the experiments

The seconds-scale P3/P4 suite is:

```powershell
python -m experiments.smoother `
  --config experiments/configs/p3_p4_smoke.json `
  --output experiments/results/p3_p4_smoke-rerun.json `
  --section all
```

The 135-completion hidden-budget/weight/cut audit is:

```powershell
python -m experiments.smoother.parameter_sweep `
  --config experiments/configs/p4_parameter_sweep.json `
  --output experiments/results/p4_parameter_sweep-rerun.json
```

The fixed 54-problem grid, 420-fit nested sampling study, and connected-path
tightness witnesses are reproduced with:

```powershell
python -m experiments.smoother.core_evidence `
  --config experiments/configs/p4_core_evidence.json `
  --output experiments/results/p4_core_evidence-rerun.json

python -m experiments.smoother.solver_crosscheck `
  --config experiments/configs/p4_solver_crosscheck.json `
  --output experiments/results/p4_solver_crosscheck-rerun.json
```

The second command independently reruns every primary-solver
`optimal_inaccurate` cell with high-accuracy SCS and compares full-scenario
re-evaluated radii.  It is a numerical cross-solver check, not a rigorous dual
lower bound.

The external scale and GPU commands are documented in
`experiments/README.md` and `experiments/scalable/README.md`.  Use new output
filenames so the checked records remain auditable.  Do not expect timing fields
to be bitwise reproducible; certificate values, scenario counts, residual
checks, and deterministic CSV transformations are the stable comparisons.

## External data provenance

- PGLib-OPF v23.07 case files are used only to obtain graph topology and
  positive branch weights.  The repository retains the upstream CC BY 4.0 data
  notice and per-case attribution metadata.
- SuiteSparse matrices are converted to undirected sparsity graphs.  The
  repository retains the CC BY 4.0 text and source metadata; this transformation
  is stated in the results and paper.
- The benchmark then constructs its own unit-grounded smoother `I + L`.  It does
  not reproduce or make claims about either dataset's original estimator.

The project code and manuscript must receive an explicit project-wide license
before an archival release. Third-party data are not covered by that future
license and retain the terms recorded in their dataset directories.
