# Reproducibility guide

The repository separates theorem-matched hidden-completion experiments from
fixed-known-graph stress tests.  The former produce exact finite-scenario or
unbounded certificates.  The latter measure graph-smoother signal error and
resource use; they are not general WLS or Kalman MSE experiments.

## Checked environments

The committed CPU results were produced with Python 3.13.5 on Windows 11
(build 26100).  Install the exact numerical stack in an isolated environment:

```powershell
python -m venv .venv-repro
.venv-repro/Scripts/python -m pip install --upgrade pip
.venv-repro/Scripts/python -m pip install -r experiments/requirements-repro-cpu.txt
```

The paper figures were checked with MiKTeX 25.12 (XeTeX 4.16), Latexmk 4.87,
and Poppler 24.04.0.  The manuscript itself uses pdfLaTeX; each standalone TikZ
figure uses XeLaTeX and is included as vector PDF.

The committed accelerator crossover was measured on an NVIDIA GeForce RTX
5080 with Python 3.13.5, CuPy 14.2.0, CUDA runtime 12.9 (locally installed
toolkit 12.8), and NVIDIA driver API 13.02.  Install its isolated stack with:

```powershell
python -m venv experiments/.venv-gpu
experiments/.venv-gpu/Scripts/python -m pip install --upgrade pip
experiments/.venv-gpu/Scripts/python -m pip install -r experiments/scalable/requirements-repro-gpu.txt
```

GPU wall-clock values are hardware-specific.  The recorded comparison fixes
OpenMP and OpenBLAS to one CPU thread, uses the same float64 CSR operator on
both devices, and reports both resident-kernel and transfer-inclusive timings.

## One-command paper verification

From the repository root, the following command runs the 35 regression tests,
regenerates every paper CSV and the ablation table from committed JSON, builds
the seven standalone TikZ figures and PNG previews, and compiles the manuscript:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/build_paper.ps1
```

Use `-SkipTests` for a formatting-only build or `-SkipData` to leave committed
CSV/table artifacts untouched.  The final PDF is `paper/build/main.pdf`.

## Re-run the experiments

The seconds-scale P3/P4 suite is:

```powershell
.venv-repro/Scripts/python -m experiments.smoother `
  --config experiments/configs/p3_p4_smoke.json `
  --output experiments/results/p3_p4_smoke-rerun.json `
  --section all
```

The 135-completion hidden-budget/weight/cut audit is:

```powershell
.venv-repro/Scripts/python -m experiments.smoother.parameter_sweep `
  --config experiments/configs/p4_parameter_sweep.json `
  --output experiments/results/p4_parameter_sweep-rerun.json
```

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

See `release/PUBLIC_RELEASE_AUDIT.md` before publishing any repository history.
