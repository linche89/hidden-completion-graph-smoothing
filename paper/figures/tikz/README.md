# Paper figures: TikZ/PGFPlots

This folder contains the English-first, standalone sources for the paper's
theorem and experiment figures.  Every numerical coordinate comes from a CSV
file under `data/`; `data/manifest.json` pins the SHA-256 of every source JSON.
The shared `figure_style.tex` uses the Okabe--Ito colour-blind-safe palette,
serif LaTeX typography, explicit line styles and journal-scale labels.

## Regenerate data

From the repository root:

```powershell
python -m experiments.smoother.paper_data `
  --p3-p4 experiments/results/p3_p4_smoke.json `
  --external experiments/results/p4_external_scale.json `
  --core-evidence experiments/results/p4_core_evidence.json `
  --solver-crosscheck experiments/results/p4_solver_crosscheck.json `
  --gpu experiments/results/smoother_gpu_scale.json `
        experiments/results/smoother_gpu_scale_recheck1.json `
        experiments/results/smoother_gpu_scale_recheck2.json
```

## Compile and render

The required tools are XeLaTeX and Poppler's `pdftoppm`.  From this folder:

```powershell
$figures = @(
  'theorem_structure',
  'budget_pareto',
  'scenario_scaling',
  'core_design_benchmark',
  'sampling_budget_curve',
  'path_tightness',
  'local_method_tradeoff',
  'external_stress',
  'cpu_gpu_crossover'
)
foreach ($figure in $figures) {
  xelatex -interaction=nonstopmode -halt-on-error "$figure.tex"
  pdftoppm -png -r 200 -singlefile "$figure.pdf" $figure
}
```

From the repository root, `scripts/build_tsp.ps1` regenerates the data, runs
the regression suite, compiles the nine figures used by the TSP manuscript and
supplement, and builds both PDFs in one command. Compilation logs and `.aux`
files are intentionally ignored; editable `.tex`, CSV, vector PDF, and PNG
previews are retained.
