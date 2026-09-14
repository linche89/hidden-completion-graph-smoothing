# Paper build

`main.tex` is the scope-frozen formal proof draft.  It is organized around one
geometric main theorem and three application corollaries:

1. full-input Schatten reductions and the sharp `p=2` boundary;
2. the necessary-and-sufficient spectral SDP;
3. finite-budget convergence of certified risk.

Build from this directory with

```powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

The generated files live under `paper/build/` and are not tracked.  The paper
uses the established random-forest/DPP mean-projection identity as a cited
input and explicitly avoids claiming that identity as original.

## TikZ/PGFPlots figures

The publication figures live in `figures/tikz/`.  Labels, axes, legends and
captions are English-first.  Each figure is a standalone, editable `.tex`
source; numerical plots read generated CSV files rather than copied numbers.
The checked PDF and 200-dpi PNG preview are stored beside each source.

Regenerate the data from the repository root:

```powershell
python -m experiments.smoother.paper_data `
  --p3-p4 experiments/results/p3_p4_smoke.json `
  --external experiments/results/p4_external_scale.json `
  --gpu experiments/results/smoother_gpu_scale.json `
        experiments/results/smoother_gpu_scale_recheck1.json `
        experiments/results/smoother_gpu_scale_recheck2.json
```

Compile one figure from `paper/figures/tikz/` with:

```powershell
xelatex -interaction=nonstopmode -halt-on-error theorem_structure.tex
pdftoppm -png -r 200 -singlefile theorem_structure.pdf theorem_structure
```

See `figures/tikz/README.md`, `figures/tikz/CAPTIONS.md` and
`figures/tikz/QA.md` for the complete batch command and review record.
