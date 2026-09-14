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
