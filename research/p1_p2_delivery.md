# P1/P2 delivery status (scope frozen)

> Status date: 2026-09-14
> Baseline retained: `a97e62b` / `theorem-baseline-2026-09-14`
> Theory scope: unchanged. No directed, signed, block-valued, dense-SPD, or
> nonlinear extension was added.

## P1: formal mathematical paper draft

- [x] Created `paper/main.tex` and `paper/references.bib`.
- [x] Compressed the result into one geometric main theorem and three
  corollaries.
- [x] Defined completion classes, common local rules, closure, and completion
  suprema before the main theorem.
- [x] Wrote the forest-mixture upper bound as a cited known input.
- [x] Wrote the connected maximum-degree-two path converse.
- [x] Proved the complete finite-`h` exposed-vertex characterization.
- [x] Proved the finite full-input all-Schatten reduction with the canonical
  hidden-input block `C(X-X^2)^(1/2)`.
- [x] Stated and proved the sharp unbounded Schatten `p=2` boundary and its
  two-vertex counterexample.
- [x] Derived finite and unbounded necessary-and-sufficient spectral SDPs.
- [x] Separated graph-smoother error from generic WLS/Kalman MSE claims.
- [x] Consolidated references and explicitly attributed the matrix-forest,
  random-forest, mean-projection, trace-convexity, and spectral-LMI tools.
- [x] Compiled a 10-page PDF with no undefined references, layout warnings,
  overfull boxes, or underfull boxes; visually checked every rendered page.

Build command:

```powershell
Set-Location paper
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

## P2: exact scenario solver

- [x] Canonical set-partition enumeration.
- [x] Weak hidden-budget allocation enumeration.
- [x] Exact pruning: retain `k=0` and `sum(k)=h` for finite spectral risk.
- [x] Construction of `X_(pi,k)` and the canonical full hidden-input error.
- [x] Lossless spectral-norm SDP in CVXPY.
- [x] `Q_loc` constraints for explicit zero patterns, `r`-hop access masks,
  shared coefficients, fixed entries, generic affine equalities, and row sums.
- [x] JSON output of every active worst-case scenario and its exact
  post-solve risk.
- [x] Small `(q,h)` exhaustive count/uniqueness regression tests.
- [x] Analytic regression for the finite `q=2,h=1` optimum and the unbounded
  shared-opcode optimum.
- [x] Machine-readable CLI plus installable package metadata.

Commands:

```powershell
python -m pip install -e .
python -m unittest discover -v
python -m hidden_completion solve experiments/configs/exact_finite_q2_h1.json
python -m hidden_completion solve experiments/configs/exact_unbounded_shared_opcode.json
```

The SDP is an exact finite reformulation of the topology supremum. Returned
decimal values are numerical conic-solver approximations; the implementation
re-evaluates the exact scenario formula at the returned `Q` and reports both
that radius and the solver objective.

## Verification snapshot

- Full discovered unit-test suite: 24 tests passed.
- Frozen finite-budget theorem checks: 6 tests passed.
- Hidden partial-partition stress checks: 5 checks passed.
- Sharp-boundary/counterexample checks: 6 checks passed.
- Finite analytic radius: `0.47140452079113604`.
- Shared-opcode analytic radius: `1.088303688353231`.

At the time of this P1/P2 delivery, P3--P5 were intentionally outside its
scope.  Their subsequently completed initial implementation is tracked in
`research/p3_p5_delivery.md`.
