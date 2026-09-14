# P3--P5 delivery status (theory scope frozen)

> Status date: 2026-09-14
>
> Theory baseline: `a97e62b` / `theorem-baseline-2026-09-14`
>
> Initial experimental release line: `1157898` through `6d1e183`

No directed, signed, block-valued, dense-SPD, or nonlinear extension was added.
Theorem-matched certificates and fixed-known-graph stress tests remain
different experiment classes throughout the code, records, figures, and paper.

## P3: primary application model

- [x] Fixed the application as unit-grounded graph-regularized state smoothing.
- [x] Defined the `r`-round zero pattern on the visible communication graph.
- [x] Implemented rounds, directed messages, scalar-hop, byte-hop, FLOPs, and
  sparse coefficient-memory accounting for an explicit schedule.
- [x] Separated offline scenario/model/solver cost from online rule execution.
- [x] Included shared-coefficient/shared-opcode rules as the constrained second
  application.
- [x] Kept graph-signal smoother MSE distinct from generic WLS/Kalman MSE.

The frozen semantics and accounting contract are in
`research/p3_experiment_contract.md`.

## P4: numerical experiments

### Methods

- [x] Global exact smoother on a fixed known graph.
- [x] Ordinary induced-neighborhood truncation.
- [x] Neumann and Chebyshev local polynomial iterations.
- [x] Random physical-completion sampling with exact post-evaluation.
- [x] Safe but coarse `0 <= X <= I` contraction relaxation.
- [x] Exact finite-`h` scenario certificate.
- [x] Exact unbounded partial-partition certificate.

### Data and parameter coverage

- [x] Random cyclic, grid, geometric, and small-world synthetic graphs.
- [x] A 135-completion factorial audit over hidden budgets
  `{0,1,2,4,8}`, weak/moderate/strong edge weights, and three graph densities.
  The realized positive-hidden port cuts span 2--15 edges and total cut weight
  `0.0344`--`508.44`.
- [x] Six PGLib-OPF v23.07 topology transformations and three SuiteSparse
  sparsity-graph transformations, up to 129,164 vertices.
- [x] Explicit theorem-matched versus external-stress labels and provenance.

### Metrics and checks

- [x] Worst-case certified operator error and active scenarios.
- [x] Empirical signal error/MSE.
- [x] Communication rounds and byte-hop.
- [x] CPU/GPU runtime and transfer-inclusive timing.
- [x] Sampled RSS/peak-memory fields.
- [x] Scenario generation, conic modeling, SDP, and verification time.
- [x] Vertex-pruning counts, timing, and radius-equivalence check.
- [x] Exact `q/(q+h)` budget-convergence data and figure.

Verification snapshot:

- 35 discovered regression tests pass.
- The factorial sweep has 135/135 expected rows, zero certificate violation,
  and multiple realized port-cut sizes.
- Canonical sweep result:
  `experiments/results/p4_parameter_sweep.json`.
- Canonical sweep SHA-256:
  `4847a787c0c47f9f48c4f072a57e1033d74e86d5de9769d7318a938d8faab115`.

## P5: paper figures and release

- [x] English-first TikZ main-theorem structure figure.
- [x] Hidden-budget/error/scenario-count Pareto figure.
- [x] Exact-certificate versus random-sampling failure figure.
- [x] CPU/GPU crossover figure.
- [x] Ablation/NO-GO table and additional scaling/stress figures.
- [x] Editable TikZ/PGFPlots sources, CSV provenance manifest, vector PDFs,
  PNG previews, and visual-QA record.
- [x] Exact CPU/GPU environment files and one-command paper build.
- [x] Public-release redistribution audit.
- [ ] Replace the Git publication email and amend/export the intended public
  baseline after the author supplies the address.
- [ ] Select source-code and manuscript licenses.
- [ ] Create the `github.com/linche89/...` repository after the author supplies
  the repository name and confirms publication.
- [ ] Re-check the target venue's preprint/PDF policy once the venue is fixed.

The current private history must not be pushed directly: it contains 69
third-party literature PDFs.  The safe release route is a clean-history export
described in `release/PUBLIC_RELEASE_AUDIT.md`.
