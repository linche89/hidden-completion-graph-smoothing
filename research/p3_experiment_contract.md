# P3 experiment contract: graph-regularized state smoothing

> Frozen on: 2026-09-14
> Scope: experimental semantics only; the mathematical completion class is unchanged.

## 1. Application and two experiment tracks

The primary application is the scalar, unit-grounded graph smoother

\[
  x_G(y)=(I+L_G)^{-1}y
  =\arg\min_x\left\{
  \frac12\lVert x-y\rVert_2^2+
  \frac12\sum_{\{i,j\}\in E}w_{ij}(x_i-x_j)^2
  \right\}.
\]

Every reported result belongs to exactly one of the following tracks.

1. `theorem_matching`: the first `q` vertices are visible ports, at most `h`
   unit-grounded vertices are hidden, and the common rule is
   `Q E_Gamma y`.  Certified error means the full-input operator norm
   \(\sup_G\lVert C E_\Gamma(I+L_G)^{-1}-Q E_\Gamma\rVert_2\).
2. `external_application_stress`: one fully specified graph is known and the
   target is its exact smoother.  Errors are fixed-instance operator error or
   empirical graph-signal error.  PGLib contributes topology and positive
   branch weights; SuiteSparse contributes a sparsity graph.  Neither source
   is represented as a hidden-completion theorem instance.

In particular, no smoother error in this repository is called a generic WLS
or Kalman measurement-noise MSE.  The older DC-WLS benchmarks remain separate.

## 2. `r`-round local rule

Let `H=(Gamma,E_H)` be the undirected communication graph on the visible
ports, and attach output row `a` to port `o(a)`.  A matrix is `r`-round
implementable when

\[
  Q_{aj}=0\quad\text{if}\quad d_H(o(a),j)>r.
\]

This is a support/causality condition.  It does not assert that every such
matrix has a unique distributed implementation.  The reference online
accounting uses deterministic lexicographic shortest-path unicast: each
nonzero dependency sends one scalar from `j` to `o(a)`.  Consequently:

- rounds are the largest used distance;
- message count equals scalar-hop count for scalar packets;
- byte-hop is scalar-hop count times the configured bytes per scalar;
- arithmetic for row `a` with `k_a` nonzeros is `k_a` multiplies and
  `max(k_a-1,0)` additions;
- coefficient storage is `nnz(Q)` scalars plus sparse indices.

This explicit schedule is a reproducible upper bound, not a communication
lower bound.  Shared-opcode constraints are a second application: selected
entries of `Q` use one common coefficient while retaining the same support
and communication accounting.

## 3. Iterative and truncation protocols

For a fixed known graph, a degree-`r` polynomial in `I+L` uses `r` neighbor
SpMVs and hence at most `r` communication hops.  Each SpMV exchanges one
scalar in each direction on every undirected edge for each right-hand side.
The benchmark reports directed messages, byte-hop, SpMV count and an explicit
FLOP estimate.

The ordinary truncation baseline gathers the induced radius-`r` subgraph and
solves its unit-grounded smoother.  Its model-collection payload and local
dense factor/solve estimates are reported separately.  This baseline may use
more local computation than a degree-`r` polynomial even though both have the
same causal radius.

## 4. Offline versus online cost

Offline design includes scenario enumeration, vertex pruning, conic model
construction, SDP solution, polynomial coefficient construction and reusable
local factorizations.  Online cost includes input communication, application
of the already designed rule, iterative SpMVs and per-signal solves.  Wall
time, peak sampled RSS, matrix storage and arithmetic counts are never merged
into one number.

The exact scenario solver is boundary-size fixed-parameter: its offline cost
depends on `q` and `h`, not on a subsequently realized hidden graph.  Online
application of `Q` depends on `nnz(Q)` and the communication graph, not on the
number of enumerated scenarios.

## 5. Baselines and certificate meanings

- `exact_finite_h`: lossless finite-`h` scenarios and exact vertex pruning.
- `exact_unbounded`: lossless partial-partition certificate.
- `random_completion_sampling`: minimax design over finitely sampled physical
  graphs; its training radius is not a certificate.  Every returned rule is
  re-evaluated against the exact scenario sets.
- `psd_contraction_relaxation`: replace the terminal hull by
  \(0\preceq X\preceq I\).  The spectral robust problem is solved losslessly
  through the equivalent norm-bounded full-block LMI.  It is a certified but
  generally conservative outer relaxation.
- `endpoint_only`: checks only the zero and identity projections.  It is a
  deliberate NO-GO ablation and is never labeled safe.
- `global_exact`, `induced_local_truncation`, `neumann` and `chebyshev` are
  fixed-known-graph methods in the external stress track.

## 6. Required result fields

Machine-readable results record the normalized configuration, RNG seeds,
code commit and dirty flag, platform and package versions.  Robust-design
records include training error, exact finite and unbounded certified radii,
active worst scenarios, scenario/model/solver time, pruning counts, sampled
RSS, empirical graph-signal MSE, hidden count, edge-weight range and port cut.

Fixed-graph records include method access class, rounds, messages, byte-hop,
FLOPs, runtime, stored bytes, empirical MSE, relative signal error and—when
the instance is small enough—the exact fixed-instance operator norm.  External
datasets carry their source SHA-256 and transformation description.
