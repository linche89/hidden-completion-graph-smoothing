# Suggested English captions

## `theorem_structure`

Logical structure of the result.  The forest-mixture identity supplies the
upper inclusion, while connected maximum-degree-two path completions expose
the listed vertices and give the converse.  The exact finite and unbounded
completion hulls then imply the full-input Schatten reduction, the lossless
spectral SDP, and the sharp finite-budget rate.

## `budget_pareto`

Finite hidden-budget hierarchy for the two-port fixed-rule audit.  (a) The
observed gap between finite-budget and unbounded squared risks agrees with the
exact operator-Hausdorff distance (q/(q+h)) for every tested budget.  (b) A
larger hidden budget increases both the certified radius and the offline
scenario count; exact vertex pruning removes redundant LMIs without changing
the certificate.

## `certificate_sampling_failure`

Certificate comparison for (q=3), (h=2), a one-hop path support and shared
diagonal coefficients.  Physical training and held-out maxima do not certify
the completion class.  In particular, the rule designed from 12 random
training completions understates its exact finite-budget worst case by 0.295.
The PSD-contraction result is a safe outer certificate (reported separately
from the bars), whereas the endpoint-only construction is a deliberate NO-GO
ablation.

## `local_method_tradeoff`

Fixed-known-graph accuracy and communication tradeoff, averaged over random
cyclic, grid, geometric and small-world graphs with 36 nodes and 16 signal
draws per graph.  Induced local truncation, Neumann iteration and Chebyshev
iteration share a causal radius but have different communication and local
computation models.  The vertical quantity is empirical graph-signal MSE, not
general WLS or Kalman measurement-noise MSE.

## `scenario_scaling`

Exact scenario counts at hidden budget (h=2).  Finite candidates, retained
finite-hull vertices and unbounded partial partitions all grow
combinatorially with boundary size (q).  The method is fixed-parameter in
the visible boundary, not polynomial in (q).

## `cpu_gpu_crossover`

CPU/GPU crossover for 40 rounds of the fixed local Richardson smoother on
unit-grounded two-dimensional grids.  Curves show medians over three process
runs; the blue band spans the observed resident-kernel range.  On the tested
single-thread CPU and RTX 5080, the GPU becomes faster between 10,000 and
50,176 states, and vector-transfer overhead does not change the crossover.

## `external_stress`

External graph-smoother stress tests.  (a) Relative signal error after 16
Chebyshev rounds and (b) sparse global-reference solve time.  PGLib cases
contribute topology and positive branch weights, while SuiteSparse cases
contribute symmetric sparsity graphs; neither transformation is presented as
the original dataset's WLS model or as a theorem-matched hidden completion.
