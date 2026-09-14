# Research TODO

## Current TSP paper: frozen scope

- [x] Target exact robust local graph-smoother design for a modest observed
  port boundary `q`.
- [x] Allow the hidden network to be arbitrarily large in the unbounded
  completion model; small `q` does **not** mean a small physical graph.
- [x] Keep fixed-known-graph solvers and hardware timing as secondary stress
  tests rather than presenting them as the paper's new algorithm.
- [ ] Complete author metadata, EDICS selection, bibliography verification,
  final licenses, and submission-system checks.
- [ ] Freeze and tag the submitted TSP source, result records, and anonymous
  review artifact before starting the project below.

Scope guardrail: do not expand the current manuscript into a general
large-boundary or full-graph solver paper. Its exact contribution is
fixed-parameter in the observed boundary size and independent of the eventual
hidden-network size.

## Next project: scalable robust design for large observed boundaries

Working title: **Scalable Certified Graph-Smoother Design Under Hidden
Completions**.

### Core problem

Develop a method for large `q` that avoids explicit enumeration in

```text
minimize_Q  max_{P in partial partitions of [q]} || C P - Q ||_2
```

and, subsequently, in the finite-budget problem over `(partition,
hidden-allocation)` scenarios. The output should include a local rule, a valid
worst-case upper certificate, a computable lower bound, and an optimality gap.

Important distinction:

- large hidden-node count with modest `q`: already covered exactly by the TSP
  paper;
- large `q`: open scalability problem;
- a known full graph with all-node output: mature sparse-solver/graph-filter
  problem and not the intended contribution.

### P0: complexity and separation audit

- [ ] Formulate the exact worst-partition separation oracle for a fixed `Q`.
- [ ] Determine whether the general separation problem is NP-hard or admits a
  polynomial formulation. Exponentially many exposed vertices alone is not a
  hardness proof.
- [ ] Classify tractable cases: path/tree protocol overlays, bounded treewidth,
  diagonal or low-rank `C`, shared coefficients, fixed output rank, and
  symmetry classes.
- [ ] Separate the difficulty of choosing a partition from the finite-`h`
  hidden-budget allocation subproblem.

### P1: exact delayed-constraint prototype

- [ ] Implement a master SDP containing only discovered scenarios.
- [ ] Implement an exact small-`q` oracle by reusing the current enumerator.
- [ ] Add branch-and-bound, symmetry/orbit reduction, duplicate elimination,
  warm starts, and active-scenario caching.
- [ ] Stop only with a certified separation bound; do not label a heuristic
  search as exact.
- [ ] Cross-check every result against complete enumeration for `q <= 8`.

### P2: scalable certified approximation

- [ ] Develop a relaxation of the partition oracle with a rigorous upper
  bound and a feasible rounding/lower-bound construction.
- [ ] Report the resulting robust-design optimality gap, not only empirical
  test error.
- [ ] Explore column generation, cutting planes, mixed-integer conic models,
  dynamic programming on structured overlays, and randomized rounding.
- [ ] Preserve zero-pattern, `r`-hop, row-sum, and shared-coefficient local
  rule support.

### P3: experiments

- [ ] Use the current exact solver as ground truth for small `q`.
- [ ] Scale through `q = 10, 20, 50, 100` where feasible.
- [ ] Measure oracle calls, generated/active scenarios, certified gap, runtime,
  peak memory, and warm-start benefit.
- [ ] Compare full enumeration, random completion sampling, the generic
  positive-contraction relaxation, delayed constraints, and the proposed
  certified approximation.
- [ ] Keep theorem-matched hidden-completion evidence separate from known-graph
  PGLib/SuiteSparse stress tests.

### Launch criterion

Start this as a new repository or a branch created from the frozen TSP tag.
The first milestone is not a large experiment: it is a defensible complexity
result or a certified separation oracle. Until that milestone is reached, do
not promise unrestricted exact large-`q` scalability.
