# Exact hidden-completion scenario solver

This package implements the lossless finite-scenario reductions stated in
`paper/tsp/main.tex`. It solves

\[
\min_{Q\in\mathcal Q_{\rm loc}}\ \sup_G
\left\|C[I_q\ 0](I+L_G)^{-1}-Q[I_q\ 0]\right\|_2.
\]

Two completion classes are supported:

- `finite`: at most `h` unit-grounded hidden vertices;
- `unbounded`: any finite number of unit-grounded hidden vertices, with the
  supremum reduced to partial-partition projections.

The numerical optimization uses CVXPY, but the formulation itself is exact:
there is no topology sampling or uncertainty relaxation.  Floating-point
answers retain the tolerance of the selected conic solver.

## Install and test

```powershell
python -m pip install -r hidden_completion/requirements.txt
# Or install the package and `hidden-completion` console command:
python -m pip install -e .
python -m unittest experiments.test_exact_scenario_solver -v
```

## Solve a configured problem

```powershell
python -m hidden_completion solve experiments/configs/exact_finite_q2_h1.json
python -m hidden_completion solve experiments/configs/exact_unbounded_shared_opcode.json
```

The JSON result contains the optimal `Q`, the directly re-evaluated exact
scenario radius, scenario counts, solver timing, and every active worst-case
scenario.  Enumerate scenarios without solving an SDP via

```powershell
python -m hidden_completion enumerate --q 3 --h 2 --vertices-only
python -m hidden_completion enumerate --q 3
```

The first command returns the 21 exact vertices of the 31-candidate
finite-`h` set; the second returns the 15 partial partitions (`B_4`).

## Local-rule schema

`local_rule` is the intersection of any supplied affine restrictions:

- `allowed_mask`: Boolean matrix; `false` forces the corresponding entry of
  `Q` to zero;
- `zeros`: list of `[row, column]` entries forced to zero;
- `r_hop`: an `adjacency`, one port index per `output_nodes` row, and a
  nonnegative `radius`; entries outside the accessible balls are zero;
- `shared_coefficients`: lists of matrix entries that must have one common
  coefficient;
- `fixed_entries`: `[row, column, value]` triples (object form is also
  accepted);
- `linear_equalities`: generic affine equalities with
  `[row, column, coefficient]` terms and an `rhs`;
- `row_sums`: one prescribed value per output row, or one shared scalar.

For example, the following describes `Q = alpha I_2`:

```json
{
  "zeros": [[0, 1], [1, 0]],
  "shared_coefficients": [[[0, 0], [1, 1]]]
}
```

Finite spectral optimization automatically uses the proved exact pruning
rule: keep `k=0` and allocations with total `h`.  Set
`"prune_vertices": false` only to audit the redundant formulation.  Do not
reuse that pruning rule for finite Schatten exponents below two: the terminal
hull remains exact, but the full-input risk need not be convex in `X`.
