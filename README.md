# Exact Robust Local Graph-Smoother Design Under Hidden-Topology Completions

Python reference implementation and reproducibility artifact for the
manuscript **“Exact Robust Local Graph-Smoother Design Under Hidden-Topology
Completions.”** The manuscript is being prepared for submission to *IEEE
Transactions on Signal Processing*.

The project designs one port-local linear approximation to the graph
Tikhonov smoother

$$
    (I+L_G)^{-1}
$$

when the topology and inputs beyond a labeled visible port set are unknown.
For a finite hidden-node budget, the infinite completion family is reduced
exactly to scenarios indexed by set partitions and hidden-budget allocations.
The unbounded model reduces to partial partitions. Spectral-norm design under
semidefinite-representable locality rules is then solved by a finite SDP.

## Scope

The exact method is fixed-parameter in the number of visible ports `q`; the
hidden physical network itself may be arbitrarily large. It is intended for a
modest observed boundary and is not presented as a general fast smoother for
unrestricted large `q`.

The implementation supports:

- finite-`h` and unbounded hidden-completion scenario enumeration;
- exact finite-budget vertex pruning;
- spectral-norm SDP design with zero patterns, `r`-hop masks, affine
  equalities, fixed entries, row sums, and shared coefficients;
- active worst-case scenario reporting;
- theorem-matched numerical studies and fixed-known-graph stress tests.

## Quick start

Python 3.11 or newer is required.

```powershell
python -m pip install -e ".[experiments]"
python -m unittest discover -s experiments -p "test_*.py" -v
```

Solve representative finite-budget and unbounded problems:

```powershell
python -m hidden_completion solve experiments/configs/exact_finite_q2_h1.json
python -m hidden_completion solve experiments/configs/exact_unbounded_shared_opcode.json
```

Each result reports the optimized local operator, its re-evaluated exact
radius, scenario counts, solver diagnostics, and active worst-case topologies.
See [hidden_completion/README.md](hidden_completion/README.md) for the complete
configuration schema.

## Reproduce the TSP manuscript

Install the pinned CPU dependencies, then run the checked build:

```powershell
python -m pip install -r experiments/requirements-repro-cpu.txt
powershell -ExecutionPolicy Bypass -File scripts/build_tsp.ps1
```

The build runs the regression suite, regenerates the paper data, compiles the
standalone TikZ figures, and writes:

- `output/pdf/tsp_manuscript.pdf`;
- `output/pdf/tsp_supplement.pdf`.

MiKTeX/XeLaTeX, Latexmk, and Poppler are needed for the full paper build. The
Python experiments can be run without a TeX installation. Detailed commands
and checked environments are in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Repository layout

- `hidden_completion/`: exact scenario generator and SDP solver;
- `experiments/`: regression tests, experiment runners, configurations, and
  checked numerical records;
- `paper/tsp/`: TSP manuscript and supplement;
- `paper/figures/tikz/`: editable TikZ/PGFPlots sources and generated data;
- `datasets/raw/`: licensed third-party stress-test data and provenance;
- `scripts/`: reproducible paper-build and artifact-packaging commands.

## External benchmark data

The `.m` files under `datasets/raw/pglib_opf_v23.07/` are upstream PGLib-OPF
cases in MATPOWER data format. They are **not MATLAB implementation code** and
are parsed directly by the Python experiment loader; MATLAB and Octave are not
required. PGLib-OPF and SuiteSparse data remain under the licenses and
attribution notices stored in their respective dataset directories.

## Author

Che Lin — <linmo891104@gmail.com>
