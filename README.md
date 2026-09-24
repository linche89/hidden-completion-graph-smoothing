# Hidden-Completion Graph Smoothing

Python implementation for **Exact Robust Local Graph-Smoother Design under
Hidden-Topology Completions**, by Che Lin.

[Manuscript](paper/tsp/manuscript.pdf) ·
[Supplement](paper/tsp/supplement.pdf) ·
[Solver guide](hidden_completion/README.md) ·
[Reproducibility](REPRODUCIBILITY.md)

Design a local linear approximation to the graph Tikhonov smoother
$(I+L_G)^{-1}$ when the remote topology and input signals are unknown.
The method reduces worst-case design to a finite semidefinite program,
with either a bounded or unbounded number of hidden nodes.

The model assumes undirected graphs, nonnegative edge weights, and unit
grounding. Exact scenario enumeration is intended for small visible port sets.

## Quick start

Requires **Python 3.11+**. The example runs on CPU.

```bash
git clone https://github.com/linche89/hidden-completion-graph-smoothing.git
cd hidden-completion-graph-smoothing
python -m pip install -e .
python -m hidden_completion solve experiments/configs/exact_finite_q2_h1.json
```

Returns JSON with the local operator `Q`, worst-case error radius
(approximately `0.471405` for this example), and solver diagnostics.
See the [solver guide](hidden_completion/README.md) for unbounded completions
and local-rule configuration.

## Reproducibility

Install the experiment dependencies and run the tests:

```bash
python -m pip install -e ".[experiments]"
python -m unittest discover -s experiments -p "test_*.py" -v
```

The [reproducibility guide](REPRODUCIBILITY.md) covers pinned environments,
experiment reruns, and paper builds. PDF builds additionally require LaTeX,
Latexmk, and Poppler.

## Repository

| Directory | Contents |
| --- | --- |
| [`hidden_completion/`](hidden_completion/) | Scenario enumeration and SDP solver |
| [`experiments/`](experiments/) | Tests, configurations, and recorded results |
| [`paper/`](paper/) | Manuscript, supplement, and figure sources |
| [`datasets/raw/`](datasets/raw/) | Benchmark data, provenance, and licenses |

## License and contact

Code: [MIT](LICENSE). The manuscript and original figures are excluded;
benchmark data retain their upstream licenses.

Che Lin · [linmo891104@gmail.com](mailto:linmo891104@gmail.com)
