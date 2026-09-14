# Anonymous TSP software artifact

This artifact contains the exact-scenario solver, theorem-matched experiments,
checked numerical records, LaTeX sources, and vector-figure inputs used by the
anonymous manuscript. It intentionally excludes the Git repository/history,
author metadata,
internal research notes, and unrelated literature files.

## Environment and checks

From the artifact root on Windows PowerShell:

```powershell
python -m pip install -r experiments/requirements-repro-cpu.txt
python -m unittest discover -s experiments -p "test_*.py" -v
```

The primary evidence and independent cross-solver audit can be rerun with new
output filenames:

```powershell
python -m experiments.smoother.core_evidence `
  --config experiments/configs/p4_core_evidence.json `
  --output experiments/results/p4_core_evidence-rerun.json

python -m experiments.smoother.solver_crosscheck `
  --config experiments/configs/p4_solver_crosscheck.json `
  --output experiments/results/p4_solver_crosscheck-rerun.json
```

The cross-solver audit is numerical corroboration; it does not claim a rigorous
dual lower bound or uniqueness of the returned optimizer. To regenerate data,
English TikZ figures, and both IEEEtran PDFs from the archived records:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/build_tsp.ps1
```

`MANIFEST.sha256` records every packaged file and is generated before the ZIP
archive is written.
