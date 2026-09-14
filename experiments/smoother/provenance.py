"""Environment and Git provenance for machine-readable experiment records."""

from __future__ import annotations

import os
import platform
from pathlib import Path
import subprocess
import sys
from typing import Any

import clarabel
import cvxpy
import networkx
import numpy
import psutil
import scipy


def _git(arguments: list[str], cwd: Path) -> str | None:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def experiment_provenance(repository_root: str | Path) -> dict[str, Any]:
    root = Path(repository_root).resolve()
    status = _git(["status", "--porcelain"], root)
    return {
        "repository_root": str(root),
        "git_commit": _git(["rev-parse", "HEAD"], root),
        "git_dirty": bool(status) if status is not None else None,
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "packages": {
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "networkx": networkx.__version__,
            "cvxpy": cvxpy.__version__,
            "clarabel": clarabel.__version__,
            "psutil": psutil.__version__,
        },
        "logical_cpus": os.cpu_count(),
        "thread_environment": {
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
        },
    }
