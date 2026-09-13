"""Same-tolerance CPU/CuPy sparse SpMV and CG on one PGLib DC-WLS matrix."""

from __future__ import annotations

import argparse
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import scipy
import scipy.sparse.linalg as scipy_spla

from .pglib_dc import load_dc_wls


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tol", type=float, default=1e-8)
    parser.add_argument("--maxiter", type=int, default=20000)
    parser.add_argument("--spmv-repeats", type=int, default=200)
    parser.add_argument("--preconditioner", choices=["none", "jacobi"], default="jacobi")
    args = parser.parse_args()

    import cupy as cp
    import cupyx
    import cupyx.scipy.sparse as cupy_sparse
    import cupyx.scipy.sparse.linalg as cupy_spla

    instance, preprocessing = load_dc_wls(
        args.case,
        seed=20260913,
        flow_fraction=1.0,
        injection_fraction=0.25,
        angle_fraction=0.05,
        flow_holder="from",
    )
    j = instance.j
    rhs = instance.rhs
    probe = np.random.default_rng(7).standard_normal(j.shape[0])

    cpu_spmv_start = time.perf_counter()
    cpu_value = probe
    for _ in range(args.spmv_repeats):
        cpu_value = j @ probe
    cpu_spmv_seconds = time.perf_counter() - cpu_spmv_start

    cpu_iterations = 0
    def cpu_callback(_: np.ndarray) -> None:
        nonlocal cpu_iterations
        cpu_iterations += 1

    cpu_preconditioner = None
    diagonal_inverse = None
    if args.preconditioner == "jacobi":
        diagonal = j.diagonal()
        if np.any(diagonal <= 0):
            raise RuntimeError("Jacobi preconditioner requires a positive diagonal")
        diagonal_inverse = 1.0 / diagonal
        cpu_preconditioner = scipy_spla.LinearOperator(
            j.shape, matvec=lambda vector: diagonal_inverse * vector, dtype=j.dtype
        )
    cpu_cg_start = time.perf_counter()
    cpu_solution, cpu_info = scipy_spla.cg(
        j,
        rhs,
        rtol=args.tol,
        atol=0.0,
        maxiter=args.maxiter,
        M=cpu_preconditioner,
        callback=cpu_callback,
    )
    cpu_cg_seconds = time.perf_counter() - cpu_cg_start
    cpu_residual = float(np.linalg.norm(j @ cpu_solution - rhs) / np.linalg.norm(rhs))

    transfer_start = time.perf_counter()
    j_gpu = cupy_sparse.csr_matrix(j)
    rhs_gpu = cp.asarray(rhs)
    probe_gpu = cp.asarray(probe)
    gpu_preconditioner = None
    if diagonal_inverse is not None:
        diagonal_inverse_gpu = cp.asarray(diagonal_inverse)
        gpu_preconditioner = cupy_spla.LinearOperator(
            j_gpu.shape,
            matvec=lambda vector: diagonal_inverse_gpu * vector,
            dtype=j_gpu.dtype,
        )
    cp.cuda.Stream.null.synchronize()
    host_to_device_seconds = time.perf_counter() - transfer_start

    # Warmup is separated from measured kernel time.
    warmup_start = time.perf_counter()
    _ = j_gpu @ probe_gpu
    cp.cuda.Stream.null.synchronize()
    warmup_seconds = time.perf_counter() - warmup_start

    gpu_spmv_start = time.perf_counter()
    for _ in range(args.spmv_repeats):
        gpu_value = j_gpu @ probe_gpu
    cp.cuda.Stream.null.synchronize()
    gpu_spmv_seconds = time.perf_counter() - gpu_spmv_start

    gpu_iterations = 0
    def gpu_callback(_: cp.ndarray) -> None:
        nonlocal gpu_iterations
        gpu_iterations += 1

    gpu_cg_start = time.perf_counter()
    gpu_solution, gpu_info = cupy_spla.cg(
        j_gpu,
        rhs_gpu,
        rtol=args.tol,
        atol=0.0,
        maxiter=args.maxiter,
        M=gpu_preconditioner,
        callback=gpu_callback,
    )
    cp.cuda.Stream.null.synchronize()
    gpu_cg_seconds = time.perf_counter() - gpu_cg_start

    device_to_host_start = time.perf_counter()
    gpu_solution_host = cp.asnumpy(gpu_solution)
    gpu_value_host = cp.asnumpy(gpu_value)
    cp.cuda.Stream.null.synchronize()
    device_to_host_seconds = time.perf_counter() - device_to_host_start
    gpu_residual = float(np.linalg.norm(j @ gpu_solution_host - rhs) / np.linalg.norm(rhs))
    spmv_relative_difference = float(
        np.linalg.norm(gpu_value_host - cpu_value) / max(np.linalg.norm(cpu_value), 1e-300)
    )

    payload = {
        "schema_version": 1,
        "case": args.case.name,
        "matrix": {"shape": list(j.shape), "nnz": int(j.nnz)},
        "tolerance": args.tol,
        "maxiter": args.maxiter,
        "preconditioner": args.preconditioner,
        "spmv_repeats": args.spmv_repeats,
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cupy": cp.__version__,
            "cupyx": getattr(cupyx, "__version__", None),
            "cuda_runtime": int(cp.cuda.runtime.runtimeGetVersion()),
            "cuda_driver": int(cp.cuda.runtime.driverGetVersion()),
            "gpu": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
            "mpi": False,
        },
        "preprocessing_seconds": preprocessing,
        "transfer_seconds": {
            "host_to_device_matrix_and_vectors": host_to_device_seconds,
            "device_to_host_solution_and_spmv": device_to_host_seconds,
        },
        "warmup_seconds": warmup_seconds,
        "algorithm_seconds": {
            "cpu_spmv_total": cpu_spmv_seconds,
            "gpu_spmv_total": gpu_spmv_seconds,
            "cpu_cg": cpu_cg_seconds,
            "gpu_cg": gpu_cg_seconds,
        },
        "correctness": {
            "spmv_relative_cpu_gpu_difference": spmv_relative_difference,
            "cpu_cg_info": int(cpu_info),
            "gpu_cg_info": int(gpu_info),
            "cpu_cg_iterations": cpu_iterations,
            "gpu_cg_iterations": gpu_iterations,
            "cpu_relative_residual": cpu_residual,
            "gpu_relative_residual": gpu_residual,
        },
        "interpretation": (
            f"Only same-matrix SpMV and {args.preconditioner}-preconditioned CG are compared. "
            "Transfer, warmup, and preprocessing are separate. No MPI and no sparse direct "
            "GPU solve are claimed."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
