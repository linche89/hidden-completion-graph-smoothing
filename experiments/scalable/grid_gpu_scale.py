"""CPU/GPU scale benchmark for sparse 2-D grid SPD operators.

No NetworkX, dense system matrix, or inverse is used.  The stable family is a grid
graph Laplacian plus positive mass.  The stress family is the principal submatrix
obtained by grounding one corner node, so only one reference/anchor removes the
global gauge.  Starting from zero, k Richardson updates produce a polynomial of
degree at most k-1 in the matrix.  The literal loop schedules k SpMVs, although its
first multiplication is by the known zero vector and can be omitted in a message
passing implementation.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import time
from pathlib import Path
from typing import Any

import numpy as np
import psutil
import scipy
import scipy.sparse as sp


def grid_laplacian(side: int) -> sp.csr_matrix:
    """Direct sparse 4-neighbor graph Laplacian for a side-by-side grid."""
    if side < 2:
        raise ValueError("side must be at least two")
    n = side * side
    degree = np.full(n, 4.0, dtype=np.float64)
    degree[:side] -= 1.0
    degree[-side:] -= 1.0
    degree[::side] -= 1.0
    degree[side - 1 :: side] -= 1.0
    horizontal = -np.ones(n - 1, dtype=np.float64)
    horizontal[np.arange(1, n) % side == 0] = 0.0
    vertical = -np.ones(n - side, dtype=np.float64)
    matrix = sp.diags(
        (vertical, horizontal, degree, horizontal, vertical),
        (-side, -1, 0, 1, side),
        shape=(n, n),
        format="csr",
        dtype=np.float64,
    )
    matrix.eliminate_zeros()
    matrix.sort_indices()
    return matrix


def build_family(side: int, family: str, mass: float) -> tuple[sp.csr_matrix, dict[str, Any]]:
    laplacian = grid_laplacian(side)
    if family == "stable_mass":
        if mass <= 0:
            raise ValueError("stable_mass requires mass > 0")
        matrix = laplacian + mass * sp.eye(laplacian.shape[0], format="csr")
        lower_bound = mass
        upper_bound = 8.0 + mass
        step_size = 2.0 / (lower_bound + upper_bound)
        metadata = {
            "family": family,
            "mass": mass,
            "grounded_nodes": 0,
            "spectral_bounds_used": [lower_bound, upper_bound],
        }
    elif family == "single_ground":
        # Removing row/column zero fixes one reference angle.  The resulting
        # principal submatrix is SPD for the connected grid, but its smallest
        # eigenvalue shrinks with side.
        matrix = laplacian[1:, 1:].tocsr()
        upper_bound = 8.0
        step_size = 1.0 / upper_bound
        metadata = {
            "family": family,
            "mass": 0.0,
            "grounded_nodes": 1,
            "spectral_bounds_used": [None, upper_bound],
        }
    else:
        raise ValueError(f"unknown family {family}")
    matrix.eliminate_zeros()
    matrix.sort_indices()
    metadata["richardson_step_size"] = step_size
    return matrix, metadata


def csr_bytes(matrix: sp.csr_matrix) -> int:
    return int(matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes)


def cpu_fixed_richardson(matrix: sp.csr_matrix, rhs: np.ndarray, steps: int,
                         step_size: float) -> tuple[np.ndarray, float]:
    started = time.perf_counter()
    value = np.zeros_like(rhs)
    for _ in range(steps):
        value += step_size * (rhs - matrix @ value)
    return value, time.perf_counter() - started


def one_case(config: dict[str, Any], cp: Any, cupy_sparse: Any) -> dict[str, Any]:
    process = psutil.Process(os.getpid())
    host_rss_checkpoints = {"start": process.memory_info().rss}
    side = int(config["side"])
    family = str(config["family"])
    mass = float(config.get("mass", 0.2))
    batches = [int(x) for x in config["batches"]]
    spmv_repeats = int(config["spmv_repeats"])
    spmm_repeats = int(config["spmm_repeats"])
    local_steps = int(config["local_steps"])
    seed = int(config["seed"])

    build_started = time.perf_counter()
    matrix, family_metadata = build_family(side, family, mass)
    build_seconds = time.perf_counter() - build_started
    host_rss_checkpoints["after_matrix_build"] = process.memory_info().rss
    n = matrix.shape[0]
    rng = np.random.default_rng(seed)
    vector = np.ascontiguousarray(rng.standard_normal(n), dtype=np.float64)
    truth = np.ascontiguousarray(rng.standard_normal(n), dtype=np.float64)
    rhs = np.ascontiguousarray(matrix @ truth, dtype=np.float64)
    assert vector.flags.c_contiguous and truth.flags.c_contiguous and rhs.flags.c_contiguous
    host_rss_checkpoints["after_vectors"] = process.memory_info().rss

    # CPU SpMV kernel timing: output is consumed by a checksum after the loop.
    _ = matrix @ vector
    cpu_spmv_started = time.perf_counter()
    for _ in range(spmv_repeats):
        cpu_spmv_result = matrix @ vector
    cpu_spmv_seconds = time.perf_counter() - cpu_spmv_started
    cpu_spmv_checksum = float(np.sum(cpu_spmv_result))

    step_size = float(family_metadata["richardson_step_size"])
    _warm_cpu, cpu_iteration_warmup = cpu_fixed_richardson(matrix, rhs, 2, step_size)
    cpu_local, cpu_local_seconds = cpu_fixed_richardson(
        matrix, rhs, local_steps, step_size
    )
    cpu_local_error = float(np.linalg.norm(cpu_local - truth) / np.linalg.norm(truth))
    cpu_local_residual = float(
        np.linalg.norm(matrix @ cpu_local - rhs) / max(np.linalg.norm(rhs), 1e-300)
    )

    # Context/allocator initialization is reported separately from all kernels.
    context_started = time.perf_counter()
    cp.cuda.Device(0).use()
    _ = cp.empty(1, dtype=cp.float64)
    cp.cuda.Stream.null.synchronize()
    context_seconds = time.perf_counter() - context_started
    free_before, total_device = cp.cuda.runtime.memGetInfo()

    matrix_transfer_started = time.perf_counter()
    matrix_gpu = cupy_sparse.csr_matrix(matrix)
    cp.cuda.Stream.null.synchronize()
    matrix_h2d_seconds = time.perf_counter() - matrix_transfer_started
    free_after_matrix, _ = cp.cuda.runtime.memGetInfo()

    vector_transfer_started = time.perf_counter()
    vector_gpu = cp.asarray(vector)
    cp.cuda.Stream.null.synchronize()
    vector_h2d_seconds = time.perf_counter() - vector_transfer_started

    spmv_warmup_started = time.perf_counter()
    _ = matrix_gpu @ vector_gpu
    cp.cuda.Stream.null.synchronize()
    gpu_spmv_warmup = time.perf_counter() - spmv_warmup_started
    gpu_spmv_started = time.perf_counter()
    for _ in range(spmv_repeats):
        gpu_spmv_result = matrix_gpu @ vector_gpu
    cp.cuda.Stream.null.synchronize()
    gpu_spmv_seconds = time.perf_counter() - gpu_spmv_started
    spmv_d2h_started = time.perf_counter()
    gpu_spmv_host = cp.asnumpy(gpu_spmv_result)
    cp.cuda.Stream.null.synchronize()
    spmv_d2h_seconds = time.perf_counter() - spmv_d2h_started
    spmv_difference = float(
        np.linalg.norm(gpu_spmv_host - cpu_spmv_result)
        / max(np.linalg.norm(cpu_spmv_result), 1e-300)
    )

    # Warm up the complete local recurrence separately from SpMV and SpMM.
    rhs_transfer_started = time.perf_counter()
    rhs_gpu = cp.asarray(rhs)
    cp.cuda.Stream.null.synchronize()
    rhs_h2d_seconds = time.perf_counter() - rhs_transfer_started
    iteration_warmup_started = time.perf_counter()
    local_gpu = cp.zeros_like(rhs_gpu)
    for _ in range(2):
        local_gpu += step_size * (rhs_gpu - matrix_gpu @ local_gpu)
    cp.cuda.Stream.null.synchronize()
    gpu_iteration_warmup = time.perf_counter() - iteration_warmup_started
    gpu_local_started = time.perf_counter()
    local_gpu = cp.zeros_like(rhs_gpu)
    for _ in range(local_steps):
        local_gpu += step_size * (rhs_gpu - matrix_gpu @ local_gpu)
    cp.cuda.Stream.null.synchronize()
    gpu_local_seconds = time.perf_counter() - gpu_local_started
    local_d2h_started = time.perf_counter()
    gpu_local_host = cp.asnumpy(local_gpu)
    cp.cuda.Stream.null.synchronize()
    local_d2h_seconds = time.perf_counter() - local_d2h_started
    gpu_local_error = float(np.linalg.norm(gpu_local_host - truth) / np.linalg.norm(truth))
    gpu_local_residual = float(
        np.linalg.norm(matrix @ gpu_local_host - rhs) / max(np.linalg.norm(rhs), 1e-300)
    )

    batch_results: list[dict[str, Any]] = []
    for batch in batches:
        dense_cpu = np.ascontiguousarray(rng.standard_normal((n, batch)), dtype=np.float64)
        assert dense_cpu.flags.c_contiguous
        _ = matrix @ dense_cpu
        cpu_started = time.perf_counter()
        for _ in range(spmm_repeats):
            cpu_dense_result = matrix @ dense_cpu
        cpu_seconds = time.perf_counter() - cpu_started

        transfer_started = time.perf_counter()
        dense_gpu = cp.asarray(dense_cpu, order="C")
        cp.cuda.Stream.null.synchronize()
        h2d_seconds = time.perf_counter() - transfer_started
        if not dense_gpu.flags.c_contiguous:
            raise RuntimeError("GPU dense batch is not C contiguous")

        warmup_started = time.perf_counter()
        _ = matrix_gpu @ dense_gpu
        cp.cuda.Stream.null.synchronize()
        warmup_seconds = time.perf_counter() - warmup_started
        gpu_started = time.perf_counter()
        for _ in range(spmm_repeats):
            gpu_dense_result = matrix_gpu @ dense_gpu
        cp.cuda.Stream.null.synchronize()
        gpu_seconds = time.perf_counter() - gpu_started
        d2h_started = time.perf_counter()
        gpu_dense_host = cp.asnumpy(gpu_dense_result)
        cp.cuda.Stream.null.synchronize()
        d2h_seconds = time.perf_counter() - d2h_started
        difference = float(
            np.linalg.norm(gpu_dense_host - cpu_dense_result)
            / max(np.linalg.norm(cpu_dense_result), 1e-300)
        )
        batch_results.append(
            {
                "batch": batch,
                "layout": "C",
                "repeats": spmm_repeats,
                "dense_input_bytes": int(dense_cpu.nbytes),
                "cpu_kernel_seconds": cpu_seconds,
                "gpu_kernel_seconds": gpu_seconds,
                "kernel_speedup_cpu_over_gpu": cpu_seconds / gpu_seconds,
                "h2d_seconds": h2d_seconds,
                "warmup_seconds": warmup_seconds,
                "d2h_seconds": d2h_seconds,
                "gpu_end_to_end_excluding_matrix_preprocess": h2d_seconds
                + gpu_seconds
                + d2h_seconds,
                "relative_cpu_gpu_difference": difference,
            }
        )
        host_rss_checkpoints[f"batch_{batch}_before_release"] = process.memory_info().rss
        del dense_cpu, dense_gpu, cpu_dense_result, gpu_dense_result, gpu_dense_host
        cp.get_default_memory_pool().free_all_blocks()

    free_final, _ = cp.cuda.runtime.memGetInfo()
    matrix_device_bytes = int(free_before - free_after_matrix)
    result = {
        "case": config,
        "family_metadata": family_metadata,
        "matrix": {
            "shape": list(matrix.shape),
            "nnz": int(matrix.nnz),
            "dtype": str(matrix.dtype),
            "host_csr_bytes": csr_bytes(matrix),
            "measured_device_matrix_allocation_bytes": matrix_device_bytes,
        },
        "construction_seconds": build_seconds,
        "gpu_context_allocator_init_seconds": context_seconds,
        "matrix_h2d_preprocess_seconds": matrix_h2d_seconds,
        "memory": {
            "host_rss_checkpoints": host_rss_checkpoints,
            "host_rss_checkpoint_max": max(host_rss_checkpoints.values()),
            "gpu_total_bytes": int(total_device),
            "gpu_free_before_matrix": int(free_before),
            "gpu_free_after_matrix": int(free_after_matrix),
            "gpu_free_final": int(free_final),
        },
        "spmv": {
            "repeats": spmv_repeats,
            "vector_h2d_seconds": vector_h2d_seconds,
            "cpu_kernel_seconds": cpu_spmv_seconds,
            "gpu_warmup_seconds": gpu_spmv_warmup,
            "gpu_kernel_seconds": gpu_spmv_seconds,
            "gpu_d2h_seconds": spmv_d2h_seconds,
            "gpu_end_to_end_excluding_matrix_preprocess": vector_h2d_seconds
            + gpu_spmv_seconds
            + spmv_d2h_seconds,
            "kernel_speedup_cpu_over_gpu": cpu_spmv_seconds / gpu_spmv_seconds,
            "relative_cpu_gpu_difference": spmv_difference,
            "cpu_checksum": cpu_spmv_checksum,
        },
        "spmm": batch_results,
        "fixed_local_richardson": {
            "steps_rounds": local_steps,
            "step_size": step_size,
            "rhs_h2d_seconds": rhs_h2d_seconds,
            "cpu_warmup_seconds": cpu_iteration_warmup,
            "gpu_warmup_seconds": gpu_iteration_warmup,
            "cpu_kernel_seconds": cpu_local_seconds,
            "gpu_kernel_seconds": gpu_local_seconds,
            "gpu_d2h_seconds": local_d2h_seconds,
            "gpu_end_to_end_excluding_matrix_preprocess": rhs_h2d_seconds
            + gpu_local_seconds
            + local_d2h_seconds,
            "kernel_speedup_cpu_over_gpu": cpu_local_seconds / gpu_local_seconds,
            "cpu_relative_state_error": cpu_local_error,
            "gpu_relative_state_error": gpu_local_error,
            "cpu_relative_residual": cpu_local_residual,
            "gpu_relative_residual": gpu_local_residual,
            "relative_cpu_gpu_solution_difference": float(
                np.linalg.norm(cpu_local - gpu_local_host)
                / max(np.linalg.norm(cpu_local), 1e-300)
            ),
        },
    }
    del matrix_gpu, vector_gpu, rhs_gpu, local_gpu
    cp.get_default_memory_pool().free_all_blocks()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))

    import cupy as cp
    import cupyx.scipy.sparse as cupy_sparse

    run_started = time.perf_counter()
    global_warmup_started = time.perf_counter()
    tiny_cpu = sp.eye(32, format="csr", dtype=np.float64)
    tiny_gpu = cupy_sparse.csr_matrix(tiny_cpu)
    tiny_vector = cp.ones(32, dtype=cp.float64)
    _ = tiny_gpu @ tiny_vector
    cp.cuda.Stream.null.synchronize()
    global_cuda_cusparse_warmup_seconds = time.perf_counter() - global_warmup_started
    del tiny_gpu, tiny_vector
    cp.get_default_memory_pool().free_all_blocks()
    results = []
    for case in config["cases"]:
        case_started = time.perf_counter()
        result = one_case(case, cp, cupy_sparse)
        result["case_total_seconds"] = time.perf_counter() - case_started
        results.append(result)
    payload = {
        "schema_version": 1,
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cupy": cp.__version__,
            "gpu": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
            "cuda_runtime": int(cp.cuda.runtime.runtimeGetVersion()),
            "cuda_driver": int(cp.cuda.runtime.driverGetVersion()),
            "logical_cpus": os.cpu_count(),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
            "mpi": False,
        },
        "semantics": {
            "matrix": "same float64 CSR values and indices on CPU/GPU",
            "dense_layout": "C contiguous on CPU and GPU",
            "fixed_iteration": (
                "from x0=0, k Richardson updates give degree at most k-1 and "
                "depend on at most k-1 hops; the literal benchmark schedules k SpMVs"
            ),
            "timing": "construction, context, matrix H2D, per-operation warmup, kernels, D2H separated",
            "claim": "performance experiment only; no theorem or MPI claim",
        },
        "global_cuda_cusparse_lazy_init_warmup_seconds": global_cuda_cusparse_warmup_seconds,
        "total_seconds": time.perf_counter() - run_started,
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
