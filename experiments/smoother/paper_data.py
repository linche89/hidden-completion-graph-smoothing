"""Export checked experiment results for TikZ/PGFPlots paper figures.

The exporter deliberately performs no plotting.  It converts the canonical
JSON experiment records into small CSV files consumed by the standalone TikZ
sources under ``paper/figures/tikz`` and emits the LaTeX ablation table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics
from typing import Any, Iterable, Mapping, Sequence


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _number(value: Any) -> str:
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return format(value, ".16g")
    return str(value)


def _write_csv(
    destination: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _number(row[key]) for key in fieldnames})


def _export_budget(robust: Mapping[str, Any], data_directory: Path) -> None:
    rows = robust["budget_sweep"]["rows"]
    fields = (
        "h",
        "hausdorff",
        "risk_gap",
        "candidates",
        "vertices",
        "certified_radius",
    )
    _write_csv(
        data_directory / "budget_pareto.csv",
        fields,
        (
            {
                "h": row["h"],
                "hausdorff": row["exact_hausdorff_distance"],
                "risk_gap": row["fixed_Q_squared_risk_gap"],
                "candidates": row["candidate_scenarios"],
                "vertices": row["retained_vertices"],
                "certified_radius": row["optimized_finite_radius"],
            }
            for row in rows
        ),
    )


def _export_certificate_failure(robust: Mapping[str, Any], data_directory: Path) -> None:
    selected = (
        ("Exact finite-h", "exact_finite_h"),
        ("Random sampling", "random_completion_sampling"),
        ("PSD contraction", "psd_contraction_relaxation"),
        ("Endpoints only", "endpoint_only_no_go"),
    )
    rows = []
    for position, (label, key) in enumerate(selected):
        method = robust["methods"][key]
        rows.append(
            {
                "x": position,
                "method": label,
                "training": method["training_completion_evaluation"][
                    "maximum_operator_2_norm_error"
                ],
                "held_out": method["held_out_completion_evaluation"][
                    "maximum_operator_2_norm_error"
                ],
                "exact": method["exact_finite_h_evaluation"]["certified_radius"],
            }
        )
    _write_csv(
        data_directory / "certificate_failure.csv",
        ("x", "method", "training", "held_out", "exact"),
        rows,
    )


def _resource_byte_hop(payload: Mapping[str, Any]) -> float:
    value = payload.get("total_byte_hop", payload.get("byte_hop", 0.0))
    return float(value)


def _export_local_tradeoff(stress: Mapping[str, Any], data_directory: Path) -> None:
    cases = [
        case
        for case in stress["cases"]
        if case["classification"] == "fixed_graph_synthetic_stress"
    ]
    if not cases:
        raise ValueError("no fixed_graph_synthetic_stress cases found")
    round_values = [row["rounds"] for row in cases[0]["round_sweep"]]
    methods = {
        "trunc": "induced_local_truncation",
        "neumann": "neumann",
        "chebyshev": "chebyshev",
    }
    rows = []
    for index, rounds in enumerate(round_values):
        row: dict[str, Any] = {"rounds": rounds}
        for short_name, method_name in methods.items():
            method_rows = [
                case["round_sweep"][index]["methods"][method_name] for case in cases
            ]
            row[f"{short_name}_mse"] = statistics.fmean(
                method["signal_metrics"]["empirical_mse_per_node"]
                for method in method_rows
            )
            row[f"{short_name}_byte_hop"] = statistics.fmean(
                _resource_byte_hop(method["resources"]) for method in method_rows
            )
        rows.append(row)
    fields = (
        "rounds",
        "trunc_mse",
        "neumann_mse",
        "chebyshev_mse",
        "trunc_byte_hop",
        "neumann_byte_hop",
        "chebyshev_byte_hop",
    )
    _write_csv(data_directory / "local_method_tradeoff.csv", fields, rows)


def _export_scenario_scaling(robust: Mapping[str, Any], data_directory: Path) -> None:
    rows = (
        {
            "q": row["q"],
            "finite_candidates": row["finite_candidates"],
            "finite_vertices": row["finite_vertices"],
            "unbounded_partial_partitions": row["unbounded_partial_partitions"],
        }
        for row in robust["scenario_scaling"]
    )
    _write_csv(
        data_directory / "scenario_scaling.csv",
        (
            "q",
            "finite_candidates",
            "finite_vertices",
            "unbounded_partial_partitions",
        ),
        rows,
    )


def _export_gpu_crossover(
    gpu_payloads: Sequence[Mapping[str, Any]], data_directory: Path
) -> None:
    if not gpu_payloads:
        raise ValueError("at least one GPU result is required")
    node_counts = [row["matrix"]["shape"][0] for row in gpu_payloads[0]["results"]]
    kernel_samples: list[list[float]] = [[] for _ in node_counts]
    endpoint_samples: list[list[float]] = [[] for _ in node_counts]
    for payload in gpu_payloads:
        observed = [row["matrix"]["shape"][0] for row in payload["results"]]
        if observed != node_counts:
            raise ValueError("GPU repetitions do not use the same matrix sizes")
        for index, row in enumerate(payload["results"]):
            iteration = row["fixed_local_richardson"]
            kernel_samples[index].append(
                float(iteration["kernel_speedup_cpu_over_gpu"])
            )
            endpoint_samples[index].append(
                float(iteration["cpu_kernel_seconds"])
                / float(iteration["gpu_end_to_end_excluding_matrix_preprocess"])
            )
    rows = []
    for index, nodes in enumerate(node_counts):
        rows.append(
            {
                "nodes": nodes,
                "kernel_median": statistics.median(kernel_samples[index]),
                "kernel_min": min(kernel_samples[index]),
                "kernel_max": max(kernel_samples[index]),
                "end_to_end_median": statistics.median(endpoint_samples[index]),
            }
        )
    _write_csv(
        data_directory / "cpu_gpu_crossover.csv",
        ("nodes", "kernel_median", "kernel_min", "kernel_max", "end_to_end_median"),
        rows,
    )


def _safe_label(case: Mapping[str, Any], index: int) -> str:
    source = case["source"]
    raw = source.get("name") or source.get("dataset") or source.get("path") or f"case_{index}"
    return Path(str(raw)).stem.replace(" ", "_")


def _export_external_stress(stress: Mapping[str, Any], data_directory: Path) -> None:
    groups = {
        "PGLib topology": "external_pglib.csv",
        "SuiteSparse sparsity graph": "external_suitesparse.csv",
    }
    for source_kind, filename in groups.items():
        rows = []
        for index, case in enumerate(stress["cases"]):
            if case["source"]["source_kind"] != source_kind:
                continue
            final_row = case["round_sweep"][-1]
            rows.append(
                {
                    "nodes": case["graph"]["nodes"],
                    "relative_error": final_row["methods"]["chebyshev"][
                        "signal_metrics"
                    ]["empirical_relative_l2_error"],
                    "solve_seconds": case["global_exact"]["online"]["solve_seconds"],
                    "label": _safe_label(case, index),
                }
            )
        _write_csv(
            data_directory / filename,
            ("nodes", "relative_error", "solve_seconds", "label"),
            rows,
        )


def _write_ablation_table(robust: Mapping[str, Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    methods = robust["methods"]
    sampling_gap = methods["random_completion_sampling"][
        "finite_certificate_gap_over_training"
    ]
    endpoint = methods["endpoint_only_no_go"]
    endpoint_gap = (
        endpoint["exact_finite_h_evaluation"]["certified_radius"]
        - endpoint["training_radius_recomputed"]
    )
    coarse = methods["psd_contraction_relaxation"]
    coarse_gap = (
        coarse["certified_radius"]
        - coarse["exact_finite_h_evaluation"]["certified_radius"]
    )
    pruning = robust["vertex_pruning_ablation"]
    text = rf"""% Generated by experiments.smoother.paper_data; do not edit by hand.
\begin{{tabular}}{{p{{0.28\linewidth}}p{{0.22\linewidth}}p{{0.40\linewidth}}}}
\toprule
Ablation & Observed effect & Safe interpretation \\
\midrule
Random completion sampling & certificate gap {sampling_gap:.3f} & Training maximum is not a worst-case certificate. \\
Endpoint-only scenarios & missed radius {endpoint_gap:.3f} & Connected partition directions cannot be deleted. \\
$0\preceq X\preceq I$ relaxation & conservatism {coarse_gap:.3f} & Safe outer certificate, but it loses graph-terminal geometry. \\
Exact vertex pruning & {pruning['pruned_lmis']} vs. {pruning['unpruned_lmis']} LMIs & Same radius to {pruning['absolute_radius_difference']:.1e}; fewer offline constraints. \\
Raw WLS/Kalman relabeling & not performed & Reported MSE is graph-signal smoother error only. \\
\bottomrule
\end{{tabular}}
"""
    destination.write_text(text, encoding="utf-8")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_manifest(
    sources: Sequence[Path], data_directory: Path, destination: Path
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "generator": "experiments.smoother.paper_data",
        "semantics": {
            "robust": "theorem-matched hidden-completion experiment",
            "synthetic_stress": "fixed-graph graph-smoother stress test",
            "external_stress": "topology/sparsity transformation stress test; not raw WLS",
            "gpu": "fixed 40-round I+L Richardson kernel benchmark",
        },
        "sources": [
            {"path": path.as_posix(), "sha256": _sha256(path)} for path in sources
        ],
        "outputs": sorted(path.name for path in data_directory.glob("*.csv")),
    }
    destination.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p3-p4", type=Path, required=True)
    parser.add_argument("--external", type=Path, required=True)
    parser.add_argument("--gpu", type=Path, nargs="+", required=True)
    parser.add_argument(
        "--data-dir", type=Path, default=Path("paper/figures/tikz/data")
    )
    parser.add_argument(
        "--table", type=Path, default=Path("paper/tables/p4_ablation.tex")
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("paper/figures/tikz/data/manifest.json"),
    )
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    p3_p4 = _load(arguments.p3_p4)
    external = _load(arguments.external)
    gpu = [_load(path) for path in arguments.gpu]
    robust = p3_p4["robust"]
    _export_budget(robust, arguments.data_dir)
    _export_certificate_failure(robust, arguments.data_dir)
    _export_local_tradeoff(p3_p4["stress"], arguments.data_dir)
    _export_scenario_scaling(robust, arguments.data_dir)
    _export_gpu_crossover(gpu, arguments.data_dir)
    _export_external_stress(external["stress"], arguments.data_dir)
    _write_ablation_table(robust, arguments.table)
    sources = [arguments.p3_p4, arguments.external, *arguments.gpu]
    _write_manifest(sources, arguments.data_dir, arguments.manifest)
    print(arguments.data_dir)


if __name__ == "__main__":
    main()
