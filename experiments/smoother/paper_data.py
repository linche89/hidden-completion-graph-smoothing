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


def _latex_scientific(value: float, digits: int = 2) -> str:
    """Format a finite scalar as LaTeX scientific notation."""

    mantissa, exponent = f"{value:.{digits}e}".split("e")
    return rf"{mantissa}\times 10^{{{int(exponent)}}}"


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


def _linear_quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("cannot take a quantile of an empty sequence")
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _export_core_round_tradeoff(core: Mapping[str, Any], data_directory: Path) -> None:
    source_rows = core["multi_instance"]["round_tradeoffs"]
    categories = (
        (0, "zero_hop_radius", "0-hop"),
        (1, "one_hop_radius", "1-hop"),
        (2, "full_port_radius", "Full port"),
    )
    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for x, field, label in categories:
        values = []
        for row in source_rows:
            normalized = float(row[field]) / float(row["zero_hop_radius"])
            values.append(normalized)
            raw_rows.append(
                {
                    "q": row["q"],
                    "h": row["h"],
                    "x": x,
                    "rule": label,
                    "normalized_radius": normalized,
                }
            )
        summary_rows.append(
            {
                "x": x,
                "rule": label,
                "minimum": min(values),
                "q25": _linear_quantile(values, 0.25),
                "median": statistics.median(values),
                "q75": _linear_quantile(values, 0.75),
                "maximum": max(values),
            }
        )
    _write_csv(
        data_directory / "core_round_tradeoff_raw.csv",
        ("q", "h", "x", "rule", "normalized_radius"),
        raw_rows,
    )
    _write_csv(
        data_directory / "core_round_tradeoff_summary.csv",
        ("x", "rule", "minimum", "q25", "median", "q75", "maximum"),
        summary_rows,
    )


def _write_ecdf(destination: Path, values: Sequence[float]) -> None:
    ordered = sorted(max(0.0, 100.0 * float(value)) for value in values)
    _write_csv(
        destination,
        ("overhead_percent", "fraction"),
        (
            {
                "overhead_percent": value,
                "fraction": (index + 1) / len(ordered),
            }
            for index, value in enumerate(ordered)
        ),
    )


def _export_core_baseline_cdfs(core: Mapping[str, Any], data_directory: Path) -> None:
    rows = core["multi_instance"]["rows"]
    sampling = [
        replicate["relative_true_risk_over_exact_design"]
        for row in rows
        for replicate in row["random_sampling"]["replicates"]
    ]
    contraction = [
        row["psd_contraction"]["relative_true_risk_over_exact_design"]
        for row in rows
    ]
    _write_ecdf(data_directory / "core_sampling_overhead_cdf.csv", sampling)
    _write_ecdf(data_directory / "core_psd_overhead_cdf.csv", contraction)


def _export_sampling_budget_curve(core: Mapping[str, Any], data_directory: Path) -> None:
    rows = []
    rows_by_case: dict[str, list[dict[str, Any]]] = {}
    for case_index, case in enumerate(core["sampling_curve"]["cases"]):
        case_rows: list[dict[str, Any]] = []
        for summary in case["summary_by_sample_size"]:
            under = summary["relative_underreport"]
            regret = summary["relative_true_risk_over_exact_design"]
            exported = {
                    "case_index": case_index,
                    "case_id": case["case_id"],
                    "sample_size": summary["sample_size"],
                    "under_q25_percent": 100.0 * under["q25"],
                    "under_median_percent": 100.0 * under["median"],
                    "under_q75_percent": 100.0 * under["q75"],
                    "regret_q25_percent": 100.0 * regret["q25"],
                    "regret_median_percent": 100.0 * regret["median"],
                    "regret_q75_percent": 100.0 * regret["q75"],
                    "failure_rate_percent": 100.0 * summary["failure_rate"],
                    "wilson_low_percent": 100.0
                    * summary["failure_rate_wilson95_low"],
                    "wilson_high_percent": 100.0
                    * summary["failure_rate_wilson95_high"],
                }
            rows.append(exported)
            case_rows.append(exported)
        rows_by_case[str(case["case_id"])] = case_rows
    _write_csv(
        data_directory / "core_sampling_budget.csv",
        (
            "case_index",
            "case_id",
            "sample_size",
            "under_q25_percent",
            "under_median_percent",
            "under_q75_percent",
            "regret_q25_percent",
            "regret_median_percent",
            "regret_q75_percent",
            "failure_rate_percent",
            "wilson_low_percent",
            "wilson_high_percent",
        ),
        rows,
    )
    fields = (
        "case_index",
        "case_id",
        "sample_size",
        "under_q25_percent",
        "under_median_percent",
        "under_q75_percent",
        "regret_q25_percent",
        "regret_median_percent",
        "regret_q75_percent",
        "failure_rate_percent",
        "wilson_low_percent",
        "wilson_high_percent",
    )
    for case_id, case_rows in rows_by_case.items():
        _write_csv(
            data_directory / f"core_sampling_budget_{case_id}.csv",
            fields,
            case_rows,
        )


def _export_path_tightness(core: Mapping[str, Any], data_directory: Path) -> None:
    rows = []
    rows_by_witness: dict[int, list[dict[str, Any]]] = {}
    for witness in core["path_tightness"]["witnesses"]:
        index = witness["witness_index"]
        block_count = len(witness["scenario"]["partition"])
        hidden_used = witness["scenario"]["hidden_used"]
        short_label = f"{block_count} blocks / {hidden_used} hidden"
        for row in witness["rows"]:
            exported = {
                    "witness": index,
                    "label": short_label,
                    "strong_conductance": row["strong_conductance"],
                    "physical_over_canonical": row["physical_over_canonical"],
                    "relative_radius_gap": row["relative_radius_gap"],
                    "port_block_operator_distance": row[
                        "port_block_operator_distance"
                    ],
                }
            rows.append(exported)
            rows_by_witness.setdefault(int(index), []).append(exported)
    _write_csv(
        data_directory / "core_path_tightness.csv",
        (
            "witness",
            "label",
            "strong_conductance",
            "physical_over_canonical",
            "relative_radius_gap",
            "port_block_operator_distance",
        ),
        rows,
    )
    fields = (
        "witness",
        "label",
        "strong_conductance",
        "physical_over_canonical",
        "relative_radius_gap",
        "port_block_operator_distance",
    )
    for index, witness_rows in rows_by_witness.items():
        _write_csv(
            data_directory / f"core_path_witness_{index}.csv",
            fields,
            witness_rows,
        )
    grouped: dict[float, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(float(row["strong_conductance"]), []).append(row)
    envelope = [
        {
            "strong_conductance": strong,
            "maximum_relative_radius_gap": max(
                float(row["relative_radius_gap"]) for row in selected
            ),
            "maximum_port_block_distance": max(
                float(row["port_block_operator_distance"]) for row in selected
            ),
        }
        for strong, selected in sorted(grouped.items())
    ]
    _write_csv(
        data_directory / "core_path_envelope.csv",
        (
            "strong_conductance",
            "maximum_relative_radius_gap",
            "maximum_port_block_distance",
        ),
        envelope,
    )


def _write_core_number_macros(
    core: Mapping[str, Any],
    crosscheck: Mapping[str, Any],
    destination: Path,
) -> None:
    rows = core["multi_instance"]["rows"]
    principal = next(
        row
        for row in rows
        if row["q"] == 4
        and row["h"] == 1
        and row["family"] == "constant_preserving_state"
        and row["rule_name"] == "one_hop"
    )
    sampling_case = next(
        case
        for case in core["sampling_curve"]["cases"]
        if case["case_id"] == "constant_preserving_state"
    )
    sampling_32 = next(
        row for row in sampling_case["summary_by_sample_size"] if row["sample_size"] == 32
    )
    one_hop_rows = [
        row
        for row in rows
        if row["family"] == "constant_preserving_state"
        and row["rule_name"] == "one_hop"
    ]
    psd_one_hop = sorted(
        float(row["psd_contraction"]["relative_true_risk_over_exact_design"])
        for row in one_hop_rows
    )
    final_path_gaps = [
        float(witness["rows"][-1]["relative_radius_gap"])
        for witness in core["path_tightness"]["witnesses"]
    ]
    final_block_gaps = [
        float(witness["rows"][-1]["port_block_operator_distance"])
        for witness in core["path_tightness"]["witnesses"]
    ]
    total_sampling_runs = sum(
        len(case["records"]) for case in core["sampling_curve"]["cases"]
    )
    total_sampling_failures = sum(
        bool(row["underreport_failure"])
        for case in core["sampling_curve"]["cases"]
        for row in case["records"]
    )
    resources = principal["online_resources"]
    crosscheck_summary = crosscheck["summary"]
    primary_statuses = [str(row["exact_finite"]["status"]) for row in rows]
    primary_epigraph_gaps = [
        abs(
            float(row["exact_finite"]["optimization_squared_radius"])
            - float(row["exact_finite"]["raw_verified_radius"]) ** 2
        )
        for row in rows
    ]
    exact_total_seconds = [
        float(row["exact_finite"]["timing_seconds"]["total"])
        for row in rows
    ]
    slowest_exact_row = max(
        rows,
        key=lambda row: float(row["exact_finite"]["timing_seconds"]["total"]),
    )
    values = {
        "CoreProblemCount": core["multi_instance"]["grid"]["problem_count"],
        "CoreSamplingDesignCount": sum(
            len(row["random_sampling"]["replicates"]) for row in rows
        ),
        "PrincipalExactRadius": f"{principal['exact_finite']['radius']:.6f}",
        "PrincipalUnboundedRadius": f"{principal['exact_unbounded']['radius']:.6f}",
        "PrincipalFiniteGainPct": f"{100.0 * principal['finite_budget_radius_reduction_from_unbounded']:.1f}",
        "PrincipalPSDTrueRadius": f"{principal['psd_contraction']['exact_finite_radius']:.6f}",
        "PrincipalPSDOverheadPct": f"{100.0 * principal['psd_contraction']['relative_true_risk_over_exact_design']:.1f}",
        "PrincipalSampleUnderMedianPct": f"{100.0 * sampling_32['relative_underreport']['median']:.1f}",
        "PrincipalSampleRegretMedianPct": f"{100.0 * sampling_32['relative_true_risk_over_exact_design']['median']:.1f}",
        "OneHopMedianGainPct": f"{100.0 * core['multi_instance']['summary']['one_hop_radius_reduction_from_zero']['median']:.1f}",
        "FullPortMedianGainPct": f"{100.0 * core['multi_instance']['summary']['full_port_radius_reduction_from_zero']['median']:.1f}",
        "OneHopPSDMedianOverheadPct": f"{100.0 * statistics.median(psd_one_hop):.1f}",
        "SamplingCurveRunCount": total_sampling_runs,
        "SamplingCurveFailureCount": total_sampling_failures,
        "SamplingWilsonLowPct": f"{100.0 * sampling_32['failure_rate_wilson95_low']:.1f}",
        "PathWitnessCount": len(core["path_tightness"]["witnesses"]),
        "PathMaxFinalRiskGap": _latex_scientific(max(final_path_gaps)),
        "PathMaxFinalBlockGap": _latex_scientific(max(final_block_gaps)),
        "PrincipalRounds": resources["rounds"],
        "PrincipalScalarHops": resources["scalar_hops"],
        "PrincipalOnlineFlops": resources["online_flops"],
        "PrincipalStorageBytes": resources["sparse_coefficient_storage_bytes"],
        "CrosscheckCount": crosscheck_summary["selected_cells"],
        "CrosscheckMaxRelativeRadiusDelta": _latex_scientific(
            crosscheck_summary["maximum_absolute_relative_raw_radius_delta"]
        ),
        "CrosscheckMaxDeploymentRadiusDelta": _latex_scientific(
            crosscheck_summary["maximum_absolute_relative_deployed_radius_delta"]
        ),
        "CrosscheckMaxEpigraphGap": _latex_scientific(
            crosscheck_summary["maximum_absolute_epigraph_squared_risk_gap"]
        ),
        "CrosscheckMaxLocalResidual": _latex_scientific(
            crosscheck_summary["maximum_local_rule_residual"]
        ),
        "CorePrimaryOptimalCount": primary_statuses.count("optimal"),
        "CorePrimaryInaccurateCount": primary_statuses.count(
            "optimal_inaccurate"
        ),
        "CoreExactMedianTotalSeconds": f"{statistics.median(exact_total_seconds):.3f}",
        "CoreExactMaxTotalSeconds": f"{max(exact_total_seconds):.2f}",
        "CoreExactMaxQ": slowest_exact_row["q"],
        "CoreExactMaxH": slowest_exact_row["h"],
        "CoreMaxEpigraphGap": _latex_scientific(max(primary_epigraph_gaps)),
        "CoreMaxDeploymentRadiusChange": _latex_scientific(
            max(
                abs(float(row["exact_finite"]["deployment_radius_change"]))
                for row in rows
            )
        ),
        "CoreMinDeployedCertificateEigenvalue": _latex_scientific(
            min(
                float(row["exact_finite"]["certificate_min_eigenvalue"])
                for row in rows
            )
        ),
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines = ["% Generated by experiments.smoother.paper_data; do not edit by hand."]
    lines.extend(
        rf"\newcommand{{\{name}}}{{{value}}}" for name, value in values.items()
    )
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
            "core_evidence": "fixed theorem-matched breadth, sampling, and path-tightness audit",
            "solver_crosscheck": "independent numerical SCS check of flagged Clarabel cells",
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
    parser.add_argument("--core-evidence", type=Path, required=True)
    parser.add_argument("--solver-crosscheck", type=Path, required=True)
    parser.add_argument(
        "--data-dir", type=Path, default=Path("paper/figures/tikz/data")
    )
    parser.add_argument(
        "--core-table",
        type=Path,
        default=Path("paper/tsp/tables/core_numbers.tex"),
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
    core = _load(arguments.core_evidence)
    crosscheck = _load(arguments.solver_crosscheck)
    robust = p3_p4["robust"]
    _export_budget(robust, arguments.data_dir)
    _export_local_tradeoff(p3_p4["stress"], arguments.data_dir)
    _export_scenario_scaling(robust, arguments.data_dir)
    _export_gpu_crossover(gpu, arguments.data_dir)
    _export_external_stress(external["stress"], arguments.data_dir)
    _export_core_round_tradeoff(core, arguments.data_dir)
    _export_core_baseline_cdfs(core, arguments.data_dir)
    _export_sampling_budget_curve(core, arguments.data_dir)
    _export_path_tightness(core, arguments.data_dir)
    _write_core_number_macros(core, crosscheck, arguments.core_table)
    sources = [
        arguments.p3_p4,
        arguments.external,
        arguments.core_evidence,
        arguments.solver_crosscheck,
        *arguments.gpu,
    ]
    _write_manifest(sources, arguments.data_dir, arguments.manifest)
    print(arguments.data_dir)


if __name__ == "__main__":
    main()
