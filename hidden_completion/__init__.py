"""Exact scenario reduction for hidden graph completions.

The public API deliberately separates combinatorial scenario generation from
the numerical SDP solve.  This keeps the theorem-level objects inspectable and
makes it possible to regression-test the enumeration without an SDP solver.
"""

from .local_rules import LocalRuleSpec, r_hop_allowed_mask
from .metrics import canonical_error, full_input_schatten_norm, scenario_error_gram
from .scenarios import (
    Scenario,
    finite_scenarios,
    partial_partition_scenarios,
    scenario_count,
    vertex_count,
)
from .solver import ExactSolveResult, solve_exact_spectral

__all__ = [
    "ExactSolveResult",
    "LocalRuleSpec",
    "Scenario",
    "canonical_error",
    "finite_scenarios",
    "full_input_schatten_norm",
    "partial_partition_scenarios",
    "r_hop_allowed_mask",
    "scenario_count",
    "scenario_error_gram",
    "solve_exact_spectral",
    "vertex_count",
]
