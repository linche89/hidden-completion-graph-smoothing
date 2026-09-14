"""Affine descriptions of implementable common local rules Q."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray


BoolArray = NDArray[np.bool_]
SUPPORTED_FIELDS = frozenset(
    {
        "allowed_mask",
        "zeros",
        "r_hop",
        "shared_coefficients",
        "fixed_entries",
        "linear_equalities",
        "row_sums",
    }
)


def _entry(entry: object, rows: int, columns: int) -> tuple[int, int]:
    if not isinstance(entry, (list, tuple)) or len(entry) != 2:
        raise ValueError(f"matrix entry must be [row, column], got {entry!r}")
    row, column = entry
    if (
        not isinstance(row, int)
        or isinstance(row, bool)
        or not isinstance(column, int)
        or isinstance(column, bool)
    ):
        raise ValueError(f"matrix entry indices must be integers, got {entry!r}")
    if not (0 <= row < rows and 0 <= column < columns):
        raise ValueError(f"matrix entry {entry!r} is outside shape {(rows, columns)}")
    return row, column


def r_hop_allowed_mask(
    adjacency: NDArray[np.float64] | list[list[float]],
    output_nodes: list[int] | tuple[int, ...],
    radius: int,
) -> BoolArray:
    """Return the zero-pattern mask realizable in ``radius`` synchronous hops."""

    graph = np.asarray(adjacency)
    if graph.ndim != 2 or graph.shape[0] != graph.shape[1]:
        raise ValueError("adjacency must be square")
    if not np.issubdtype(graph.dtype, np.number) or not np.all(np.isfinite(graph)):
        raise ValueError("adjacency must contain only finite numbers")
    if not isinstance(radius, int) or isinstance(radius, bool) or radius < 0:
        raise ValueError("radius must be a nonnegative integer")
    q = graph.shape[0]
    undirected = np.logical_or(graph != 0, graph.T != 0)
    np.fill_diagonal(undirected, False)
    mask = np.zeros((len(output_nodes), q), dtype=bool)
    for row, source in enumerate(output_nodes):
        if not isinstance(source, int) or isinstance(source, bool) or not 0 <= source < q:
            raise ValueError(f"invalid output node {source!r}")
        distance = [-1] * q
        distance[source] = 0
        queue: deque[int] = deque([source])
        while queue:
            vertex = queue.popleft()
            if distance[vertex] == radius:
                continue
            for neighbor in np.flatnonzero(undirected[vertex]):
                neighbor_int = int(neighbor)
                if distance[neighbor_int] < 0:
                    distance[neighbor_int] = distance[vertex] + 1
                    queue.append(neighbor_int)
        mask[row] = np.asarray([0 <= value <= radius for value in distance])
    return mask


@dataclass(slots=True)
class LocalRuleSpec:
    """Convex affine constraints defining ``Q_loc``.

    Supported JSON fields are ``allowed_mask``, ``zeros``, ``r_hop``,
    ``shared_coefficients``, ``fixed_entries``, ``linear_equalities``, and
    ``row_sums``.  Multiple fields are intersected.
    """

    data: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any] | None) -> "LocalRuleSpec":
        if data is None:
            return cls()
        if not isinstance(data, dict):
            raise ValueError("local_rule must be a JSON object")
        unknown = set(data) - SUPPORTED_FIELDS
        if unknown:
            raise ValueError(f"unsupported local_rule fields: {sorted(unknown)}")
        return cls(dict(data))

    def _check_supported_fields(self) -> None:
        unknown = set(self.data) - SUPPORTED_FIELDS
        if unknown:
            raise ValueError(f"unsupported local_rule fields: {sorted(unknown)}")

    def allowed_mask(self, rows: int, columns: int) -> BoolArray:
        self._check_supported_fields()
        mask = np.ones((rows, columns), dtype=bool)
        if "allowed_mask" in self.data:
            supplied = np.asarray(self.data["allowed_mask"], dtype=bool)
            if supplied.shape != (rows, columns):
                raise ValueError(
                    f"allowed_mask has shape {supplied.shape}, expected {(rows, columns)}"
                )
            mask &= supplied

        if "r_hop" in self.data:
            specification = self.data["r_hop"]
            if not isinstance(specification, dict):
                raise ValueError("r_hop must be an object")
            missing = {"adjacency", "output_nodes", "radius"} - set(specification)
            unknown = set(specification) - {"adjacency", "output_nodes", "radius"}
            if missing or unknown:
                raise ValueError(
                    f"r_hop fields invalid; missing={sorted(missing)}, unknown={sorted(unknown)}"
                )
            hop_mask = r_hop_allowed_mask(
                specification["adjacency"],
                specification["output_nodes"],
                specification["radius"],
            )
            if hop_mask.shape != (rows, columns):
                raise ValueError(
                    "r_hop output_nodes length and adjacency size must match the shape of Q"
                )
            mask &= hop_mask

        for entry in self.data.get("zeros", []):
            mask[_entry(entry, rows, columns)] = False
        return mask

    def cvxpy_constraints(self, variable: Any) -> tuple[list[Any], BoolArray]:
        """Build CVXPY constraints for a matrix variable."""

        rows, columns = variable.shape
        mask = self.allowed_mask(rows, columns)
        constraints: list[Any] = []
        for row, column in zip(*np.nonzero(~mask), strict=True):
            constraints.append(variable[int(row), int(column)] == 0)

        for group in self.data.get("shared_coefficients", []):
            if not isinstance(group, list) or len(group) < 2:
                raise ValueError("each shared_coefficients group must contain at least two entries")
            entries = [_entry(item, rows, columns) for item in group]
            anchor = entries[0]
            constraints.extend(variable[item] == variable[anchor] for item in entries[1:])

        for item in self.data.get("fixed_entries", []):
            if isinstance(item, dict):
                row, column = _entry([item.get("row"), item.get("column")], rows, columns)
                value = item.get("value")
            elif isinstance(item, (list, tuple)) and len(item) == 3:
                row, column = _entry(item[:2], rows, columns)
                value = item[2]
            else:
                raise ValueError(f"invalid fixed entry {item!r}")
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not np.isfinite(value):
                raise ValueError(f"fixed entry value must be finite, got {value!r}")
            constraints.append(variable[row, column] == float(value))

        for equality in self.data.get("linear_equalities", []):
            if not isinstance(equality, dict) or "terms" not in equality or "rhs" not in equality:
                raise ValueError("linear_equalities entries need terms and rhs")
            expression: Any = 0.0
            for term in equality["terms"]:
                if not isinstance(term, (list, tuple)) or len(term) != 3:
                    raise ValueError("a linear equality term must be [row, column, coefficient]")
                row, column = _entry(term[:2], rows, columns)
                coefficient = float(term[2])
                if not np.isfinite(coefficient):
                    raise ValueError("linear equality coefficients must be finite")
                expression += coefficient * variable[row, column]
            right_hand_side = float(equality["rhs"])
            if not np.isfinite(right_hand_side):
                raise ValueError("linear equality right-hand sides must be finite")
            constraints.append(expression == right_hand_side)

        if "row_sums" in self.data:
            values = self.data["row_sums"]
            if isinstance(values, (int, float)):
                values = [float(values)] * rows
            if not isinstance(values, list) or len(values) != rows:
                raise ValueError("row_sums must be a scalar or one value per output row")
            numeric_values = [float(value) for value in values]
            if not np.all(np.isfinite(numeric_values)):
                raise ValueError("row_sums values must be finite")
            constraints.extend(
                variable[row, :] @ np.ones(columns) == numeric_values[row]
                for row in range(rows)
            )

        return constraints, mask

    def summary(self, rows: int, columns: int) -> dict[str, object]:
        mask = self.allowed_mask(rows, columns)
        return {
            "free_entries_before_equalities": int(mask.sum()),
            "forced_zero_entries": int(mask.size - mask.sum()),
            "shared_groups": len(self.data.get("shared_coefficients", [])),
            "fixed_entries": len(self.data.get("fixed_entries", [])),
            "linear_equalities": len(self.data.get("linear_equalities", [])),
            "has_row_sum_constraints": "row_sums" in self.data,
        }
