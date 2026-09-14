"""Finite and unbounded hidden-completion scenario enumeration."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Iterator, Literal

import numpy as np
from numpy.typing import NDArray


Array = NDArray[np.float64]
Block = tuple[int, ...]
Partition = tuple[Block, ...]


def _validate_q(q: int) -> None:
    if not isinstance(q, int) or isinstance(q, bool) or q < 1:
        raise ValueError(f"q must be a positive integer, got {q!r}")


def set_partitions(items: tuple[int, ...]) -> Iterator[Partition]:
    """Yield every set partition once in canonical block order.

    Blocks and the blocks within a partition are ordered by their least
    element.  The routine is intended for the boundary-size fixed-parameter
    regime; its Bell-number output size is unavoidable.
    """

    if not items:
        yield tuple()
        return
    if tuple(sorted(set(items))) != items:
        raise ValueError("items must be strictly increasing and duplicate-free")

    first, rest = items[0], items[1:]
    for tail in set_partitions(rest):
        yield ((first,),) + tail
        for index in range(len(tail)):
            blocks = list(tail)
            blocks[index] = (first,) + blocks[index]
            yield tuple(sorted(blocks, key=lambda block: block[0]))


def _weak_compositions(total: int, parts: int) -> Iterator[tuple[int, ...]]:
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _weak_compositions(total - first, parts - 1):
            yield (first,) + rest


def hidden_allocations(h: int, block_count: int) -> Iterator[tuple[int, ...]]:
    """Yield all nonnegative allocations with total at most ``h``."""

    if not isinstance(h, int) or isinstance(h, bool) or h < 0:
        raise ValueError(f"h must be a nonnegative integer, got {h!r}")
    if block_count < 1:
        raise ValueError("block_count must be positive")
    for total in range(h + 1):
        yield from _weak_compositions(total, block_count)


def finite_scenario_matrix(q: int, partition: Partition, allocation: tuple[int, ...]) -> Array:
    """Construct X_(pi,k) for a full partition and hidden allocation."""

    _validate_q(q)
    if len(partition) != len(allocation):
        raise ValueError("partition and allocation must have the same length")
    flattened = sorted(vertex for block in partition for vertex in block)
    if flattened != list(range(q)):
        raise ValueError("partition must be a full partition of range(q)")
    if any(
        not isinstance(value, int) or isinstance(value, bool) or value < 0
        for value in allocation
    ):
        raise ValueError("allocation entries must be nonnegative integers")

    matrix = np.zeros((q, q), dtype=float)
    for block, hidden_count in zip(partition, allocation, strict=True):
        indices = np.asarray(block, dtype=int)
        matrix[np.ix_(indices, indices)] = 1.0 / (len(block) + hidden_count)
    return matrix


def partial_projection_matrix(q: int, partition: Partition) -> Array:
    """Construct a partial-partition orthogonal projection."""

    _validate_q(q)
    if any(not block for block in partition):
        raise ValueError("partial-partition blocks must be nonempty")
    flattened = [vertex for block in partition for vertex in block]
    if len(flattened) != len(set(flattened)) or any(not 0 <= vertex < q for vertex in flattened):
        raise ValueError("partial partition contains duplicate or invalid port indices")
    matrix = np.zeros((q, q), dtype=float)
    for block in partition:
        indices = np.asarray(block, dtype=int)
        matrix[np.ix_(indices, indices)] = 1.0 / len(block)
    return matrix


@dataclass(frozen=True, slots=True)
class Scenario:
    """One exact terminal scenario.

    ``allocation`` is present only for finite-budget scenarios.  In the
    unbounded model, ports outside the union of ``partition`` are inactive.
    """

    kind: Literal["finite", "partial"]
    q: int
    partition: Partition
    matrix: Array
    allocation: tuple[int, ...] | None = None
    is_vertex: bool = True

    @property
    def active_ports(self) -> tuple[int, ...]:
        return tuple(sorted(vertex for block in self.partition for vertex in block))

    @property
    def inactive_ports(self) -> tuple[int, ...]:
        active = set(self.active_ports)
        return tuple(vertex for vertex in range(self.q) if vertex not in active)

    @property
    def label(self) -> str:
        blocks = "|".join("{" + ",".join(map(str, block)) + "}" for block in self.partition)
        if self.kind == "finite":
            allocation = ",".join(map(str, self.allocation or ()))
            return f"finite:{blocks};k=({allocation})"
        inactive = ",".join(map(str, self.inactive_ports))
        return f"partial:{blocks or 'empty'};inactive=({inactive})"

    def to_dict(self, *, include_matrix: bool = True) -> dict[str, object]:
        record: dict[str, object] = {
            "kind": self.kind,
            "label": self.label,
            "partition": [list(block) for block in self.partition],
            "active_ports": list(self.active_ports),
            "inactive_ports": list(self.inactive_ports),
            "is_vertex": self.is_vertex,
        }
        if self.allocation is not None:
            record["allocation"] = list(self.allocation)
            record["hidden_used"] = sum(self.allocation)
        if include_matrix:
            record["matrix"] = self.matrix.tolist()
        return record


def finite_scenarios(q: int, h: int, *, vertices_only: bool = False) -> list[Scenario]:
    """Enumerate the exact finite-``h`` candidate set or its hull vertices."""

    _validate_q(q)
    if not isinstance(h, int) or isinstance(h, bool) or h < 0:
        raise ValueError(f"h must be a nonnegative integer, got {h!r}")

    answer: list[Scenario] = []
    for partition in set_partitions(tuple(range(q))):
        for allocation in hidden_allocations(h, len(partition)):
            is_vertex = h == 0 or not any(allocation) or sum(allocation) == h
            if vertices_only and not is_vertex:
                continue
            answer.append(
                Scenario(
                    kind="finite",
                    q=q,
                    partition=partition,
                    allocation=allocation,
                    matrix=finite_scenario_matrix(q, partition, allocation),
                    is_vertex=is_vertex,
                )
            )
    return answer


def partial_partition_scenarios(q: int) -> list[Scenario]:
    """Enumerate all B_(q+1) unbounded partial-partition projections."""

    _validate_q(q)
    ports = tuple(range(q))
    answer: list[Scenario] = []
    for active_count in range(q + 1):
        for active in combinations(ports, active_count):
            for partition in set_partitions(active):
                answer.append(
                    Scenario(
                        kind="partial",
                        q=q,
                        partition=partition,
                        matrix=partial_projection_matrix(q, partition),
                    )
                )
    return answer


def stirling_second_kind(n: int, k: int) -> int:
    if n == k == 0:
        return 1
    if n <= 0 or k <= 0 or k > n:
        return 0
    table = [[0] * (k + 1) for _ in range(n + 1)]
    table[0][0] = 1
    for row in range(1, n + 1):
        for column in range(1, min(row, k) + 1):
            table[row][column] = table[row - 1][column - 1] + column * table[row - 1][column]
    return table[n][k]


def bell_number(n: int) -> int:
    if n < 0:
        raise ValueError("n must be nonnegative")
    return sum(stirling_second_kind(n, k) for k in range(n + 1))


def scenario_count(q: int, h: int | None = None) -> int:
    """Return N(q,h), or B_(q+1) when ``h`` is ``None``."""

    _validate_q(q)
    if h is None:
        return bell_number(q + 1)
    if not isinstance(h, int) or isinstance(h, bool) or h < 0:
        raise ValueError("h must be nonnegative or None")
    return sum(stirling_second_kind(q, blocks) * comb(h + blocks, blocks) for blocks in range(1, q + 1))


def vertex_count(q: int, h: int) -> int:
    """Return the exact number of vertices of the finite terminal hull."""

    _validate_q(q)
    if not isinstance(h, int) or isinstance(h, bool) or h < 0:
        raise ValueError("h must be nonnegative")
    if h == 0:
        return bell_number(q)
    return sum(
        stirling_second_kind(q, blocks) * (1 + comb(h + blocks - 1, blocks - 1))
        for blocks in range(1, q + 1)
    )
