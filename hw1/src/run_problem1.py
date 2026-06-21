"""Exercise 1: max-flow, cuts, and capacity augmentation.

Run from the hw1 directory:

    python src/run_problem1.py
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


NODES = ["o", "a", "b", "c", "d"]
SOURCE = "o"
SINK = "d"


@dataclass(frozen=True)
class Edge:
    """Directed edge with the label and capacity used in Exercise 1."""

    label: str
    tail: str
    head: str
    capacity: int


BASE_EDGES = [
    Edge("e1", "o", "a", 3),
    Edge("e2", "a", "d", 3),
    Edge("e3", "o", "b", 3),
    Edge("e4", "b", "c", 3),
    Edge("e5", "c", "d", 2),
    Edge("e6", "a", "b", 1),
]

ALLOCATION_CACHE: dict[tuple[int, int], np.ndarray] = {}


def ensure_output_dirs() -> None:
    """Create output folders for CSV tables and figures."""
    Path("results").mkdir(exist_ok=True)
    Path("figures").mkdir(exist_ok=True)


def powerset(items: list[str]) -> Iterable[tuple[str, ...]]:
    """Yield all subsets of the non-terminal nodes when enumerating cuts."""
    for size in range(len(items) + 1):
        yield from combinations(items, size)


def enumerate_cuts(nodes: list[str], source: str, sink: str) -> list[frozenset[str]]:
    """All source-sink cuts S with source in S and sink not in S."""
    middle_nodes = [node for node in nodes if node not in {source, sink}]
    cuts = []
    for subset in powerset(middle_nodes):
        cuts.append(frozenset((source, *subset)))
    return cuts


def crossing_edge_indices(edges: list[Edge], cut: frozenset[str]) -> list[int]:
    """Return indices of directed edges crossing from S to V minus S."""
    return [
        idx
        for idx, edge in enumerate(edges)
        if edge.tail in cut and edge.head not in cut
    ]


def cut_capacity(edges: list[Edge], cut: frozenset[str]) -> int:
    """Compute the total capacity of a source-sink cut."""
    return sum(edges[idx].capacity for idx in crossing_edge_indices(edges, cut))


def cut_table(edges: list[Edge], nodes: list[str]) -> pd.DataFrame:
    """Build the table of all source-sink cuts and their capacities."""
    rows = []
    universe = set(nodes)
    for cut in enumerate_cuts(nodes, SOURCE, SINK):
        crossing = crossing_edge_indices(edges, cut)
        rows.append(
            {
                "S": "{" + ",".join(sorted(cut, key=nodes.index)) + "}",
                "V_minus_S": "{" + ",".join(sorted(universe - cut, key=nodes.index)) + "}",
                "crossing_edges": ", ".join(edges[idx].label for idx in crossing),
                "capacity": sum(edges[idx].capacity for idx in crossing),
            }
        )
    return pd.DataFrame(rows)


def build_networkx_graph(edges: list[Edge]) -> nx.DiGraph:
    """Build the directed capacitated graph used by NetworkX max-flow."""
    graph = nx.DiGraph()
    graph.add_nodes_from(NODES)
    for edge in edges:
        graph.add_edge(edge.tail, edge.head, capacity=edge.capacity, label=edge.label)
    return graph


def max_flow_value(edges: list[Edge]) -> int:
    """Compute the maximum o-d flow value for the current edge set."""
    graph = build_networkx_graph(edges)
    value, _ = nx.maximum_flow(graph, SOURCE, SINK, capacity="capacity")
    return int(value)


def integer_compositions(total: int, parts: int) -> Iterable[tuple[int, ...]]:
    """Yield all nonnegative integer vectors of length parts summing exactly to total."""
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in integer_compositions(total - first, parts - 1):
            yield (first, *rest)


def allocations_up_to(max_extra: int, parts: int) -> np.ndarray:
    """All nonnegative integer allocations with total budget at most max_extra.

    The array is cached because part 1(c) solves the same allocation problem for
    every candidate added link.
    """
    key = (max_extra, parts)
    if key not in ALLOCATION_CACHE:
        rows = [
            allocation
            for total in range(max_extra + 1)
            for allocation in integer_compositions(total, parts)
        ]
        ALLOCATION_CACHE[key] = np.array(rows, dtype=np.int16)
    return ALLOCATION_CACHE[key]


def cut_masks(edges: list[Edge], nodes: list[str]) -> tuple[list[int], list[list[int]]]:
    """Return base cut capacities and 0/1 masks for crossing edges."""
    base_capacities = []
    masks = []
    for cut in enumerate_cuts(nodes, SOURCE, SINK):
        crossing = set(crossing_edge_indices(edges, cut))
        base_capacities.append(sum(edges[idx].capacity for idx in crossing))
        masks.append([1 if idx in crossing else 0 for idx in range(len(edges))])
    return base_capacities, masks


def throughput_from_allocation(
    base_cut_capacities: list[int],
    masks: list[list[int]],
    allocation: tuple[int, ...],
) -> int:
    """By max-flow/min-cut, throughput equals the minimum augmented cut."""
    return min(
        base + sum(mask[idx] * allocation[idx] for idx in range(len(allocation)))
        for base, mask in zip(base_cut_capacities, masks)
    )


def format_allocation(edges: list[Edge], allocation: tuple[int, ...]) -> str:
    """Format an integer capacity allocation for console and CSV output."""
    return "; ".join(
        f"{edge.label}({edge.tail}->{edge.head})={units}"
        for edge, units in zip(edges, allocation)
    )


def optimize_extra_capacity(
    edges: list[Edge],
    nodes: list[str],
    max_extra: int,
) -> pd.DataFrame:
    """Exact search over integer allocations of x units across the given edges."""
    base_cut_capacities, masks = cut_masks(edges, nodes)
    base_array = np.array(base_cut_capacities, dtype=np.int16)
    mask_array = np.array(masks, dtype=np.int16)
    allocations = allocations_up_to(max_extra, len(edges))
    budgets = allocations.sum(axis=1)
    throughputs = (allocations @ mask_array.T + base_array).min(axis=1)
    rows = []

    for x in range(max_extra + 1):
        candidate_indices = np.flatnonzero(budgets == x)
        candidate_throughputs = throughputs[candidate_indices]
        best_local_index = int(candidate_indices[int(candidate_throughputs.argmax())])
        best_throughput = int(throughputs[best_local_index])
        best_allocation = tuple(int(value) for value in allocations[best_local_index])
        row = {
            "x": x,
            "best_throughput": best_throughput,
            "allocation": format_allocation(edges, best_allocation),
        }
        for edge, units in zip(edges, best_allocation):
            row[f"extra_{edge.label}"] = units
        rows.append(row)

    return pd.DataFrame(rows)


def candidate_added_links(nodes: list[str], edges: list[Edge]) -> list[tuple[str, str]]:
    """Return missing directed links that can be tested in Exercise 1(c)."""
    existing = {(edge.tail, edge.head) for edge in edges}
    return [
        (tail, head)
        for tail, head in product(nodes, nodes)
        if tail != head and (tail, head) not in existing
    ]


def optimize_added_link_case(max_extra: int) -> pd.DataFrame:
    """Test every missing directed link with base capacity 1 and keep the best."""
    best_by_x: dict[int, dict[str, object]] = {}

    for tail, head in candidate_added_links(NODES, BASE_EDGES):
        new_label = f"new_{tail}_{head}"
        edges = [*BASE_EDGES, Edge(new_label, tail, head, 1)]
        candidate_results = optimize_extra_capacity(edges, NODES, max_extra)

        for _, candidate_row in candidate_results.iterrows():
            x = int(candidate_row["x"])
            throughput = int(candidate_row["best_throughput"])
            current_best = best_by_x.get(x)
            link_name = f"{tail}->{head}"

            if (
                current_best is None
                or throughput > int(current_best["best_throughput"])
                or (
                    throughput == int(current_best["best_throughput"])
                    and link_name < str(current_best["added_link"])
                )
            ):
                row = candidate_row.to_dict()
                row["added_link"] = link_name
                row["added_link_capacity"] = 1
                best_by_x[x] = row

    rows = [best_by_x[x] for x in range(max_extra + 1)]
    front_columns = [
        "x",
        "added_link",
        "added_link_capacity",
        "best_throughput",
        "allocation",
    ]
    remaining_columns = [
        column for column in rows[0].keys() if column not in front_columns
    ]
    return pd.DataFrame(rows)[front_columns + remaining_columns]


def plot_throughput(
    dataframe: pd.DataFrame,
    path: Path,
    title: str,
    comparison: pd.DataFrame | None = None,
) -> None:
    """Plot best achievable throughput as the extra-capacity budget changes."""
    plt.figure(figsize=(8, 5))
    plt.plot(
        dataframe["x"],
        dataframe["best_throughput"],
        marker="o",
        label="best with added link" if comparison is not None else "best throughput",
    )
    if comparison is not None:
        plt.plot(
            comparison["x"],
            comparison["best_throughput"],
            marker="s",
            linestyle="--",
            label="existing links only",
        )
    plt.xlabel("extra capacity units x")
    plt.ylabel("maximum o-d throughput")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def print_section(title: str) -> None:
    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))


def main() -> None:
    ensure_output_dirs()

    print_section("Exercise 1(a): cuts")
    cuts = cut_table(BASE_EDGES, NODES)
    cuts.to_csv("results/problem1_cuts.csv", index=False)
    print(cuts.to_string(index=False))

    min_cut_capacity = int(cuts["capacity"].min())
    max_flow = max_flow_value(BASE_EDGES)
    print(f"\nMinimum cut capacity: {min_cut_capacity}")
    print(f"NetworkX maximum flow value: {max_flow}")
    print(f"Max-flow/min-cut theorem verified: {max_flow == min_cut_capacity}")

    print_section("Exercise 1(b): existing-link capacity augmentation")
    b_results = optimize_extra_capacity(BASE_EDGES, NODES, max_extra=20)
    b_results.to_csv("results/problem1_b_results.csv", index=False)
    print(b_results[["x", "best_throughput", "allocation"]].to_string(index=False))
    plot_throughput(
        b_results,
        Path("figures/problem1_b_throughput.png"),
        "Exercise 1(b): best throughput with existing links",
    )

    print_section("Exercise 1(c): one added link plus capacity augmentation")
    c_results = optimize_added_link_case(max_extra=20)
    c_results.to_csv("results/problem1_c_results.csv", index=False)
    print(
        c_results[
            ["x", "added_link", "best_throughput", "allocation"]
        ].to_string(index=False)
    )
    plot_throughput(
        c_results,
        Path("figures/problem1_c_added_link_throughput.png"),
        "Exercise 1(c): best throughput with one added link",
        comparison=b_results,
    )

    comparison = c_results[["x", "best_throughput"]].merge(
        b_results[["x", "best_throughput"]],
        on="x",
        suffixes=("_with_added_link", "_existing_only"),
    )
    comparison["improvement"] = (
        comparison["best_throughput_with_added_link"]
        - comparison["best_throughput_existing_only"]
    )
    print("\nComparison with 1(b):")
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
