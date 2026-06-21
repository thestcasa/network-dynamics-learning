"""Exercise 2: Katz centrality and distributed PageRank.

Run from the hw1 directory:

    python src/run_problem2.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


NODES = [f"n{i}" for i in range(1, 16)]

# Figure 2 reconstruction.
# The dense right-hand component is the complete graph on n1,...,n6.
# No visually ambiguous edges were identified in the rendered PDF crop.
EDGES = [
    # Clique on n1,...,n6
    ("n1", "n2"),
    ("n1", "n3"),
    ("n1", "n4"),
    ("n1", "n5"),
    ("n1", "n6"),
    ("n2", "n3"),
    ("n2", "n4"),
    ("n2", "n5"),
    ("n2", "n6"),
    ("n3", "n4"),
    ("n3", "n5"),
    ("n3", "n6"),
    ("n4", "n5"),
    ("n4", "n6"),
    ("n5", "n6"),
    # Attachments and path to the left side
    ("n6", "n15"),
    ("n6", "n7"),
    ("n7", "n8"),
    ("n8", "n9"),
    # Star centered at n9
    ("n9", "n10"),
    ("n9", "n11"),
    ("n9", "n12"),
    ("n9", "n13"),
    ("n9", "n14"),
]

AMBIGUOUS_EDGES = []

KATZ_BETA = 0.15
PAGERANK_BETA = 0.15
PAGERANK_TOL = 1e-12
PAGERANK_MAX_ITER = 100_000

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"


def build_graph() -> nx.Graph:
    """Build the undirected graph reconstructed from the homework figure."""
    graph = nx.Graph()
    graph.add_nodes_from(NODES)
    graph.add_edges_from(EDGES)
    return graph


def adjacency_matrix(graph: nx.Graph) -> np.ndarray:
    """Return the adjacency matrix in the fixed node order used in the tables."""
    return nx.to_numpy_array(graph, nodelist=NODES, dtype=float)


def uniform_mu(n: int) -> np.ndarray:
    """Return the uniform preference vector used by Katz and PageRank."""
    return np.ones(n, dtype=float) / n


def katz_centrality(
    W: np.ndarray, beta: float, mu: np.ndarray
) -> Tuple[np.ndarray, float, bool]:
    """Solve the Katz linear system and report the spectral admissibility check."""
    eigenvalues = np.linalg.eigvalsh(W)
    lambda_max = float(np.max(eigenvalues))
    spectral_ok = beta < 1.0 / lambda_max
    x = np.linalg.solve(np.eye(W.shape[0]) - beta * W, mu)
    return x, lambda_max, spectral_ok


def transition_matrix_for_undirected_graph(graph: nx.Graph) -> np.ndarray:
    """Column-stochastic random-walk matrix P with P[i,j] = 1/deg(j) if j -> i."""
    W = adjacency_matrix(graph)
    degrees = W.sum(axis=0)
    if np.any(degrees == 0):
        raise ValueError("PageRank transition matrix needs non-isolated nodes.")
    return W / degrees


def distributed_pagerank(
    P: np.ndarray,
    beta: float,
    mu: np.ndarray,
    tol: float = PAGERANK_TOL,
    max_iter: int = PAGERANK_MAX_ITER,
) -> Tuple[np.ndarray, int, float]:
    """Iterate x(t+1) = beta P x(t) + (1-beta) mu until L1 difference is small."""
    x = mu.copy()
    teleport_weight = 1.0 - beta
    for iteration in range(1, max_iter + 1):
        x_next = beta * (P @ x) + teleport_weight * mu
        diff = float(np.linalg.norm(x_next - x, ord=1))
        if diff < tol:
            return x_next, iteration, diff
        x = x_next
    raise RuntimeError(f"PageRank did not converge in {max_iter} iterations.")


def sorted_centrality_table(values: np.ndarray, column_name: str) -> pd.DataFrame:
    """Create a ranked node table with centrality values and node degrees."""
    table = pd.DataFrame({"node": NODES, column_name: values})
    table["rank"] = (
        table[column_name].round(12).rank(ascending=False, method="min").astype(int)
    )
    table["degree"] = [int(build_graph().degree(node)) for node in NODES]
    return table.sort_values(["rank", "node"]).reset_index(drop=True)


def save_bar_plot(
    table: pd.DataFrame, value_column: str, title: str, output_path: Path
) -> None:
    """Save a node-ordered bar plot for a centrality vector."""
    ordered = table.sort_values(
        "node", key=lambda s: s.str.extract(r"(\d+)").astype(int)[0]
    )
    plt.figure(figsize=(9, 4.8))
    plt.bar(ordered["node"], ordered[value_column], color="#2f6f73")
    plt.xlabel("Node")
    plt.ylabel(value_column)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def print_edge_list() -> None:
    print("\nStep 1 -- Figure 2 edge list reconstruction")
    for edge in EDGES:
        print(f"  {edge[0]} -- {edge[1]}")
    if AMBIGUOUS_EDGES:
        print("\nWARNING: visually ambiguous edges were marked:")
        for edge, note in AMBIGUOUS_EDGES:
            print(f"  {edge}: {note}")
    else:
        print("\nNo ambiguous edges marked.")


def print_table(title: str, table: pd.DataFrame, value_column: str) -> None:
    print(f"\n{title}")
    print(
        table[["rank", "node", "degree", value_column]].to_string(
            index=False, float_format="%.10f"
        )
    )


def compare_n6_n9(graph: nx.Graph, katz: np.ndarray, pagerank: np.ndarray) -> None:
    idx = {node: i for i, node in enumerate(NODES)}
    print("\n2(c) -- Comparison of n6 and n9")
    for node in ["n6", "n9"]:
        neighbors = sorted(graph.neighbors(node), key=lambda name: int(name[1:]))
        print(f"{node}:")
        print(f"  degree = {graph.degree(node)}")
        print(f"  neighbors = {', '.join(neighbors)}")
        print(f"  Katz = {katz[idx[node]]:.10f}")
        print(f"  PageRank = {pagerank[idx[node]]:.10f}")
    print("Structural comments:")
    print("  n6 connects the six-node clique to the left-side path and to n15.")
    print("  n9 is the center of a leaf star, but its neighbors mostly have low degree.")
    print("  Katz rewards raw walks, so n6 benefits strongly from the dense clique.")
    print(
        "  PageRank degree-normalizes outgoing mass, "
        "so degree-one leaves send all their mass to n9."
    )


def beta_sensitivity(P: np.ndarray, mu: np.ndarray) -> pd.DataFrame:
    """Evaluate PageRank for the requested beta values and compare n6 to n9."""
    rows = []
    idx_n6 = NODES.index("n6")
    idx_n9 = NODES.index("n9")
    for beta in [0.0, 0.25, 0.5, 0.75, 1.0]:
        values, iterations, diff = distributed_pagerank(P, beta, mu)
        rows.append(
            {
                "beta": beta,
                "PR(n6)": values[idx_n6],
                "PR(n9)": values[idx_n9],
                "PR(n6)-PR(n9)": values[idx_n6] - values[idx_n9],
                "iterations": iterations,
                "final_l1_diff": diff,
            }
        )
    return pd.DataFrame(rows)


def is_monotone(values: pd.Series, tol: float = 1e-12) -> str:
    """Classify the tested beta-sensitivity sequence up to a numerical tolerance."""
    diffs = np.diff(values.to_numpy())
    nondecreasing = np.all(diffs >= -tol)
    nonincreasing = np.all(diffs <= tol)
    if nondecreasing and nonincreasing:
        return "constant"
    if nondecreasing:
        return "nondecreasing"
    if nonincreasing:
        return "nonincreasing"
    return "not monotone"


def save_beta_difference_plot(table: pd.DataFrame, output_path: Path) -> None:
    """Plot PR(n6) - PR(n9) as beta varies."""
    plt.figure(figsize=(7, 4.5))
    plt.plot(table["beta"], table["PR(n6)-PR(n9)"], marker="o", color="#8b3a3a")
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.xlabel("beta")
    plt.ylabel("PR(n6) - PR(n9)")
    plt.title("PageRank difference between n6 and n9")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)

    graph = build_graph()
    print_edge_list()
    print(f"\nGraph summary: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")

    W = adjacency_matrix(graph)
    mu = uniform_mu(len(NODES))

    katz, lambda_max, spectral_ok = katz_centrality(W, KATZ_BETA, mu)
    katz_table = sorted_centrality_table(katz, "katz")
    katz_table.to_csv(RESULTS_DIR / "problem2_katz.csv", index=False)
    save_bar_plot(
        katz_table,
        "katz",
        "Katz centrality, beta = 0.15",
        FIGURES_DIR / "problem2_katz.png",
    )

    print("\n2(a) -- Katz centrality")
    print(f"lambda_max(W) = {lambda_max:.10f}")
    print(f"1/lambda_max(W) = {1.0 / lambda_max:.10f}")
    print(f"beta = {KATZ_BETA:.2f}")
    print(f"spectral condition beta < 1/lambda_max(W): {spectral_ok}")
    print_table("Katz values", katz_table, "katz")

    P = transition_matrix_for_undirected_graph(graph)
    pagerank, iterations, final_diff = distributed_pagerank(P, PAGERANK_BETA, mu)
    pagerank_table = sorted_centrality_table(pagerank, "pagerank")
    pagerank_table.to_csv(RESULTS_DIR / "problem2_pagerank_beta015.csv", index=False)
    save_bar_plot(
        pagerank_table,
        "pagerank",
        "Distributed PageRank, beta = 0.15",
        FIGURES_DIR / "problem2_pagerank_beta015.png",
    )

    print("\n2(b) -- Distributed PageRank")
    print("Update rule: x_i(t+1) = (1-beta) mu_i + beta * sum_{j in N(i)} x_j(t)/deg(j)")
    print(f"beta = {PAGERANK_BETA:.2f}")
    print(f"tolerance = {PAGERANK_TOL:.1e}")
    print(f"iterations = {iterations}")
    print(f"final L1 difference = {final_diff:.3e}")
    print_table("PageRank values", pagerank_table, "pagerank")

    compare_n6_n9(graph, katz, pagerank)

    sensitivity = beta_sensitivity(P, mu)
    sensitivity.to_csv(RESULTS_DIR / "problem2_pagerank_beta_sensitivity.csv", index=False)
    save_beta_difference_plot(sensitivity, FIGURES_DIR / "problem2_pagerank_difference_beta.png")
    monotonicity = is_monotone(sensitivity["PR(n6)-PR(n9)"])

    print("\n2(d) -- PageRank beta sensitivity")
    print(sensitivity.to_string(index=False, float_format="%.10f"))
    print(f"Monotonicity of PR(n6)-PR(n9) over tested beta values: {monotonicity}")
    print("\nSaved outputs:")
    print(f"  {RESULTS_DIR / 'problem2_katz.csv'}")
    print(f"  {RESULTS_DIR / 'problem2_pagerank_beta015.csv'}")
    print(f"  {RESULTS_DIR / 'problem2_pagerank_beta_sensitivity.csv'}")
    print(f"  {FIGURES_DIR / 'problem2_katz.png'}")
    print(f"  {FIGURES_DIR / 'problem2_pagerank_beta015.png'}")
    print(f"  {FIGURES_DIR / 'problem2_pagerank_difference_beta.png'}")


if __name__ == "__main__":
    main()
