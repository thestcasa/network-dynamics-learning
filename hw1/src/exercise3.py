"""Exercise 3: traffic assignment.

Run from the project root:

    python src/exercise3.py
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import warnings

import cvxpy as cp
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.io


SOURCE = 1
SINK = 17
CAPACITY_FACTOR = 0.999
RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures")


@dataclass(frozen=True)
class TrafficData:
    f_given: np.ndarray
    capacities: np.ndarray
    incidence: np.ndarray
    free_flow_time: np.ndarray


def ensure_output_dirs() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)


def load_data() -> TrafficData:
    f_given = scipy.io.loadmat("data/flow.mat")["flow"].reshape(28,)
    capacities = scipy.io.loadmat("data/capacities.mat")["capacities"].reshape(28,)
    incidence = scipy.io.loadmat("data/traffic.mat")["traffic"]
    free_flow_time = scipy.io.loadmat("data/traveltime.mat")["traveltime"].reshape(28,)
    return TrafficData(
        f_given=f_given.astype(float),
        capacities=capacities.astype(float),
        incidence=incidence.astype(float),
        free_flow_time=free_flow_time.astype(float),
    )


def reconstruct_edges(B: np.ndarray) -> pd.DataFrame:
    rows = []
    for edge_idx in range(B.shape[1]):
        tail_candidates = np.flatnonzero(np.isclose(B[:, edge_idx], 1.0))
        head_candidates = np.flatnonzero(np.isclose(B[:, edge_idx], -1.0))
        if len(tail_candidates) != 1 or len(head_candidates) != 1:
            raise ValueError(f"Column {edge_idx + 1} is not a valid directed incidence column.")
        rows.append(
            {
                "edge": edge_idx + 1,
                "tail": int(tail_candidates[0] + 1),
                "head": int(head_candidates[0] + 1),
            }
        )
    return pd.DataFrame(rows)


def build_graph(edges: pd.DataFrame, l: np.ndarray, C: np.ndarray) -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.add_nodes_from(range(1, int(edges[["tail", "head"]].to_numpy().max()) + 1))
    for row in edges.itertuples(index=False):
        idx = int(row.edge) - 1
        graph.add_edge(
            int(row.tail),
            int(row.head),
            edge=int(row.edge),
            weight=float(l[idx]),
            capacity=float(C[idx]),
        )
    return graph


def shortest_path_results(graph: nx.DiGraph) -> tuple[pd.DataFrame, list[int], list[int], float]:
    path_nodes = nx.shortest_path(graph, SOURCE, SINK, weight="weight")
    path_edges = []
    total_time = 0.0
    rows = []
    for tail, head in zip(path_nodes[:-1], path_nodes[1:]):
        data = graph[tail][head]
        path_edges.append(int(data["edge"]))
        total_time += float(data["weight"])
        rows.append(
            {
                "order": len(rows) + 1,
                "tail": tail,
                "head": head,
                "edge": int(data["edge"]),
                "free_flow_time": float(data["weight"]),
            }
        )
    table = pd.DataFrame(rows)
    table["path_nodes"] = " -> ".join(str(node) for node in path_nodes)
    table["path_edges"] = " -> ".join(str(edge) for edge in path_edges)
    table["total_free_flow_time"] = total_time
    return table, path_nodes, path_edges, total_time


def maxflow_results(graph: nx.DiGraph, edges: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    value, flow_dict = nx.maximum_flow(graph, SOURCE, SINK, capacity="capacity")
    rows = []
    for row in edges.itertuples(index=False):
        flow_value = float(flow_dict[int(row.tail)].get(int(row.head), 0.0))
        rows.append(
            {
                "edge": int(row.edge),
                "tail": int(row.tail),
                "head": int(row.head),
                "capacity": float(graph[int(row.tail)][int(row.head)]["capacity"]),
                "maxflow_edge_flow": flow_value,
                "positive_flow": flow_value > 1e-8,
                "maxflow_value": float(value),
            }
        )
    return pd.DataFrame(rows), float(value)


def exogenous_inflows(B: np.ndarray, f_given: np.ndarray) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    nu = B @ f_given
    nu_new = np.zeros_like(nu)
    nu_new[0] = nu[0]
    nu_new[16] = -nu[0]
    table = pd.DataFrame(
        {
            "node": np.arange(1, B.shape[0] + 1),
            "nu_original": nu,
            "nu_new": nu_new,
        }
    )
    return nu, nu_new, table


def delay(f: np.ndarray, l: np.ndarray, C: np.ndarray) -> np.ndarray:
    return l / (1.0 - f / C)


def tau_prime(f: np.ndarray, l: np.ndarray, C: np.ndarray) -> np.ndarray:
    return (l / C) / np.square(1.0 - f / C)


def total_travel_time_value(f: np.ndarray, l: np.ndarray, C: np.ndarray) -> float:
    return float(np.sum(f * delay(f, l, C)))


def total_travel_time_expression(f: cp.Variable, l: np.ndarray, C: np.ndarray) -> cp.Expression:
    return cp.sum(cp.multiply(l * C, cp.inv_pos(1.0 - cp.multiply(1.0 / C, f))) - l * C)


def beckmann_expression(f: cp.Variable, l: np.ndarray, C: np.ndarray) -> cp.Expression:
    return cp.sum(cp.multiply(-l * C, cp.log(1.0 - cp.multiply(1.0 / C, f))))


def solve_convex_flow(
    objective: cp.Expression,
    f: cp.Variable,
    B: np.ndarray,
    nu_new: np.ndarray,
    C: np.ndarray,
    label: str,
    prefer_ecos: bool = False,
) -> tuple[str, float, np.ndarray]:
    problem = cp.Problem(
        cp.Minimize(objective),
        [
            B @ f == nu_new,
            f >= 0.0,
            f <= CAPACITY_FACTOR * C,
        ],
    )
    if prefer_ecos:
        try:
            value = solve_with_ecos(problem)
        except cp.error.SolverError:
            value = solve_with_clarabel(problem)
    else:
        value = solve_with_clarabel(problem)
    if problem.status == cp.OPTIMAL_INACCURATE:
        try:
            value = solve_with_ecos(problem)
        except cp.error.SolverError:
            pass

    print(f"{label} solver status: {problem.status}")
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE} or f.value is None:
        raise RuntimeError(f"{label} did not solve successfully: {problem.status}")

    flow = np.asarray(f.value, dtype=float).reshape(-1)
    flow[np.abs(flow) < 1e-7] = 0.0
    return str(problem.status), float(value), flow


def solve_with_clarabel(problem: cp.Problem) -> float:
    try:
        return float(
            problem.solve(
                solver=cp.CLARABEL,
                verbose=False,
                tol_gap_abs=1e-10,
                tol_gap_rel=1e-10,
                tol_feas=1e-10,
                max_iter=1000,
            )
        )
    except cp.error.SolverError:
        return float(problem.solve(verbose=False))


def solve_with_ecos(problem: cp.Problem) -> float:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        return float(
            problem.solve(
                solver=cp.ECOS,
                verbose=False,
                abstol=1e-10,
                reltol=1e-10,
                feastol=1e-10,
                max_iters=1000,
            )
        )


def solve_social_optimum(data: TrafficData, nu_new: np.ndarray) -> tuple[str, float, np.ndarray, pd.DataFrame]:
    f = cp.Variable(data.capacities.size)
    objective = total_travel_time_expression(f, data.free_flow_time, data.capacities)
    status, value, flow = solve_convex_flow(
        objective,
        f,
        data.incidence,
        nu_new,
        data.capacities,
        "Social optimum",
    )
    delays = delay(flow, data.free_flow_time, data.capacities)
    table = edge_flow_table(flow, data, status, value)
    table["delay"] = delays
    table["edge_total_travel_time"] = flow * delays
    return status, value, flow, table


def solve_wardrop(
    data: TrafficData,
    nu_new: np.ndarray,
    tolls: np.ndarray | None = None,
    label: str = "Wardrop",
) -> tuple[str, float, np.ndarray]:
    f = cp.Variable(data.capacities.size)
    objective = beckmann_expression(f, data.free_flow_time, data.capacities)
    if tolls is not None:
        objective += tolls @ f
    return solve_convex_flow(
        objective,
        f,
        data.incidence,
        nu_new,
        data.capacities,
        label,
        prefer_ecos=tolls is not None,
    )


def solve_additional_delay(data: TrafficData, nu_new: np.ndarray) -> tuple[str, float, np.ndarray]:
    f = cp.Variable(data.capacities.size)
    total_expression = total_travel_time_expression(f, data.free_flow_time, data.capacities)
    objective = total_expression - data.free_flow_time @ f
    return solve_convex_flow(
        objective,
        f,
        data.incidence,
        nu_new,
        data.capacities,
        "Additional-delay optimum",
    )


def edge_flow_table(
    flow: np.ndarray,
    data: TrafficData,
    status: str,
    objective_value: float,
    edges: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if edges is None:
        edges = reconstruct_edges(data.incidence)
    table = edges.copy()
    table["flow"] = flow
    table["capacity"] = data.capacities
    table["free_flow_time"] = data.free_flow_time
    table["capacity_utilization"] = flow / data.capacities
    table["solver_status"] = status
    table["objective_value"] = objective_value
    return table


def save_comparison_plot(
    social: np.ndarray,
    wardrop: np.ndarray,
    tolled: np.ndarray,
    output_path: Path,
) -> None:
    edge_numbers = np.arange(1, social.size + 1)
    plt.figure(figsize=(10, 5.5))
    plt.plot(edge_numbers, social, marker="o", label="social optimum")
    plt.plot(edge_numbers, wardrop, marker="s", label="Wardrop")
    plt.plot(edge_numbers, tolled, marker="^", label="tolled Wardrop")
    plt.xlabel("edge")
    plt.ylabel("flow")
    plt.title("Exercise 3 flow comparison")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def markdown_table(dataframe: pd.DataFrame, columns: list[str] | None = None, digits: int = 6) -> str:
    if columns is not None:
        dataframe = dataframe[columns]

    formatted = dataframe.copy()
    for column in formatted.columns:
        if pd.api.types.is_float_dtype(formatted[column]):
            formatted[column] = formatted[column].map(lambda value: f"{value:.{digits}f}")

    headers = [str(column) for column in formatted.columns]
    rows = [[str(value) for value in row] for row in formatted.itertuples(index=False, name=None)]
    widths = [max(len(headers[idx]), *(len(row[idx]) for row in rows)) for idx in range(len(headers))]
    header_line = "| " + " | ".join(header.ljust(widths[idx]) for idx, header in enumerate(headers)) + " |"
    separator = "| " + " | ".join("-" * width for width in widths) + " |"
    body = [
        "| " + " | ".join(row[idx].ljust(widths[idx]) for idx in range(len(headers))) + " |"
        for row in rows
    ]
    return "\n".join([header_line, separator, *body])


def append_results_log(
    path_nodes: list[int],
    path_edges: list[int],
    path_time: float,
    maxflow_value: float,
    nu_table: pd.DataFrame,
    social_table: pd.DataFrame,
    wardrop_table: pd.DataFrame,
    tolls_table: pd.DataFrame,
    additional_table: pd.DataFrame,
    poa: float,
    tolled_norm: float,
    tolled_add_norm: float,
) -> None:
    log_path = Path("RESULTS_LOG.md")
    existing = log_path.read_text(encoding="utf-8") if log_path.exists() else "# Results Log\n"
    exercise3_start = existing.find("\n## Exercise 3")
    if exercise3_start != -1:
        existing = existing[:exercise3_start].rstrip() + "\n"
    else:
        existing = existing.rstrip() + "\n"

    social_summary = social_table[["edge", "tail", "head", "flow", "delay", "edge_total_travel_time"]]
    wardrop_summary = wardrop_table[
        ["edge", "tail", "head", "wardrop_flow", "wardrop_delay", "wardrop_total_travel_time"]
    ]
    tolls_summary = tolls_table[["edge", "omega", "f_tolled", "f_star", "flow_difference"]]
    additional_summary = additional_table[
        ["edge", "f_star_add", "omega_add", "f_tolled_add", "flow_difference"]
    ]

    lines = [
        existing,
        "## Exercise 3",
        "",
        "### 3(a): shortest path",
        "",
        f"Path nodes: {' -> '.join(str(node) for node in path_nodes)}",
        f"Path edges: {' -> '.join(str(edge) for edge in path_edges)}",
        f"Total free-flow travel time: {path_time:.6f}",
        "",
        "### 3(b): maximum flow",
        "",
        f"Maximum flow from node 1 to node 17: {maxflow_value:.6f}",
        "",
        "### 3(c): exogenous inflow",
        "",
        markdown_table(nu_table, ["node", "nu_original", "nu_new"], digits=3),
        "",
        "### 3(d): social optimum",
        "",
        f"Solver status: {social_table['solver_status'].iloc[0]}",
        f"Objective value / total travel time: {social_table['objective_value'].iloc[0]:.6f}",
        "",
        markdown_table(social_summary, digits=6),
        "",
        "### 3(e): Wardrop equilibrium",
        "",
        f"Solver status: {wardrop_table['solver_status'].iloc[0]}",
        f"Total travel time: {wardrop_table['wardrop_total_travel_time'].sum():.6f}",
        f"Price of anarchy: {poa:.6f}",
        "",
        markdown_table(wardrop_summary, digits=6),
        "",
        "### 3(f): marginal-cost tolls",
        "",
        f"Solver status: {tolls_table['solver_status'].iloc[0]}",
        f"||f_tolled - f_star||_2: {tolled_norm:.6e}",
        "",
        markdown_table(tolls_summary, digits=6),
        "",
        "### 3(g): additional-delay objective",
        "",
        f"Optimum solver status: {additional_table['optimum_status'].iloc[0]}",
        f"Tolled Wardrop solver status: {additional_table['tolled_status'].iloc[0]}",
        f"||f_tolled_add - f_star_add||_2: {tolled_add_norm:.6e}",
        "",
        markdown_table(additional_summary, digits=6),
        "",
        "Generated files:",
        "",
        "- `results/exercise3_shortest_path.csv`",
        "- `results/exercise3_maxflow.csv`",
        "- `results/exercise3_nu.csv`",
        "- `results/exercise3_social_optimum.csv`",
        "- `results/exercise3_wardrop.csv`",
        "- `results/exercise3_tolls_total_travel_time.csv`",
        "- `results/exercise3_additional_delay.csv`",
        "- `figures/exercise3_flow_comparison.png`",
        "",
    ]
    log_path.write_text("\n".join(lines), encoding="utf-8")


def print_section(title: str) -> None:
    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))


def main() -> None:
    ensure_output_dirs()
    data = load_data()
    edges = reconstruct_edges(data.incidence)
    graph = build_graph(edges, data.free_flow_time, data.capacities)

    print_section("3(a): shortest path")
    shortest_table, path_nodes, path_edges, path_time = shortest_path_results(graph)
    shortest_table.to_csv(RESULTS_DIR / "exercise3_shortest_path.csv", index=False)
    print(f"path nodes: {path_nodes}")
    print(f"path edges: {path_edges}")
    print(f"total free-flow travel time: {path_time:.6f}")

    print_section("3(b): maximum flow")
    maxflow_table, maxflow_value = maxflow_results(graph, edges)
    maxflow_table.to_csv(RESULTS_DIR / "exercise3_maxflow.csv", index=False)
    print(f"max-flow value: {maxflow_value:.6f}")
    print(maxflow_table[maxflow_table["positive_flow"]].to_string(index=False, float_format="%.6f"))

    print_section("3(c): exogenous inflow")
    nu, nu_new, nu_table = exogenous_inflows(data.incidence, data.f_given)
    nu_table.to_csv(RESULTS_DIR / "exercise3_nu.csv", index=False)
    print(f"nu = {nu}")
    print(f"nu_new = {nu_new}")

    print_section("3(d): social optimum")
    social_status, social_value, f_star, social_table = solve_social_optimum(data, nu_new)
    social_table.to_csv(RESULTS_DIR / "exercise3_social_optimum.csv", index=False)
    print(f"objective value: {social_value:.6f}")
    print(social_table[["edge", "tail", "head", "flow", "delay"]].to_string(index=False, float_format="%.6f"))

    print_section("3(e): Wardrop equilibrium")
    wardrop_status, _, f_wardrop = solve_wardrop(data, nu_new, label="Wardrop")
    wardrop_delays = delay(f_wardrop, data.free_flow_time, data.capacities)
    wardrop_total = total_travel_time_value(f_wardrop, data.free_flow_time, data.capacities)
    social_total = total_travel_time_value(f_star, data.free_flow_time, data.capacities)
    poa = wardrop_total / social_total if social_total > 0 else np.nan
    wardrop_table = edges.copy()
    wardrop_table["wardrop_flow"] = f_wardrop
    wardrop_table["wardrop_delay"] = wardrop_delays
    wardrop_table["wardrop_total_travel_time"] = f_wardrop * wardrop_delays
    wardrop_table["social_flow"] = f_star
    wardrop_table["social_total_travel_time"] = social_table["edge_total_travel_time"]
    wardrop_table["total_travel_time_social"] = social_total
    wardrop_table["total_travel_time_wardrop"] = wardrop_total
    wardrop_table["price_of_anarchy"] = poa
    wardrop_table["solver_status"] = wardrop_status
    wardrop_table.to_csv(RESULTS_DIR / "exercise3_wardrop.csv", index=False)
    print(f"total travel time at Wardrop: {wardrop_total:.6f}")
    print(f"total travel time at social optimum: {social_total:.6f}")
    print(f"price of anarchy: {poa:.6f}")

    print_section("3(f): marginal-cost tolls")
    omega = f_star * tau_prime(f_star, data.free_flow_time, data.capacities)
    toll_status, _, f_tolled = solve_wardrop(data, nu_new, tolls=omega, label="Tolled Wardrop")
    tolled_norm = float(np.linalg.norm(f_tolled - f_star, ord=2))
    tolls_table = edges.copy()
    tolls_table["omega"] = omega
    tolls_table["f_star"] = f_star
    tolls_table["f_tolled"] = f_tolled
    tolls_table["flow_difference"] = f_tolled - f_star
    tolls_table["norm_difference"] = tolled_norm
    tolls_table["solver_status"] = toll_status
    tolls_table.to_csv(RESULTS_DIR / "exercise3_tolls_total_travel_time.csv", index=False)
    print(f"||f_tolled - f_star||_2 = {tolled_norm:.6e}")

    print_section("3(g): additional-delay objective")
    add_status, add_value, f_star_add = solve_additional_delay(data, nu_new)
    omega_add = -data.free_flow_time + f_star_add * tau_prime(f_star_add, data.free_flow_time, data.capacities)
    tolled_add_status, _, f_tolled_add = solve_wardrop(
        data,
        nu_new,
        tolls=omega_add,
        label="Additional-delay tolled Wardrop",
    )
    tolled_add_norm = float(np.linalg.norm(f_tolled_add - f_star_add, ord=2))
    additional_table = edges.copy()
    additional_table["f_star_add"] = f_star_add
    additional_table["additional_delay_objective"] = add_value
    additional_table["omega_add"] = omega_add
    additional_table["f_tolled_add"] = f_tolled_add
    additional_table["flow_difference"] = f_tolled_add - f_star_add
    additional_table["norm_difference"] = tolled_add_norm
    additional_table["optimum_status"] = add_status
    additional_table["tolled_status"] = tolled_add_status
    additional_table.to_csv(RESULTS_DIR / "exercise3_additional_delay.csv", index=False)
    print(f"additional-delay optimum objective: {add_value:.6f}")
    print(f"||f_tolled_add - f_star_add||_2 = {tolled_add_norm:.6e}")

    save_comparison_plot(
        f_star,
        f_wardrop,
        f_tolled,
        FIGURES_DIR / "exercise3_flow_comparison.png",
    )
    append_results_log(
        path_nodes,
        path_edges,
        path_time,
        maxflow_value,
        nu_table,
        social_table,
        wardrop_table,
        tolls_table,
        additional_table,
        poa,
        tolled_norm,
        tolled_add_norm,
    )

    print("\nSaved Exercise 3 CSV files and figures/exercise3_flow_comparison.png")


if __name__ == "__main__":
    main()
