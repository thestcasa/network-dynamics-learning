"""Exercise 3: traffic assignment.

Run from the hw1 directory:

    python src/run_problem3.py
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
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures")


@dataclass(frozen=True)
class TrafficData:
    """Traffic network arrays loaded from the homework MAT files."""

    f_given: np.ndarray
    capacities: np.ndarray
    incidence: np.ndarray
    free_flow_time: np.ndarray


def ensure_output_dirs() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)


def load_data() -> TrafficData:
    """Load flow, capacity, incidence, and free-flow travel-time arrays."""
    f_given = scipy.io.loadmat(DATA_DIR / "flow.mat")["flow"].reshape(28,)
    capacities = scipy.io.loadmat(DATA_DIR / "capacities.mat")["capacities"].reshape(28,)
    incidence = scipy.io.loadmat(DATA_DIR / "traffic.mat")["traffic"]
    free_flow_time = scipy.io.loadmat(DATA_DIR / "traveltime.mat")["traveltime"].reshape(28,)
    return TrafficData(
        f_given=f_given.astype(float),
        capacities=capacities.astype(float),
        incidence=incidence.astype(float),
        free_flow_time=free_flow_time.astype(float),
    )


def reconstruct_edges(B: np.ndarray) -> pd.DataFrame:
    """Recover tail and head nodes from each column of the incidence matrix."""
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
    """Build the directed traffic graph with travel times and capacities."""
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


def shortest_path_results(
    graph: nx.DiGraph,
) -> tuple[pd.DataFrame, list[int], list[int], float]:
    """Compute the shortest 1-to-17 path using free-flow travel times."""
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
    """Compute the maximum feasible 1-to-17 flow under edge capacities."""
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


def exogenous_inflows(
    B: np.ndarray, f_given: np.ndarray
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Build the original and reassigned net inflow vectors for the network."""
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
    """Evaluate edge delays l_e / (1 - f_e / C_e)."""
    return l / (1.0 - f / C)


def tau_prime(f: np.ndarray, l: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Derivative of the delay function, used for marginal-cost tolls."""
    return (l / C) / np.square(1.0 - f / C)


def total_travel_time_value(f: np.ndarray, l: np.ndarray, C: np.ndarray) -> float:
    """Return sum_e f_e tau_e(f_e) for a concrete flow vector."""
    return float(np.sum(f * delay(f, l, C)))


def total_travel_time_expression(
    f: cp.Variable, l: np.ndarray, C: np.ndarray
) -> cp.Expression:
    """CVXPY expression for total travel time."""
    return cp.sum(cp.multiply(l * C, cp.inv_pos(1.0 - cp.multiply(1.0 / C, f))) - l * C)


def beckmann_expression(f: cp.Variable, l: np.ndarray, C: np.ndarray) -> cp.Expression:
    """CVXPY Beckmann potential whose minimizer gives a Wardrop equilibrium."""
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
    """Solve a convex traffic-assignment problem with B f = nu constraints."""
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
    """Solve a CVXPY problem with CLARABEL, falling back to CVXPY defaults."""
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
    """Solve a CVXPY problem with ECOS while suppressing solver warnings."""
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


def solve_social_optimum(
    data: TrafficData, nu_new: np.ndarray
) -> tuple[str, float, np.ndarray, pd.DataFrame]:
    """Minimize total travel time over feasible edge flows."""
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
    """Solve the Wardrop equilibrium, optionally with edge tolls."""
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
    """Minimize the additional-delay objective from part 3(g)."""
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
    """Combine edge metadata, flow values, and solver information in one table."""
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
    """Plot the three main flow vectors for visual comparison."""
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
    shortest_table.to_csv(RESULTS_DIR / "problem3_shortest_path.csv", index=False)
    print(f"path nodes: {path_nodes}")
    print(f"path edges: {path_edges}")
    print(f"total free-flow travel time: {path_time:.6f}")

    print_section("3(b): maximum flow")
    maxflow_table, maxflow_value = maxflow_results(graph, edges)
    maxflow_table.to_csv(RESULTS_DIR / "problem3_maxflow.csv", index=False)
    print(f"max-flow value: {maxflow_value:.6f}")
    print(
        maxflow_table[maxflow_table["positive_flow"]].to_string(
            index=False, float_format="%.6f"
        )
    )

    print_section("3(c): exogenous inflow")
    nu, nu_new, nu_table = exogenous_inflows(data.incidence, data.f_given)
    nu_table.to_csv(RESULTS_DIR / "problem3_nu.csv", index=False)
    print(f"nu = {nu}")
    print(f"nu_new = {nu_new}")

    print_section("3(d): social optimum")
    social_status, social_value, f_star, social_table = solve_social_optimum(data, nu_new)
    social_table.to_csv(RESULTS_DIR / "problem3_social_optimum.csv", index=False)
    print(f"objective value: {social_value:.6f}")
    print(
        social_table[["edge", "tail", "head", "flow", "delay"]].to_string(
            index=False, float_format="%.6f"
        )
    )

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
    wardrop_table.to_csv(RESULTS_DIR / "problem3_wardrop.csv", index=False)
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
    tolls_table.to_csv(RESULTS_DIR / "problem3_tolls_total_travel_time.csv", index=False)
    print(f"||f_tolled - f_star||_2 = {tolled_norm:.6e}")

    print_section("3(g): additional-delay objective")
    add_status, add_value, f_star_add = solve_additional_delay(data, nu_new)
    omega_add = -data.free_flow_time + f_star_add * tau_prime(
        f_star_add, data.free_flow_time, data.capacities
    )
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
    additional_table.to_csv(RESULTS_DIR / "problem3_additional_delay.csv", index=False)
    print(f"additional-delay optimum objective: {add_value:.6f}")
    print(f"||f_tolled_add - f_star_add||_2 = {tolled_add_norm:.6e}")

    save_comparison_plot(
        f_star,
        f_wardrop,
        f_tolled,
        FIGURES_DIR / "problem3_flow_comparison.png",
    )

    print("\nSaved Exercise 3 CSV files and figures/problem3_flow_comparison.png")


if __name__ == "__main__":
    main()
