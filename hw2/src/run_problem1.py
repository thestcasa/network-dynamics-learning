"""Run all Problem 1 computations and output CSV/figure files."""

from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from constants import DEFAULT_SEED, FIGURES_DIR, LAMBDA, NODE_TO_INDEX, NODES, RESULTS_DIR
from plotting import save_figure
from problem1_opinion import (
    asymptotic_state,
    graph_summary,
    laplacian,
    left_consensus_vector,
    remove_edges,
    simulate_fdg,
)
from problem1_single_ctmc import (
    generator,
    simulate_hitting_times,
    simulate_return_times,
    stationary_distribution,
    theoretical_hitting_times,
    theoretical_return_time,
)
from utils import confidence_interval_95, write_rows_csv, write_single_row_csv

RETURN_RUNS = 100_000
HITTING_RUNS = 100_000
CONSENSUS_MC_RUNS = 50_000
FDG_TIME_MAX = 60.0
FDG_TIME_POINTS = 600


def format_float(value):
    return f"{float(value):.12g}"


def plot_trajectory(times, trajectory, title, path):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for index, node in enumerate(NODES):
        ax.plot(times, trajectory[:, index], label=node, linewidth=1.8)
    ax.set_xlabel("time")
    ax.set_ylabel("opinion value")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(title="node", ncol=len(NODES), loc="best")
    save_figure(fig, path)


def solve_return_time():
    rng = np.random.default_rng(DEFAULT_SEED + 11)
    start = NODE_TO_INDEX["a"]
    samples = simulate_return_times(LAMBDA, start, RETURN_RUNS, rng)
    mean, std, stderr, ci_low, ci_high = confidence_interval_95(samples)
    theoretical = theoretical_return_time(LAMBDA, start)
    row = {
        "problem": "1a_1b",
        "quantity": "return_time_to_a_including_initial_holding_time",
        "start_node": "a",
        "target_node": "a",
        "initial_holding_time_included": "yes",
        "seed": DEFAULT_SEED + 11,
        "num_runs": RETURN_RUNS,
        "simulation_mean": format_float(mean),
        "simulation_std": format_float(std),
        "simulation_standard_error": format_float(stderr),
        "simulation_ci95_low": format_float(ci_low),
        "simulation_ci95_high": format_float(ci_high),
        "theoretical_value": format_float(theoretical),
        "simulation_minus_theory": format_float(mean - theoretical),
        "relative_error": format_float((mean - theoretical) / theoretical),
        "theory_method": "stationary_cycle_formula_1_over_pi_i_omega_i",
    }
    output = RESULTS_DIR / "problem1_return_time_a.csv"
    write_single_row_csv(output, row)
    return row


def solve_hitting_time():
    rng = np.random.default_rng(DEFAULT_SEED + 12)
    start = NODE_TO_INDEX["o"]
    target = NODE_TO_INDEX["d"]
    samples = simulate_hitting_times(LAMBDA, start, target, HITTING_RUNS, rng)
    mean, std, stderr, ci_low, ci_high = confidence_interval_95(samples)
    h = theoretical_hitting_times(LAMBDA, target)
    theoretical = h[start]
    row = {
        "problem": "1c_1d",
        "quantity": "hitting_time_o_to_d",
        "start_node": "o",
        "target_node": "d",
        "seed": DEFAULT_SEED + 12,
        "num_runs": HITTING_RUNS,
        "simulation_mean": format_float(mean),
        "simulation_std": format_float(std),
        "simulation_standard_error": format_float(stderr),
        "simulation_ci95_low": format_float(ci_low),
        "simulation_ci95_high": format_float(ci_high),
        "theoretical_value": format_float(theoretical),
        "simulation_minus_theory": format_float(mean - theoretical),
        "relative_error": format_float((mean - theoretical) / theoretical),
        "theory_method": "linear_system_Q_h_equals_minus_one_with_h_d_zero",
    }
    output = RESULTS_DIR / "problem1_hitting_o_to_d.csv"
    write_single_row_csv(output, row)
    return row, h


def solve_fdg_original():
    times = np.linspace(0.0, FDG_TIME_MAX, FDG_TIME_POINTS)
    x0 = np.array([1.0, -0.5, 2.0, 0.25, 1.5], dtype=float)
    trajectory = simulate_fdg(LAMBDA, x0, times)
    pi = left_consensus_vector(LAMBDA)
    consensus_value = float(pi @ x0)
    summary = graph_summary(LAMBDA, NODES)
    final_values = trajectory[-1]
    output_fig = FIGURES_DIR / "problem1_fdg_original_trajectory.png"
    plot_trajectory(times, trajectory, "Problem 1(e): original French-DeGroot trajectory", output_fig)
    row = {
        "problem": "1e",
        "fdg_convention": "continuous_time_dxdt_equals_minus_Lx_L_diag_Lambda1_minus_Lambda",
        "initial_condition": json.dumps(dict(zip(NODES, map(float, x0)))),
        "sccs": json.dumps(summary["sccs"]),
        "sink_components": json.dumps(summary["sink_components"]),
        "num_sink_components": summary["num_sink_components"],
        "converges_to_consensus_for_every_initial_condition": "yes" if summary["single_sink_reachable_from_all"] else "no",
        "consensus_vector_pi": json.dumps(dict(zip(NODES, map(float, pi)))),
        "consensus_value_for_initial_condition": format_float(consensus_value),
        "max_abs_final_minus_consensus": format_float(np.max(np.abs(final_values - consensus_value))),
        "figure": str(output_fig.relative_to(ROOT)),
    }
    write_single_row_csv(RESULTS_DIR / "problem1_fdg_original_consensus.csv", row)
    return row, pi


def solve_consensus_variance(pi):
    rng = np.random.default_rng(DEFAULT_SEED + 13)
    variances = np.array([1.0, 2.0, 2.0, 2.0, 1.0], dtype=float)
    stds = np.sqrt(variances)
    samples = rng.normal(loc=0.0, scale=stds, size=(CONSENSUS_MC_RUNS, len(NODES)))
    consensus_samples = samples @ pi
    sample_variance = float(np.var(consensus_samples, ddof=1))
    theoretical_variance = float(np.sum(pi * pi * variances))
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.hist(consensus_samples, bins=60, density=True, color="#4c78a8", edgecolor="white", alpha=0.85)
    ax.axvline(0.0, color="black", linewidth=1.2, label="theoretical mean")
    ax.set_xlabel("consensus value")
    ax.set_ylabel("density")
    ax.set_title("Problem 1(f): Monte Carlo consensus value distribution")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best")
    output_fig = FIGURES_DIR / "problem1_consensus_value_histogram.png"
    save_figure(fig, output_fig)
    row = {
        "problem": "1f",
        "seed": DEFAULT_SEED + 13,
        "num_runs": CONSENSUS_MC_RUNS,
        "initial_mean": 0.0,
        "variances_by_node": json.dumps(dict(zip(NODES, map(float, variances)))),
        "consensus_vector_pi": json.dumps(dict(zip(NODES, map(float, pi)))),
        "theoretical_variance": format_float(theoretical_variance),
        "monte_carlo_variance": format_float(sample_variance),
        "monte_carlo_mean": format_float(np.mean(consensus_samples)),
        "monte_carlo_minus_theory": format_float(sample_variance - theoretical_variance),
        "relative_error": format_float((sample_variance - theoretical_variance) / theoretical_variance),
        "figure": str(output_fig.relative_to(ROOT)),
    }
    write_single_row_csv(RESULTS_DIR / "problem1_consensus_variance.csv", row)
    return row


def solve_removed_g():
    edges = [("d", "a"), ("d", "c"), ("a", "c"), ("b", "c")]
    modified = remove_edges(LAMBDA, NODE_TO_INDEX, edges)
    summary = graph_summary(modified, NODES)
    times = np.linspace(0.0, 120.0, FDG_TIME_POINTS)
    x0 = np.array([0.0, 3.0, -1.0, 2.0, 5.0], dtype=float)
    trajectory = simulate_fdg(modified, x0, times)
    limit = asymptotic_state(modified, x0)
    output_fig = FIGURES_DIR / "problem1_fdg_removed_edges_g.png"
    plot_trajectory(times, trajectory, "Problem 1(g): French-DeGroot after edge removals", output_fig)
    row = {
        "problem": "1g",
        "removed_edges": json.dumps(edges),
        "initial_condition": json.dumps(dict(zip(NODES, map(float, x0)))),
        "sccs": json.dumps(summary["sccs"]),
        "sink_components": json.dumps(summary["sink_components"]),
        "num_sink_components": summary["num_sink_components"],
        "converges_to_global_consensus_for_every_initial_condition": "yes" if summary["single_sink_reachable_from_all"] else "no",
        "asymptotic_state_for_initial_condition": json.dumps(dict(zip(NODES, map(float, limit)))),
        "interpretation": "two_sink_components_so_limit_depends_on_sink_component_initial_values_and_transient_routing",
        "figure": str(output_fig.relative_to(ROOT)),
    }
    write_single_row_csv(RESULTS_DIR / "problem1_removed_edges_g_summary.csv", row)
    return row


def solve_removed_h():
    edges = [("b", "o"), ("d", "a")]
    modified = remove_edges(LAMBDA, NODE_TO_INDEX, edges)
    summary = graph_summary(modified, NODES)
    times = np.linspace(0.0, 120.0, FDG_TIME_POINTS)
    cases = [
        ("case1", np.array([0.0, 3.0, -1.0, 2.0, 5.0], dtype=float)),
        ("case2", np.array([10.0, -4.0, 1.0, -2.0, 0.5], dtype=float)),
    ]
    rows = []
    for label, x0 in cases:
        trajectory = simulate_fdg(modified, x0, times)
        limit = asymptotic_state(modified, x0)
        output_fig = FIGURES_DIR / f"problem1_fdg_removed_edges_h_{label}.png"
        plot_trajectory(times, trajectory, f"Problem 1(h): French-DeGroot {label}", output_fig)
        rows.append(
            {
                "problem": "1h",
                "case": label,
                "removed_edges": json.dumps(edges),
                "initial_condition": json.dumps(dict(zip(NODES, map(float, x0)))),
                "sccs": json.dumps(summary["sccs"]),
                "sink_components": json.dumps(summary["sink_components"]),
                "num_sink_components": summary["num_sink_components"],
                "converges_to_global_consensus_for_every_initial_condition": "yes" if summary["single_sink_reachable_from_all"] else "no",
                "asymptotic_state_for_initial_condition": json.dumps(dict(zip(NODES, map(float, limit)))),
                "max_minus_min_asymptotic_state": format_float(np.max(limit) - np.min(limit)),
                "interpretation": "single_sink_component_reachable_from_all_so_global_consensus_depends_on_sink_component_initial_values",
                "figure": str(output_fig.relative_to(ROOT)),
            }
        )
    write_rows_csv(RESULTS_DIR / "problem1_removed_edges_h_summary.csv", rows)
    return rows


def update_results_log(outputs):
    log_path = ROOT / "RESULTS_LOG.md"
    content = f"""# RESULTS_LOG.md - Homework 2

This file tracks all numerical results, generated figures, assumptions, and report-writing status for HW2.

Do not insert invented values. Every numerical result must come from a script output or a checked theoretical computation.

---

## Initialization - 2026-06-04

- Folder structure created for `data/`, `src/`, `figures/`, `results/`, and `report/`.
- Placeholder source files created in `src/`.
- Minimal LaTeX report skeleton created at `report/report.tex`.
- Placeholder `report/report.pdf` created because `pdflatex` is not available in this environment.
- No simulations were run during initialization.

---

## Problem 1 run - 2026-06-04

### Global methods and assumptions

- Work performed only inside `hw2/`.
- Node order: `{NODES}`.
- Matrix convention: `Lambda[i,j] = rate/weight from source node i to destination node j`.
- CTMC generator: `Q = Lambda - diag(Lambda @ 1)`.
- Return-time convention for P1(a)-P1(b): initial holding time in node `a` is included. The simulated time starts immediately at node `a`, samples the first holding time in `a`, and stops at the next entrance to `a` after leaving.
- French-DeGroot convention for P1(e)-P1(h): continuous-time `dx/dt = -Lx`, with `L = diag(Lambda @ 1) - Lambda`, because lecture-specific notes were not available in this workspace.
- Python version: `{platform.python_version()}`.
- OS: `{platform.system()} {platform.release()}`.
- Main libraries: numpy, matplotlib, Python standard library.

### P1(a)-P1(b): return time to node a

- Method: Monte Carlo CTMC simulation for P1(a); stationary CTMC cycle formula `E_a[T_a^+] = 1 / (pi_a omega_a)` for P1(b).
- Seed: `{outputs['return']['seed']}`.
- Number of Monte Carlo runs: `{outputs['return']['num_runs']}`.
- Result file: `results/problem1_return_time_a.csv`.
- Simulation mean: `{outputs['return']['simulation_mean']}`.
- Simulation standard deviation: `{outputs['return']['simulation_std']}`.
- Simulation standard error: `{outputs['return']['simulation_standard_error']}`.
- 95 percent CI: [`{outputs['return']['simulation_ci95_low']}`, `{outputs['return']['simulation_ci95_high']}`].
- Theoretical value: `{outputs['return']['theoretical_value']}`.
- Simulation minus theory: `{outputs['return']['simulation_minus_theory']}`.
- Relative error: `{outputs['return']['relative_error']}`.

### P1(c)-P1(d): hitting time from node o to node d

- Method: Monte Carlo CTMC simulation for P1(c); linear system `sum_j Q[i,j] h_j = -1` for `i != d`, with `h_d = 0`, for P1(d).
- Seed: `{outputs['hitting']['seed']}`.
- Number of Monte Carlo runs: `{outputs['hitting']['num_runs']}`.
- Result file: `results/problem1_hitting_o_to_d.csv`.
- Simulation mean: `{outputs['hitting']['simulation_mean']}`.
- Simulation standard deviation: `{outputs['hitting']['simulation_std']}`.
- Simulation standard error: `{outputs['hitting']['simulation_standard_error']}`.
- 95 percent CI: [`{outputs['hitting']['simulation_ci95_low']}`, `{outputs['hitting']['simulation_ci95_high']}`].
- Theoretical value: `{outputs['hitting']['theoretical_value']}`.
- Simulation minus theory: `{outputs['hitting']['simulation_minus_theory']}`.
- Relative error: `{outputs['hitting']['relative_error']}`.

### P1(e): original French-DeGroot dynamics

- Method: exact matrix exponential via eigendecomposition for `exp(-Lt)x(0)`; SCC/sink-component graph check.
- Result file: `results/problem1_fdg_original_consensus.csv`.
- Figure: `figures/problem1_fdg_original_trajectory.png`.
- SCCs: `{outputs['fdg']['sccs']}`.
- Sink components: `{outputs['fdg']['sink_components']}`.
- Converges to consensus for every initial condition: `{outputs['fdg']['converges_to_consensus_for_every_initial_condition']}`.
- Consensus vector pi: `{outputs['fdg']['consensus_vector_pi']}`.
- Consensus value for saved initial condition: `{outputs['fdg']['consensus_value_for_initial_condition']}`.
- Max absolute final-minus-consensus on plotted horizon: `{outputs['fdg']['max_abs_final_minus_consensus']}`.

### P1(f): variance of consensus value

- Method: theoretical variance `sum_i pi_i^2 Var(x_i(0))`; Monte Carlo sampling of independent zero-mean Gaussian initial states with requested variances.
- Seed: `{outputs['variance']['seed']}`.
- Number of Monte Carlo runs: `{outputs['variance']['num_runs']}`.
- Result file: `results/problem1_consensus_variance.csv`.
- Figure: `figures/problem1_consensus_value_histogram.png`.
- Initial variances by node: `{outputs['variance']['variances_by_node']}`.
- Theoretical variance: `{outputs['variance']['theoretical_variance']}`.
- Monte Carlo variance: `{outputs['variance']['monte_carlo_variance']}`.
- Monte Carlo mean: `{outputs['variance']['monte_carlo_mean']}`.
- Monte Carlo minus theory: `{outputs['variance']['monte_carlo_minus_theory']}`.
- Relative error: `{outputs['variance']['relative_error']}`.

### P1(g): removed edges `(d,a), (d,c), (a,c), (b,c)`

- Method: set requested directed entries of Lambda to zero; compute SCCs and sink SCCs; simulate a representative trajectory.
- Result file: `results/problem1_removed_edges_g_summary.csv`.
- Figure: `figures/problem1_fdg_removed_edges_g.png`.
- SCCs: `{outputs['g']['sccs']}`.
- Sink components: `{outputs['g']['sink_components']}`.
- Converges to global consensus for every initial condition: `{outputs['g']['converges_to_global_consensus_for_every_initial_condition']}`.
- Asymptotic state for saved initial condition: `{outputs['g']['asymptotic_state_for_initial_condition']}`.
- Interpretation note: `{outputs['g']['interpretation']}`.

### P1(h): removed edges `(b,o), (d,a)`

- Method: set requested directed entries of Lambda to zero; compute SCCs and sink SCCs; simulate two representative trajectories.
- Result file: `results/problem1_removed_edges_h_summary.csv`.
- Figures: `figures/problem1_fdg_removed_edges_h_case1.png`, `figures/problem1_fdg_removed_edges_h_case2.png`.
- Case 1 asymptotic state: `{outputs['h'][0]['asymptotic_state_for_initial_condition']}`.
- Case 2 asymptotic state: `{outputs['h'][1]['asymptotic_state_for_initial_condition']}`.
- SCCs: `{outputs['h'][0]['sccs']}`.
- Sink components: `{outputs['h'][0]['sink_components']}`.
- Converges to global consensus for every initial condition: `{outputs['h'][0]['converges_to_global_consensus_for_every_initial_condition']}`.
- Interpretation note: `{outputs['h'][0]['interpretation']}`.

### Report-writing notes and ambiguities

- Confirm whether course lecture notes use the same continuous-time French-DeGroot convention.
- State explicitly in the report that return times include the initial holding time in node `a`.
- These factual notes were used as source material for the final report prose.

---

# Problem 2 - Many particles

Not solved in this step.

---

# Problem 3 - Open network

Not solved in this step.
"""
    log_path.write_text(content, encoding="utf-8")


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    outputs = {}
    outputs["return"] = solve_return_time()
    outputs["hitting"], _ = solve_hitting_time()
    outputs["fdg"], pi = solve_fdg_original()
    outputs["variance"] = solve_consensus_variance(pi)
    outputs["g"] = solve_removed_g()
    outputs["h"] = solve_removed_h()
    update_results_log(outputs)

    print("Problem 1 completed.")
    print("Generated result files in results/ and figures in figures/.")


if __name__ == "__main__":
    main()
