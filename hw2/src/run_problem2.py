"""Run Problem 2 simulations and save requested artifacts."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from constants import FIGURES_DIR, LAMBDA, NODE_TO_INDEX, NODES, RESULTS_DIR
from plotting import save_figure
from problem2_many_particles import (
    ctmc_components,
    simulate_node_count_runs,
    simulate_particle_return_runs,
    stationary_distribution,
)
from utils import confidence_interval_95, write_rows_csv, write_single_row_csv


NUM_PARTICLES = 100
PARTICLE_RETURN_RUNS = 2000
NODE_COUNT_RUNS = 5000
HORIZON = 60.0
TIME_STEP = 0.1
PARTICLE_SEED = 20260624
NODE_SEED = 20260625


def fmt(value):
    """Format numerical results compactly without hiding precision."""
    return f"{float(value):.12g}"


def read_problem1_return_result(path):
    """Read the saved Problem 1 return-time result for comparison."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    row = rows[0]
    return {
        "simulation_mean": float(row["simulation_mean"]),
        "simulation_standard_error": float(row["simulation_standard_error"]),
        "theoretical_value": float(row["theoretical_value"]),
        "convention": row["initial_holding_time_included"],
        "num_runs": int(row["num_runs"]),
        "seed": int(row["seed"]),
    }


def write_particle_return_csv(run_means, p1_result):
    """Save P2(a) particle-perspective return-time summary."""
    mean, std, stderr, ci_low, ci_high = confidence_interval_95(run_means)
    theory = p1_result["theoretical_value"]
    row = {
        "problem": "2a",
        "quantity": "average_return_time_to_a_across_100_particles",
        "start_node": "a",
        "target_node": "a",
        "initial_holding_time_included": "yes",
        "seed": PARTICLE_SEED,
        "num_particles": NUM_PARTICLES,
        "num_monte_carlo_runs": PARTICLE_RETURN_RUNS,
        "total_particle_return_times": NUM_PARTICLES * PARTICLE_RETURN_RUNS,
        "mean_of_run_particle_averages": fmt(mean),
        "std_of_run_particle_averages": fmt(std),
        "standard_error_of_run_particle_averages": fmt(stderr),
        "ci95_low": fmt(ci_low),
        "ci95_high": fmt(ci_high),
        "problem1_simulation_mean": fmt(p1_result["simulation_mean"]),
        "problem1_simulation_standard_error": fmt(
            p1_result["simulation_standard_error"]
        ),
        "problem1_theoretical_value": fmt(theory),
        "difference_from_problem1_simulation_mean": fmt(
            mean - p1_result["simulation_mean"]
        ),
        "difference_from_problem1_theory": fmt(mean - theory),
        "relative_error_vs_problem1_theory": fmt((mean - theory) / theory),
    }
    output_path = RESULTS_DIR / "problem2_particle_return_times.csv"
    write_single_row_csv(output_path, row)
    return output_path, row


def write_time_series_csv(sample_times, mean_counts, std_counts):
    """Save average node-count trajectories on the sampling grid."""
    rows = []
    for time_index, time_value in enumerate(sample_times):
        for node_index, node in enumerate(NODES):
            std = std_counts[time_index, node_index]
            rows.append(
                {
                    "time": fmt(time_value),
                    "node": node,
                    "mean_count": fmt(mean_counts[time_index, node_index]),
                    "std_count": fmt(std),
                    "standard_error": fmt(std / np.sqrt(NODE_COUNT_RUNS)),
                    "seed": NODE_SEED,
                    "num_particles": NUM_PARTICLES,
                    "num_monte_carlo_runs": NODE_COUNT_RUNS,
                }
            )
    output_path = RESULTS_DIR / "problem2_node_time_series.csv"
    write_rows_csv(output_path, rows)
    return output_path


def write_final_counts_csv(final_counts, pi):
    """Save final count summary and stationary expected-count comparison."""
    rows = []
    mean = np.mean(final_counts, axis=0)
    std = np.std(final_counts, axis=0, ddof=1)
    expected = NUM_PARTICLES * pi
    for node_index, node in enumerate(NODES):
        rows.append(
            {
                "node": node,
                "time": fmt(HORIZON),
                "simulated_average_final_count": fmt(mean[node_index]),
                "simulated_final_count_std": fmt(std[node_index]),
                "simulated_final_count_standard_error": fmt(
                    std[node_index] / np.sqrt(NODE_COUNT_RUNS)
                ),
                "stationary_probability": fmt(pi[node_index]),
                "theoretical_stationary_expected_count": fmt(expected[node_index]),
                "simulated_minus_theoretical": fmt(mean[node_index] - expected[node_index]),
                "relative_error": fmt((mean[node_index] - expected[node_index]) / expected[node_index]),
                "seed": NODE_SEED,
                "num_particles": NUM_PARTICLES,
                "num_monte_carlo_runs": NODE_COUNT_RUNS,
            }
        )
    output_path = RESULTS_DIR / "problem2_node_final_counts.csv"
    write_rows_csv(output_path, rows)
    return output_path, rows


def plot_node_count_trajectories(sample_times, mean_counts):
    """Plot mean node-count trajectories."""
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    colors = ["#334155", "#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
    for node_index, node in enumerate(NODES):
        ax.plot(sample_times, mean_counts[:, node_index], label=node, linewidth=1.8, color=colors[node_index])
    ax.set_xlabel("time")
    ax.set_ylabel("mean number of particles")
    ax.set_title("Problem 2 node-count trajectories")
    ax.grid(True, alpha=0.25)
    ax.legend(title="node", ncol=len(NODES), fontsize=9)
    output_path = FIGURES_DIR / "problem2_node_counts_over_time.png"
    save_figure(fig, output_path)
    return output_path


def plot_stationary_comparison(final_rows):
    """Plot simulated final counts against stationary expected counts."""
    simulated = np.array(
        [float(row["simulated_average_final_count"]) for row in final_rows]
    )
    theoretical = np.array(
        [float(row["theoretical_stationary_expected_count"]) for row in final_rows]
    )
    x_values = np.arange(len(NODES))
    width = 0.36

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.bar(x_values - width / 2, simulated, width, label="simulation", color="#1f77b4")
    ax.bar(x_values + width / 2, theoretical, width, label="N*pi", color="#ff7f0e")
    ax.set_xticks(x_values)
    ax.set_xticklabels(NODES)
    ax.set_xlabel("node")
    ax.set_ylabel("number of particles")
    ax.set_title("Final counts at T = 60 vs stationary expectation")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    output_path = FIGURES_DIR / "problem2_stationary_distribution_comparison.png"
    save_figure(fig, output_path)
    return output_path


def update_results_log(particle_row, final_rows, stationary_pi, artifact_paths):
    """Replace the Problem 2 section in RESULTS_LOG.md with factual notes."""
    log_path = Path(__file__).resolve().parents[1] / "RESULTS_LOG.md"
    text = log_path.read_text(encoding="utf-8")
    start_candidates = ["## Problem 2 - Many particles", "# Problem 2 - Many particles"]
    end_candidates = ["## Problem 3 - Open network", "# Problem 3 - Open network"]
    start_positions = [text.find(marker) for marker in start_candidates]
    end_positions = [text.find(marker) for marker in end_candidates]
    start = min(position for position in start_positions if position >= 0)
    end = min(position for position in end_positions if position >= 0)

    final_count_lines = "\n".join(
        [
            "- Node `{node}`: simulated average final count `{simulated_average_final_count}`, "
            "stationary expected count `{theoretical_stationary_expected_count}`, "
            "difference `{simulated_minus_theoretical}`.".format(**row)
            for row in final_rows
        ]
    )
    pi_text = ", ".join(
        f"{node}: {fmt(stationary_pi[index])}" for index, node in enumerate(NODES)
    )

    new_section = f"""## Problem 2 - Many particles

Solved on 2026-06-04.

### Methods and assumptions

- Work performed only inside `hw2/`.
- Same closed network as Problem 1; node order `['o', 'a', 'b', 'c', 'd']`.
- Matrix convention: `Lambda[i,j]` is the rate from source node `i` to destination node `j`.
- CTMC definitions: `omega = Lambda @ 1`, `Q = Lambda - diag(omega)`, `P_jump = diag(omega)^(-1) Lambda`.
- `P_jump` was used only as the embedded jump-chain transition matrix after a departure event.
- No uniformized transition matrix was used.
- P2(a) uses the same return-time convention as Problem 1: the initial holding time in node `a` is included.
- P2(b) uses a direct event-based node-count CTMC: node `i` has departure rate `n_i(t) omega_i`.

### P2(a): particle perspective

- Seed: `{PARTICLE_SEED}`.
- Number of particles per Monte Carlo run: `{NUM_PARTICLES}`.
- Number of Monte Carlo runs: `{PARTICLE_RETURN_RUNS}`.
- Total particle return times simulated: `{NUM_PARTICLES * PARTICLE_RETURN_RUNS}`.
- Result file: `results/problem2_particle_return_times.csv`.
- Mean of run particle averages: `{particle_row['mean_of_run_particle_averages']}`.
- Standard deviation of run particle averages: `{particle_row['std_of_run_particle_averages']}`.
- Standard error of run particle averages: `{particle_row['standard_error_of_run_particle_averages']}`.
- 95 percent CI: [`{particle_row['ci95_low']}`, `{particle_row['ci95_high']}`].
- Problem 1 simulation mean: `{particle_row['problem1_simulation_mean']}`.
- Problem 1 theoretical value: `{particle_row['problem1_theoretical_value']}`.
- Difference from Problem 1 simulation mean: `{particle_row['difference_from_problem1_simulation_mean']}`.
- Difference from Problem 1 theory: `{particle_row['difference_from_problem1_theory']}`.

### P2(b): node perspective

- Seed: `{NODE_SEED}`.
- Number of particles: `{NUM_PARTICLES}`.
- Initial condition: all particles in node `a`.
- Simulation horizon: `{fmt(HORIZON)}`.
- Number of Monte Carlo runs: `{NODE_COUNT_RUNS}`.
- Time-series sampling step: `{fmt(TIME_STEP)}`.
- Result files:
  - `results/problem2_node_final_counts.csv`
  - `results/problem2_node_time_series.csv`
- Figures:
  - `figures/problem2_node_counts_over_time.png`
  - `figures/problem2_stationary_distribution_comparison.png`
- Stationary distribution pi: `{pi_text}`.
- Final count comparison:
{final_count_lines}

### Generated artifacts

- `{artifact_paths['particle']}`
- `{artifact_paths['final_counts']}`
- `{artifact_paths['time_series']}`
- `{artifact_paths['trajectory_figure']}`
- `{artifact_paths['stationary_figure']}`

### Report-writing notes

- Explain that the particle-perspective average agrees with Problem 1 because each particle follows the same single-particle CTMC independently.
- Explain that average final node counts at `T = 60` are close to `N*pi`, with deviations attributable to finite time and Monte Carlo error.
- These factual notes were used as source material for the final Problem 2 report prose.

---

"""
    log_path.write_text(text[:start] + new_section + text[end:], encoding="utf-8")


def main():
    omega, q_matrix, p_jump = ctmc_components(LAMBDA)
    pi = stationary_distribution(q_matrix)
    start_index = NODE_TO_INDEX["a"]

    p1_result = read_problem1_return_result(RESULTS_DIR / "problem1_return_time_a.csv")
    return_run_means = simulate_particle_return_runs(
        NUM_PARTICLES,
        PARTICLE_RETURN_RUNS,
        start_index,
        omega,
        p_jump,
        PARTICLE_SEED,
    )
    particle_path, particle_row = write_particle_return_csv(return_run_means, p1_result)

    sample_times = np.round(np.arange(0.0, HORIZON + TIME_STEP / 2.0, TIME_STEP), 10)
    mean_counts, std_counts, final_counts = simulate_node_count_runs(
        NUM_PARTICLES,
        HORIZON,
        sample_times,
        NODE_COUNT_RUNS,
        omega,
        p_jump,
        NODE_SEED,
    )
    time_series_path = write_time_series_csv(sample_times, mean_counts, std_counts)
    final_counts_path, final_rows = write_final_counts_csv(final_counts, pi)
    trajectory_figure_path = plot_node_count_trajectories(sample_times, mean_counts)
    stationary_figure_path = plot_stationary_comparison(final_rows)

    artifact_paths = {
        "particle": particle_path.as_posix(),
        "final_counts": final_counts_path.as_posix(),
        "time_series": time_series_path.as_posix(),
        "trajectory_figure": trajectory_figure_path.as_posix(),
        "stationary_figure": stationary_figure_path.as_posix(),
    }
    update_results_log(particle_row, final_rows, pi, artifact_paths)

    print("Problem 2 complete.")
    print(f"Particle return mean: {particle_row['mean_of_run_particle_averages']}")
    print(f"Stationary distribution: {[fmt(value) for value in pi]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
