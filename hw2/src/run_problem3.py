"""Run Problem 3 open-network simulations and save requested artifacts."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from constants import DEFAULT_SEED, FIGURES_DIR, LAMBDA_OPEN, NODES, RESULTS_DIR
from plotting import save_figure
from problem3_open_network import (
    fixed_rate_capacity,
    mean_absorption_times,
    open_network_components,
    simulate_open_network,
    simulate_replicates,
    theoretical_loads,
    total_population_slope,
)
from utils import write_rows_csv


PROPORTIONAL_SEED = DEFAULT_SEED + 300
FIXED_SEED = DEFAULT_SEED + 301
PROPORTIONAL_SCAN_SEED = DEFAULT_SEED + 302
FIXED_SCAN_SEED = DEFAULT_SEED + 303

PROPORTIONAL_MAIN_LAMBDA = 100.0
PROPORTIONAL_MAIN_HORIZON = 60.0
PROPORTIONAL_SCAN_HORIZON = 120.0
PROPORTIONAL_SCAN_REPLICATES = 2
PROPORTIONAL_SCAN_LAMBDAS = [10.0, 50.0, 100.0, 200.0, 500.0]

FIXED_MAIN_LAMBDA = 2.0
FIXED_MAIN_HORIZON = 6000.0
FIXED_SCAN_HORIZON = 6000.0
FIXED_SCAN_REPLICATES = 3
FIXED_SCAN_LAMBDAS = [0.5, 1.0, 1.2, 1.3, 4.0 / 3.0, 1.35, 1.5, 2.0]

PROPORTIONAL_DRIFT_FRACTION = 0.05
FIXED_DRIFT_THRESHOLD = 0.005
FIXED_FINAL_TOTAL_THRESHOLD = 1000.0


def fmt(value):
    """Format numerical results compactly without hiding precision."""
    return f"{float(value):.12g}"


def sample_grid(horizon, step):
    """Return an inclusive sampling grid for the requested horizon."""
    return np.round(np.arange(0.0, horizon + step / 2.0, step), 10)


def write_time_series_csv(path, result, arrival_rate, service_mode, seed):
    """Save one sampled trajectory in long format."""
    rows = []
    for time_index, time_value in enumerate(result.sample_times):
        total = float(np.sum(result.counts[time_index]))
        for node_index, node in enumerate(NODES):
            rows.append(
                {
                    "time": fmt(time_value),
                    "node": node,
                    "count": fmt(result.counts[time_index, node_index]),
                    "total_count": fmt(total),
                    "lambda": fmt(arrival_rate),
                    "service_mode": service_mode,
                    "seed": seed,
                }
            )
    write_rows_csv(path, rows)


def plot_counts(path, result, title):
    """Plot node counts over time for one trajectory."""
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    colors = ["#334155", "#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
    for node_index, node in enumerate(NODES):
        ax.plot(
            result.sample_times,
            result.counts[:, node_index],
            label=node,
            linewidth=1.4,
            color=colors[node_index],
        )
    ax.set_xlabel("time")
    ax.set_ylabel("number of particles")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(title="node", ncol=len(NODES), fontsize=9)
    save_figure(fig, path)


def summarize_scan_replicates(
    lambda_value,
    service_mode,
    horizon,
    replicate_results,
    omega,
    proportional_drift_limit,
    fixed_capacity,
):
    """Summarize replicated trajectories for one lambda value."""
    final_totals = np.array(
        [np.sum(result.counts[-1]) for result in replicate_results], dtype=float
    )
    slopes = np.array(
        [
            total_population_slope(result.sample_times, result.counts)
            for result in replicate_results
        ],
        dtype=float,
    )
    mean_slope = float(np.mean(slopes))
    mean_final_total = float(np.mean(final_totals))

    if service_mode == "proportional":
        drift_limit = proportional_drift_limit * lambda_value
        simulation_classification = (
            "stable_by_scan" if mean_slope <= drift_limit else "blow_up_by_scan"
        )
        theoretical_classification = "stable_for_any_finite_lambda"
        capacity = "unbounded"
        bottleneck_node = ""
    else:
        drift_limit = FIXED_DRIFT_THRESHOLD
        simulation_classification = (
            "stable_by_scan"
            if mean_slope <= drift_limit
            and mean_final_total <= FIXED_FINAL_TOTAL_THRESHOLD
            else "blow_up_by_scan"
        )
        theoretical_classification = (
            "stable_below_bottleneck"
            if lambda_value < fixed_capacity
            else "not_stable_at_or_above_bottleneck"
        )
        capacity = fmt(fixed_capacity)
        _, bottleneck_index, _ = fixed_rate_capacity(omega)
        bottleneck_node = NODES[bottleneck_index]

    row = {
        "service_mode": service_mode,
        "lambda": fmt(lambda_value),
        "horizon": fmt(horizon),
        "num_replicates": len(replicate_results),
        "mean_final_total": fmt(mean_final_total),
        "max_final_total": fmt(np.max(final_totals)),
        "mean_second_half_total_slope": fmt(mean_slope),
        "max_second_half_total_slope": fmt(np.max(slopes)),
        "drift_threshold": fmt(drift_limit),
        "final_total_threshold": (
            "" if service_mode == "proportional" else fmt(FIXED_FINAL_TOTAL_THRESHOLD)
        ),
        "simulation_classification": simulation_classification,
        "theoretical_classification": theoretical_classification,
        "theoretical_capacity_lambda": capacity,
        "theoretical_bottleneck_node": bottleneck_node,
    }
    return row


def run_scan(
    lambdas,
    service_mode,
    horizon,
    sample_step,
    num_replicates,
    seed,
    omega,
    p_jump_cdf,
    fixed_capacity,
):
    """Run the stability scan for one service mode."""
    rows = []
    sample_times = sample_grid(horizon, sample_step)
    for lambda_index, lambda_value in enumerate(lambdas):
        replicate_results = simulate_replicates(
            lambda_value,
            horizon,
            sample_times,
            omega,
            p_jump_cdf,
            service_mode,
            num_replicates,
            seed + 1000 * lambda_index,
        )
        rows.append(
            summarize_scan_replicates(
                lambda_value,
                service_mode,
                horizon,
                replicate_results,
                omega,
                PROPORTIONAL_DRIFT_FRACTION,
                fixed_capacity,
            )
        )
    return rows


def plot_scan(path, scan_rows, title):
    """Plot scan drift and final totals against lambda."""
    lambdas = np.array([float(row["lambda"]) for row in scan_rows], dtype=float)
    slopes = np.array(
        [float(row["mean_second_half_total_slope"]) for row in scan_rows], dtype=float
    )
    final_totals = np.array(
        [float(row["mean_final_total"]) for row in scan_rows], dtype=float
    )

    fig, axes = plt.subplots(2, 1, figsize=(7.6, 6.4), sharex=True)
    axes[0].plot(lambdas, slopes, marker="o", color="#1f77b4")
    axes[0].axhline(0.0, color="#334155", linewidth=1.0)
    axes[0].set_ylabel("mean second-half drift")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(lambdas, final_totals, marker="o", color="#d62728")
    axes[1].set_xlabel("lambda")
    axes[1].set_ylabel("mean final total")
    axes[1].grid(True, alpha=0.25)
    fig.suptitle(title)
    save_figure(fig, path)


def update_results_log(
    omega,
    absorption_times,
    fixed_capacity,
    fixed_bottleneck_index,
    fixed_capacity_by_node,
    proportional_scan_rows,
    fixed_scan_rows,
):
    """Replace the Problem 3 section in RESULTS_LOG.md with factual notes."""
    log_path = Path(__file__).resolve().parents[1] / "RESULTS_LOG.md"
    text = log_path.read_text(encoding="utf-8")
    marker = "# Problem 3 - Open network"
    start = text.find(marker)
    if start < 0:
        marker = "## Problem 3 - Open network"
        start = text.find(marker)
    if start < 0:
        raise ValueError("Could not find Problem 3 section in RESULTS_LOG.md")

    proportional_lambdas = ", ".join(
        f"`{row['lambda']}`" for row in proportional_scan_rows
    )
    fixed_lambdas = ", ".join(f"`{row['lambda']}`" for row in fixed_scan_rows)
    omega_text = ", ".join(
        f"{node}: {fmt(omega[index])}" for index, node in enumerate(NODES)
    )
    absorption_text = ", ".join(
        f"{node}: {fmt(absorption_times[index])}" for index, node in enumerate(NODES)
    )
    fixed_capacity_text = ", ".join(
        f"{node}: {fmt(fixed_capacity_by_node[index])}"
        for index, node in enumerate(NODES)
    )

    proportional_scan_summary = "\n".join(
        "- Lambda `{lambda}`: mean final total `{mean_final_total}`, "
        "mean second-half drift `{mean_second_half_total_slope}`, "
        "classification `{simulation_classification}`.".format(**row)
        for row in proportional_scan_rows
    )
    fixed_scan_summary = "\n".join(
        "- Lambda `{lambda}`: mean final total `{mean_final_total}`, "
        "mean second-half drift `{mean_second_half_total_slope}`, "
        "classification `{simulation_classification}`, "
        "theoretical classification `{theoretical_classification}`.".format(**row)
        for row in fixed_scan_rows
    )

    new_section = f"""# Problem 3 - Open network

Solved on 2026-06-04.

### Methods and assumptions

- Work performed only inside `hw2/`.
- Node order: `['o', 'a', 'b', 'c', 'd']`.
- Matrix convention: `Lambda_open[i,j]` is the rate from source node `i` to destination node `j`.
- External arrivals enter node `o` according to a Poisson process with rate `lambda`.
- Service clock rates use `omega = Lambda_open @ 1`, with the assignment override `omega_d = 7/4`.
- Omega values used: `{omega_text}`.
- Event simulator supports service modes `proportional` and `fixed`.
- In proportional mode, node `i` has service rate `omega_i N_i(t)`.
- In fixed mode, node `i` has clock rate `omega_i`; empty-node clock ticks do nothing.
- For nodes `o,a,b,c`, a nonempty service event forwards one particle using normalized outgoing rates from `Lambda_open`.
- For node `d`, a nonempty service event removes one particle from the system.

### P3(a): proportional-rate scenario

- Main simulation seed: `{PROPORTIONAL_SEED}`.
- Main simulation lambda: `{fmt(PROPORTIONAL_MAIN_LAMBDA)}`.
- Main simulation horizon: `{fmt(PROPORTIONAL_MAIN_HORIZON)}`.
- Result file: `results/problem3_proportional_timeseries_lambda100.csv`.
- Figure: `figures/problem3_proportional_lambda100_counts.png`.
- Stability scan seed base: `{PROPORTIONAL_SCAN_SEED}`.
- Stability scan lambdas tested: {proportional_lambdas}.
- Stability scan horizon: `{fmt(PROPORTIONAL_SCAN_HORIZON)}`.
- Stability scan replicates per lambda: `{PROPORTIONAL_SCAN_REPLICATES}`.
- Blow-up criterion for scan: classify as blow-up if the mean fitted slope of total population over the second half of the horizon exceeds `{fmt(PROPORTIONAL_DRIFT_FRACTION)}` times the tested lambda.
- Result file: `results/problem3_proportional_stability_scan.csv`.
- Figure: `figures/problem3_proportional_stability_scan.png`.
- Theoretical motivation note: proportional service makes the open network a linear infinite-server network; the expected single-particle remaining lifetimes are `{absorption_text}`, so finite lambda gives finite expected population.
- Proportional scan summary:
{proportional_scan_summary}

### P3(b): fixed-rate scenario

- Main simulation seed: `{FIXED_SEED}`.
- Main simulation lambda: `{fmt(FIXED_MAIN_LAMBDA)}`.
- Main simulation horizon: `{fmt(FIXED_MAIN_HORIZON)}`.
- Result file: `results/problem3_fixed_timeseries_lambda2.csv`.
- Figure: `figures/problem3_fixed_lambda2_counts.png`.
- Stability scan seed base: `{FIXED_SCAN_SEED}`.
- Stability scan lambdas tested: {fixed_lambdas}.
- Stability scan horizon: `{fmt(FIXED_SCAN_HORIZON)}`.
- Stability scan replicates per lambda: `{FIXED_SCAN_REPLICATES}`.
- Blow-up criterion for scan: classify as blow-up if the mean fitted slope of total population over the second half of the horizon exceeds `{fmt(FIXED_DRIFT_THRESHOLD)}` or if the mean final total exceeds `{fmt(FIXED_FINAL_TOTAL_THRESHOLD)}`.
- Result file: `results/problem3_fixed_stability_scan.csv`.
- Figure: `figures/problem3_fixed_stability_scan.png`.
- Theoretical routing loads per unit lambda are saved in the run code as `[1, 1/2, 5/8, 3/4, 1]`.
- Fixed-rate capacity by node: `{fixed_capacity_text}`.
- Theoretical fixed-rate bottleneck node: `{NODES[fixed_bottleneck_index]}`.
- Theoretical largest stable input rate: `{fmt(fixed_capacity)}`.
- Fixed scan summary:
{fixed_scan_summary}

### Generated artifacts

- `results/problem3_proportional_timeseries_lambda100.csv`
- `results/problem3_proportional_stability_scan.csv`
- `results/problem3_fixed_timeseries_lambda2.csv`
- `results/problem3_fixed_stability_scan.csv`
- `figures/problem3_proportional_lambda100_counts.png`
- `figures/problem3_proportional_stability_scan.png`
- `figures/problem3_fixed_lambda2_counts.png`
- `figures/problem3_fixed_stability_scan.png`

### Report-writing notes

- These factual notes were used as source material for the final Problem 3 report prose.
- State explicitly that the proportional-rate scenario has no finite theoretical upper bound for finite input rates.
- State explicitly that the fixed-rate scenario has bottleneck capacity at node `{NODES[fixed_bottleneck_index]}`.

"""
    log_path.write_text(text[:start] + new_section, encoding="utf-8")


def main():
    omega, _, p_jump_cdf = open_network_components(LAMBDA_OPEN)
    fixed_capacity, fixed_bottleneck_index, fixed_capacity_by_node = (
        fixed_rate_capacity(omega)
    )
    absorption_times = mean_absorption_times(LAMBDA_OPEN, omega)

    proportional_sample_times = sample_grid(PROPORTIONAL_MAIN_HORIZON, 0.05)
    proportional_rng = np.random.default_rng(PROPORTIONAL_SEED)
    proportional_result = simulate_open_network(
        PROPORTIONAL_MAIN_LAMBDA,
        PROPORTIONAL_MAIN_HORIZON,
        proportional_sample_times,
        omega,
        p_jump_cdf,
        "proportional",
        proportional_rng,
    )
    write_time_series_csv(
        RESULTS_DIR / "problem3_proportional_timeseries_lambda100.csv",
        proportional_result,
        PROPORTIONAL_MAIN_LAMBDA,
        "proportional",
        PROPORTIONAL_SEED,
    )
    plot_counts(
        FIGURES_DIR / "problem3_proportional_lambda100_counts.png",
        proportional_result,
        "Problem 3 proportional service, lambda = 100",
    )

    fixed_sample_times = sample_grid(FIXED_MAIN_HORIZON, 2.0)
    fixed_rng = np.random.default_rng(FIXED_SEED)
    fixed_result = simulate_open_network(
        FIXED_MAIN_LAMBDA,
        FIXED_MAIN_HORIZON,
        fixed_sample_times,
        omega,
        p_jump_cdf,
        "fixed",
        fixed_rng,
    )
    write_time_series_csv(
        RESULTS_DIR / "problem3_fixed_timeseries_lambda2.csv",
        fixed_result,
        FIXED_MAIN_LAMBDA,
        "fixed",
        FIXED_SEED,
    )
    plot_counts(
        FIGURES_DIR / "problem3_fixed_lambda2_counts.png",
        fixed_result,
        "Problem 3 fixed service, lambda = 2",
    )

    proportional_scan_rows = run_scan(
        PROPORTIONAL_SCAN_LAMBDAS,
        "proportional",
        PROPORTIONAL_SCAN_HORIZON,
        2.0,
        PROPORTIONAL_SCAN_REPLICATES,
        PROPORTIONAL_SCAN_SEED,
        omega,
        p_jump_cdf,
        fixed_capacity,
    )
    write_rows_csv(
        RESULTS_DIR / "problem3_proportional_stability_scan.csv",
        proportional_scan_rows,
    )
    plot_scan(
        FIGURES_DIR / "problem3_proportional_stability_scan.png",
        proportional_scan_rows,
        "Proportional-rate stability scan",
    )

    fixed_scan_rows = run_scan(
        FIXED_SCAN_LAMBDAS,
        "fixed",
        FIXED_SCAN_HORIZON,
        20.0,
        FIXED_SCAN_REPLICATES,
        FIXED_SCAN_SEED,
        omega,
        p_jump_cdf,
        fixed_capacity,
    )
    write_rows_csv(
        RESULTS_DIR / "problem3_fixed_stability_scan.csv",
        fixed_scan_rows,
    )
    plot_scan(
        FIGURES_DIR / "problem3_fixed_stability_scan.png",
        fixed_scan_rows,
        "Fixed-rate stability scan",
    )

    update_results_log(
        omega,
        absorption_times,
        fixed_capacity,
        fixed_bottleneck_index,
        fixed_capacity_by_node,
        proportional_scan_rows,
        fixed_scan_rows,
    )

    print("Problem 3 complete.")
    print(f"Proportional scan lambdas: {PROPORTIONAL_SCAN_LAMBDAS}")
    print(f"Fixed-rate theoretical capacity: {fmt(fixed_capacity)}")
    print(f"Fixed-rate bottleneck node: {NODES[fixed_bottleneck_index]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
