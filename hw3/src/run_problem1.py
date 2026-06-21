"""Problem 1.1 (epidemic on a known k-regular graph) and 1.2 (PA graph validation).

Run from the hw3/ directory:

    python src/run_problem1.py
"""

from __future__ import annotations

import numpy as np

from constants import (
    BETA,
    DEFAULT_SEED,
    FIGURES_DIR,
    K_KREGULAR,
    K_PA,
    N_INIT_INFECTED,
    N_KREGULAR,
    N_PA,
    N_RUNS,
    RESULTS_DIR,
    RHO,
    WEEKS,
)
from graphs import degree_stats, k_regular, preferential_attachment
from plotting import plot_histogram, plot_lines
from sir import run_many
from utils import format_float, write_single_row_csv, write_weekly_csv


def run_kregular_epidemic():
    """Problem 1.1: SIR on the symmetric k-regular graph, averaged over N runs."""
    adjacency = k_regular(N_KREGULAR, K_KREGULAR)
    mean_degree, degrees = degree_stats(adjacency)
    assert np.all(degrees == K_KREGULAR), "k-regular graph must have constant degree k"
    print(
        f"[1.1] k-regular graph: n={N_KREGULAR}, k={K_KREGULAR}, "
        f"mean degree={mean_degree:.3f}"
    )

    summary = run_many(
        adjacency,
        BETA,
        RHO,
        N_INIT_INFECTED,
        WEEKS,
        N_RUNS,
        base_seed=DEFAULT_SEED + 1000,
    )
    weeks = summary["weeks"]

    write_weekly_csv(
        RESULTS_DIR / "problem1_kregular_newly_infected.csv",
        weeks,
        {
            "newly_infected_mean": summary["newly_infected_mean"],
            "newly_infected_std": summary["newly_infected_std"],
        },
    )
    write_weekly_csv(
        RESULTS_DIR / "problem1_kregular_sir_totals.csv",
        weeks,
        {
            "S_mean": summary["S_mean"],
            "I_mean": summary["I_mean"],
            "R_mean": summary["R_mean"],
            "S_std": summary["S_std"],
            "I_std": summary["I_std"],
            "R_std": summary["R_std"],
        },
    )

    plot_lines(
        weeks,
        {"newly infected": summary["newly_infected_mean"]},
        "week",
        "average number of newly infected",
        f"Problem 1.1: newly infected per week (k-regular, k={K_KREGULAR}, N={N_RUNS})",
        FIGURES_DIR / "problem1_kregular_newly_infected.png",
        bands={"newly infected": summary["newly_infected_std"]},
    )
    plot_lines(
        weeks,
        {
            "susceptible": summary["S_mean"],
            "infected": summary["I_mean"],
            "recovered": summary["R_mean"],
        },
        "week",
        "average number of individuals",
        f"Problem 1.1: S/I/R totals per week (k-regular, N={N_RUNS})",
        FIGURES_DIR / "problem1_kregular_sir_totals.png",
    )
    peak_week = int(np.argmax(summary["newly_infected_mean"]))
    print(
        f"[1.1] peak newly-infected at week {peak_week} "
        f"(~{summary['newly_infected_mean'][peak_week]:.1f}); "
        f"final recovered ~{summary['R_mean'][-1]:.1f}"
    )


def validate_pa_graph():
    """Problem 1.2: build a preferential-attachment graph and check its average degree."""
    rng = np.random.default_rng(DEFAULT_SEED + 1100)
    adjacency = preferential_attachment(N_PA, K_PA, rng)
    mean_degree, degrees = degree_stats(adjacency)
    print(
        f"[1.2] PA graph: n={N_PA}, target k={K_PA}, "
        f"realized mean degree={mean_degree:.3f}, max degree={degrees.max()}"
    )

    write_single_row_csv(
        RESULTS_DIR / "problem1_pa_degree_summary.csv",
        {
            "n": N_PA,
            "target_k": K_PA,
            "realized_mean_degree": format_float(mean_degree),
            "min_degree": int(degrees.min()),
            "max_degree": int(degrees.max()),
            "seed": DEFAULT_SEED + 1100,
        },
    )
    plot_histogram(
        degrees, "degree", "number of nodes",
        f"Problem 1.2: degree distribution of the PA graph (n={N_PA}, k={K_PA})",
        FIGURES_DIR / "problem1_pa_degree_distribution.png",
        bins=range(int(degrees.min()), int(degrees.max()) + 2),
    )


def main():
    run_kregular_epidemic()
    validate_pa_graph()


if __name__ == "__main__":
    main()
