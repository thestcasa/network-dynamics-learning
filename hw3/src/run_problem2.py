"""Problem 2: SIR epidemic on a preferential-attachment graph, no vaccination.

Run from the hw3/ directory:

    python src/run_problem2.py
"""

from __future__ import annotations

import numpy as np

from constants import (
    BETA,
    DEFAULT_SEED,
    FIGURES_DIR,
    K_PA,
    N_INIT_INFECTED,
    N_PA,
    N_RUNS,
    RESULTS_DIR,
    RHO,
    WEEKS,
)
from graphs import degree_stats, preferential_attachment
from plotting import plot_lines
from sir import run_many
from utils import write_weekly_csv


def main():
    rng = np.random.default_rng(DEFAULT_SEED + 2100)
    adjacency = preferential_attachment(N_PA, K_PA, rng)
    mean_degree, degrees = degree_stats(adjacency)
    print(
        f"[2] PA graph: n={N_PA}, target k={K_PA}, "
        f"realized mean degree={mean_degree:.3f}, max degree={degrees.max()}"
    )

    summary = run_many(
        adjacency,
        BETA,
        RHO,
        N_INIT_INFECTED,
        WEEKS,
        N_RUNS,
        base_seed=DEFAULT_SEED + 2000,
    )
    weeks = summary["weeks"]

    write_weekly_csv(
        RESULTS_DIR / "problem2_newly_infected.csv",
        weeks,
        {
            "newly_infected_mean": summary["newly_infected_mean"],
            "newly_infected_std": summary["newly_infected_std"],
        },
    )
    write_weekly_csv(
        RESULTS_DIR / "problem2_sir_totals.csv",
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
        f"Problem 2: newly infected per week (PA graph, k={K_PA}, N={N_RUNS})",
        FIGURES_DIR / "problem2_newly_infected.png",
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
        f"Problem 2: S/I/R totals per week (PA graph, N={N_RUNS})",
        FIGURES_DIR / "problem2_sir_totals.png",
    )
    peak_week = int(np.argmax(summary["newly_infected_mean"]))
    print(
        f"[2] peak newly-infected at week {peak_week} "
        f"(~{summary['newly_infected_mean'][peak_week]:.1f}); "
        f"final recovered ~{summary['R_mean'][-1]:.1f}"
    )


if __name__ == "__main__":
    main()
