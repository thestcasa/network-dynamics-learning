"""Problem 4: estimate (k, beta, rho) for the H1N1 Sweden 2009 pandemic.

Coordinate (local) search over a 3x3x3 neighbourhood in (k, beta, rho), minimizing the
RMSE between the simulated average newly-infected curve and the real scaled data I0(t).
Each evaluation builds a preferential-attachment graph (cached per k, since the graph
depends only on k) and simulates the vaccinated epidemic of Section 1.3.

The search is run twice, side by side: once at N_RUNS_SEARCH (=100; the single-seed RMSE
objective is noisy, so we test 100 for a more stable, result; see constants.py)
and once at the assignment-specified N_RUNS_SEARCH_SPEC (=10), whose artifacts are
suffixed "_n10". Both tracks share the same base seeds, so the N=10 runs are literally
the first 10 of the N=100 runs.

Run from the hw3/ directory:

    python src/run_problem4.py
"""

from __future__ import annotations

import numpy as np

from constants import (
    BETA0_SEARCH,
    DEFAULT_SEED,
    DELTA_BETA,
    DELTA_K,
    DELTA_RHO,
    FIGURES_DIR,
    I0_SWEDEN,
    K0_SEARCH,
    N_RUNS_SEARCH,
    N_RUNS_SEARCH_SPEC,
    N_SWEDEN,
    RESULTS_DIR,
    RHO0_SEARCH,
    VACC_SWEDEN,
    WEEKS,
)
from graphs import preferential_attachment
from plotting import plot_lines
from sir import run_many
from utils import format_float, write_rows_csv, write_single_row_csv, write_weekly_csv

I0 = np.array(I0_SWEDEN, dtype=float)
N_INIT_SWEDEN = int(I0_SWEDEN[0])  # initial infected = first entry of I0(t) = 1
MAX_REFINEMENTS = 2  # halve (d_beta, d_rho) at most twice after convergence


def get_graph(k, graph_cache):
    """Return a cached PA graph for average degree k, built once per k value."""
    if k not in graph_cache:
        rng = np.random.default_rng(DEFAULT_SEED + 4000 + k)
        graph_cache[k] = preferential_attachment(N_SWEDEN, k, rng)
    return graph_cache[k]


def rmse_to_real(newly_infected_mean):
    """RMSE over weeks 1..15 between simulated I(t) and the real I0(t)."""
    simulated = newly_infected_mean[1:WEEKS + 1]
    target = I0[1:WEEKS + 1]
    return float(np.sqrt(np.mean((simulated - target) ** 2)))


def evaluate(k, beta, rho, graph_cache, rmse_cache, base_seed, n_runs):
    """Evaluate one (k, beta, rho); memoized. Returns (rmse, newly_infected_mean)."""
    key = (int(k), round(float(beta), 4), round(float(rho), 4))
    if key in rmse_cache:
        return rmse_cache[key]
    adjacency = get_graph(int(k), graph_cache)
    summary = run_many(
        adjacency,
        beta,
        rho,
        N_INIT_SWEDEN,
        WEEKS,
        n_runs,
        base_seed=base_seed,
        vacc_schedule=VACC_SWEDEN,
    )
    rmse = rmse_to_real(summary["newly_infected_mean"])
    rmse_cache[key] = (rmse, summary["newly_infected_mean"])
    return rmse_cache[key]


def neighbourhood(center, deltas):
    """Yield the 3x3x3 (k, beta, rho) grid around center with clamping."""
    kc, bc, rc = center
    dk, db, dr = deltas
    k_values = sorted({max(2, int(kc + s * dk)) for s in (-1, 0, 1)})
    beta_values = sorted(
        {min(0.95, max(0.05, round(bc + s * db, 4))) for s in (-1, 0, 1)}
    )
    rho_values = sorted({min(0.95, max(0.05, round(rc + s * dr, 4))) for s in (-1, 0, 1)})
    for k in k_values:
        for beta in beta_values:
            for rho in rho_values:
                yield (k, beta, rho)


def coordinate_search(n_runs):
    """Run the coordinate search and return the best params plus the full log."""
    graph_cache = {}
    rmse_cache = {}
    log_rows = []
    base_seed = DEFAULT_SEED + 4100

    center = (K0_SEARCH, round(BETA0_SEARCH, 4), round(RHO0_SEARCH, 4))
    deltas = [DELTA_K, DELTA_BETA, DELTA_RHO]
    iteration = 0
    refinements = 0

    while True:
        # Coordinate descent to convergence at the current resolution.
        while True:
            iteration += 1
            best_point = None
            best_rmse = np.inf
            for point in neighbourhood(center, deltas):
                rmse, _ = evaluate(*point, graph_cache, rmse_cache, base_seed, n_runs)
                log_rows.append(
                    {
                        "iteration": iteration,
                        "k": point[0],
                        "beta": format_float(point[1]),
                        "rho": format_float(point[2]),
                        "rmse": format_float(rmse),
                    }
                )
                if rmse < best_rmse:
                    best_rmse = rmse
                    best_point = point
            print(
                f"[4] iter {iteration:>2} (deltas={tuple(deltas)}): "
                f"best={best_point} RMSE={best_rmse:.4f}"
            )
            if best_point == center:
                break
            center = best_point

        if refinements >= MAX_REFINEMENTS:
            break
        deltas = [deltas[0], round(deltas[1] / 2, 4), round(deltas[2] / 2, 4)]
        refinements += 1
        print(f"[4] refining: halving beta/rho steps to {deltas[1]}, {deltas[2]}")

    best_rmse, best_curve = evaluate(*center, graph_cache, rmse_cache, base_seed, n_runs)
    return center, best_rmse, best_curve, graph_cache, log_rows


def run_track(n_runs_search, n_runs_final, suffix):
    """Run one search-and-fit track and write its artifacts.

    'n_runs_search' is the Monte Carlo sample size used during the coordinate search,
    'n_runs_final' the size of the final presentation run, and 'suffix' is appended to
    every output filename (empty for the N=100 headline track, "_n10" for the
    assignment-specified track). Returns the best (k, beta, rho), the search RMSE, and
    the final RMSE.
    """
    center, best_rmse, _, graph_cache, log_rows = coordinate_search(n_runs_search)
    best_k, best_beta, best_rho = center
    print(
        f"[4] BEST (N={n_runs_search}): k={best_k}, beta={best_beta}, rho={best_rho}, "
        f"search RMSE={best_rmse:.4f}"
    )

    write_rows_csv(RESULTS_DIR / f"problem4_search_log{suffix}.csv", log_rows)

    # Final presentation run with the best parameters.
    adjacency = graph_cache[best_k]
    summary = run_many(
        adjacency,
        best_beta,
        best_rho,
        N_INIT_SWEDEN,
        WEEKS,
        n_runs_final,
        base_seed=DEFAULT_SEED + 4200,
        vacc_schedule=VACC_SWEDEN,
    )
    weeks = summary["weeks"]
    final_rmse = rmse_to_real(summary["newly_infected_mean"])

    write_single_row_csv(
        RESULTS_DIR / f"problem4_best_params{suffix}.csv",
        {
            "n": N_SWEDEN,
            "best_k": best_k,
            "best_beta": format_float(best_beta),
            "best_rho": format_float(best_rho),
            "n_runs_search": n_runs_search,
            "search_rmse": format_float(best_rmse),
            "n_runs_final": n_runs_final,
            "final_rmse": format_float(final_rmse),
            "n_init_infected": N_INIT_SWEDEN,
            "seed": DEFAULT_SEED + 4200,
        },
    )
    write_weekly_csv(
        RESULTS_DIR / f"problem4_fit_vs_real{suffix}.csv",
        weeks,
        {
            "model_newly_infected_mean": summary["newly_infected_mean"],
            "real_newly_infected": I0,
        },
    )
    write_weekly_csv(
        RESULTS_DIR / f"problem4_sirv_totals{suffix}.csv",
        weeks,
        {
            "S_mean": summary["S_mean"],
            "I_mean": summary["I_mean"],
            "R_mean": summary["R_mean"],
            "V_mean": summary["V_mean"],
        },
    )

    plot_lines(
        weeks,
        {
            "model (best params)": summary["newly_infected_mean"],
            "real data $I_0(t)$": I0,
        },
        "week",
        "average number of newly infected",
        f"Problem 4: fitted vs real newly infected "
        f"(k={best_k}, $\\beta$={best_beta}, $\\rho$={best_rho}, N={n_runs_search})",
        FIGURES_DIR / f"problem4_fit_vs_real{suffix}.png",
    )
    plot_lines(
        weeks,
        {
            "susceptible": summary["S_mean"],
            "infected": summary["I_mean"],
            "recovered": summary["R_mean"],
            "vaccinated": summary["V_mean"],
        },
        "week",
        "average number of individuals",
        f"Problem 4: S/I/R/V totals under the best model (N={n_runs_final})",
        FIGURES_DIR / f"problem4_sirv_totals{suffix}.png",
    )
    print(f"[4] final RMSE (N={n_runs_final}) at best params = {final_rmse:.4f}")
    return center, best_rmse, final_rmse


def main():
    # Headline track: stable N=100 search (artifacts keep their original names).
    run_track(N_RUNS_SEARCH, 100, "")
    # Assignment-specified track: search and run at N=10, artifacts suffixed "_n10".
    print(f"\n[4] re-running the search at the assignment-specified N={N_RUNS_SEARCH_SPEC}")
    run_track(N_RUNS_SEARCH_SPEC, N_RUNS_SEARCH_SPEC, "_n10")


if __name__ == "__main__":
    main()
