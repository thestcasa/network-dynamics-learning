"""Part 2: network coordination / anti-coordination game (n = 3).

Computes the pure Nash equilibria for n1 in {3, 2, 1, 0}; draws the asynchronous
best-response transition graph for n1 = 3 and n1 = 0; and reports the limiting
distribution conditioned on X(0) = (+1, -1, +1) for both the noiseless best-response
dynamics and the noisy (logit) best-response dynamics in the vanishing-noise limit.

Run from the hw3/ directory:

    python src/run_games.py
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402

from constants import FIGURES_DIR, RESULTS_DIR  # noqa: E402
from games import (all_profiles, best_response_transition, limit_distribution_from_initial,  # noqa: E402
                   logit_transition, nash_equilibria, profile_index_map,
                   stationary_distribution, utility, vanishing_noise_limit)
from plotting import save_figure  # noqa: E402
from utils import format_float, write_rows_csv  # noqa: E402

N = 3
X0 = (1, -1, 1)  # initial configuration for the limit distributions
EPSILONS = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]  # vanishing-noise sequence

# Lattice positions for the 8 states: x = number of +1 actions, y spread within a level.
LAYOUT_Y = {0: [0.0], 1: [-1.0, 0.0, 1.0], 2: [-1.0, 0.0, 1.0], 3: [0.0]}


def label(profile):
    """Compact string label for a profile, e.g. (+1,-1,+1) -> '+-+'."""
    return "".join("+" if action == 1 else "-" for action in profile)


def layout_positions():
    """Return a dict profile -> (x, y) arranged by number of +1 actions."""
    positions = {}
    buckets = {0: [], 1: [], 2: [], 3: []}
    for profile in all_profiles(N):
        buckets[sum(1 for a in profile if a == 1)].append(profile)
    for level, profiles in buckets.items():
        for profile, y in zip(profiles, LAYOUT_Y[level]):
            positions[profile] = (float(level), y)
    return positions


def write_nash_tables():
    """Write per-profile utilities + Nash flags for each n1, and a summary of NE sets."""
    profiles = all_profiles(N)
    summary_rows = []
    for n1 in (3, 2, 1, 0):
        equilibria = set(nash_equilibria(N, n1))
        rows = []
        for profile in profiles:
            row = {"profile": label(profile)}
            for player in range(N):
                row[f"u{player + 1}"] = format_float(utility(player, profile, n1))
            row["is_nash"] = int(profile in equilibria)
            rows.append(row)
        write_rows_csv(RESULTS_DIR / f"games_nash_n1_{n1}.csv", rows)
        ne_labels = [label(p) for p in profiles if p in equilibria]
        summary_rows.append({"n1": n1, "num_nash": len(ne_labels),
                             "nash_equilibria": " ".join(ne_labels)})
        print(f"[games] n1={n1}: {len(ne_labels)} pure NE -> {ne_labels}")
    write_rows_csv(RESULTS_DIR / "games_nash_summary.csv", summary_rows)


def draw_transition_graph(n1, path):
    """Draw the asynchronous best-response transition graph for the given n1."""
    transition = best_response_transition(N, n1)
    index = profile_index_map(N)
    profiles = all_profiles(N)
    equilibria = set(nash_equilibria(N, n1))
    positions = layout_positions()

    graph = nx.DiGraph()
    for profile in profiles:
        graph.add_node(label(profile))
    self_loop_prob = {}
    edge_labels = {}
    for source in profiles:
        for target in profiles:
            prob = transition[index[source], index[target]]
            if prob <= 0:
                continue
            if source == target:
                self_loop_prob[label(source)] = prob
            else:
                graph.add_edge(label(source), label(target))
                edge_labels[(label(source), label(target))] = _fraction(prob)

    pos = {label(p): positions[p] for p in profiles}
    node_colors = ["#d96459" if p in equilibria else "#9bc1d4" for p in profiles]

    fig, ax = plt.subplots(figsize=(9, 6))
    nx.draw_networkx_nodes(graph, pos, node_color=node_colors, node_size=1500,
                           edgecolors="black", ax=ax)
    nx.draw_networkx_labels(graph, pos, font_size=12, font_family="monospace", ax=ax)
    nx.draw_networkx_edges(graph, pos, ax=ax, node_size=1500, arrowsize=18,
                           connectionstyle="arc3,rad=0.12", min_target_margin=18)
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, ax=ax,
                                 font_size=9, label_pos=0.5,
                                 connectionstyle="arc3,rad=0.12")
    # Annotate self-loop (staying) probabilities next to absorbing / indifferent states.
    for node, prob in self_loop_prob.items():
        x, y = pos[node]
        ax.annotate(f"stay {_fraction(prob)}", (x, y), xytext=(0, 22),
                    textcoords="offset points", ha="center", fontsize=8, color="#444444")
    ax.set_title(f"Asynchronous best-response transition graph (n={N}, n1={n1})\n"
                 f"red = Nash equilibrium (absorbing); edge labels are probabilities")
    ax.axis("off")
    save_figure(fig, path)


def _fraction(prob):
    """Render a probability that is a multiple of 1/3 as a compact fraction."""
    numerator = int(round(prob * N))
    if numerator == N:
        return "1"
    return f"{numerator}/{N}"


def limit_distributions(n1):
    """Return (br_limit, noisy_limit, epsilon_table) conditioned on X(0) for this n1."""
    index = profile_index_map(N)
    x0_index = index[X0]

    br_transition = best_response_transition(N, n1)
    br_limit = limit_distribution_from_initial(br_transition, x0_index)

    # Noisy (logit) chain is ergodic for any eps > 0, so its limit is the stationary
    # distribution, independent of X(0). We report the exact vanishing-noise limit
    # (uniform over the potential maximizers) and track the numerical logit stationary
    # as eps -> 0 to confirm it.
    epsilon_table = [(epsilon, stationary_distribution(logit_transition(N, n1, epsilon)))
                     for epsilon in EPSILONS]
    noisy_limit = vanishing_noise_limit(N, n1)
    return br_limit, noisy_limit, epsilon_table


def write_and_plot_limits():
    """Compute, save, and plot the limit distributions for n1 = 3 and n1 = 0."""
    profiles = all_profiles(N)
    rows = []
    results = {}
    for n1 in (3, 0):
        br_limit, noisy_limit, epsilon_table = limit_distributions(n1)
        results[n1] = (br_limit, noisy_limit)
        for profile in profiles:
            i = profiles.index(profile)
            rows.append({
                "n1": n1, "state": label(profile),
                "br_limit": format_float(br_limit[i]),
                "noisy_vanishing_noise_limit": format_float(noisy_limit[i]),
            })
        # Console summary
        br_support = {label(p): round(br_limit[profiles.index(p)], 3)
                      for p in profiles if br_limit[profiles.index(p)] > 1e-6}
        noisy_support = {label(p): round(noisy_limit[profiles.index(p)], 3)
                         for p in profiles if noisy_limit[profiles.index(p)] > 1e-6}
        print(f"[games] n1={n1} | BR limit from {label(X0)}: {br_support}")
        print(f"[games] n1={n1} | noisy vanishing-noise limit: {noisy_support}")
        _plot_epsilon_convergence(n1, epsilon_table)

    write_rows_csv(RESULTS_DIR / "games_limit_distributions.csv", rows,
                   fieldnames=["n1", "state", "br_limit",
                               "noisy_vanishing_noise_limit"])
    _plot_limit_bars(results)


def _plot_limit_bars(results):
    """Grouped bar chart of BR vs noisy limit distributions for n1 = 3 and n1 = 0."""
    profiles = all_profiles(N)
    labels = [label(p) for p in profiles]
    x = np.arange(len(labels))
    width = 0.38
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, n1 in zip(axes, (3, 0)):
        br_limit, noisy_limit = results[n1]
        ax.bar(x - width / 2, br_limit, width, label="noiseless BR", color="#d96459")
        ax.bar(x + width / 2, noisy_limit, width, label="noisy BR (eps->0)",
               color="#5b8c9b")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontfamily="monospace")
        ax.set_xlabel("state")
        ax.set_title(f"n1 = {n1}")
        ax.grid(True, alpha=0.25, axis="y")
    axes[0].set_ylabel(f"limit probability  (X(0) = {label(X0)})")
    axes[0].legend(loc="upper center")
    save_figure(fig, FIGURES_DIR / "games_limit_distributions.png")


def _plot_epsilon_convergence(n1, epsilon_table):
    """Plot the logit stationary distribution as epsilon -> 0 (vanishing noise)."""
    profiles = all_profiles(N)
    epsilons = [e for e, _ in epsilon_table]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for i, profile in enumerate(profiles):
        values = [stationary[i] for _, stationary in epsilon_table]
        if max(values) > 1e-3:
            ax.plot(epsilons, values, marker="o", markersize=4, label=label(profile))
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("noise level epsilon (log scale, decreasing ->)")
    ax.set_ylabel("stationary probability")
    ax.set_title(f"Logit stationary distribution as epsilon -> 0 (n1 = {n1})")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    save_figure(fig, FIGURES_DIR / f"games_logit_convergence_n1_{n1}.png")


def main():
    write_nash_tables()
    draw_transition_graph(3, FIGURES_DIR / "games_transition_n1_3.png")
    draw_transition_graph(0, FIGURES_DIR / "games_transition_n1_0.png")
    write_and_plot_limits()


if __name__ == "__main__":
    main()
