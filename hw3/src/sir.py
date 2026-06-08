"""Discrete-time SIR epidemic engine for Homework 3.

State codes per node: 0 = S (susceptible), 1 = I (infected), 2 = R (recovered),
3 = V (vaccinated). One time step is one week.

Update convention (synchronous, from the start-of-week state X(t), as dictated by the
assignment's transition equations which both condition on X(t) to yield X(t+1)):

  * a susceptible node with 'm' infected neighbours becomes infected with
    probability '1 - (1 - beta) ** m';
  * an infected node recovers with probability 'rho';
  * both transitions are computed from the same start-of-week infected set and applied
    together, so a node infected this week cannot recover until the following week.

When a vaccination schedule is given, vaccination is applied at the *start* of each
week (before infection spreads): the cumulative number of vaccinated individuals is
brought up to 'round(frac * n)', with the newcomers chosen uniformly at random among
those not yet vaccinated (any S/I/R state). Vaccinated nodes can neither infect nor be
infected.

Weekly outputs are length 'weeks + 1' and indexed by week 0..weeks; week 0 is the
initial configuration, with 'newly_infected[0]' equal to the number of seeds.
"""

from __future__ import annotations

import numpy as np

S, I, R, V = 0, 1, 2, 3


def simulate_sir(adjacency, beta, rho, n_init_infected, weeks, rng, vacc_schedule=None):
    """Simulate one SIR run and return per-week time series.

    'vacc_schedule' (optional) is a sequence of length 'weeks' of cumulative
    vaccinated *percentages*; entry 'w-1' is the target applied at the start of week
    'w' (so week 0 carries no vaccination). Returns a dict of length 'weeks + 1'
    arrays: 'newly_infected', 'newly_vaccinated', 'S', 'I', 'R', 'V'.
    """
    n = adjacency.shape[0]
    state = np.zeros(n, dtype=np.int8)

    # Initial configuration: place the infected seeds at random.
    seeds = rng.choice(n, size=n_init_infected, replace=False)
    state[seeds] = I

    newly_infected = np.zeros(weeks + 1)
    newly_vaccinated = np.zeros(weeks + 1)
    totals = {label: np.zeros(weeks + 1) for label in ("S", "I", "R", "V")}

    newly_infected[0] = n_init_infected
    _record_totals(totals, state, 0)

    for week in range(1, weeks + 1):
        # 1) Vaccination at the start of the week (before infection spreads).
        if vacc_schedule is not None:
            target = int(round(vacc_schedule[week - 1] / 100.0 * n))
            already = int(np.count_nonzero(state == V))
            increment = max(0, target - already)
            if increment > 0:
                eligible = np.flatnonzero(state != V)
                increment = min(increment, eligible.size)
                chosen = rng.choice(eligible, size=increment, replace=False)
                state[chosen] = V
                newly_vaccinated[week] = increment

        # 2) Synchronous infection + recovery from the start-of-week infected set.
        infected_mask = (state == I).astype(np.float64)
        m = adjacency.dot(infected_mask)  # number of infected neighbours per node

        susceptible = state == S
        infection_prob = 1.0 - (1.0 - beta) ** m
        new_infections = susceptible & (rng.random(n) < infection_prob)

        infectious = state == I
        recoveries = infectious & (rng.random(n) < rho)

        state[new_infections] = I
        state[recoveries] = R

        newly_infected[week] = int(np.count_nonzero(new_infections))
        _record_totals(totals, state, week)

    return {
        "newly_infected": newly_infected,
        "newly_vaccinated": newly_vaccinated,
        "S": totals["S"],
        "I": totals["I"],
        "R": totals["R"],
        "V": totals["V"],
    }


def _record_totals(totals, state, week):
    """Record the S/I/R/V counts for the given week, asserting conservation."""
    counts = {
        "S": int(np.count_nonzero(state == S)),
        "I": int(np.count_nonzero(state == I)),
        "R": int(np.count_nonzero(state == R)),
        "V": int(np.count_nonzero(state == V)),
    }
    assert sum(counts.values()) == state.size, "S+I+R+V must equal the population"
    for label, value in counts.items():
        totals[label][week] = value


def run_many(adjacency, beta, rho, n_init_infected, weeks, n_runs, base_seed,
             vacc_schedule=None):
    """Average 'simulate_sir' over 'n_runs' independent runs.

    Each run uses 'np.random.default_rng(base_seed + run)' for reproducibility.
    Returns a dict with 'weeks' (the array 0..weeks) and, for every quantity,
    '<name>_mean' and '<name>_std' arrays.
    """
    quantities = ("newly_infected", "newly_vaccinated", "S", "I", "R", "V")
    collected = {name: np.zeros((n_runs, weeks + 1)) for name in quantities}

    for run in range(n_runs):
        rng = np.random.default_rng(base_seed + run)
        result = simulate_sir(adjacency, beta, rho, n_init_infected, weeks, rng,
                              vacc_schedule=vacc_schedule)
        for name in quantities:
            collected[name][run] = result[name]

    summary = {"weeks": np.arange(weeks + 1)}
    for name in quantities:
        summary[f"{name}_mean"] = collected[name].mean(axis=0)
        summary[f"{name}_std"] = collected[name].std(axis=0)
    return summary
