"""Problem 3 open-network event simulator.

Matrix convention: Lambda_open[i, j] is the rate from source node i to
destination node j. External arrivals enter node o. The assignment gives node d
a service clock of 7/4 even though row d of Lambda_open is zero.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


SERVICE_MODES = {"proportional", "fixed"}


@dataclass(frozen=True)
class SimulationResult:
    """Sampled output from one open-network simulation run."""

    sample_times: np.ndarray
    counts: np.ndarray
    event_count: int


def open_network_components(lambda_open):
    """Return omega and embedded destination probabilities for the open network."""
    lambda_open = np.asarray(lambda_open, dtype=float)
    omega = lambda_open @ np.ones(lambda_open.shape[0])
    omega = omega.astype(float)
    omega[-1] = 7.0 / 4.0

    p_jump = np.zeros_like(lambda_open, dtype=float)
    row_sums = lambda_open @ np.ones(lambda_open.shape[0])
    positive = row_sums > 0.0
    p_jump[positive, :] = lambda_open[positive, :] / row_sums[positive, None]
    p_jump_cdf = np.cumsum(p_jump, axis=1)
    p_jump_cdf[positive, -1] = 1.0
    return omega, p_jump, p_jump_cdf


def service_rates(counts, omega, service_mode):
    """Return node clock rates under the selected service rule."""
    if service_mode == "proportional":
        return omega * counts
    if service_mode == "fixed":
        return omega.copy()
    raise ValueError(f"Unknown service mode: {service_mode!r}")


def simulate_open_network(
    arrival_rate,
    horizon,
    sample_times,
    omega,
    p_jump_cdf,
    service_mode,
    rng,
):
    """Simulate one event-based open-network trajectory on a sampling grid.

    Events are external arrivals to node o and node clock ticks. For nodes
    o,a,b,c, a nonempty service event forwards one particle according to the
    normalized outgoing rates. For node d, a nonempty service event removes one
    particle from the system.
    """
    if service_mode not in SERVICE_MODES:
        raise ValueError(f"service_mode must be one of {sorted(SERVICE_MODES)}")
    if arrival_rate < 0.0:
        raise ValueError("arrival_rate must be nonnegative")

    n_nodes = len(omega)
    destination_node = n_nodes - 1
    counts = np.zeros(n_nodes, dtype=int)
    sampled_counts = np.empty((len(sample_times), n_nodes), dtype=float)

    time = 0.0
    event_count = 0
    sample_index = 0
    while sample_index < len(sample_times) and sample_times[sample_index] <= time:
        sampled_counts[sample_index] = counts
        sample_index += 1

    while time < horizon:
        node_rates = service_rates(counts, omega, service_mode)
        total_rate = float(arrival_rate + np.sum(node_rates))
        if total_rate <= 0.0:
            break

        next_time = time + rng.exponential(1.0 / total_rate)
        while (
            sample_index < len(sample_times)
            and sample_times[sample_index] <= next_time
        ):
            sampled_counts[sample_index] = counts
            sample_index += 1

        if next_time > horizon:
            time = next_time
            break

        event_threshold = rng.random() * total_rate
        if event_threshold < arrival_rate:
            counts[0] += 1
        else:
            service_threshold = event_threshold - arrival_rate
            service_cdf = np.cumsum(node_rates)
            source = int(np.searchsorted(service_cdf, service_threshold, side="right"))
            if counts[source] > 0:
                counts[source] -= 1
                if source != destination_node:
                    destination = int(
                        np.searchsorted(
                            p_jump_cdf[source], rng.random(), side="right"
                        )
                    )
                    counts[destination] += 1
        event_count += 1
        time = next_time

    while sample_index < len(sample_times):
        sampled_counts[sample_index] = counts
        sample_index += 1

    return SimulationResult(sample_times, sampled_counts, event_count)


def simulate_replicates(
    arrival_rate,
    horizon,
    sample_times,
    omega,
    p_jump_cdf,
    service_mode,
    num_replicates,
    seed,
):
    """Run independent replicates and return the list of sampled trajectories."""
    rng = np.random.default_rng(seed)
    return [
        simulate_open_network(
            arrival_rate,
            horizon,
            sample_times,
            omega,
            p_jump_cdf,
            service_mode,
            rng,
        )
        for _ in range(num_replicates)
    ]


def total_population_slope(sample_times, counts, start_fraction=0.5):
    """Fit a line to total population over the final part of a trajectory."""
    totals = np.sum(counts, axis=1)
    start_index = int(np.floor(len(sample_times) * start_fraction))
    x_values = sample_times[start_index:]
    y_values = totals[start_index:]
    if len(x_values) < 2:
        return 0.0
    slope, _ = np.polyfit(x_values, y_values, deg=1)
    return float(slope)


def theoretical_loads(lambda_value):
    """Return fixed-rate traffic loads implied by the acyclic routing graph."""
    return np.array(
        [
            lambda_value,
            lambda_value / 2.0,
            5.0 * lambda_value / 8.0,
            3.0 * lambda_value / 4.0,
            lambda_value,
        ],
        dtype=float,
    )


def fixed_rate_capacity(omega):
    """Return the bottleneck capacity for the fixed-rate network."""
    load_per_unit_lambda = theoretical_loads(1.0)
    ratios = omega / load_per_unit_lambda
    return float(np.min(ratios)), int(np.argmin(ratios)), ratios


def mean_absorption_times(lambda_open, omega):
    """Expected remaining lifetime of one particle from each node.

    This is used for proportional-rate stability reasoning. It solves
    tau_i = 1/omega_i + sum_j P_ij tau_j for i != d and tau_d = 1/omega_d.
    """
    _, p_jump, _ = open_network_components(lambda_open)
    n_nodes = len(omega)
    system = np.eye(n_nodes)
    rhs = np.zeros(n_nodes, dtype=float)
    for i in range(n_nodes):
        rhs[i] = 1.0 / omega[i]
        if i != n_nodes - 1:
            system[i, :] -= p_jump[i, :]
    return np.linalg.solve(system, rhs)
