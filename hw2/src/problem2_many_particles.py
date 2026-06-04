"""Problem 2 many-particle continuous-time random walk routines.

Matrix convention: Lambda[i, j] is the transition rate from source node i to
destination node j. The simulations below use direct CTMC events. The matrix
P_jump is therefore the embedded jump-chain transition matrix after a particle
leaves a node, not a uniformized transition matrix with self-loops.
"""

from __future__ import annotations

import numpy as np


def ctmc_components(lambda_matrix):
    """Return omega, Q, and the embedded jump-chain matrix P_jump."""
    lambda_matrix = np.asarray(lambda_matrix, dtype=float)
    omega = lambda_matrix @ np.ones(lambda_matrix.shape[0])
    q_matrix = lambda_matrix - np.diag(omega)
    p_jump = np.zeros_like(lambda_matrix, dtype=float)
    positive = omega > 0
    p_jump[positive, :] = lambda_matrix[positive, :] / omega[positive, None]
    return omega, q_matrix, p_jump


def stationary_distribution(q_matrix):
    """Solve pi Q = 0 with sum_i pi_i = 1."""
    q_matrix = np.asarray(q_matrix, dtype=float)
    n_nodes = q_matrix.shape[0]
    system = q_matrix.T.copy()
    rhs = np.zeros(n_nodes, dtype=float)
    system[-1, :] = 1.0
    rhs[-1] = 1.0
    pi = np.linalg.solve(system, rhs)
    return pi / np.sum(pi)


def simulate_return_time(start_index, omega, p_jump_cdf, rng):
    """Simulate T_start^+, including the initial holding time at start_index."""
    current = start_index
    total_time = 0.0
    has_left_start = False

    while True:
        total_time += rng.exponential(1.0 / omega[current])
        next_node = int(np.searchsorted(p_jump_cdf[current], rng.random(), side="right"))
        if current == start_index and next_node != start_index:
            has_left_start = True
        current = next_node
        if has_left_start and current == start_index:
            return total_time


def simulate_particle_return_run(num_particles, start_index, omega, p_jump_cdf, rng):
    """Return the average return time across num_particles independent particles."""
    values = [
        simulate_return_time(start_index, omega, p_jump_cdf, rng)
        for _ in range(num_particles)
    ]
    return float(np.mean(values))


def simulate_particle_return_runs(
    num_particles, num_runs, start_index, omega, p_jump, seed
):
    """Simulate repeated batches of independent particle return times."""
    rng = np.random.default_rng(seed)
    p_jump_cdf = np.cumsum(p_jump, axis=1)
    p_jump_cdf[:, -1] = 1.0
    run_means = np.empty(num_runs, dtype=float)
    for run_index in range(num_runs):
        run_means[run_index] = simulate_particle_return_run(
            num_particles, start_index, omega, p_jump_cdf, rng
        )
    return run_means


def simulate_node_count_run(num_particles, horizon, sample_times, omega, p_jump_cdf, rng):
    """Simulate the many-particle node-count CTMC on a fixed sampling grid."""
    n_nodes = len(omega)
    counts = np.zeros(n_nodes, dtype=int)
    counts[1] = num_particles
    sampled_counts = np.empty((len(sample_times), n_nodes), dtype=float)

    time = 0.0
    sample_index = 0
    while sample_index < len(sample_times) and sample_times[sample_index] <= time:
        sampled_counts[sample_index] = counts
        sample_index += 1

    while time < horizon:
        node_rates = counts * omega
        total_rate = float(np.sum(node_rates))
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

        source_threshold = rng.random() * total_rate
        source = int(np.searchsorted(np.cumsum(node_rates), source_threshold, side="right"))
        destination = int(np.searchsorted(p_jump_cdf[source], rng.random(), side="right"))
        counts[source] -= 1
        counts[destination] += 1
        time = next_time

    while sample_index < len(sample_times):
        sampled_counts[sample_index] = counts
        sample_index += 1

    return sampled_counts


def simulate_node_count_runs(
    num_particles, horizon, sample_times, num_runs, omega, p_jump, seed
):
    """Return sample mean, sample std, and final counts for node-count runs."""
    rng = np.random.default_rng(seed)
    p_jump_cdf = np.cumsum(p_jump, axis=1)
    p_jump_cdf[:, -1] = 1.0
    n_times = len(sample_times)
    n_nodes = len(omega)
    count_sum = np.zeros((n_times, n_nodes), dtype=float)
    count_sumsq = np.zeros((n_times, n_nodes), dtype=float)
    final_counts = np.empty((num_runs, n_nodes), dtype=float)

    for run_index in range(num_runs):
        sampled_counts = simulate_node_count_run(
            num_particles, horizon, sample_times, omega, p_jump_cdf, rng
        )
        count_sum += sampled_counts
        count_sumsq += sampled_counts * sampled_counts
        final_counts[run_index] = sampled_counts[-1]

    mean = count_sum / num_runs
    if num_runs > 1:
        variance = (count_sumsq - num_runs * mean * mean) / (num_runs - 1)
        std = np.sqrt(np.maximum(variance, 0.0))
    else:
        std = np.zeros_like(mean)
    return mean, std, final_counts
