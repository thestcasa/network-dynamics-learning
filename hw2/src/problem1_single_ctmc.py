"""Problem 1 single-particle continuous-time Markov chain routines."""

from __future__ import annotations

import numpy as np


def exit_rates(lambda_matrix):
    """Return omega_i = sum_j Lambda[i, j]."""
    return np.asarray(lambda_matrix, dtype=float).sum(axis=1)


def generator(lambda_matrix):
    """Return the CTMC generator Q = Lambda - diag(Lambda @ 1)."""
    lam = np.asarray(lambda_matrix, dtype=float)
    return lam - np.diag(exit_rates(lam))


def jump_matrix(lambda_matrix):
    """Return embedded jump-chain probabilities Lambda[i, j] / omega_i."""
    lam = np.asarray(lambda_matrix, dtype=float)
    omega = exit_rates(lam)
    p_jump = np.zeros_like(lam, dtype=float)
    active = omega > 0
    p_jump[active] = lam[active] / omega[active, None]
    return p_jump


def stationary_distribution(lambda_matrix):
    """Compute the stationary distribution pi satisfying pi^T Q = 0."""
    q = generator(lambda_matrix)
    n = q.shape[0]
    system = q.T.copy()
    rhs = np.zeros(n)
    system[-1, :] = 1.0
    rhs[-1] = 1.0
    pi = np.linalg.solve(system, rhs)
    pi[np.abs(pi) < 1e-14] = 0.0
    return pi / pi.sum()


def simulate_return_times(lambda_matrix, start_index, runs, rng):
    """Simulate first positive return times including the initial holding time.

    Starting in `start_index`, the clock in that state is sampled first. The
    returned time is the elapsed time until the chain next enters start_index.
    """
    omega = exit_rates(lambda_matrix)
    p_jump = jump_matrix(lambda_matrix)
    samples = np.empty(runs, dtype=float)
    for run in range(runs):
        state = start_index
        elapsed = 0.0
        left_start = False
        while True:
            elapsed += rng.exponential(1.0 / omega[state])
            state = int(rng.choice(p_jump.shape[0], p=p_jump[state]))
            if state != start_index:
                left_start = True
            if left_start and state == start_index:
                samples[run] = elapsed
                break
    return samples


def theoretical_return_time(lambda_matrix, start_index):
    """Compute E_i[T_i^+] including the initial holding time.

    The calculation uses the CTMC cycle identity E_i[T_i^+] = 1 / (pi_i omega_i),
    where pi is the stationary distribution and 1 / omega_i is the mean holding
    time per visit to state i.
    """
    pi = stationary_distribution(lambda_matrix)
    omega = exit_rates(lambda_matrix)
    return float(1.0 / (pi[start_index] * omega[start_index]))


def simulate_hitting_times(lambda_matrix, start_index, target_index, runs, rng):
    """Simulate CTMC hitting times from start_index to target_index."""
    omega = exit_rates(lambda_matrix)
    p_jump = jump_matrix(lambda_matrix)
    samples = np.empty(runs, dtype=float)
    for run in range(runs):
        state = start_index
        elapsed = 0.0
        while state != target_index:
            elapsed += rng.exponential(1.0 / omega[state])
            state = int(rng.choice(p_jump.shape[0], p=p_jump[state]))
        samples[run] = elapsed
    return samples


def theoretical_hitting_times(lambda_matrix, target_index):
    """Solve sum_j Q[i, j] h_j = -1 for i != target, with h_target = 0."""
    q = generator(lambda_matrix)
    n = q.shape[0]
    unknown = [idx for idx in range(n) if idx != target_index]
    system = q[np.ix_(unknown, unknown)]
    rhs = -np.ones(len(unknown))
    solution = np.linalg.solve(system, rhs)
    h = np.zeros(n, dtype=float)
    h[unknown] = solution
    return h
