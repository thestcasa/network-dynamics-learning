"""Network coordination / anti-coordination game for Homework 3, Part 2.

Players 'V = {0, ..., n-1}' (0-indexed) with actions 'A = {-1, +1}' interact on the
complete graph (each utility sums over all other players). Players '0..n1-1' are
coordinators, 'n1..n-1' are anti-coordinators:

    u_i = (1/2) * sum_{j != i} |x_i + x_j|     if i is a coordinator
    u_i = (1/2) * sum_{j != i} |x_i - x_j|     if i is an anti-coordinator

For binary +/-1 actions these reduce to "number of neighbours matching x_i" and
"number of neighbours differing from x_i" respectively.

We study the discrete-time asynchronous dynamics on the 2^n profiles: the (noiseless)
best-response chain and the noisy best-response (logit) chain, the latter in the
vanishing-noise limit.
"""

from __future__ import annotations

from itertools import product

import numpy as np
from scipy.special import expit  # numerically stable logistic sigmoid


def all_profiles(n):
    """Return all 2^n action profiles as tuples of +/-1, in a fixed canonical order."""
    return list(product((-1, 1), repeat=n))


def profile_index_map(n):
    """Return a dict mapping each profile tuple to its index in 'all_profiles(n)'."""
    return {profile: index for index, profile in enumerate(all_profiles(n))}


def utility(player, profile, n1):
    """Return the utility of 'player' under 'profile' (computed from the formula)."""
    xi = profile[player]
    others = (profile[j] for j in range(len(profile)) if j != player)
    if player < n1:  # coordinator
        return sum(abs(xi + xj) for xj in others) / 2.0
    return sum(abs(xi - xj) for xj in others) / 2.0  # anti-coordinator


def _flip(profile, player):
    """Return 'profile' with the action of 'player' negated."""
    flipped = list(profile)
    flipped[player] = -flipped[player]
    return tuple(flipped)


def _with_action(profile, player, action):
    """Return 'profile' with 'player' set to 'action' (+/-1)."""
    updated = list(profile)
    updated[player] = action
    return tuple(updated)


def nash_equilibria(n, n1):
    """Return the list of pure Nash equilibria (no strictly profitable deviation)."""
    equilibria = []
    for profile in all_profiles(n):
        is_nash = True
        for player in range(n):
            if utility(player, _flip(profile, player), n1) > utility(player, profile, n1):
                is_nash = False
                break
        if is_nash:
            equilibria.append(profile)
    return equilibria


def best_response_transition(n, n1):
    """Return the 2^n x 2^n transition matrix of the asynchronous best-response chain.

    At each step a player is chosen uniformly at random (probability 1/n) and switches
    to its unique best response only if that strictly increases its utility (ties keep
    the current action). Absorbing states are exactly the pure Nash equilibria.
    """
    profiles = all_profiles(n)
    index = profile_index_map(n)
    size = len(profiles)
    transition = np.zeros((size, size))
    for source, profile in enumerate(profiles):
        for player in range(n):
            flipped = _flip(profile, player)
            if utility(player, flipped, n1) > utility(player, profile, n1):
                target = flipped
            else:
                target = profile
            transition[source, index[target]] += 1.0 / n
    return transition


def logit_transition(n, n1, epsilon):
    """Return the 2^n x 2^n transition matrix of the asynchronous logit (noisy BR) chain.

    A player chosen uniformly at random sets x_i = +1 with probability
    'exp(u_i(+1)/eps) / (exp(u_i(+1)/eps) + exp(u_i(-1)/eps))'. The chain is ergodic
    for any epsilon > 0, so its limit distribution is the unique stationary one.
    """
    profiles = all_profiles(n)
    index = profile_index_map(n)
    size = len(profiles)
    transition = np.zeros((size, size))
    for source, profile in enumerate(profiles):
        for player in range(n):
            u_plus = utility(player, _with_action(profile, player, 1), n1)
            u_minus = utility(player, _with_action(profile, player, -1), n1)
            prob_plus = float(expit((u_plus - u_minus) / epsilon))
            target_plus = index[_with_action(profile, player, 1)]
            target_minus = index[_with_action(profile, player, -1)]
            transition[source, target_plus] += (1.0 / n) * prob_plus
            transition[source, target_minus] += (1.0 / n) * (1.0 - prob_plus)
    return transition


def stationary_distribution(transition):
    """Return the stationary distribution of an (ergodic) row-stochastic matrix.

    Solved as the least-squares solution of 'pi (P - I) = 0' with 'sum(pi) = 1'.
    This stays well-conditioned (and symmetry-preserving) for small noise levels, where
    an eigenvector solve becomes unreliable because the chain is nearly reducible.
    """
    size = transition.shape[0]
    matrix = np.vstack([transition.T - np.eye(size), np.ones(size)])
    rhs = np.zeros(size + 1)
    rhs[-1] = 1.0
    solution, *_ = np.linalg.lstsq(matrix, rhs, rcond=None)
    solution = np.clip(solution, 0.0, None)
    return solution / solution.sum()


def potential(profile, n1):
    """Exact potential of the pure game (n1 == n coordination, or n1 == 0 anti-).

    Coordination: number of agreeing pairs; anti-coordination: number of disagreeing
    pairs. Defined only for the fully homogeneous cases used for the limit (n1 in
    {0, n}); the mixed games are not needed for the requested limit distributions.
    """
    n = len(profile)
    if n1 == n:  # all coordinators
        return sum(
            1
            for i in range(n)
            for j in range(i + 1, n)
            if profile[i] == profile[j]
        )
    if n1 == 0:  # all anti-coordinators
        return sum(
            1
            for i in range(n)
            for j in range(i + 1, n)
            if profile[i] != profile[j]
        )
    raise ValueError("exact potential is only defined for n1 in {0, n}")


def vanishing_noise_limit(n, n1):
    """Return the exact vanishing-noise (stochastically stable) limit distribution.

    For log-linear (logit) learning on an exact potential game the stationary
    distribution is the Gibbs measure 'pi_eps(x) ~ exp(potential(x)/eps)'; as
    'eps -> 0' it concentrates uniformly on the potential maximizers.
    """
    profiles = all_profiles(n)
    potentials = np.array([potential(p, n1) for p in profiles], dtype=float)
    maximizers = potentials == potentials.max()
    return maximizers / maximizers.sum()


def limit_distribution_from_initial(transition, initial_index, max_iters=100000,
                                    tol=1e-15):
    """Return 'lim_{t->inf} P(X(t) = . | X(0) = initial_index)' by power iteration.

    Suitable for absorbing chains (the noiseless best-response dynamics): the row
    distribution converges to the absorption probabilities over the Nash equilibria.
    """
    size = transition.shape[0]
    distribution = np.zeros(size)
    distribution[initial_index] = 1.0
    for _ in range(max_iters):
        nxt = distribution @ transition
        if np.max(np.abs(nxt - distribution)) < tol:
            distribution = nxt
            break
        distribution = nxt
    return distribution
