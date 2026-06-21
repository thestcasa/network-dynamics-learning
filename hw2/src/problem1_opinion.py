"""Problem 1 French-DeGroot opinion dynamics routines."""

from __future__ import annotations

import numpy as np

from problem1_single_ctmc import exit_rates


def laplacian(lambda_matrix):
    """Return L = diag(Lambda @ 1) - Lambda for dx/dt = -L x."""
    lam = np.asarray(lambda_matrix, dtype=float)
    return np.diag(exit_rates(lam)) - lam


def left_consensus_vector(lambda_matrix, component=None):
    """Return normalized nonnegative left nullvector pi^T L = 0.

    If component is given, the vector is computed on that submatrix and expanded
    to the full node set with zeros outside the component.
    """
    lam = np.asarray(lambda_matrix, dtype=float)
    n = lam.shape[0]
    if component is None:
        indices = list(range(n))
    else:
        indices = list(component)
    sub_l = laplacian(lam[np.ix_(indices, indices)])
    m = len(indices)
    system = sub_l.T.copy()
    rhs = np.zeros(m)
    system[-1, :] = 1.0
    rhs[-1] = 1.0
    weights = np.linalg.solve(system, rhs)
    weights[np.abs(weights) < 1e-12] = 0.0
    weights = weights / weights.sum()
    full = np.zeros(n, dtype=float)
    for local, global_index in enumerate(indices):
        full[global_index] = weights[local]
    return full


def matrix_exponential_action(matrix, vector, times):
    """Compute exp(matrix * t) @ vector for each t using eigendecomposition."""
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    coeffs = np.linalg.solve(eigenvectors, vector)
    trajectories = []
    for time in times:
        value = eigenvectors @ (np.exp(eigenvalues * time) * coeffs)
        trajectories.append(np.real_if_close(value, tol=1000).real)
    return np.vstack(trajectories)


def simulate_fdg(lambda_matrix, x0, times):
    """Simulate continuous-time French-DeGroot dx/dt = -L x."""
    return matrix_exponential_action(-laplacian(lambda_matrix), np.asarray(x0, dtype=float), times)


def remove_edges(lambda_matrix, node_to_index, edges):
    """Return a copy of Lambda with the directed edge weights set to zero."""
    modified = np.array(lambda_matrix, dtype=float, copy=True)
    for source, destination in edges:
        modified[node_to_index[source], node_to_index[destination]] = 0.0
    return modified


def adjacency(lambda_matrix, tol=1e-12):
    """Return adjacency lists for positive Lambda[i, j]."""
    lam = np.asarray(lambda_matrix, dtype=float)
    return [[j for j, value in enumerate(row) if value > tol] for row in lam]


def strongly_connected_components(lambda_matrix):
    """Compute SCCs with Tarjan's algorithm."""
    graph = adjacency(lambda_matrix)
    index = 0
    stack = []
    on_stack = set()
    indices = {}
    lowlinks = {}
    components = []

    def visit(node):
        nonlocal index
        indices[node] = index
        lowlinks[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)

        for neighbor in graph[node]:
            if neighbor not in indices:
                visit(neighbor)
                lowlinks[node] = min(lowlinks[node], lowlinks[neighbor])
            elif neighbor in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[neighbor])

        if lowlinks[node] == indices[node]:
            component = []
            while True:
                item = stack.pop()
                on_stack.remove(item)
                component.append(item)
                if item == node:
                    break
            components.append(sorted(component))

    for node in range(len(graph)):
        if node not in indices:
            visit(node)
    return sorted(components, key=lambda comp: (comp[0], len(comp)))


def sink_components(lambda_matrix):
    """Return SCCs with no outgoing edges to another SCC."""
    components = strongly_connected_components(lambda_matrix)
    component_of = {}
    for idx, comp in enumerate(components):
        for node in comp:
            component_of[node] = idx
    graph = adjacency(lambda_matrix)
    sinks = []
    for idx, comp in enumerate(components):
        outgoing = False
        for node in comp:
            for neighbor in graph[node]:
                if component_of[neighbor] != idx:
                    outgoing = True
                    break
            if outgoing:
                break
        if not outgoing:
            sinks.append(comp)
    return sinks


def reaches_component(lambda_matrix, start, component):
    """Return whether start can reach any node in component."""
    targets = set(component)
    graph = adjacency(lambda_matrix)
    seen = {start}
    stack = [start]
    while stack:
        node = stack.pop()
        if node in targets:
            return True
        for neighbor in graph[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return False


def graph_summary(lambda_matrix, node_names):
    """Summarize SCCs, sink SCCs, and consensus condition."""
    sccs = strongly_connected_components(lambda_matrix)
    sinks = sink_components(lambda_matrix)
    sink_reachable_from_all = []
    for sink in sinks:
        sink_reachable_from_all.append(
            all(
                reaches_component(lambda_matrix, node, sink)
                for node in range(len(node_names))
            )
        )
    return {
        "sccs": [[node_names[i] for i in comp] for comp in sccs],
        "sink_components": [[node_names[i] for i in comp] for comp in sinks],
        "num_sink_components": len(sinks),
        "single_sink_reachable_from_all": len(sinks) == 1 and bool(sink_reachable_from_all[0]),
    }


def asymptotic_state(lambda_matrix, x0):
    """Approximate the asymptotic state by evaluating exp(-Lt)x0 at large t."""
    return simulate_fdg(lambda_matrix, x0, np.array([250.0]))[0]
