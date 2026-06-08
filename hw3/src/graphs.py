"""Graph generators for Homework 3.

Two undirected graph families, returned as symmetric 'scipy.sparse' adjacency
matrices (CSR) with 0/1 entries and a zero diagonal:

* 'k_regular(n, k)'        -- the deterministic symmetric k-regular ring of 1.1.1;
* 'preferential_attachment(n, k, rng)' -- the random graph of 1.1.2.

Both are parameterized only by '(n, k)' so that changing the average degree means
changing a single argument (required for reuse in Problem 4).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def k_regular(n, k):
    """Return the symmetric k-regular ring adjacency on 'n' nodes.

    Each node 'i' is linked to the 'k' nodes whose index is closest modulo 'n',
    i.e. 'i +/- 1, ..., i +/- k/2'. 'k' must be even (so the graph is regular).
    """
    if k % 2 != 0:
        raise ValueError("k_regular requires an even k so every node has degree k")
    half = k // 2
    rows = []
    cols = []
    for offset in range(1, half + 1):
        neighbours = (np.arange(n) + offset) % n
        rows.extend(range(n))
        cols.extend(neighbours.tolist())
    # add both directions so the matrix is symmetric
    data = np.ones(2 * len(rows), dtype=np.int8)
    all_rows = np.array(rows + cols)
    all_cols = np.array(cols + rows)
    adjacency = sp.coo_matrix((data, (all_rows, all_cols)), shape=(n, n))
    adjacency = adjacency.tocsr()
    adjacency.data[:] = 1  # collapse any accidental duplicates to 0/1
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()
    return adjacency


def preferential_attachment(n, k, rng):
    """Return a preferential-attachment graph with 'n' nodes, avg degree ~ 'k'.

    Start from a complete graph on 'k0 = k + 1' nodes. Each later node attaches to
    'c existing nodes chosen without replacement, with probability proportional to
    their current degree. For odd 'k' we alternate 'floor(k/2)' and 'ceil(k/2)'
    new links per added node so the average degree still converges to 'k'.
    """
    k0 = k + 1
    if n < k0:
        raise ValueError(f"need at least k0 = k + 1 = {k0} nodes, got n = {n}")

    floor_c = k // 2
    ceil_c = k - floor_c  # equals floor_c for even k, floor_c + 1 for odd k

    # Adjacency as a list of sets, plus a running degree array.
    neighbours = [set() for _ in range(n)]
    degree = np.zeros(n, dtype=np.int64)

    # t = 1: complete graph on the first k0 nodes.
    for i in range(k0):
        for j in range(i + 1, k0):
            neighbours[i].add(j)
            neighbours[j].add(i)
    degree[:k0] = k0 - 1

    # t >= 2: add the remaining nodes one at a time.
    add_count = 0  # counts added nodes, to alternate floor/ceil for odd k
    for new_node in range(k0, n):
        c = floor_c if (add_count % 2 == 0) else ceil_c
        add_count += 1
        c = min(c, new_node)  # cannot link to more nodes than currently exist

        weights = degree[:new_node].astype(float)
        total = weights.sum()
        probabilities = weights / total if total > 0 else None
        targets = rng.choice(new_node, size=c, replace=False, p=probabilities)

        for target in targets:
            target = int(target)
            neighbours[new_node].add(target)
            neighbours[target].add(new_node)
        degree[new_node] += len(targets)
        degree[targets] += 1

    # Build the symmetric CSR adjacency from the neighbour sets.
    rows = []
    cols = []
    for i in range(n):
        for j in neighbours[i]:
            rows.append(i)
            cols.append(j)
    data = np.ones(len(rows), dtype=np.int8)
    adjacency = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()
    return adjacency


def degree_stats(adjacency):
    """Return '(mean_degree, degree_array)' for a symmetric 0/1 adjacency."""
    degrees = np.asarray(adjacency.sum(axis=1)).ravel().astype(int)
    return float(degrees.mean()), degrees
