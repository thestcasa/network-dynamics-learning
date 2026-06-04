# Homework 1 — Network Dynamics and Learning

This file is the working context for solving Homework 1.

The goal is to solve the full homework in Python and produce:

- well-commented Python code;
- numerical results for every sub-question;
- plots where requested;
- a final LaTeX/PDF-style lab report readable independently from the code;
- a README explaining how to reproduce the results.

Do not invent numerical values in the report. Use only values produced by the scripts.

---

## Available files

The project folder contains:

```text
hw1_en_2025-26.pdf
traffic.mat
capacities.mat
traveltime.mat
flow.mat
```

The `.mat` files are needed for Exercise 3.

---

## Suggested project structure

```text
hw1/
├── HOMEWORK1_CONTEXT.md
├── data/
│   ├── traffic.mat
│   ├── capacities.mat
│   ├── traveltime.mat
│   └── flow.mat
├── src/
│   ├── exercise1.py
│   ├── exercise2.py
│   ├── exercise3.py
│   └── utils.py
├── figures/
├── results/
├── report/
│   ├── report.tex
│   └── report.pdf
└── README.md
```

All scripts must run from the project root.

Recommended dependencies:

```bash
pip install numpy scipy networkx cvxpy matplotlib pandas
```

---

# Exercise 1 — Max-flow and cuts

## Network

Nodes:

```text
o, a, b, c, d
```

Directed edges and capacities:

```text
e1: o -> a, capacity 3
e2: a -> d, capacity 3
e3: o -> b, capacity 3
e4: b -> c, capacity 3
e5: c -> d, capacity 2
e6: a -> b, capacity 1
```

## Exercise 1(a)

Compute all `(o,d)` cuts.

A valid cut is a set `S` such that:

```text
o in S
d not in S
```

The intermediate nodes are `a,b,c`, so there are:

```text
2^3 = 8 cuts
```

For each cut:

- list crossing edges from `S` to `V \ S`;
- compute the cut capacity;
- identify the minimum cut;
- verify the result with the max-flow/min-cut theorem.

Output:

- table of all cuts;
- minimum cut capacity;
- maximum flow value.

---

## Exercise 1(b)

Given `x >= 0` extra integer units of capacity, distribute them over the existing six edges to maximize throughput from `o` to `d`.

Use:

```text
x = 0, 1, ..., 20
```

For each `x`:

- search over integer capacity allocations;
- compute max-flow;
- save best allocation;
- save best throughput.

Output:

- table with `x`, best throughput, best allocation;
- plot: maximum throughput vs `x`.

Important:

- The graph is small, so brute force is acceptable.
- Code should still be clean and general.

---

## Exercise 1(c)

Now add one new directed link with base capacity `1`, then distribute `x >= 0` extra integer capacity units.

Consider all possible directed links between distinct nodes that are not already present.

For each `x = 0,...,20`:

- test every candidate added link;
- add base capacity `1`;
- distribute the extra `x` capacity units;
- compute max-flow;
- report the best added link and best throughput.

Output:

- table with `x`, best added link, best throughput, best allocation;
- plot: maximum throughput vs `x`;
- comparison with Exercise 1(b).

---

# Exercise 2 — Katz centrality and PageRank

## Graph

Simple undirected graph with nodes:

```text
n1, n2, ..., n15
```

The graph must be reconstructed from Figure 2 of the homework PDF.

Important:

- Write the assumed edge list explicitly in `src/exercise2.py`.
- Print the edge list.
- Mention in the report if any edge is visually ambiguous.
- Make the edge list easy to modify.

---

## Exercise 2(a)

Compute Katz centrality with:

```text
beta = 0.15
uniform intrinsic centrality mu
```

Use adjacency matrix `W`.

Formula:

```text
x = (I - beta W)^(-1) mu
```

Check the spectral condition:

```text
beta < 1 / lambda_max(W)
```

Output:

- Katz centrality table;
- bar plot;
- comment on the most central nodes.

---

## Exercise 2(b)

Write and implement a distributed PageRank algorithm with:

```text
beta = 0.15
uniform intrinsic centrality mu
```

Use an iterative update rule.

Output:

- description of the distributed update;
- PageRank table;
- convergence tolerance and number of iterations;
- bar plot.

---

## Exercise 2(c)

Compare nodes:

```text
n6 and n9
```

Discuss:

- Katz centrality of `n6` vs `n9`;
- PageRank centrality of `n6` vs `n9`;
- degree;
- neighbors’ importance;
- bridge/hub role;
- difference between Katz and PageRank normalization.

---

## Exercise 2(d)

Compute PageRank for:

```text
beta in {0, 1/4, 1/2, 3/4, 1}
```

For each beta:

- compute PageRank;
- report `PR(n6)`;
- report `PR(n9)`;
- compute `PR(n6) - PR(n9)`.

Output:

- table;
- plot of `PR(n6) - PR(n9)` vs beta;
- state whether the difference is monotone in beta.

Explain:

- beta = 0: only intrinsic centrality matters, so all nodes should be equal;
- beta = 1: centrality depends only on graph structure / random-walk behavior.

---

# Exercise 3 — Traffic assignment

Use the `.mat` files.

Load data with:

```python
import scipy.io

f = scipy.io.loadmat("data/flow.mat")["flow"].reshape(28,)
C = scipy.io.loadmat("data/capacities.mat")["capacities"].reshape(28,)
B = scipy.io.loadmat("data/traffic.mat")["traffic"]
l = scipy.io.loadmat("data/traveltime.mat")["traveltime"].reshape(28,)
```

The incidence matrix `B` has:

- rows = nodes;
- columns = links;
- `+1` at the tail node;
- `-1` at the head node.

There are:

- 17 nodes;
- 28 directed links.

Delay function:

```text
tau_e(f_e) = l_e / (1 - f_e / C_e),    0 <= f_e < C_e
```

For numerical stability, use:

```text
0 <= f_e <= 0.999 * C_e
```

---

## Exercise 3(a)

Find the shortest path from node `1` to node `17`.

Use:

- graph reconstructed from incidence matrix `B`;
- edge weights `l`.

Output:

- path nodes;
- path edges;
- total free-flow travel time.

---

## Exercise 3(b)

Find the maximum flow from node `1` to node `17`.

Use:

- capacities `C`.

Output:

- max-flow value;
- flow on edges if useful;
- optional minimum cut.

---

## Exercise 3(c)

Given `flow.mat`, compute:

```text
nu = B @ f
```

Then define the exogenous inflow vector for the optimization problems:

```text
nu_new[0] = nu[0]
nu_new[16] = -nu[0]
all other entries = 0
```

Remember:

- Python indexing starts from 0;
- node 1 is index 0;
- node 17 is index 16.

Output:

- original `nu`;
- reduced exogenous inflow vector `nu_new`.

---

## Exercise 3(d)

Find the social optimum minimizing total travel time:

```text
minimize sum_e f_e * tau_e(f_e)
```

Equivalent convex expression:

```text
sum_e (l_e * C_e / (1 - f_e / C_e) - l_e * C_e)
```

Subject to:

```text
B @ f = nu_new
0 <= f <= 0.999 * C
```

Use CVXPY.

Output:

- optimal flow `f_star`;
- objective value;
- edge delays;
- total travel time.

---

## Exercise 3(e)

Find Wardrop equilibrium by minimizing Beckmann potential:

```text
sum_e integral_0^{f_e} tau_e(s) ds
```

For the given delay:

```text
-l_e * C_e * log(1 - f_e / C_e)
```

Subject to the same constraints.

Output:

- Wardrop flow `f_wardrop`;
- total travel time at Wardrop;
- comparison with social optimum;
- optional price of anarchy:

```text
PoA = total_cost(Wardrop) / total_cost(Social optimum)
```

---

## Exercise 3(f)

Introduce marginal-cost tolls so that Wardrop equilibrium coincides with the social optimum.

Use:

```text
omega_e = f_star_e * tau'_e(f_star_e)
```

where:

```text
tau'_e(f_e) = (l_e / C_e) / (1 - f_e / C_e)^2
```

Tolled Wardrop objective:

```text
Beckmann potential + sum_e omega_e * f_e
```

Output:

- toll vector `omega`;
- tolled Wardrop flow `f_tolled`;
- norm difference:

```text
||f_tolled - f_star||_2
```

Expected result:

- the norm should be close to zero up to solver tolerance.

---

## Exercise 3(g)

Now the system cost is total additional travel time compared to free flow:

```text
psi_e(f_e) = f_e * (tau_e(f_e) - l_e)
```

Solve:

```text
minimize sum_e psi_e(f_e)
```

Subject to the same flow constraints.

Derivative:

```text
psi'_e(f_e) = tau_e(f_e) - l_e + f_e * tau'_e(f_e)
```

Toll that decentralizes this optimum:

```text
omega_e = psi'_e(f_star_add_e) - tau_e(f_star_add_e)
        = -l_e + f_star_add_e * tau'_e(f_star_add_e)
```

Output:

- additional-delay optimum `f_star_add`;
- toll vector `omega_add`;
- tolled Wardrop flow `f_tolled_add`;
- norm difference:

```text
||f_tolled_add - f_star_add||_2
```

Expected result:

- the norm should be close to zero up to solver tolerance.

---

# Report structure

The final report should have:

```text
1. Introduction
2. Exercise 1 — Max-flow and cuts
3. Exercise 2 — Katz centrality and PageRank
4. Exercise 3 — Traffic assignment
5. Conclusion
6. Reproducibility appendix
```

The report must:

- explain the method;
- present numerical results;
- include required plots;
- avoid huge raw arrays unless formatted as tables;
- be readable without opening the code.

---

# Figures to generate

At minimum:

```text
figures/exercise1_b_throughput.png
figures/exercise1_c_added_link_throughput.png
figures/exercise2_katz.png
figures/exercise2_pagerank_beta015.png
figures/exercise2_pagerank_difference_beta.png
```

Exercise 3 may also generate comparison plots/tables if useful.

---

# Final checklist

Before finalizing, verify:

```text
[ ] Exercise 1(a): all 8 cuts listed
[ ] Exercise 1(b): capacity augmentation solved and plotted
[ ] Exercise 1(c): added-link case solved and plotted
[ ] Exercise 2(a): Katz centrality computed
[ ] Exercise 2(b): distributed PageRank algorithm described and implemented
[ ] Exercise 2(c): n6 and n9 compared
[ ] Exercise 2(d): beta sensitivity analyzed
[ ] Exercise 3(a): shortest path computed
[ ] Exercise 3(b): max-flow computed
[ ] Exercise 3(c): nu computed
[ ] Exercise 3(d): social optimum solved
[ ] Exercise 3(e): Wardrop equilibrium solved
[ ] Exercise 3(f): marginal-cost tolls verified
[ ] Exercise 3(g): additional-delay tolls verified
[ ] All figures saved
[ ] All scripts run from project root
[ ] No absolute paths
[ ] Report compiles
[ ] README completed
```