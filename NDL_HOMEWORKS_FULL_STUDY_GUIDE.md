# Network Dynamics and Learning — Complete Homework Study Guide

---

## Table of Contents

- [0. How to Use This Guide](#0-how-to-use-this-guide)
- [1. Repository Map](#1-repository-map)
- [2. Mathematical Background Needed for All Homeworks](#2-mathematical-background-needed-for-all-homeworks)
  - [2.1 Graphs and Directed Graphs](#21-graphs-and-directed-graphs)
  - [2.2 Network Flows](#22-network-flows)
  - [2.3 Centrality Measures](#23-centrality-measures)
  - [2.4 Markov Chains and Random Walks](#24-markov-chains-and-random-walks)
  - [2.5 Opinion Dynamics](#25-opinion-dynamics)
  - [2.6 Interacting Particle Systems](#26-interacting-particle-systems)
  - [2.7 Epidemic Models](#27-epidemic-models)
  - [2.8 Network Games](#28-network-games)
- [3. Homework 1 — Full Study Guide](#3-homework-1--full-study-guide)
- [4. Homework 2 — Full Study Guide](#4-homework-2--full-study-guide)
- [5. Homework 3 — Study Guide and Implementation Roadmap](#5-homework-3--study-guide-and-implementation-roadmap)
- [6. Cross-Homework Connections](#6-cross-homework-connections)
- [7. Master Formula Sheet](#7-master-formula-sheet)
- [8. Master Oral Exam Questions](#8-master-oral-exam-questions)
- [9. Common Mistakes and How to Avoid Them](#9-common-mistakes-and-how-to-avoid-them)
- [10. What I Should Be Able to Explain in 5 Minutes](#10-what-i-should-be-able-to-explain-in-5-minutes)

---

## 0. How to Use This Guide

This document is a **complete study guide and oral defense manual** for all three homework assignments of the Network Dynamics and Learning course. It is designed so that you can:

1. **Study from zero**: every concept is explained starting from the mathematical foundations.
2. **Prepare for the oral exam**: every exercise is mapped to theory, code, results, and likely exam questions.
3. **Understand the codebase**: every source file is explained in terms of what it computes, why, and how.
4. **Verify results**: all numerical values cited here come from actual CSV/result files in the repository. Nothing is invented.

**Reading strategy**:
- If you need to review the theory first, start with **Section 2**.
- If you need to prepare a specific homework, jump to **Sections 3, 4, or 5**.
- If you need a quick formula lookup, use **Section 7**.
- If you need to practice oral defense questions, use **Section 8**.
- If you are looking for common traps, see **Section 9**.

**Notation conventions used throughout**:
- Matrices are written as capital letters: $W$, $B$, $\Lambda$, $Q$, $P$, $L$.
- Vectors are lowercase: $f$, $\nu$, $\pi$, $\omega$, $\mu$, $x$, $h$.
- Scalars: $\beta$, $\lambda$, $\rho$, $n$, $k$, $T$.
- Node sets: $V$; edge sets: $E$.
- "Source" means where flow originates; "sink" means where flow is absorbed.

---

## 1. Repository Map

The repository root contains three homework folders:

```
network-dynamics-learning/
├── hw1/          ← Homework 1: max-flow, centrality, traffic assignment
├── hw2/          ← Homework 2: CTMC, opinion dynamics, particle systems
├── hw3/          ← Homework 3: epidemics, games (currently empty)
└── NDL_HOMEWORKS_FULL_STUDY_GUIDE.md   ← This file
```

### hw1/ Structure

| Path | Type | Purpose |
|------|------|---------|
| `hw1/hw1_en_2025-26.pdf` | PDF | Official assignment text |
| `hw1/HOMEWORK1_CONTEXT.md` | Markdown | Working context: problem statements, conventions, data formats |
| `hw1/ORAL_DEFENSE_HW1.md` | Markdown | Comprehensive oral defense guide with Q&A |
| `hw1/RESULTS_LOG.md` | Markdown | All numerical results with tables |
| `hw1/README.md` | Markdown | How to run, dependencies, output files |
| `hw1/src/exercise1.py` | Python | Max-flow, cuts, capacity augmentation |
| `hw1/src/exercise2.py` | Python | Katz centrality and distributed PageRank |
| `hw1/src/exercise3.py` | Python | Traffic assignment, social optimum, Wardrop, tolls |
| `hw1/src/utils.py` | Python | Shared utilities |
| `hw1/data/` | MATLAB `.mat` | Input data for Exercise 3 (flow, capacities, traffic, traveltime) |
| `hw1/results/` | CSV | Generated numerical results |
| `hw1/figures/` | PNG | Generated plots |
| `hw1/report/` | LaTeX | Report source and compiled PDF |

### hw2/ Structure

| Path | Type | Purpose |
|------|------|---------|
| `hw2/data/hw2_en_2025-26.pdf` | PDF | Official assignment text |
| `hw2/HW2_CONTEXT.md` | Markdown | Problem statements and conventions |
| `hw2/CODEBASE_CONTEXT.md` | Markdown | Codebase documentation |
| `hw2/HW2_STUDY_GUIDE.md` | Markdown | Existing study guide for HW2 |
| `hw2/RESULTS_LOG.md` | Markdown | All numerical results with seeds and CIs |
| `hw2/README.md` | Markdown | How to run |
| `hw2/src/constants.py` | Python | Shared constants: node order, Lambda matrices, paths |
| `hw2/src/problem1_single_ctmc.py` | Python | Single-particle CTMC routines |
| `hw2/src/problem1_opinion.py` | Python | French-DeGroot dynamics |
| `hw2/src/problem2_many_particles.py` | Python | Many-particle simulations |
| `hw2/src/problem3_open_network.py` | Python | Open network simulator |
| `hw2/src/run_problem1.py` | Python | Main script for Problem 1 |
| `hw2/src/run_problem2.py` | Python | Main script for Problem 2 |
| `hw2/src/run_problem3.py` | Python | Main script for Problem 3 |
| `hw2/results/` | CSV | 13 result files |
| `hw2/figures/` | PNG | 11 figure files |
| `hw2/report/` | LaTeX | Report source |

### hw3/ Structure

**HW3 is currently empty.** No source code, results, figures, or PDF have been added yet. Section 5 of this guide provides a theoretical and implementation roadmap based on the expected HW3 topics (epidemic models and network games).

---

## 2. Mathematical Background Needed for All Homeworks

### 2.1 Graphs and Directed Graphs

**Nodes and edges.** A graph $G = (V, E)$ consists of a set of nodes $V$ and a set of edges $E$. In an undirected graph, each edge is an unordered pair $\{u, v\}$. In a directed graph (digraph), each edge (arc) is an ordered pair $(u, v)$, where $u$ is the **tail** and $v$ is the **head**.

**Adjacency matrix.** For an undirected graph with $n$ nodes, the adjacency matrix $W \in \{0,1\}^{n \times n}$ has $W_{ij} = 1$ if nodes $i$ and $j$ are connected. $W$ is symmetric. For a weighted graph, $W_{ij}$ can be any nonneg real number.

**Incidence matrix.** For a directed graph with $n$ nodes and $m$ edges, the incidence matrix $B \in \{-1,0,+1\}^{n \times m}$ has, for each column (edge $e$ from $u$ to $v$):
- $B_{u,e} = +1$ (tail)
- $B_{v,e} = -1$ (head)
- All other entries zero.

The incidence matrix encodes flow conservation: if $f$ is a flow vector, then $Bf = \nu$, where $\nu_i$ is the net exogenous inflow at node $i$.

**Strongly connected components (SCCs).** In a directed graph, an SCC is a maximal subset of nodes where every node can reach every other node via directed paths. SCCs determine convergence properties of dynamics on directed networks.

**Sink components.** A sink SCC is an SCC with no outgoing edges to other SCCs. Sink components act as absorbing regions for dynamics.

> **What to remember**: The direction of edges matters fundamentally. In HW1 Exercise 1, only edges *leaving* $S$ count toward cut capacity. In HW2, the row-source convention $\Lambda_{ij}$ = rate from $i$ to $j$ must be used consistently.

### 2.2 Network Flows

**Feasible flow.** A flow $f$ on a directed capacitated network is feasible if:
1. $0 \leq f_e \leq c_e$ for every edge $e$ (capacity constraint).
2. At every intermediate node, flow is conserved: inflow = outflow (or equivalently, $Bf = \nu$).

**Source and sink.** The source node injects flow into the network; the sink node absorbs it. In HW1 Exercise 1, the source is $o$ and the sink is $d$. In HW1 Exercise 3, the source is node 1 and the sink is node 17.

**Source-sink cut.** A cut $(S, V \setminus S)$ is a partition of the node set with the source in $S$ and the sink not in $S$. The **directed cut capacity** is the sum of capacities of edges whose tail is in $S$ and whose head is in $V \setminus S$. Edges entering $S$ do NOT count.

$$C(S) = \sum_{\substack{e = (u,v): \\ u \in S,\; v \notin S}} c_e$$

**Max-flow / min-cut theorem.** The maximum feasible flow from source to sink equals the minimum cut capacity:

$$\max_{\text{feasible } f} \text{value}(f) = \min_{S:\; s \in S,\; t \notin S} C(S)$$

This is the central theorem connecting flow optimization to graph structure.

**Shortest path.** The path from source to sink minimizing total edge weight. In HW1 Exercise 3, free-flow travel times $l_e$ are used as weights.

**Traffic assignment.** Unlike max-flow, traffic assignment considers congestion: the travel time on a link increases with the flow on that link.

**Social optimum.** The flow assignment minimizing total system travel time:

$$\min_f \sum_e f_e \cdot \tau_e(f_e) \quad \text{s.t. } Bf = \nu,\; 0 \leq f \leq C$$

**Wardrop equilibrium.** The selfish-routing equilibrium where no driver can reduce travel time by switching routes. All used routes between the same OD pair have equal travel time. It is found by minimizing the **Beckmann potential**:

$$\min_f \sum_e \int_0^{f_e} \tau_e(s)\, ds \quad \text{s.t. } Bf = \nu,\; 0 \leq f \leq C$$

**Marginal-cost tolls.** To make selfish users replicate the social optimum, charge each link a toll equal to the congestion externality:

$$\omega_e = f_e^* \cdot \tau'_e(f_e^*)$$

This makes the perceived cost (delay + toll) equal to the marginal social cost.

> **What to remember**: Social optimum minimizes total cost; Wardrop minimizes the Beckmann potential; tolls bridge the gap. The price of anarchy measures how much worse selfish routing is compared to the social optimum.

### 2.3 Centrality Measures

**Katz centrality.** Measures node importance by counting discounted walks of all lengths:

$$x = (I - \beta W)^{-1} \mu$$

where $W$ is the adjacency matrix, $\beta$ is a discount factor, and $\mu$ is the intrinsic centrality vector. The walk expansion is:

$$x = \sum_{k=0}^{\infty} \beta^k W^k \mu$$

This converges if and only if $\beta < 1/\lambda_{\max}(W)$, where $\lambda_{\max}(W)$ is the spectral radius (largest eigenvalue) of $W$.

**Interpretation of Katz.** A node is central if many walks of all lengths reach it. Dense subgraphs (cliques) create many walks and boost Katz centrality for their members.

**PageRank.** A centrality measure based on a random walk with teleportation:

$$x_i^{(t+1)} = (1-\beta)\mu_i + \beta \sum_{j \in N(i)} \frac{x_j^{(t)}}{\deg(j)}$$

In matrix form: $x^{(t+1)} = \beta P x^{(t)} + (1-\beta)\mu$, where $P$ is the column-stochastic transition matrix.

**Key difference between Katz and PageRank.** Katz uses raw adjacency weights (so high-degree nodes contribute more total mass). PageRank normalizes each node's contribution by its degree (so a degree-1 leaf sends all its mass to its only neighbor).

**Role of $\beta$.** 
- At $\beta = 0$: centrality equals intrinsic $\mu$, so all nodes are equal for uniform $\mu$.
- At $\beta = 1$ (PageRank): centrality becomes the stationary distribution of the random walk. For connected, non-bipartite undirected graphs, $\pi_i = \deg(i) / (2|E|)$.
- For intermediate $\beta$: a blend of intrinsic importance and network structure.

> **What to remember**: Katz favors nodes involved in many walks (dense neighborhoods). PageRank favors nodes that receive concentrated normalized mass from their neighbors (e.g., star centers with degree-1 leaves).

### 2.4 Markov Chains and Random Walks

**Discrete-time Markov chain (DTMC).** A stochastic process $\{X_n\}$ where the next state depends only on the current state: $\Pr(X_{n+1} = j \mid X_n = i) = P_{ij}$.

**Continuous-time Markov chain (CTMC).** A stochastic process $\{X(t)\}$ in continuous time. The dynamics are governed by:

1. **Transition rate matrix** $\Lambda$: $\Lambda_{ij} \geq 0$ is the rate of transitions from $i$ to $j$ (for $i \neq j$).
2. **Exit rates**: $\omega_i = \sum_j \Lambda_{ij}$ is the total rate of leaving node $i$.
3. **Generator**: $Q = \Lambda - \text{diag}(\omega)$, so $Q_{ij} = \Lambda_{ij}$ for $i \neq j$ and $Q_{ii} = -\omega_i$. Rows of $Q$ sum to zero.
4. **Embedded jump chain**: $P^{\text{jump}}_{ij} = \Lambda_{ij} / \omega_i$ gives the probability of jumping to $j$ given a departure from $i$.

**Simulation algorithm for CTMC:**
1. If at node $i$, draw holding time $\Delta t \sim \text{Exp}(\omega_i)$.
2. Jump to node $j$ with probability $\Lambda_{ij}/\omega_i$.
3. Update time $t \leftarrow t + \Delta t$ and position $i \leftarrow j$.

**Exponential waiting times.** The memoryless property of the exponential distribution is what makes CTMCs "Markov" in continuous time. If multiple independent exponential clocks are ticking (one per outgoing edge with rate $\Lambda_{ij}$), the minimum waiting time is exponential with rate $\omega_i = \sum_j \Lambda_{ij}$, and the clock that rings first identifies the next jump destination.

**Stationary distribution.** A probability vector $\pi$ satisfying $\pi^T Q = 0$ (equivalently, $\pi^T \Lambda = \pi^T \text{diag}(\omega)$). For an irreducible CTMC, $\pi$ is unique and gives the long-run fraction of time spent in each state.

**Return time.** The expected time to return to state $i$, starting from $i$:

$$\mathbb{E}_i[T_i^+] = \frac{1}{\pi_i \cdot \omega_i}$$

This formula connects the return time to the stationary distribution: the chain returns to $i$ more often if $\pi_i$ is large (state is visited frequently) and $\omega_i$ is large (departures happen quickly, so cycles are short).

**Hitting time.** The expected time to reach state $d$ for the first time, starting from state $i$. Let $h_i = \mathbb{E}_i[T_d]$. Then $h_d = 0$ and for $i \neq d$:

$$\sum_j Q_{ij} h_j = -1 \quad \Leftrightarrow \quad h_i = \frac{1}{\omega_i} + \sum_j P^{\text{jump}}_{ij} h_j$$

This is a linear system solvable by standard methods.

> **What to remember**: $\Lambda$ contains nonneg transition rates; $Q$ has negative diagonal and rows sum to zero; $P^{\text{jump}}$ is the embedded jump matrix. Never confuse these three.

### 2.5 Opinion Dynamics

**French-DeGroot model.** Agents update their opinions toward weighted averages of their neighbors' opinions. In continuous time:

$$\frac{dx}{dt} = -Lx, \quad L = \text{diag}(\Lambda \mathbf{1}) - \Lambda$$

where $L$ is the graph Laplacian. Each node $i$ updates as:

$$\frac{dx_i}{dt} = \sum_j \Lambda_{ij}(x_j - x_i)$$

The solution is $x(t) = e^{-Lt} x(0)$.

**Consensus.** All opinions converge to the same value: $\lim_{t \to \infty} x_i(t) = x^*$ for all $i$. This happens if and only if the underlying directed graph has a **unique sink SCC reachable from all nodes**.

**Consensus value.**

$$x^* = \sum_i \pi_i \cdot x_i(0)$$

where $\pi$ is the normalized left eigenvector of $L$ associated with eigenvalue 0: $\pi^T L = 0$, $\sum_i \pi_i = 1$.

**Influence vector.** The vector $\pi$ determines each node's influence on the consensus value. If $\pi_i$ is large, node $i$'s initial opinion matters more.

**Variance of consensus value.** If initial opinions are independent random variables:

$$\text{Var}(x^*) = \sum_i \pi_i^2 \cdot \text{Var}(x_i(0))$$

Cross terms vanish by independence.

**Effect of edge removal.** Removing directed edges changes the SCC structure. If edge removal creates multiple sink SCCs, global consensus is lost: different sink components converge to different values, and transient nodes converge to weighted combinations of sink component values.

> **What to remember**: Consensus depends on graph structure (unique reachable sink SCC). The consensus value is a weighted average of initial opinions, with weights given by the left eigenvector of the Laplacian.

### 2.6 Interacting Particle Systems

**Particle perspective.** Track each particle individually as an independent CTMC. If all particles follow the same transition rates, their average return time matches the single-particle theoretical value. Using $N$ particles reduces Monte Carlo variance but does not change the expected return time.

**Node perspective.** Track the number of particles at each node: $N(t) = (N_1(t), \ldots, N_n(t))$. In a closed network with proportional service:
- Node $i$ forwards particles at total rate $N_i(t) \cdot \omega_i$.
- Destination chosen by embedded jump probabilities $P^{\text{jump}}_{ij}$.

At stationarity, the expected number of particles at node $i$ is $N \cdot \pi_i$.

**Open networks.** Particles arrive from outside (Poisson process with rate $\lambda$) and eventually exit. Two service rules:

1. **Proportional rate**: service rate at node $i$ is $\omega_i \cdot N_i(t)$. Capacity grows with population. Every finite arrival rate gives finite expected population (stable for all $\lambda$).

2. **Fixed rate**: service clock at node $i$ has rate $\omega_i$ regardless of queue size. Empty ticks do nothing. Bottleneck capacity limits the maximum sustainable input rate. If the routing load at node $i$ per unit input is $\rho_i$, the maximum stable $\lambda$ is:

$$\lambda_{\max} = \min_i \frac{\omega_i}{\rho_i}$$

> **What to remember**: Proportional service scales with demand (no blow-up). Fixed service has bottlenecks (blow-up above critical $\lambda$). The particle perspective and node perspective are two views of the same system.

### 2.7 Epidemic Models

**SIR model.** Each node is in one of three states:
- **S** (Susceptible): can be infected.
- **I** (Infected): can spread the disease and can recover.
- **R** (Recovered): immune and cannot infect or be infected.

**Infection rule.** If node $i$ is susceptible and has $m$ infected neighbors, the probability of becoming infected in one time step is:

$$p_{\text{infection}} = 1 - (1-\beta)^m$$

where $\beta$ is the per-neighbor infection probability. This comes from assuming independent Bernoulli trials for each infected neighbor.

**Recovery rule.** An infected node recovers with probability $\rho$ per time step, independently of neighbors.

**Preferential attachment graphs.** Start with a complete graph of $k+1$ nodes. Add nodes one at a time, each connecting to existing nodes with probability proportional to their current degree. This creates scale-free networks with hubs (high-degree nodes) that are important for epidemic spread.

**Parameter estimation (RMSE).** Given observed epidemic data, find parameters $(k, \beta, \rho)$ that minimize the root mean squared error between simulated and observed newly infected counts:

$$\text{RMSE} = \sqrt{\frac{1}{T}\sum_{t=1}^T \left(\hat{I}(t) - I_0(t)\right)^2}$$

where $\hat{I}(t)$ is the simulated average and $I_0(t)$ is the observed data.

**Vaccination.** Vaccinated nodes cannot infect or be infected. Random vaccination among non-vaccinated nodes reduces epidemic spread by removing potential transmission paths.

### 2.8 Network Games

**Setup.** A set of $n$ players connected by a graph. Each player $i$ chooses an action $a_i \in \{-1, +1\}$.

**Coordination games.** A coordination player wants to match the actions of their neighbors. Utility increases when $a_i = a_j$ for neighbors $j$.

**Anti-coordination games.** An anti-coordination player wants to differ from neighbors. Utility increases when $a_i \neq a_j$.

**Mixed networks.** Some players are coordinators, others are anti-coordinators. The number of each type ($n_1$ coordinators vs $n - n_1$ anti-coordinators) determines the Nash equilibria.

**Nash equilibrium.** An action profile where no player can unilaterally improve their utility by changing their action.

**Best-response dynamics.** Players take turns updating to their best response given current neighbors' actions:
- **Asynchronous**: one randomly chosen player updates at each step.
- **Noisy best response**: players choose the better action with higher probability but not certainty. A noise parameter $\epsilon$ controls randomness.

**Vanishing noise limit.** As $\epsilon \to 0$, the noisy best-response dynamics selects among Nash equilibria. The transition graph over action profiles determines which equilibria are absorbing.

**Transition graph.** The $2^n$ possible action profiles are nodes. Edges represent transitions under best-response dynamics. Analysis of this graph reveals which equilibria are reachable and stable.

> **What to remember**: Nash equilibria depend on the mix of coordination/anti-coordination players. Best-response dynamics can converge to different equilibria depending on initial conditions and noise.

---

## 3. Homework 1 — Full Study Guide

### 3.1 What Homework 1 Is About

Homework 1 covers three main topics:

1. **Exercise 1**: Max-flow, cuts, and capacity augmentation on a small directed network.
2. **Exercise 2**: Katz centrality and distributed PageRank on an undirected graph.
3. **Exercise 3**: Traffic assignment with social optimum, Wardrop equilibrium, and tolls on a 17-node traffic network.

All code is in `hw1/src/`. Results are in `hw1/results/`. Figures are in `hw1/figures/`.

### 3.2 HW1 Exercise 1 — Network Cuts, Capacity Expansion, Extra Link

#### 3.2.1 Network Setup

The graph is directed with 5 nodes: $V = \{o, a, b, c, d\}$. Source is $o$, sink is $d$.

| Edge | Arc | Capacity |
|------|-----|----------|
| $e_1$ | $o \to a$ | 3 |
| $e_2$ | $a \to d$ | 3 |
| $e_3$ | $o \to b$ | 3 |
| $e_4$ | $b \to c$ | 3 |
| $e_5$ | $c \to d$ | 2 |
| $e_6$ | $a \to b$ | 1 |

**Code reference**: `hw1/src/exercise1.py`, lines 21–41 define `NODES`, `SOURCE`, `SINK`, and `BASE_EDGES`.

#### 3.2.2 Part (a): All Cuts and Minimum Cut

**Problem**: Enumerate all $o$-$d$ cuts and find the minimum cut.

**Theory**: A cut is a set $S$ with $o \in S$ and $d \notin S$. The three intermediate nodes $\{a, b, c\}$ can each be in or out of $S$, giving $2^3 = 8$ cuts. The directed cut capacity counts only edges leaving $S$.

**Results** (from `hw1/results/exercise1_cuts.csv`):

| $S$ | $V \setminus S$ | Crossing edges | Capacity |
|-----|-----------------|----------------|----------|
| $\{o\}$ | $\{a,b,c,d\}$ | $e_1, e_3$ | 6 |
| $\{o,a\}$ | $\{b,c,d\}$ | $e_2, e_3, e_6$ | 7 |
| $\{o,b\}$ | $\{a,c,d\}$ | $e_1, e_4$ | 6 |
| $\{o,c\}$ | $\{a,b,d\}$ | $e_1, e_3, e_5$ | 8 |
| $\{o,a,b\}$ | $\{c,d\}$ | $e_2, e_4$ | 6 |
| $\{o,a,c\}$ | $\{b,d\}$ | $e_2, e_3, e_5, e_6$ | 9 |
| $\{o,b,c\}$ | $\{a,d\}$ | $e_1, e_5$ | **5** |
| $\{o,a,b,c\}$ | $\{d\}$ | $e_2, e_5$ | **5** |

**Minimum cut = 5**, achieved by $\{o,b,c\}$ and $\{o,a,b,c\}$.

**Maximum flow = 5** (computed by NetworkX), confirming the max-flow/min-cut theorem.

**Code walkthrough**: The function `enumerate_cuts()` generates all subsets of intermediate nodes, prepending $o$. The function `crossing_edge_indices()` finds edges leaving $S$, and `cut_capacity()` sums their capacities. The function `max_flow_value()` uses NetworkX's `maximum_flow()`.

**Why the result makes sense**: The bottleneck is edge $e_5$ ($c \to d$, capacity 2). The two minimum cuts both include $c$ on the source side, so the only way to reach $d$ is through $e_5$ (capacity 2) and $e_2$ (capacity 3), totaling 5.

> **Common mistakes**:
> - Counting edges *entering* $S$ in the cut capacity. Only edges leaving $S$ count.
> - Forgetting that $o$ must be in $S$ and $d$ must not be in $S$.

> **Oral defense questions**:
> 1. *Why are there exactly 8 cuts?* — Three intermediate nodes, each in/out: $2^3 = 8$.
> 2. *Why only edges leaving S?* — Because only those carry flow toward the sink.
> 3. *What does the max-flow/min-cut theorem say?* — Maximum feasible flow equals minimum cut capacity.

#### 3.2.3 Part (b): Capacity Augmentation on Existing Links

**Problem**: Distribute $x$ extra integer capacity units over the 6 existing links to maximize throughput, for $x = 0, \ldots, 20$.

**Method**: Exhaustive enumeration over all integer compositions of $x$ into 6 parts. For each allocation, compute max-flow via the min-cut formula (throughput = minimum augmented cut capacity). This is exact because the graph is small.

**Code reference**: `hw1/src/exercise1.py`, functions `integer_compositions()`, `allocations_up_to()`, `cut_masks()`, `throughput_from_allocation()`, and `optimize_extra_capacity()`.

**Key insight**: The function `throughput_from_allocation()` (line 145) computes throughput as the minimum over all cuts of the augmented cut capacity. This avoids rerunning NetworkX for every allocation — it directly applies the max-flow/min-cut theorem.

**Results** (from `hw1/results/exercise1_b_results.csv`):

The optimal throughput follows:

$$T_b(x) = 5 + \lceil x/2 \rceil$$

For example: $T_b(0) = 5$, $T_b(1) = 6$, $T_b(2) = 6$, $T_b(3) = 7$, $\ldots$, $T_b(20) = 15$.

**Why $\lceil x/2 \rceil$?** Extra capacity must often be added to *complementary bottlenecks* before throughput increases. For instance, adding 1 unit to $e_5$ raises the cut $\{o,b,c\}$ from 5 to 6, but to also raise $\{o,a,b,c\}$, we need to add to $e_2$ as well. Bottlenecks come in pairs.

**Figure**: `hw1/figures/exercise1_b_throughput.png` — Staircase plot of throughput vs. $x$.

> **Oral defense questions**:
> 1. *Why doesn't every extra unit increase throughput?* — Capacity must be paired across bottleneck structures.
> 2. *Is the optimal allocation unique?* — No, ties exist.
> 3. *How is global optimality guaranteed?* — Exhaustive integer enumeration.

#### 3.2.4 Part (c): Adding One New Directed Link

**Problem**: Add one new directed link (base capacity 1) and distribute $x$ extra capacity units. Find the best link and allocation.

**Method**: Test all candidate directed links between distinct nodes not already present. For each candidate, solve the same capacity augmentation problem.

**Code reference**: `hw1/src/exercise1.py`, function `candidate_added_links()` and `optimize_added_link_case()`.

**Results** (from `hw1/results/exercise1_c_results.csv`):

- At $x = 0$: best link is $b \to d$, throughput = 6.
- For $x \geq 1$: best link is $o \to d$, throughput = $6 + x$.

The direct link $o \to d$ is optimal because it crosses *every* $o$-$d$ cut (since $o$ is always in $S$ and $d$ is always outside $S$). Every extra unit on this link directly increases throughput by 1.

$$T_c(x) = 6 + x$$

**Defensive note**: The $o \to d$ link may seem trivial. The assignment does not explicitly forbid direct source-sink arcs. If challenged: "I included it under the literal candidate set. If forbidden, the optimization should be rerun with a restricted candidate set."

**Comparison with part (b)**: The improvement from adding a link grows linearly: $T_c(x) - T_b(x) = 1 + x - \lceil x/2 \rceil$.

**Figure**: `hw1/figures/exercise1_c_added_link_throughput.png` — Both curves plotted together.

### 3.3 HW1 Exercise 2 — Katz and PageRank Centrality

#### 3.3.1 Graph Reconstruction

The graph is undirected with 15 nodes $n_1, \ldots, n_{15}$ and 24 edges, reconstructed from Figure 2 of the assignment PDF.

**Key structure**:
- $n_1, \ldots, n_6$ form a **complete graph** (clique, 15 edges).
- $n_6$ connects to $n_{15}$ and to $n_7$.
- Path: $n_6 - n_7 - n_8 - n_9$.
- $n_9$ is the center of a **star**: connected to $n_{10}, n_{11}, n_{12}, n_{13}, n_{14}$.

**Code reference**: `hw1/src/exercise2.py`, lines 24–52 define the edge list.

**Degrees**: $\deg(n_1) = \ldots = \deg(n_5) = 5$; $\deg(n_6) = 7$; $\deg(n_7) = \deg(n_8) = 2$; $\deg(n_9) = 6$; $\deg(n_{10}) = \ldots = \deg(n_{14}) = 1$; $\deg(n_{15}) = 1$. Total: $2|E| = 48$.

#### 3.3.2 Part (a): Katz Centrality

**Formula**: $x = (I - \beta W)^{-1} \mu$ with $\beta = 0.15$, $\mu_i = 1/15$.

**Spectral check**:
- $\lambda_{\max}(W) = 5.0710614673$
- Threshold: $1/\lambda_{\max} = 0.1971973731$
- $\beta = 0.15 < 0.1972$ ✓

**Code reference**: `hw1/src/exercise2.py`, function `katz_centrality()` (line 81) — computes eigenvalues via `eigvalsh`, then solves the linear system.

**Results** (from `hw1/results/exercise2_katz.csv`):

| Rank | Node | Degree | Katz |
|------|------|--------|------|
| 1 | $n_6$ | 7 | 0.3178 |
| 2 (tie) | $n_1$–$n_5$ | 5 | 0.2858 |
| 7 | $n_9$ | 6 | 0.1498 |

**Why $n_6$ wins**: Katz counts discounted walks. $n_6$ is inside the dense clique AND bridges to the rest of the graph, so an enormous number of walks pass through it.

**Why $n_1$–$n_5$ are tied**: They are symmetric within the clique (all have the same structural position).

**Figure**: `hw1/figures/exercise2_katz.png`

#### 3.3.3 Part (b): Distributed PageRank

**Update rule**:

$$x_i^{(t+1)} = (1 - \beta)\mu_i + \beta \sum_{j \in N(i)} \frac{x_j^{(t)}}{\deg(j)}$$

**Code reference**: `hw1/src/exercise2.py`, function `distributed_pagerank()` (line 98). The transition matrix $P$ is column-stochastic: $P_{ij} = W_{ij}/\deg(j)$ (neighbor $j$ distributes mass equally among its $\deg(j)$ neighbors).

**Convergence**: With $\beta = 0.15$, tolerance $10^{-12}$, converged in **15 iterations** (final L1 diff $2.011 \times 10^{-13}$).

**Results** (from `hw1/results/exercise2_pagerank_beta015.csv`):

| Rank | Node | Degree | PageRank |
|------|------|--------|----------|
| 1 | $n_9$ | 6 | 0.1060 |
| 2 | $n_6$ | 7 | 0.0801 |
| 3 (tie) | $n_1$–$n_5$ | 5 | 0.0663 |

**Why $n_9$ wins**: Its degree-1 leaves ($n_{10}$–$n_{14}$) each send ALL their PageRank mass to $n_9$ (since they have only one neighbor). Clique nodes split their mass among 5 neighbors, diluting what $n_6$ receives.

**Figure**: `hw1/figures/exercise2_pagerank_beta015.png`

#### 3.3.4 Part (c): Comparison of $n_6$ and $n_9$

| | $n_6$ | $n_9$ |
|---|---|---|
| Degree | 7 | 6 |
| Katz centrality | 0.3178 (rank 1) | 0.1498 (rank 7) |
| PageRank ($\beta=0.15$) | 0.0801 (rank 2) | 0.1060 (rank 1) |

**Core insight**: Katz rewards raw walk density → favors $n_6$ (clique bridge). PageRank rewards concentrated normalized incoming mass → favors $n_9$ (star center with degree-1 leaves).

#### 3.3.5 Part (d): Beta Sensitivity

For $\beta \in \{0, 1/4, 1/2, 3/4, 1\}$:

| $\beta$ | $\text{PR}(n_6)$ | $\text{PR}(n_9)$ | $\text{PR}(n_6) - \text{PR}(n_9)$ |
|---------|-------------------|-------------------|--------------------------------------|
| 0 | 0.0667 | 0.0667 | 0.0000 |
| 1/4 | 0.0875 | 0.1270 | $-0.0395$ |
| 1/2 | 0.1014 | 0.1667 | $-0.0653$ |
| 3/4 | 0.1106 | 0.1911 | $-0.0805$ |
| 1 | 0.1458 | 0.1250 | $+0.0208$ |

**The difference is NOT monotone**: it starts at 0, goes negative (n9 dominates), then returns positive at $\beta = 1$.

**At $\beta = 0$**: Only intrinsic centrality matters → all nodes equal.

**At $\beta = 1$**: PageRank becomes the stationary distribution of the random walk. For a connected, non-bipartite, undirected graph:

$$\pi_i = \frac{\deg(i)}{2|E|} = \frac{\deg(i)}{48}$$

So $\text{PR}(n_6) = 7/48 = 0.1458$ and $\text{PR}(n_9) = 6/48 = 0.1250$.

**Figure**: `hw1/figures/exercise2_pagerank_difference_beta.png`

> **Common mistakes**:
> - Claiming the difference is monotone.
> - Forgetting that $\beta = 1$ requires checking the graph is connected and non-bipartite for convergence.

> **Oral defense questions**:
> 1. *Why does $n_9$ beat $n_6$ at $\beta = 0.15$ but not at $\beta = 1$?* — At intermediate $\beta$, degree normalization gives $n_9$ concentrated leaf mass. At $\beta = 1$, only degree matters, and $\deg(n_6) = 7 > 6 = \deg(n_9)$.
> 2. *What does "distributed" mean?* — Each node updates using only messages from neighbors: $x_j/\deg(j)$.

### 3.4 HW1 Exercise 3 — Traffic Network

#### 3.4.1 Network and Data

The traffic network has **17 nodes** and **28 directed links**. Data loaded from `.mat` files:
- `data/flow.mat` → measured flow vector $f$ (28 edges)
- `data/capacities.mat` → capacity vector $C$ (28 edges)
- `data/traffic.mat` → incidence matrix $B$ (17×28)
- `data/traveltime.mat` → free-flow time vector $l$ (28 edges)

**Incidence convention**: $B_{v,e} = +1$ at the tail, $B_{v,e} = -1$ at the head.

**Delay function**: $\tau_e(f_e) = \frac{l_e}{1 - f_e/C_e}$, for $0 \leq f_e < C_e$.

**Numerical safeguard**: Constrain $f_e \leq 0.999 \cdot C_e$ to avoid infinite delay.

**Code reference**: `hw1/src/exercise3.py`, function `load_data()` (line 42), `delay()` (line 146), `tau_prime()` (line 150).

#### 3.4.2 Part (a): Shortest Path

**Method**: Use free-flow times $l_e$ as edge weights. Find shortest path from node 1 to node 17.

**Result**: Path: $1 \to 2 \to 3 \to 9 \to 13 \to 17$, link indices: 1, 2, 12, 9, 25, total free-flow time: **0.559833**.

**Code**: `shortest_path_results()` uses NetworkX's `shortest_path()` with `weight="weight"`.

#### 3.4.3 Part (b): Maximum Flow

**Result**: Max flow from node 1 to node 17 = **22,448**.

This ignores congestion delays entirely — it's a pure capacity problem.

#### 3.4.4 Part (c): Exogenous Inflow Vector $\nu$

Compute $\nu = Bf$ from the measured flow. Then define the reduced demand:

$$\nu_{\text{new},1} = 16806, \quad \nu_{\text{new},17} = -16806, \quad \text{all others } = 0$$

This models a single origin-destination (OD) demand from node 1 to node 17.

#### 3.4.5 Part (d): Social Optimum

**Objective**: Minimize total travel time $\sum_e f_e \cdot \tau_e(f_e)$.

**CVXPY formulation**: The equivalent convex expression is:

$$\sum_e \left(\frac{l_e C_e}{1 - f_e/C_e} - l_e C_e\right)$$

This is algebraically equal to $\sum_e f_e \tau_e(f_e)$ (multiply numerator and denominator).

**Constraints**: $Bf = \nu_{\text{new}}$, $0 \leq f \leq 0.999 C$.

**Result**: Social optimum total travel time = **26,142.669150**. Solver status: optimal.

**Code**: `solve_social_optimum()` in `exercise3.py` (line 236).

#### 3.4.6 Part (e): Wardrop Equilibrium

**Objective**: Minimize the Beckmann potential:

$$\sum_e \int_0^{f_e} \tau_e(s)\, ds = \sum_e \left(-l_e C_e \ln(1 - f_e/C_e)\right)$$

**Derivation of the integral**: $\int_0^{f_e} \frac{l_e}{1 - s/C_e}\, ds = -l_e C_e \ln(1 - f_e/C_e)$ (substitute $u = 1 - s/C_e$).

**Result**:
- Wardrop total travel time = **26,495.321759**
- Price of anarchy = $26495.32 / 26142.67 = \mathbf{1.013490}$

Selfish routing costs about **1.35%** more than the social optimum.

**Code**: `solve_wardrop()` and `beckmann_expression()` in `exercise3.py`.

#### 3.4.7 Part (f): Marginal-Cost Tolls

**Theory**: The marginal social cost of link $e$ is $\tau_e(f_e) + f_e \cdot \tau'_e(f_e)$. Users perceive only $\tau_e(f_e)$. The toll charges the externality:

$$\omega_e = f_e^* \cdot \tau'_e(f_e^*), \quad \text{where } \tau'_e(f_e) = \frac{l_e/C_e}{(1 - f_e/C_e)^2}$$

**Verification**: The tolled Wardrop flow should match the social optimum.

**Result**: $\|f_{\text{tolled}} - f^*\|_2 = 2.994 \times 10^{-1}$, relative error $1.326 \times 10^{-5}$.

The error is tiny relative to flows of thousands — the tolls work.

**Code**: Toll computation at line 503: `omega = f_star * tau_prime(f_star, ...)`.

#### 3.4.8 Part (g): Additional-Delay Objective

**Objective**: $\sum_e f_e(\tau_e(f_e) - l_e)$ — penalize only congestion delay above free-flow time.

**Result**: Additional-delay optimum = **15,350.353857**.

**Toll formula**: $\omega_e^{\text{add}} = -l_e + f_e^* \cdot \tau'_e(f_e^*)$

**Negative tolls** can occur because the formula subtracts $l_e$. This is mathematically correct — it can be interpreted as a subsidy.

**Verification**: $\|f_{\text{tolled,add}} - f^*_{\text{add}}\|_2 = 5.659 \times 10^{-1}$, relative error $2.507 \times 10^{-5}$.

**Figure**: `hw1/figures/exercise3_flow_comparison.png` — social, Wardrop, and tolled flows compared.

> **Oral defense questions for Exercise 3**:
> 1. *What is $Bf = \nu$?* — Flow conservation with exogenous inflow/outflow at nodes.
> 2. *Why use 0.999C?* — Numerical safeguard to keep delay finite.
> 3. *What is the Beckmann potential?* — Integral of delay; its minimizer gives Wardrop equilibrium.
> 4. *Why can tolls be negative?* — The additional-delay formula includes $-l_e$.
> 5. *How do you verify tolls worked?* — Compare target flow with tolled Wardrop flow via L2 norm.

### 3.5 HW1 Codebase Walkthrough

| File | Purpose | Key Functions | Theory Link |
|------|---------|---------------|-------------|
| `exercise1.py` | Max-flow, cuts, capacity augmentation | `enumerate_cuts()`, `cut_capacity()`, `max_flow_value()`, `optimize_extra_capacity()`, `optimize_added_link_case()` | Max-flow/min-cut theorem |
| `exercise2.py` | Katz centrality, PageRank | `katz_centrality()`, `distributed_pagerank()`, `beta_sensitivity()` | Resolvent matrix, random walk |
| `exercise3.py` | Traffic assignment | `load_data()`, `solve_social_optimum()`, `solve_wardrop()`, `solve_additional_delay()` | Beckmann potential, marginal-cost pricing |
| `utils.py` | Shared utilities | — | — |

### 3.6 HW1 Oral Defense Checklist

- [ ] I can enumerate all 8 cuts and explain the directed cut convention.
- [ ] I can state and explain the max-flow/min-cut theorem.
- [ ] I can explain why throughput grows as $\lceil x/2 \rceil$ in part (b).
- [ ] I can defend the $o \to d$ link in part (c).
- [ ] I can state the Katz formula and the spectral condition.
- [ ] I can explain why Katz favors $n_6$ and PageRank favors $n_9$.
- [ ] I can explain what happens at $\beta = 0$ and $\beta = 1$ for PageRank.
- [ ] I can derive the Beckmann integral for the given delay function.
- [ ] I can explain the price of anarchy and what it measures.
- [ ] I can derive and explain marginal-cost tolls.
- [ ] I can explain why additional-delay tolls can be negative.

---

## 4. Homework 2 — Full Study Guide

### 4.1 What Homework 2 Is About

Homework 2 studies stochastic dynamics on a 5-node directed weighted network:

1. **Problem 1 (a–d)**: Single-particle continuous-time random walk — return time and hitting time.
2. **Problem 1 (e–h)**: French-DeGroot opinion dynamics — consensus, variance, edge removal.
3. **Problem 2**: Many-particle random walk — particle and node perspectives.
4. **Problem 3**: Open network — proportional and fixed service rates, stability.

**Network**: $V = \{o, a, b, c, d\}$, index mapping $o \to 0, a \to 1, b \to 2, c \to 3, d \to 4$.

**Transition rate matrix** (rows are sources, columns are destinations):

$$\Lambda = \begin{pmatrix} 0 & 2/5 & 1/5 & 0 & 0 \\ 0 & 0 & 3/4 & 1/4 & 0 \\ 1/2 & 0 & 0 & 1/3 & 0 \\ 0 & 0 & 1/3 & 0 & 2/3 \\ 0 & 1/3 & 0 & 1/3 & 0 \end{pmatrix}$$

**Exit rates**: $\omega = \Lambda \mathbf{1} = (3/5, 1, 5/6, 1, 2/3)$.

**Code reference**: `hw2/src/constants.py` defines `LAMBDA` and `NODES`.

### 4.2 HW2 Problem 1 — Continuous-Time Random Walk

#### Parts (a)–(b): Return Time to Node $a$

**Simulation** (100,000 MC runs, seed 20260615):
- Mean return time: **6.69375** (95% CI: [6.663, 6.724])

**Theory**: $\mathbb{E}_a[T_a^+] = 1/(\pi_a \cdot \omega_a)$

The stationary distribution $\pi$ satisfies $\pi^T Q = 0$, $\sum \pi_i = 1$:

$$\pi = (0.21739, 0.14907, 0.26087, 0.18634, 0.18634)$$

With $\omega_a = 1$: $\mathbb{E}_a[T_a^+] = 1/(0.14907 \times 1) = \mathbf{6.70833}$

**Agreement**: Simulation within 0.22% of theory. ✓

**Convention**: The initial holding time in node $a$ IS included. The particle starts at $a$, waits an Exp($\omega_a$) time, leaves, and returns — the full cycle time is measured.

**Why $\omega_a$ appears**: The formula $1/(\pi_a \omega_a)$ measures the mean time between *entrances* to $a$. During one visit, the particle spends mean time $1/\omega_a$ in $a$. The fraction of time in $a$ at stationarity is $\pi_a = (1/\omega_a) / \mathbb{E}_a[T_a^+]$, giving the formula.

**Code**: `problem1_single_ctmc.py`, functions `simulate_return_times()` and `theoretical_return_time()`.

#### Parts (c)–(d): Hitting Time from $o$ to $d$

**Simulation** (100,000 MC runs, seed 20260616):
- Mean hitting time: **10.76084** (95% CI: [10.705, 10.816])

**Theory**: Solve $\sum_j Q_{ij} h_j = -1$ for $i \neq d$, with $h_d = 0$.

**Theoretical value**: $\mathbb{E}_o[T_d] = \mathbf{10.76667}$

**Agreement**: Within 0.054% of theory. ✓

**Code**: `problem1_single_ctmc.py`, functions `simulate_hitting_times()` and `theoretical_hitting_times()`.

### 4.3 HW2 Problem 1 — French-DeGroot Dynamics

#### Part (e): Original Graph Consensus

**Dynamics**: $dx/dt = -Lx$, $L = \text{diag}(\omega) - \Lambda$.

**Graph analysis**: The original graph has **one SCC** containing all 5 nodes. Since this single SCC is necessarily the unique sink, **global consensus holds for every initial condition**.

**Consensus vector**: $\pi = (0.21739, 0.14907, 0.26087, 0.18634, 0.18634)$

**Consensus value**: $x^* = \sum_i \pi_i x_i(0)$.

**Figure**: `hw2/figures/problem1_fdg_original_trajectory.png` — all trajectories converge to a single horizontal line.

**Code**: `problem1_opinion.py`, functions `laplacian()`, `simulate_fdg()`, `strongly_connected_components()`, `sink_components()`.

#### Part (f): Variance of Consensus Value

Given independent initial opinions with variances: $\text{Var}(x_o(0)) = \text{Var}(x_d(0)) = 1$, $\text{Var}(x_a(0)) = \text{Var}(x_b(0)) = \text{Var}(x_c(0)) = 2$.

$$\text{Var}(x^*) = \sum_i \pi_i^2 \cdot \text{Var}(x_i(0))$$

**Theoretical value**: **0.331970**

**Monte Carlo value** (50,000 runs): 0.331278 — within 0.21% of theory. ✓

**Figure**: `hw2/figures/problem1_consensus_value_histogram.png`

#### Part (g): Edge Removal — $\{(d,a), (d,c), (a,c), (b,c)\}$

After removing these 4 directed edges from $\Lambda$:

**New SCCs**: $\{o, a, b\}$, $\{c\}$, $\{d\}$.

**Sink components**: $\{o, a, b\}$ and $\{d\}$ — **two** sink components.

**Result**: Global consensus does NOT hold for every initial condition. The limit depends on which sink component each node can reach and on the initial conditions.

- Nodes $o, a, b$ converge among themselves.
- Node $d$ keeps its initial value (it's a singleton sink with no outgoing edges after removal).
- Node $c$ is transient and converges to a weighted combination: $x_c^\infty = (1/3) z_{oab} + (2/3) x_d(0) = 3.3821$ for the tested initial condition.

**Figure**: `hw2/figures/problem1_fdg_removed_edges_g.png`

#### Part (h): Edge Removal — $\{(b,o), (d,a)\}$

After removing these 2 directed edges:

**New SCCs**: $\{o\}$, $\{a\}$, $\{b, c, d\}$.

**Sink component**: $\{b, c, d\}$ — **one** sink reachable from all nodes.

**Result**: Global consensus DOES hold for every initial condition. The consensus value depends only on the sink component's influence weights: $\pi^{(h)} = (0, 0, 1/4, 1/4, 1/2)$.

- Case 1: $x(0) = (1, 2, 3, 4, 2)$ → consensus value = $0 + 0 + 3/4 + 1 + 1 = 2.75$ ✓
- Case 2: $x(0) = (1, 2, 0, 0, 0)$ → consensus value = $0$ ✓

**Figures**: `hw2/figures/problem1_fdg_removed_edges_h_case1.png`, `problem1_fdg_removed_edges_h_case2.png`

> **What to remember**: Edge removal changes SCCs. Multiple sinks → no global consensus. Unique reachable sink → global consensus, but only sink nodes have positive influence weight.

### 4.4 HW2 Problem 2 — Many Particles

#### Part (a): Particle Perspective

**Setup**: $N = 100$ particles, all starting at node $a$. Simulate each particle's independent CTMC and measure return time to $a$.

**Result** (2,000 MC batches, seed 20260624):
- Mean of run-averages: **6.70478**
- Problem 1 theoretical value: 6.70833
- Difference: $-0.00355$ ✓

**Why it matches**: Each particle follows the same CTMC, so the expected return time is the same as in Problem 1. Using $N = 100$ particles reduces the variance of the estimate but does not change the expectation.

#### Part (b): Node Perspective

**Setup**: $N = 100$ particles, simulate the joint count process for $T = 60$ time units (5,000 MC runs, seed 20260625).

**Node departure rate**: Node $i$ forwards particles at rate $N_i(t) \cdot \omega_i$.

**Final average counts vs. stationary expected counts ($N \cdot \pi$)**:

| Node | Simulated | Expected ($100\pi_i$) | Difference |
|------|-----------|----------------------|------------|
| $o$ | 21.725 | 21.739 | $-0.014$ |
| $a$ | 14.959 | 14.907 | $+0.052$ |
| $b$ | 25.991 | 26.087 | $-0.096$ |
| $c$ | 18.647 | 18.634 | $+0.013$ |
| $d$ | 18.678 | 18.634 | $+0.044$ |

All differences are small — the counts match $N\pi$ as expected. ✓

**Figures**: `hw2/figures/problem2_node_counts_over_time.png`, `problem2_stationary_distribution_comparison.png`

**Code**: `problem2_many_particles.py`, functions `simulate_node_count_run()` and `simulate_node_count_runs()`.

### 4.5 HW2 Problem 3 — Open Network

#### Open Network Setup

**Routing matrix**:

$$\Lambda_{\text{open}} = \begin{pmatrix} 0 & 1 & 1 & 0 & 0 \\ 0 & 0 & 1/4 & 1/4 & 1/2 \\ 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}$$

**Service rates**: $\omega = (2, 1, 1, 1, 7/4)$ (with assignment override $\omega_d = 7/4$).

Particles arrive at node $o$ via Poisson process with rate $\lambda$. At node $d$, a service event removes a particle.

#### Part (a): Proportional-Rate Service

**Rule**: Service rate at node $i$ = $\omega_i \cdot N_i(t)$.

**Main simulation**: $\lambda = 100$, $T = 60$, seed 20260904.

**Stability**: The proportional-rate system is stable for every finite $\lambda$. Each arriving particle has a finite expected lifetime before reaching $d$ and exiting. The network acts like a set of infinite servers — capacity scales with demand.

**Stability scan results** (lambdas: 10, 50, 100, 200, 500): All classified as stable.

**Figure**: `hw2/figures/problem3_proportional_lambda100_counts.png`, `problem3_proportional_stability_scan.png`

#### Part (b): Fixed-Rate Service

**Rule**: Service clock rate at node $i$ = $\omega_i$ (constant, not proportional to queue size).

**Theoretical analysis**: The routing loads per unit $\lambda$ are $(1, 1/2, 5/8, 3/4, 1)$.

Fixed-rate capacity at node $i$ is $\omega_i / \rho_i$:

| Node | $\omega_i$ | Load $\rho_i$ | Capacity $\omega_i/\rho_i$ |
|------|------------|---------------|---------------------------|
| $o$ | 2 | 1 | 2.000 |
| $a$ | 1 | 1/2 | 2.000 |
| $b$ | 1 | 5/8 | 1.600 |
| $c$ | 1 | 3/4 | **1.333** |
| $d$ | 7/4 | 1 | 1.750 |

**Bottleneck**: Node $c$, capacity $4/3 \approx 1.333$.

**Maximum stable input rate**: $\lambda_{\max} = 4/3$. Above this, queues grow without bound.

**Main simulation**: $\lambda = 2 > 4/3$, $T = 6000$. Final counts: $(233, 15, 1240, 1577, 6)$, total 3071. This confirms blow-up: node $c$ accumulates particles.

**Stability scan** confirms: stable below $4/3$, blow-up at and above $4/3$.

**Figures**: `hw2/figures/problem3_fixed_lambda2_counts.png`, `problem3_fixed_stability_scan.png`

> **Common mistakes**:
> - Confusing proportional and fixed service rates.
> - Claiming fixed-rate stability from one simulation run instead of theoretical load analysis.
> - Computing loads incorrectly (must account for routing probabilities).

### 4.6 HW2 Codebase Walkthrough

| File | Purpose | Key Functions |
|------|---------|---------------|
| `constants.py` | Shared constants | `LAMBDA`, `LAMBDA_OPEN`, `NODES`, paths |
| `problem1_single_ctmc.py` | CTMC primitives | `exit_rates()`, `generator()`, `jump_matrix()`, `stationary_distribution()`, `simulate_return_times()`, `theoretical_return_time()`, `simulate_hitting_times()`, `theoretical_hitting_times()` |
| `problem1_opinion.py` | French-DeGroot | `laplacian()`, `left_consensus_vector()`, `simulate_fdg()`, `remove_edges()`, `strongly_connected_components()`, `sink_components()`, `asymptotic_state()` |
| `problem2_many_particles.py` | Many-particle system | `simulate_particle_return_run()`, `simulate_node_count_run()`, `simulate_node_count_runs()` |
| `problem3_open_network.py` | Open network | `simulate_open_network()`, `simulate_replicates()`, `theoretical_loads()`, `fixed_rate_capacity()` |
| `run_problem1.py` | Problem 1 main | Calls all P1 computations, writes CSVs and figures |
| `run_problem2.py` | Problem 2 main | N=100, 2000 particle batches, 5000 node runs |
| `run_problem3.py` | Problem 3 main | Proportional and fixed simulations + scans |
| `utils.py` | CSV writers | `write_csv_row()`, `write_csv_rows()`, CI computation |
| `plotting.py` | Figure saving | `save_figure()` |

### 4.7 HW2 Oral Defense Checklist

- [ ] I can distinguish $\Lambda$, $Q$, $\omega$, and $P^{\text{jump}}$.
- [ ] I can derive the return time formula $1/(\pi_i \omega_i)$.
- [ ] I can set up and solve the hitting time linear system.
- [ ] I can explain consensus via SCCs and sink components.
- [ ] I can compute the variance of the consensus value.
- [ ] I can explain why edge removal (g) destroys consensus and (h) preserves it.
- [ ] I can explain particle vs. node perspectives.
- [ ] I can explain why proportional service doesn't blow up.
- [ ] I can derive the fixed-rate bottleneck as $\lambda_{\max} = 4/3$ at node $c$.

---

## 5. Homework 3 — Study Guide and Implementation Roadmap

> **Important**: HW3 is not yet implemented. No code, results, or figures exist in the repository. This section provides a **complete theoretical guide and implementation roadmap** based on expected HW3 topics. No numerical results are invented.

### 5.1 What Homework 3 Is About

HW3 is expected to cover:

1. **Part 1.1**: SIR epidemic simulation on a known $k$-regular graph.
2. **Part 1.2**: Preferential attachment random graph generation.
3. **Problem 2**: Pandemic simulation on preferential attachment graph (no vaccination).
4. **Problem 3**: Pandemic simulation with vaccination.
5. **Problem 4**: Parameter estimation for the H1N1 pandemic in Sweden.
6. **Optional**: Alternative random graph models.
7. **Network Games**: Coordination/anti-coordination games, Nash equilibria, best-response dynamics.

### 5.2 HW3 Part 1.1 — Epidemic on Known Graph

#### What Is Asked

Simulate an SIR epidemic on a **symmetric $k$-regular graph** with $n = 500$ nodes, $k = 4$.

**Parameters**:
- $\beta$ = infection probability per infected neighbor per time step
- $\rho$ = recovery probability per time step
- Infection probability with $m$ infected neighbors: $1 - (1-\beta)^m$
- Simulation over **15 weeks** (time steps)
- Average over $N = 100$ runs

**Expected outputs**: Plot of average number of susceptible, infected, and recovered nodes over time.

#### Mathematical Foundation

In a $k$-regular graph, every node has exactly $k$ neighbors. This creates a relatively homogeneous spreading environment.

The infection probability formula $1 - (1-\beta)^m$ comes from assuming each infected neighbor independently transmits with probability $\beta$. The probability of *not* being infected by any of $m$ infected neighbors is $(1-\beta)^m$, so the complement is the infection probability.

#### Implementation Roadmap

**TODO**: Create `hw3/src/sir_simulation.py`:
```python
def simulate_sir(graph, beta, rho, initial_infected, T, N_runs):
    """
    For each run:
    1. Initialize: one random node infected, rest susceptible.
    2. For t = 1, ..., T:
       a. For each susceptible node: count infected neighbors m,
          infect with prob 1 - (1-beta)^m.
       b. For each infected node: recover with prob rho.
       c. Record S(t), I(t), R(t).
    3. Average over N_runs.
    """
```

**TODO**: Create `hw3/src/graph_generation.py`:
```python
def k_regular_graph(n, k):
    """Generate a symmetric k-regular graph with n nodes."""
    # Use networkx.random_regular_graph(k, n)
```

### 5.3 HW3 Part 1.2 — Preferential Attachment Graph

#### What Is Asked

Generate a random graph using **preferential attachment**:
1. Start with a complete graph on $k + 1$ nodes.
2. Add nodes one at a time, each connecting to existing nodes with probability proportional to their current degree.
3. Each new node adds a specified number of edges.
4. Avoid duplicate edges.

The number of added links per new node should be chosen so the average degree is approximately $k$.

#### Mathematical Foundation

Preferential attachment creates **scale-free networks**: the degree distribution follows a power law $P(k) \sim k^{-\gamma}$, typically with $\gamma \approx 3$. This means a few **hubs** (high-degree nodes) and many low-degree nodes.

The probability of connecting to existing node $v$ is $\deg(v) / \sum_u \deg(u)$. This "rich get richer" mechanism creates hubs that are important for epidemic spread.

**Even/odd $k$ handling**: If $k$ is even, each new node adds $k/2$ edges. If $k$ is odd, alternating new nodes add $\lfloor k/2 \rfloor$ and $\lceil k/2 \rceil$ edges.

#### Implementation Roadmap

**TODO**: Create preferential attachment in `hw3/src/graph_generation.py`:
```python
def preferential_attachment_graph(n, k):
    """
    1. Start with K_{k+1} (complete graph).
    2. For each new node, attach to existing nodes with
       prob proportional to degree.
    3. Handle even/odd k for average degree.
    """
```

### 5.4 HW3 Problem 2 — Pandemic Without Vaccination

#### What Is Asked

Simulate SIR on a **preferential attachment graph** with $n = 500$, average degree $k = 6$.

**Expected differences from $k$-regular graph**:
- Hubs act as super-spreaders: they get infected early and infect many neighbors.
- Epidemics spread faster and may reach more of the population.
- Higher variance across runs due to heterogeneous degree distribution.

**TODO**: Run SIR simulation on preferential attachment graph and compare with $k$-regular results.

### 5.5 HW3 Problem 3 — Pandemic With Vaccination

#### What Is Asked

Introduce a vaccination schedule:
- At specified times, a certain number of nodes are vaccinated.
- Vaccinated nodes cannot infect or be infected.
- Vaccination is random among the non-vaccinated population.
- Infected individuals CAN be vaccinated (removes them from the epidemic).

**Expected impact**:
- Vaccination reduces the susceptible pool.
- Early vaccination is more effective than late vaccination.
- The epidemic curve (infected count over time) should flatten and peak lower.

**TODO**: Extend SIR simulation with vaccination logic.

### 5.6 HW3 Problem 4 — H1N1 Sweden Parameter Estimation

#### What Is Asked

Fit epidemic parameters to real H1N1 data from Sweden:
- Population scaled to $n = 934$.
- Given: vector of newly infected per week $I_0(t)$.
- Given: vaccination schedule.
- Search over grid of $(k, \beta, \rho)$ to minimize RMSE between simulated and observed $I_0(t)$.

$$\text{RMSE}(k, \beta, \rho) = \sqrt{\frac{1}{T}\sum_{t=1}^T \left(\hat{I}(t) - I_0(t)\right)^2}$$

**Implementation approach**:
1. Define a grid over $(k, \beta, \rho)$.
2. For each parameter triple, run $N$ simulations, compute average newly infected per week, compute RMSE.
3. Select the triple with minimum RMSE.
4. Optionally refine with smaller grid around the minimum.

**Computational cost**: This is expensive. If $k$ has 5 values, $\beta$ has 10 values, $\rho$ has 10 values, and $N = 100$ runs per triple, that's $5 \times 10 \times 10 \times 100 = 50,000$ simulations.

**TODO**: Create `hw3/src/parameter_search.py` with grid search and RMSE computation.

### 5.7 HW3 Optional Challenge

**Alternative random graph models**:

1. **Erdős–Rényi** $G(n, p)$: Each edge exists independently with probability $p$. Homogeneous, no hubs. Degree distribution is Binomial (approximately Poisson for large $n$).

2. **Configuration model**: Specify a degree sequence, then connect half-edges randomly. Can reproduce any desired degree distribution.

3. **Small-world networks** (Watts-Strogatz): Start with a regular lattice, rewire edges with probability $p$. Creates high clustering with short path lengths.

4. **Stochastic block models**: Nodes belong to communities; edge probability depends on community membership. Good for modeling social contact networks with group structure.

**Why preferential attachment may not be the best model for social contacts**: Real social networks have high clustering (friends of friends are often friends), which preferential attachment does not naturally produce. Small-world or stochastic block models may be more realistic.

### 5.8 HW3 Network Games and Dynamics

#### Setup

- $n$ players on a graph (e.g., $n = 3$ on a complete graph $K_3$).
- Action set: $\{-1, +1\}$.
- $n_1$ players are **coordinators**, $n - n_1$ are **anti-coordinators**.

#### Utility Functions

For a coordinator $i$:
$$u_i(a_i, a_{-i}) = \sum_{j \in N(i)} \mathbb{1}[a_i = a_j]$$

For an anti-coordinator $i$:
$$u_i(a_i, a_{-i}) = \sum_{j \in N(i)} \mathbb{1}[a_i \neq a_j]$$

#### Nash Equilibria

An action profile $a = (a_1, \ldots, a_n)$ is a Nash equilibrium if no player can improve by unilaterally changing their action. For $n = 3$ on $K_3$:

- **$n_1 = 3$ (all coordinators)**: NE are $(+1,+1,+1)$ and $(-1,-1,-1)$.
- **$n_1 = 0$ (all anti-coordinators)**: No pure NE on $K_3$ (each player wants to differ from both neighbors, impossible for 3 players with binary actions).
- **Mixed cases ($n_1 = 1, 2$)**: Depends on which nodes are coordinators.

**Total action profiles**: $2^n = 8$ for $n = 3$.

#### Best-Response Dynamics

**Asynchronous**: At each step, one random player updates to their best response.

**Noisy best response**: Player $i$ chooses action $+1$ with probability:
$$\frac{e^{u_i(+1) / \epsilon}}{e^{u_i(+1) / \epsilon} + e^{u_i(-1) / \epsilon}}$$

where $\epsilon > 0$ is the noise parameter.

**Vanishing noise limit**: As $\epsilon \to 0$, the distribution concentrates on specific equilibria determined by the stochastically stable states.

#### Transition Graph

The $2^n$ action profiles form nodes. Under best-response dynamics, edges connect profiles that differ by one player's action change. Analysis of this graph identifies:
- Absorbing states (Nash equilibria for deterministic best response).
- Recurrent classes.
- Transient states.

**TODO**: Implement `hw3/src/network_games.py`:
```python
def enumerate_nash_equilibria(n, n1, graph):
    """Check all 2^n profiles for NE condition."""

def best_response_dynamics(graph, n1, initial_profile, steps):
    """Asynchronous best-response simulation."""

def noisy_best_response(graph, n1, initial_profile, epsilon, steps):
    """Logit best-response simulation."""

def transition_graph(n, n1, graph):
    """Build transition graph over all 2^n profiles."""
```

### 5.9 HW3 Codebase Walkthrough

**TODO — files to create**:

| File | Purpose |
|------|---------|
| `hw3/src/graph_generation.py` | $k$-regular and preferential attachment graphs |
| `hw3/src/sir_simulation.py` | SIR simulation engine |
| `hw3/src/vaccination.py` | Vaccination logic |
| `hw3/src/parameter_search.py` | Grid search for H1N1 parameter estimation |
| `hw3/src/network_games.py` | Nash equilibria, best response, transition graphs |
| `hw3/src/plotting.py` | Plotting utilities |
| `hw3/src/run_problem1.py` | Main script for epidemic on known graph |
| `hw3/src/run_problem2.py` | Main script for pandemic without vaccination |
| `hw3/src/run_problem3.py` | Main script for pandemic with vaccination |
| `hw3/src/run_problem4.py` | Main script for parameter estimation |
| `hw3/src/run_games.py` | Main script for network games |

### 5.10 HW3 Oral Defense Checklist

- [ ] I can explain the SIR model and the infection probability formula.
- [ ] I can explain preferential attachment and why it creates hubs.
- [ ] I can explain why hubs accelerate epidemic spread.
- [ ] I can explain the vaccination mechanism and its effect on the epidemic curve.
- [ ] I can explain RMSE-based parameter estimation.
- [ ] I can define Nash equilibrium for network games.
- [ ] I can enumerate NE for small cases ($n = 3$, different $n_1$ values).
- [ ] I can explain asynchronous and noisy best-response dynamics.
- [ ] I can explain the transition graph and absorbing states.

---

## 6. Cross-Homework Connections

### Flow and Particles

**HW1 traffic flows** model vehicle routing on a network with congestion. **HW2 particle flows** model random walkers on a network. Both use directed graphs and conservation constraints ($Bf = \nu$ in HW1; departure rate proportional to occupancy in HW2). The key difference is that traffic flows are deterministic optimization problems, while particle flows are stochastic.

### PageRank and Stationary Distributions

**HW1 PageRank** at $\beta = 1$ gives the stationary distribution of a random walk on the undirected graph. **HW2 stationary distribution** $\pi$ is the long-run time fraction for a CTMC. Both describe how "mass" distributes in a network at equilibrium.

### CTMCs and Epidemic/Game Dynamics

**HW2 CTMCs** use exponential holding times and jump probabilities. **HW3 epidemic simulations** use discrete-time SIR (Bernoulli infections per step). **HW3 network games** use discrete best-response updates. The common thread is stochastic dynamics on graphs, where the network structure determines the evolution.

### Opinion Dynamics and Best-Response Dynamics

**HW2 French-DeGroot** has each node update toward a weighted average of neighbors — a form of consensus-seeking. **HW3 best-response dynamics** has each player update to optimize against neighbors. Both are local update rules on networks, but DeGroot is linear and continuous, while best response is nonlinear and discrete.

### Optimization and Parameter Estimation

**HW1 traffic assignment** solves convex optimization (CVXPY) to find social optimum and Wardrop equilibrium. **HW3 parameter estimation** solves a non-convex optimization (grid search) to minimize RMSE. Both require careful formulation of objectives and constraints.

### Graph Structure Effects

Across all homeworks, the graph structure determines outcomes:
- **Centrality** (HW1): cliques boost Katz, leaves boost PageRank.
- **Consensus** (HW2): sink SCCs determine convergence.
- **Epidemics** (HW3): hubs accelerate spread on scale-free networks.
- **Games** (HW3): graph topology determines Nash equilibria.
- **Bottlenecks** (HW1, HW2): min-cuts limit flow; fixed-rate bottlenecks limit throughput.

---

## 7. Master Formula Sheet

### Max-Flow and Cuts (HW1)

| Formula | Meaning |
|---------|---------|
| $C(S) = \sum_{e: \text{tail}(e) \in S,\, \text{head}(e) \notin S} c_e$ | Directed cut capacity |
| $\max f = \min_S C(S)$ | Max-flow / min-cut theorem |
| $0 \leq f_e \leq c_e$ | Capacity constraint |
| $T_b(x) = 5 + \lceil x/2 \rceil$ | HW1 existing-link throughput |
| $T_c(x) = 6 + x$ | HW1 direct-link throughput |

### Centrality (HW1)

| Formula | Meaning |
|---------|---------|
| $x = (I - \beta W)^{-1} \mu$ | Katz centrality |
| $x = \sum_{k=0}^\infty \beta^k W^k \mu$ | Walk expansion |
| $\beta < 1/\lambda_{\max}(W)$ | Spectral condition |
| $x_i^{(t+1)} = (1-\beta)\mu_i + \beta \sum_{j \in N(i)} x_j^{(t)} / \deg(j)$ | PageRank update |
| $\pi_i = \deg(i) / (2|E|)$ | Stationary distribution (undirected) |

### Traffic Assignment (HW1)

| Formula | Meaning |
|---------|---------|
| $Bf = \nu$ | Flow conservation |
| $\tau_e(f_e) = l_e / (1 - f_e/C_e)$ | Delay function |
| $\min \sum_e f_e \tau_e(f_e)$ | Social optimum |
| $\min \sum_e \int_0^{f_e} \tau_e(s)\, ds$ | Beckmann potential (Wardrop) |
| $\int_0^{f_e} \tau_e(s)\, ds = -l_e C_e \ln(1 - f_e/C_e)$ | Beckmann integral |
| $\text{PoA} = \text{cost}_{\text{Wardrop}} / \text{cost}_{\text{SO}}$ | Price of anarchy |
| $\omega_e = f_e^* \cdot \tau'_e(f_e^*)$ | Marginal-cost toll |
| $\tau'_e(f_e) = (l_e/C_e) / (1 - f_e/C_e)^2$ | Delay derivative |
| $\omega_e^{\text{add}} = -l_e + f_e^* \tau'_e(f_e^*)$ | Additional-delay toll |

### CTMC (HW2)

| Formula | Meaning |
|---------|---------|
| $\omega_i = \sum_j \Lambda_{ij}$ | Exit rate |
| $Q = \Lambda - \text{diag}(\omega)$ | CTMC generator |
| $P^{\text{jump}}_{ij} = \Lambda_{ij} / \omega_i$ | Embedded jump probability |
| $\Delta t \sim \text{Exp}(\omega_i)$ | Holding time distribution |
| $\pi^T Q = 0$, $\sum \pi_i = 1$ | Stationary distribution |
| $\mathbb{E}_i[T_i^+] = 1 / (\pi_i \omega_i)$ | Expected return time |
| $\sum_j Q_{ij} h_j = -1$ $(i \neq d)$, $h_d = 0$ | Hitting time equations |

### Opinion Dynamics (HW2)

| Formula | Meaning |
|---------|---------|
| $L = \text{diag}(\omega) - \Lambda$ | Graph Laplacian |
| $dx/dt = -Lx$ | French-DeGroot dynamics |
| $x(t) = e^{-Lt} x(0)$ | Solution |
| $x^* = \sum_i \pi_i x_i(0)$ | Consensus value |
| $\text{Var}(x^*) = \sum_i \pi_i^2 \text{Var}(x_i(0))$ | Variance of consensus |

### Epidemics (HW3)

| Formula | Meaning |
|---------|---------|
| $p_{\text{infect}} = 1 - (1-\beta)^m$ | Infection prob with $m$ infected neighbors |
| $p_{\text{recover}} = \rho$ | Recovery probability |
| $\text{RMSE} = \sqrt{\frac{1}{T}\sum_t (\hat{I}(t) - I_0(t))^2}$ | Parameter estimation objective |

### Network Games (HW3)

| Formula | Meaning |
|---------|---------|
| $u_i^{\text{coord}} = \sum_{j \in N(i)} \mathbb{1}[a_i = a_j]$ | Coordination utility |
| $u_i^{\text{anti}} = \sum_{j \in N(i)} \mathbb{1}[a_i \neq a_j]$ | Anti-coordination utility |
| NE: $u_i(a_i^*, a_{-i}^*) \geq u_i(a_i, a_{-i}^*)$ $\forall i, a_i$ | Nash equilibrium condition |

---

## 8. Master Oral Exam Questions

### HW1 Questions (20+)

**1. What is a directed cut and how is its capacity computed?**
A directed cut $(S, V \setminus S)$ has $o \in S$ and $d \notin S$. Its capacity is the sum of capacities of edges whose tail is in $S$ and whose head is in $V \setminus S$. Edges entering $S$ do not count because they cannot carry flow toward the sink.

**2. Why are there exactly 8 cuts in Exercise 1?**
Three intermediate nodes $\{a, b, c\}$ can each be in or out of $S$, giving $2^3 = 8$ possible source-side sets.

**3. State the max-flow/min-cut theorem.**
The maximum feasible source-sink flow equals the minimum source-sink cut capacity.

**4. Why is the minimum cut 5?**
The cuts $\{o,b,c\}$ and $\{o,a,b,c\}$ both have capacity 5 (edges $e_1 + e_5 = 3+2$ and $e_2 + e_5 = 3+2$ respectively). All other cuts have capacity $\geq 6$.

**5. Why does throughput grow as $\lceil x/2 \rceil$ in part (b)?**
Extra capacity must be paired across complementary bottlenecks. Adding one unit to one side of a limiting cut doesn't help if the other side still limits.

**6. Why is $o \to d$ optimal in part (c)?**
It crosses every source-sink cut. Every unit of capacity on this edge directly increases the minimum cut by 1.

**7. What is the spectral condition for Katz centrality?**
$\beta < 1/\lambda_{\max}(W)$. This ensures the resolvent $(I - \beta W)^{-1}$ exists and the walk expansion converges.

**8. Why does Katz favor $n_6$?**
$n_6$ participates in many walks through the dense 6-node clique and bridges to the rest of the graph. Katz counts raw walks without normalization.

**9. Why does PageRank favor $n_9$ at $\beta = 0.15$?**
Its degree-1 leaves ($n_{10}$–$n_{14}$) send ALL their mass to $n_9$. Clique nodes split their mass among many neighbors.

**10. What happens at $\beta = 0$ for PageRank?**
Centrality equals intrinsic $\mu$. For uniform $\mu$, all nodes are equal.

**11. What happens at $\beta = 1$ for PageRank on the Exercise 2 graph?**
PageRank becomes the stationary distribution of the random walk: $\pi_i = \deg(i)/(2|E|)$. $n_6$ has degree 7, $n_9$ has degree 6, so $n_6$ is higher.

**12. Is $\text{PR}(n_6) - \text{PR}(n_9)$ monotone in $\beta$?**
No. It is 0 at $\beta=0$, negative for intermediate $\beta$, and positive at $\beta=1$.

**13. What does $Bf = \nu$ represent?**
Flow conservation with exogenous inflow/outflow. Positive $\nu_i$ = source, negative = sink.

**14. What is the delay function and why does it go to infinity?**
$\tau_e(f_e) = l_e/(1 - f_e/C_e)$. As $f_e \to C_e$, the denominator approaches 0, so delay $\to \infty$.

**15. Why impose $f_e \leq 0.999 C_e$?**
Numerical safeguard to keep the solver inside the finite domain of the delay function.

**16. What is the social optimum?**
The feasible flow minimizing total travel time $\sum_e f_e \tau_e(f_e)$.

**17. What is the Wardrop equilibrium?**
A selfish-routing state where all used routes have equal travel time. No driver can improve by switching.

**18. How is the Beckmann potential derived?**
$\int_0^{f_e} l_e/(1-s/C_e)\, ds = -l_e C_e \ln(1-f_e/C_e)$. Its KKT conditions match Wardrop conditions.

**19. What is the price of anarchy and what value did you get?**
PoA = Wardrop cost / social optimum cost = 1.013490. Selfish routing costs 1.35% more.

**20. What are marginal-cost tolls?**
$\omega_e = f_e^* \tau'_e(f_e^*)$. They charge the congestion externality, making perceived cost equal marginal social cost.

**21. Why can additional-delay tolls be negative?**
The formula $\omega_e^{\text{add}} = -l_e + f_e^* \tau'_e(f_e^*)$ subtracts the free-flow time. This can yield negative values, interpreted as subsidies.

### HW2 Questions (20+)

**1. What is the difference between $\Lambda$, $Q$, and $P^{\text{jump}}$?**
$\Lambda$ has nonneg transition rates. $Q = \Lambda - \text{diag}(\omega)$ is the generator (rows sum to 0, negative diagonal). $P^{\text{jump}} = \text{diag}(\omega)^{-1}\Lambda$ is the embedded jump matrix.

**2. Why are holding times exponential?**
Outgoing edges have independent Poisson clocks. The minimum of exponentials is exponential with rate = sum of rates. The memoryless property gives the Markov property.

**3. Derive the return time formula.**
At stationarity, fraction of time in state $i$ is $\pi_i$. Mean visit length is $1/\omega_i$. So mean cycle = $1/(\pi_i \omega_i)$.

**4. How do you compute hitting times?**
Solve $\sum_j Q_{ij} h_j = -1$ for $i \neq d$, with $h_d = 0$. This is a linear system from first-step analysis.

**5. When does French-DeGroot converge to consensus?**
When there is a unique sink SCC reachable from all nodes. The graph structure (not the weights alone) determines this.

**6. What determines the consensus value?**
$x^* = \sum_i \pi_i x_i(0)$, where $\pi$ is the normalized left nullvector of $L$.

**7. Why does the consensus variance formula work?**
$x^*$ is a linear combination of independent random variables. $\text{Var}(\sum a_i X_i) = \sum a_i^2 \text{Var}(X_i)$ when $X_i$ are independent.

**8. What happens after removing edges in case (g)?**
Two sink components form: $\{o,a,b\}$ and $\{d\}$. Global consensus is destroyed. Different sinks converge to different values.

**9. What happens after removing edges in case (h)?**
One sink component $\{b,c,d\}$ remains, reachable from all nodes. Global consensus holds, but only sink nodes have positive influence weight: $\pi^{(h)} = (0, 0, 1/4, 1/4, 1/2)$.

**10. Why does Problem 2 particle perspective match Problem 1?**
Each of the 100 particles is an independent CTMC copy. The expected return time is the same; more particles reduce Monte Carlo variance.

**11. Why compare node counts with $N\pi$?**
At stationarity, each particle independently has probability $\pi_i$ of being at node $i$. Expected count = $N\pi_i$.

**12. What is proportional-rate service?**
Node $i$ serves at rate $\omega_i N_i(t)$. Capacity scales with queue size. No blow-up for finite $\lambda$.

**13. What is fixed-rate service?**
Node $i$ has clock rate $\omega_i$ regardless of queue size. If a clock ticks and the node is empty, nothing happens. Bottlenecks can cause blow-up.

**14. What is the fixed-rate bottleneck and critical $\lambda$?**
Node $c$ has the lowest capacity: $\omega_c / \rho_c = 1/(3/4) = 4/3$. The system is stable for $\lambda < 4/3$.

**15. How are routing loads computed?**
Solve the traffic equations: arrivals at node $o$ split according to routing probabilities. The load at each node is the fraction of arrivals that pass through it.

**16. What does "blow-up" mean in simulation?**
Sustained growth of total population, indicated by positive long-run drift, large final totals, and growth concentrated at bottleneck nodes.

**17. Why does proportional service not blow up?**
The system acts like infinite servers. Each particle has finite expected lifetime. By Little's law, $E[N] = \lambda \cdot E[\text{lifetime}] < \infty$.

**18. What convention is used for return times?**
The initial holding time in node $a$ IS included. The simulated time starts at $a$, samples the first holding time, and stops at the next entrance to $a$.

**19. What role does non-bipartiteness play (from HW1 but applies here too)?**
For a random walk on an undirected graph, non-bipartiteness prevents period-2 oscillations and ensures convergence to the unique stationary distribution.

**20. What are the main sources of numerical error?**
Monte Carlo sampling error, finite simulation horizon, randomness near stability thresholds, and floating-point roundoff.

### HW3 Questions (20+)

**1. Explain the SIR model.**
Nodes are Susceptible, Infected, or Recovered. Susceptibles get infected by neighbors with prob $1-(1-\beta)^m$. Infected recover with prob $\rho$. Recovered are immune.

**2. Why is the infection probability $1-(1-\beta)^m$?**
Each of $m$ infected neighbors independently transmits with prob $\beta$. The prob of avoiding ALL $m$ infections is $(1-\beta)^m$. Complement = infection prob.

**3. What is a $k$-regular graph?**
Every node has exactly $k$ neighbors. Homogeneous degree distribution.

**4. What is preferential attachment?**
New nodes connect to existing nodes with probability proportional to degree. "Rich get richer" creates hubs and scale-free degree distributions.

**5. How does epidemic spread differ on regular vs. preferential attachment graphs?**
Preferential attachment has hubs that act as super-spreaders, leading to faster, more explosive epidemics compared to the more uniform spread on regular graphs.

**6. What is vaccination in the SIR context?**
Vaccinated nodes become immune: they cannot infect or be infected. This removes transmission paths.

**7. How does random vaccination affect the epidemic?**
It reduces the susceptible pool, flattening the epidemic curve. However, random vaccination may not target hubs, reducing efficiency compared to targeted vaccination.

**8. What is RMSE and how is it used for parameter estimation?**
$\text{RMSE} = \sqrt{\frac{1}{T}\sum_t (\hat{I}(t) - I_0(t))^2}$. Parameters $(k, \beta, \rho)$ are chosen to minimize RMSE between simulated and observed newly infected counts.

**9. Why is grid search computationally expensive?**
Each parameter combination requires $N$ simulation runs (e.g., 100), each simulating an entire epidemic. The total number of simulations = grid points × $N$.

**10. What is a Nash equilibrium?**
An action profile where no player can improve their utility by unilaterally changing their action.

**11. How many action profiles exist for $n = 3$ players with binary actions?**
$2^3 = 8$.

**12. What is a coordination game?**
A player benefits from matching their neighbors' actions.

**13. What is an anti-coordination game?**
A player benefits from differing from their neighbors' actions.

**14. Why is there no pure NE for 3 anti-coordinators on $K_3$?**
Each player wants to differ from both neighbors. But with binary actions on a complete graph of 3 nodes, at least two players must have the same action.

**15. What is asynchronous best response?**
One random player updates to their best response at each step. Other players hold their actions fixed.

**16. What is noisy best response?**
Players choose actions with probability proportional to $\exp(u/\epsilon)$, where $\epsilon$ is a noise parameter. This introduces exploration.

**17. What happens as noise $\epsilon \to 0$?**
The dynamics converge to deterministic best response. The long-run behavior selects stochastically stable equilibria.

**18. What is the transition graph?**
A graph over all $2^n$ action profiles. Edges connect profiles that differ by one player's action under best response.

**19. What are absorbing states of the transition graph?**
Nash equilibria (no player wants to change), so no outgoing transitions.

**20. How would you model a real social contact network?**
Stochastic block models or small-world networks may be more realistic than preferential attachment, because real social networks have high clustering and community structure.

### Cross-Topic Questions (15)

**1. How do HW1 traffic flows and HW2 particle flows relate?**
Both model flow on directed networks. Traffic flows are deterministic optimization; particle flows are stochastic. Both use conservation constraints.

**2. How does PageRank relate to the CTMC stationary distribution?**
At $\beta = 1$, PageRank IS the stationary distribution of the random walk. In HW2, $\pi$ gives the CTMC stationary distribution. Both describe equilibrium mass distribution.

**3. How do min-cuts in HW1 relate to bottlenecks in HW2?**
Both identify structural limits. Min-cuts limit max-flow. Fixed-rate bottlenecks limit the maximum stable input rate. Both are properties of the graph topology.

**4. How does French-DeGroot consensus relate to best-response convergence?**
DeGroot is a linear averaging process converging to consensus (under SCC conditions). Best response is a nonlinear discrete process converging to Nash equilibria. Both are local update dynamics on networks.

**5. How is CVXPY optimization in HW1 different from RMSE grid search in HW3?**
HW1 uses convex optimization with guarantees of global optimality. HW3 uses brute-force search over a discrete parameter space with no convexity guarantee.

**6. Why does graph structure matter for all three homeworks?**
Centrality, consensus, epidemic spread, game equilibria, and bottleneck identification are all determined by the network topology.

**7. What role does stochasticity play across the homeworks?**
HW1 is deterministic (optimization). HW2 uses stochastic CTMCs with Monte Carlo validation. HW3 uses stochastic SIR epidemics and noisy best response.

**8. How do degree distributions affect results?**
In HW1, the clique boosts Katz centrality. In HW2, heterogeneous rates affect stationary distributions. In HW3, hubs in scale-free graphs accelerate epidemics.

**9. Compare the convergence guarantees: Wardrop equilibrium vs. DeGroot consensus vs. best-response dynamics.**
Wardrop: guaranteed by Beckmann convexity. DeGroot: guaranteed by unique reachable sink SCC. Best response: may cycle (anti-coordination) or converge to NE.

**10. What is the common theme of "local vs. global" across the homeworks?**
Users minimize local cost (Wardrop) vs. planner minimizes global cost (social optimum). Nodes update locally (DeGroot, best response) but global behavior emerges. Local graph properties (degree, clique membership) determine global centrality.

**11. How do simulation and theory complement each other?**
Simulation provides numerical evidence and validates theoretical predictions. Theory provides exact values, stability conditions, and convergence guarantees. Both are needed.

**12. What is the role of convex optimization vs. simulation-based optimization?**
HW1 Exercise 3 uses convex optimization (CVXPY) for guaranteed optimal solutions. HW3 uses simulation-based grid search because the epidemic model is not analytically tractable.

**13. Compare the role of the incidence matrix $B$ in HW1 and the rate matrix $\Lambda$ in HW2.**
Both encode graph structure. $B$ encodes edge directions for flow conservation. $\Lambda$ encodes transition rates for CTMC dynamics. Both are fundamental representations of directed networks.

**14. How do "externalities" appear across the homeworks?**
In HW1, congestion externalities justify tolls. In HW2, particles don't interact (independent CTMCs). In HW3, infection is an externality: being infected harms your neighbors.

**15. What mathematical skills are tested across all three homeworks?**
Linear algebra (Katz, hitting times, stationary distributions), convex optimization (traffic), probability (CTMCs, Monte Carlo), graph theory (SCCs, cuts, random graphs), and game theory (Nash equilibria).

---

## 9. Common Mistakes and How to Avoid Them

### Conceptual Mistakes

| Mistake | Why It's Wrong | How to Avoid |
|---------|---------------|--------------|
| Counting edges entering $S$ in a cut | Only edges leaving $S$ carry flow toward the sink | Remember: tail in $S$, head not in $S$ |
| Confusing $\Lambda$ with $Q$ | $\Lambda$ has nonneg entries; $Q$ has negative diagonal | Check: rows of $Q$ sum to 0, rows of $\Lambda$ sum to $\omega$ |
| Confusing $\Lambda$ with $P^{\text{jump}}$ | $\Lambda$ has rates; $P^{\text{jump}}$ has probabilities (rows sum to 1) | $P^{\text{jump}} = \text{diag}(\omega)^{-1} \Lambda$ |
| Using $P^{\text{jump}}$ as a uniformized matrix | $P^{\text{jump}}$ has no self-loops (usually); uniformized includes self-loops | Keep them separate |
| Confusing social optimum with Wardrop | Social optimum minimizes total cost; Wardrop is selfish equilibrium | Social = planner; Wardrop = users |
| Using wrong toll formula | Marginal-cost toll = $f_e^* \tau'_e(f_e^*)$, not $\tau_e(f_e^*)$ | Toll = externality, not total cost |
| Claiming consensus without graph proof | Consensus requires unique reachable sink SCC | Always check SCC structure |
| Assuming monotonicity of PR difference | The $\text{PR}(n_6) - \text{PR}(n_9)$ difference is NOT monotone | Compute for specific $\beta$ values |
| Interpreting simulation averages as exact | Simulations have Monte Carlo variance | Report confidence intervals and standard errors |
| Ignoring Monte Carlo error | Stochastic claims need many repetitions | Use sufficient runs; report CIs |
| Inventing HW3 results | No code or simulations exist yet | Mark clearly as TODO |

### Implementation Mistakes

| Mistake | Consequence | Fix |
|---------|-------------|-----|
| Wrong matrix orientation | All computations give wrong results | Verify: rows are sources, columns are destinations |
| Forgetting $0.999C$ cap | Solver crashes or gives infinite delay | Always enforce $f_e \leq 0.999 C_e$ |
| Using wrong holding time distribution | Wrong CTMC dynamics | Holding time is Exp($\omega_i$), not Exp(1) or uniform |
| Confusing node indices (0 vs 1) | Wrong node assignments | Python is 0-indexed; node 1 = index 0 |
| Not fixing random seeds | Irreproducible results | Always set seeds before simulations |

---

## 10. What I Should Be Able to Explain in 5 Minutes

### HW1 in 5 Minutes

Homework 1 solves three network problems. In **Exercise 1**, I enumerate all 8 source-sink cuts in a 5-node directed graph. The minimum cut capacity is 5, which equals the max flow by the max-flow/min-cut theorem. Capacity augmentation on existing links grows throughput as $\lceil x/2 \rceil$ because bottlenecks come in pairs. Adding a direct $o \to d$ link gives linear growth $6 + x$.

In **Exercise 2**, I compare Katz and PageRank centrality on a 15-node undirected graph. Katz uses the resolvent $(I - \beta W)^{-1}\mu$ and favors $n_6$ because of the dense 6-node clique. PageRank uses the distributed update $x_i = (1-\beta)\mu_i + \beta \sum x_j/\deg(j)$ and favors $n_9$ because degree-1 leaves send all their mass to it. At $\beta = 1$, PageRank becomes proportional to degree, reversing the ranking.

In **Exercise 3**, I use CVXPY to solve traffic assignment on a 17-node, 28-link network. The social optimum (total travel time 26,143) is better than the Wardrop equilibrium (26,495), giving a price of anarchy of 1.013. Marginal-cost tolls $\omega_e = f_e^* \tau'_e(f_e^*)$ decentralize the social optimum. Additional-delay tolls can be negative because they subtract free-flow time.

### HW2 in 5 Minutes

Homework 2 studies stochastic dynamics on a 5-node directed weighted network. In **Problem 1**, a single particle performs a continuous-time random walk governed by the rate matrix $\Lambda$. The exit rate at node $i$ is $\omega_i = \sum_j \Lambda_{ij}$. Holding times are exponential with rate $\omega_i$, and jump destinations follow $\Lambda_{ij}/\omega_i$. The expected return time to node $a$ is $1/(\pi_a \omega_a) = 6.708$, confirmed by simulation. Hitting times are solved via $\sum_j Q_{ij} h_j = -1$.

French-DeGroot dynamics $dx/dt = -Lx$ converge to consensus when there's a unique reachable sink SCC. The original graph has one SCC → consensus. Removing edges $\{(d,a),(d,c),(a,c),(b,c)\}$ creates two sinks → no global consensus. Removing $\{(b,o),(d,a)\}$ leaves one reachable sink $\{b,c,d\}$ → consensus, but only sink nodes influence the limit.

In **Problem 2**, 100 independent particles reproduce the same return time. Node counts at $T=60$ match $N\pi$. In **Problem 3**, the open network is stable under proportional service for all $\lambda$ (infinite-server behavior). Under fixed service, node $c$ is the bottleneck with critical rate $4/3$.

### HW3 in 5 Minutes

Homework 3 covers epidemic simulations and network games. **Part 1** simulates SIR epidemics on $k$-regular and preferential attachment graphs. The infection probability with $m$ infected neighbors is $1-(1-\beta)^m$. Preferential attachment creates hubs that accelerate spread. **Vaccination** removes transmission paths, flattening the epidemic curve. **Parameter estimation** for H1N1 uses RMSE minimization over a grid of $(k, \beta, \rho)$.

The **games** section studies coordination/anti-coordination on small graphs. Nash equilibria depend on the mix of player types. Best-response dynamics may converge to NE or cycle. Noisy best response with vanishing noise selects stochastically stable equilibria.

*Note: HW3 code and results are not yet implemented. Theory and roadmap are provided.*

### Full Course Connection in 5 Minutes

The course studies **network dynamics and learning** through three lenses: optimization, stochastic processes, and strategic behavior. **HW1** establishes graph foundations: max-flow/min-cut for capacity, centrality for importance, and traffic assignment for congestion pricing. **HW2** introduces stochastic dynamics: CTMCs for random walks, French-DeGroot for opinions, and particle systems for queueing. **HW3** applies these tools to epidemics and games.

The connecting themes are: (1) **graph structure determines outcomes** — cliques boost centrality, SCCs determine consensus, hubs accelerate epidemics, bottlenecks limit throughput; (2) **local vs. global** — selfish routing vs. social optimum, local opinion updates vs. global consensus, local best response vs. Nash equilibrium; (3) **simulation and theory complement each other** — Monte Carlo validates formulas, optimization produces exact solutions, graph analysis provides guarantees.
