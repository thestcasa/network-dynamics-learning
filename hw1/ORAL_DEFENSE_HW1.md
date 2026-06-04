# ORAL_DEFENSE_HW1.md

## 1. One-page executive summary

Homework 1 has three parts.

- Exercise 1 studies a small directed capacitated network from source `o` to sink `d`. I enumerated all `o-d` cuts, verified max-flow/min-cut, then optimized integer extra-capacity allocations on existing links and with one added directed link.
- Exercise 2 reconstructs the undirected graph from the assignment figure and compares Katz centrality with a distributed PageRank iteration. Katz ranks `n6` first; PageRank with `beta = 0.15` ranks `n9` first.
- Exercise 3 uses the provided traffic-network `.mat` files to compute a shortest path, max flow, `nu = Bf`, social optimum, Wardrop equilibrium, marginal-cost tolls, and additional-delay tolls.

Key numerical results:

- Exercise 1: min cut = max flow = `5`; existing-link throughput is `T_b(x) = 5 + ceil(x/2)` for the tested integer budgets; with direct added link `o -> d`, `T_c(x) = 6 + x`.
- Exercise 2: `lambda_max(W) = 5.0710614673`, so `1/lambda_max(W) = 0.1971973731`; `beta = 0.15` is valid for Katz. Katz top node is `n6`; PageRank top node is `n9`.
- Exercise 3: shortest path is `1 -> 2 -> 3 -> 9 -> 13 -> 17`, free-flow time `0.559833`; max flow is `22448`; reduced demand is `nu_new_1 = 16806`, `nu_new_17 = -16806`; social optimum total travel time is `26142.669150`; Wardrop total travel time is `26495.321759`; PoA is `1.013490`; additional-delay optimum is `15350.353857`.

Most important oral-defense points:

- A directed cut counts only edges leaving the source-side set.
- Katz counts raw walks, while PageRank degree-normalizes neighbor contributions.
- Wardrop equilibrium minimizes the Beckmann potential, not total travel time directly.
- Marginal-cost tolls decentralize the social optimum by charging congestion externalities.
- Negative additional-delay tolls can be valid because the additional-delay objective subtracts free-flow time.

Risky parts to be ready to justify:

- The direct link `o -> d` in Exercise 1(c) is a literal interpretation of the assignment because direct source-sink arcs were not explicitly forbidden.
- The Exercise 2 graph reconstruction must match the assignment figure.
- PageRank at `beta = 1` needs the connected, non-bipartite, undirected graph argument.
- Exercise 3 uses a numerical cap `f_e <= 0.999 C_e`; explain it as a domain safeguard.

---

## 2. Official assignment checklist

| Exercise | Official request | Where it is answered | Final result | What I must be ready to explain |
|---|---|---|---|---|
| 1(a) | Enumerate all cuts and find the minimum cut | Report Exercise 1, `results/exercise1_cuts.csv`, Section 3.3 here | Minimum cut capacity `5`; max flow `5` | Directed cut convention and max-flow/min-cut theorem |
| 1(b) | Add integer extra capacity on existing links | Report Exercise 1, `results/exercise1_b_results.csv`, Section 3.4 here | `T_b(x) = 5 + ceil(x/2)` | Why capacity must be paired across bottlenecks |
| 1(c) | Add one new directed link with base capacity `1`, then extra capacity | Report Exercise 1, `results/exercise1_c_results.csv`, Section 3.5 here | Literal direct-link interpretation gives `T_c(x) = 6 + x`; CSV chooses `b->d` at `x=0` as a tie/alternative optimum and `o->d` for `x >= 1` | Why `o -> d` is valid under the used candidate set |
| 2(a) | Katz centrality | Report Exercise 2, `results/exercise2_katz.csv`, Section 4.2 here | `n6` rank 1, Katz `0.3177933766` | Spectral condition and walk-count interpretation |
| 2(b) | Distributed PageRank algorithm | Report Exercise 2, `results/exercise2_pagerank_beta015.csv`, Section 4.3 here | `n9` rank 1, PageRank `0.1059574692`; `n6` rank 2 | What each node sends to neighbors |
| 2(c) | Interpret nodes `n6` and `n9` | Report Exercise 2, Section 4.4 here | Katz favors `n6`; PageRank favors `n9` | Degree normalization versus raw walks |
| 2(d) | PageRank sensitivity for `beta in {0, 1/4, 1/2, 3/4, 1}` | Report Exercise 2, `results/exercise2_pagerank_beta_sensitivity.csv`, Section 4.5 here | Difference is `0`, `-0.0394710887`, `-0.0652885611`, `-0.0805114506`, `0.0208333333` | Why the difference is not monotone |
| 3(a) | Shortest path | Report Exercise 3, `results/exercise3_shortest_path.csv`, Section 5.3 here | `1 -> 2 -> 3 -> 9 -> 13 -> 17`, time `0.559833` | Free-flow travel times as weights |
| 3(b) | Maximum flow | Report Exercise 3, `results/exercise3_maxflow.csv`, Section 5.4 here | `22448` | Capacity-only interpretation |
| 3(c) | Compute `nu = Bf` | Report Exercise 3, `results/exercise3_nu.csv`, Section 5.5 here | `nu_new_1 = 16806`, `nu_new_17 = -16806` | Incidence sign convention |
| 3(d) | Social optimum | Report Exercise 3, `results/exercise3_social_optimum.csv`, Section 5.7 here | Total travel time `26142.669150` | Convex objective and constraints |
| 3(e) | Wardrop equilibrium | Report Exercise 3, `results/exercise3_wardrop.csv`, Section 5.8 here | Total travel time `26495.321759`; PoA `1.013490` | Beckmann potential |
| 3(f) | Marginal-cost tolls and tolled Wardrop | Report Exercise 3, `results/exercise3_tolls_total_travel_time.csv`, Section 5.9 here | Error `2.993561e-1`, relative `1.325734e-5` | Toll equals externality |
| 3(g) | Additional-delay objective and tolls | Report Exercise 3, `results/exercise3_additional_delay.csv`, Section 5.10 here | Objective `15350.353857`; error `5.658613e-1`, relative `2.507184e-5` | Why tolls may be negative |

---

## 2.5 Network Dynamics and Learning concept glossary

Use this as the conceptual vocabulary list for the oral. These are the course concepts that appear in the report or in this defense file, phrased in a way that is useful when answering questions.

### Graph and flow concepts

| Concept | Meaning | Where it appears | Oral-defense sentence |
|---|---|---|---|
| Directed graph | A graph whose edges have orientation, so `u -> v` and `v -> u` are different arcs. | Exercises 1 and 3 | Direction matters for cuts, incidence matrices, and feasible routing. |
| Undirected graph | A graph where edges have no orientation. | Exercise 2 | The centrality graph is undirected, so random-walk stationary mass is proportional to degree at `beta = 1`. |
| Node | A vertex of the network. | All exercises | Nodes represent junctions, graph positions, or traffic intersections depending on the exercise. |
| Edge / arc / link | A connection between two nodes; in directed graphs it has a tail and a head. | All exercises | I use edge, arc, and link depending on context, but in traffic assignment links are directed columns of `B`. |
| Tail and head | Tail is where a directed link starts; head is where it ends. | Exercises 1 and 3 | The incidence convention is `+1` at the tail and `-1` at the head. |
| Capacity | Maximum allowable flow on an edge or link. | Exercises 1 and 3 | Capacity limits feasible flow; in Exercise 3 delay also diverges near capacity. |
| Source and sink | Source injects flow; sink receives flow. | Exercises 1 and 3 | In Exercise 1 they are `o` and `d`; in Exercise 3 they are nodes `1` and `17`. |
| Throughput | Total amount of flow sent from source to sink. | Exercise 1 | Throughput is the scalar value optimized by max flow and capacity augmentation. |
| Feasible flow | A nonnegative flow satisfying capacities and conservation. | Exercises 1 and 3 | Feasibility is about respecting both edge bounds and node balance. |
| Flow conservation | At intermediate nodes, incoming and outgoing flow must balance, except for exogenous supply/demand. | Exercises 1 and 3 | In matrix form for Exercise 3, conservation is `Bf = nu`. |
| Source-sink cut | A partition with the source on one side and the sink on the other. | Exercise 1 | A cut is a certificate upper-bounding any feasible source-sink flow. |
| Directed cut capacity | Sum of capacities of arcs leaving the source-side set. | Exercise 1 | Incoming arcs do not count because they do not carry flow across the cut toward the sink. |
| Minimum cut | The source-sink cut with smallest capacity. | Exercise 1 | It identifies the tightest structural bottleneck. |
| Maximum flow | Largest feasible source-sink throughput. | Exercises 1 and 3 | In Exercise 1 it equals the min cut; in Exercise 3 it ignores congestion delays. |
| Max-flow/min-cut theorem | Maximum source-sink flow equals minimum source-sink cut capacity. | Exercise 1 | This theorem explains why the computed max flow and min cut are both `5`. |
| Bottleneck | A capacity-limited part of the network that constrains total throughput. | Exercise 1 | Extra capacity must be placed at bottlenecks to change the max-flow value. |
| Capacity augmentation | Increasing edge capacities with an extra budget. | Exercise 1 | The optimization asks where extra capacity has the largest throughput effect. |
| Integer allocation | Extra capacity variables restricted to integer units. | Exercise 1 | The code can solve it exactly by enumerating integer compositions. |
| Candidate link set | The set of possible new arcs considered for addition. | Exercise 1(c) | The `o -> d` result depends on including direct source-sink arcs in this set. |
| Global optimum by enumeration | Checking every allowed discrete option rather than using a heuristic. | Exercise 1 | Exhaustive enumeration justifies global optimality for the tested budgets. |

### Centrality and random-walk concepts

| Concept | Meaning | Where it appears | Oral-defense sentence |
|---|---|---|---|
| Adjacency matrix `W` | Matrix with entries indicating whether nodes are connected. | Exercise 2 | Katz centrality uses `W` directly, so dense graph structure matters. |
| Degree | Number of neighbors of a node in an undirected graph. | Exercise 2 | Degree affects PageRank normalization and the `beta = 1` stationary distribution. |
| Clique | A set of nodes all connected to each other. | Exercise 2 | The clique on `n1,...,n6` creates many walks and boosts Katz centrality. |
| Bridge / connector node | A node linking two graph regions. | Exercise 2 | `n6` bridges the clique to the path and star part. |
| Star-like structure | A center node connected to several leaves. | Exercise 2 | `n9` receives concentrated PageRank mass from degree-one leaves. |
| Leaf node | A node of degree one. | Exercise 2 | Leaves send all their PageRank mass to their only neighbor. |
| Walk | A sequence of adjacent nodes, allowing repeated nodes. | Exercise 2 | Katz centrality counts weighted walks of all lengths. |
| Intrinsic centrality `mu` | Baseline importance assigned independently of network position. | Exercise 2 | Uniform `mu` means no node is favored before network effects. |
| Katz centrality | Centrality based on discounted counts of walks. | Exercise 2 | Katz favors `n6` because raw walks through the clique are numerous. |
| Spectral radius / largest eigenvalue | The largest eigenvalue controlling convergence of matrix walk expansions. | Exercise 2 | The Katz condition is `beta < 1/lambda_max(W)`. |
| Resolvent `(I-beta W)^(-1)` | Matrix inverse that sums discounted walks when the spectral condition holds. | Exercise 2 | This is why Katz has a closed-form linear-algebra solution. |
| PageRank | Centrality based on a random walk with teleportation when `beta < 1`. | Exercise 2 | PageRank uses normalized neighbor contributions, unlike Katz. |
| Distributed algorithm | An algorithm where each node updates using local neighbor messages. | Exercise 2 | Each node only needs incoming `x_j/deg(j)` messages. |
| Degree normalization | Dividing a node's contribution by its degree. | Exercise 2 | It explains why `n9` beats `n6` for PageRank at `beta = 0.15`. |
| Transition matrix | Matrix describing one step of a random walk. | Exercise 2 | For the undirected graph, columns distribute mass uniformly across neighbors. |
| Column-stochastic matrix | A matrix whose columns sum to one. | Exercise 2 | It preserves total PageRank mass during random-walk propagation. |
| Teleportation | The `(1-beta)mu` term in PageRank. | Exercise 2 | Teleportation pulls scores toward intrinsic centrality and helps convergence. |
| Stationary distribution | A probability vector unchanged by the transition matrix. | Exercise 2 | At `beta = 1`, PageRank is the random-walk stationary distribution. |
| Connected graph | A graph where every node can reach every other node. | Exercise 2 | Connectivity supports uniqueness of the stationary distribution. |
| Non-bipartite graph | A graph without purely two-period random-walk oscillation. | Exercise 2 | Non-bipartiteness supports convergence of the `beta = 1` iteration. |
| Monotonicity | Whether a sequence only increases or only decreases. | Exercise 2(d) | `PR(n6)-PR(n9)` is not monotone in the tested beta values. |

### Traffic-assignment and optimization concepts

| Concept | Meaning | Where it appears | Oral-defense sentence |
|---|---|---|---|
| Traffic network | Directed graph where links carry traffic and have capacities/travel times. | Exercise 3 | The network has `17` nodes and `28` directed links. |
| Incidence matrix `B` | Matrix encoding link directions through tail/head signs. | Exercise 3 | `Bf = nu` translates link flows into node imbalances. |
| Exogenous inflow vector `nu` | Net external injection/removal at each node. | Exercise 3 | Positive entries are sources; negative entries are sinks. |
| Reduced demand vector `nu_new` | Simplified OD demand keeping only node 1 and node 17 nonzero. | Exercise 3 | It turns the data into a single source-sink routing problem. |
| Origin-destination demand | Traffic that must be routed from an origin to a destination. | Exercise 3 | Here the demand is `16806` from node `1` to node `17`. |
| Free-flow travel time `l_e` | Link travel time at zero congestion. | Exercise 3 | Shortest path uses these values as edge weights. |
| Shortest path | Path with minimum total edge weight. | Exercise 3(a) | It is the fastest route in an empty network. |
| Capacity utilization | Ratio `f_e/C_e`. | Exercise 3 | It shows how close a link is to capacity. |
| Congestion delay | Increase in travel time caused by nonzero flow. | Exercise 3 | The delay function grows sharply as utilization approaches one. |
| Delay function `tau_e` | Flow-dependent travel time on link `e`. | Exercise 3 | `tau_e(f_e) = l_e/(1-f_e/C_e)`. |
| Domain constraint | Restriction needed for a function to be finite and well-defined. | Exercise 3 | The delay formula requires `f_e < C_e`. |
| Numerical safeguard | Practical solver constraint that prevents invalid boundary solutions. | Exercise 3 | The cap `0.999 C_e` keeps the solver away from infinite delay. |
| Convex optimization | Optimization where local optima are global under convex objective and constraints. | Exercise 3 | Social optimum and Beckmann problems are solved as convex programs. |
| Linear equality constraint | Constraint like `Bf = nu_new`. | Exercise 3 | It enforces conservation exactly. |
| Nonnegativity constraint | Constraint `f_e >= 0`. | Exercise 3 | Traffic cannot be negative on a link. |
| Social optimum | Flow minimizing total system travel time. | Exercise 3(d) | It is the planner's solution, not necessarily selfishly stable. |
| Wardrop equilibrium | Selfish-routing state where no driver can improve by switching routes. | Exercise 3(e) | It is a user equilibrium. |
| User equilibrium | Same idea as Wardrop: all used routes have minimal perceived cost. | Exercise 3 | Users minimize private travel time, not system cost. |
| Beckmann potential | Integral objective whose minimizer gives Wardrop equilibrium for separable costs. | Exercise 3(e) | Its derivative equals link delay. |
| Separable link costs | Each link cost depends only on its own flow. | Exercise 3 | This makes the Beckmann formulation applicable. |
| Total travel time | `sum_e f_e tau_e(f_e)`. | Exercise 3 | This is the social cost used for the social optimum and PoA. |
| Price of anarchy | Ratio between selfish equilibrium cost and social optimum cost. | Exercise 3(e) | Here it is `1.013490`, so selfish routing costs about `1.349%` more. |
| Marginal social cost | Extra total cost caused by one more unit of flow on a link. | Exercise 3(f) | It is `tau_e + f_e tau'_e`. |
| Externality | Cost a user imposes on others but does not personally perceive. | Exercise 3(f) | The marginal-cost toll charges exactly this externality. |
| Toll | Added perceived link cost used to influence selfish route choices. | Exercise 3(f,g) | Tolls decentralize a system-optimal target flow. |
| Decentralization | Making selfish individual optimization reproduce a planner's optimum. | Exercise 3(f,g) | Marginal-cost tolls decentralize the social optimum. |
| Tolled Wardrop equilibrium | Wardrop equilibrium under perceived cost `delay + toll`. | Exercise 3(f,g) | It is used to verify that tolls recover the target flow. |
| Additional-delay objective | Objective counting only delay above free-flow travel time. | Exercise 3(g) | It minimizes `sum_e f_e(tau_e-l_e)`. |
| Subsidy / negative toll | A negative added cost, interpreted as a payment or offset. | Exercise 3(g) | Negative additional-delay tolls can be mathematically valid. |
| Solver status | Diagnostic telling whether the numerical optimization solved successfully. | Exercise 3 | The reported statuses are `optimal`. |
| Solver tolerance | Numerical accuracy threshold, not exact symbolic equality. | Exercise 3 | Small flow-recovery errors are expected and acceptable. |
| Relative error | Error normalized by the scale of the target vector. | Exercise 3 | It shows toll recovery is tiny relative to flows of thousands. |

---

## 3. Exercise 1 oral-defense guide - Max-flow and cuts

### 3.1 Problem setup

The graph is directed and capacitated. The nodes are `V = {o, a, b, c, d}`. The source is `o`, the sink is `d`, and the directed edges are:

| Edge | Arc | Capacity |
|---|---|---:|
| `e1` | `o -> a` | 3 |
| `e2` | `a -> d` | 3 |
| `e3` | `o -> b` | 3 |
| `e4` | `b -> c` | 3 |
| `e5` | `c -> d` | 2 |
| `e6` | `a -> b` | 1 |

Throughput means the total amount of feasible flow that can be sent from `o` to `d`. A feasible flow respects capacity constraints on each directed edge and flow conservation at intermediate nodes.

### 3.2 Cut definition

An `o-d` cut is a set `S` such that `o in S` and `d not in S`. The three intermediate nodes `a,b,c` can independently be either in `S` or outside `S`, so there are `2^3 = 8` cuts.

The directed-cut convention is important: the cut capacity is the sum of capacities of edges whose tail is in `S` and whose head is outside `S`. Edges entering `S` do not count because they do not carry flow from the source side toward the sink side.

### 3.3 My result for part (a)

| `S` | `V \ S` | Crossing edges | Capacity |
|---|---|---|---:|
| `{o}` | `{a,b,c,d}` | `e1, e3` | 6 |
| `{o,a}` | `{b,c,d}` | `e2, e3, e6` | 7 |
| `{o,b}` | `{a,c,d}` | `e1, e4` | 6 |
| `{o,c}` | `{a,b,d}` | `e1, e3, e5` | 8 |
| `{o,a,b}` | `{c,d}` | `e2, e4` | 6 |
| `{o,a,c}` | `{b,d}` | `e2, e3, e5, e6` | 9 |
| `{o,b,c}` | `{a,d}` | `e1, e5` | 5 |
| `{o,a,b,c}` | `{d}` | `e2, e5` | 5 |

The minimum cut value is `5`. The maximum flow value is also `5`. By the max-flow/min-cut theorem, the maximum feasible throughput from source to sink equals the minimum capacity of any source-sink cut.

The "minimum capacity to remove" is equal to the min-cut value because removing all edges in a minimum cut disconnects `o` from `d`, and no smaller total removed capacity can block all feasible `o-d` flow.

### 3.4 My result for part (b)

The optimization distributes an integer budget `x` of extra capacity over the six existing links. The variables `y_i` are nonnegative integers: `y_i` is how many extra capacity units are assigned to edge `e_i`.

The optimal throughput pattern in the generated CSV is:

`5, 6, 6, 7, 7, 8, 8, ..., 15` for `x = 0,...,20`.

The closed form is:

`T_b(x) = 5 + ceil(x/2)`.

Intuitively, not every extra capacity unit immediately increases throughput. Some bottlenecks come in source-side/sink-side pairs, so one isolated extra unit may only relax one side of a limiting cut without creating a full additional route to `d`.

Important edge cases for the oral explanation:

- For `x = 0`, no capacity is added.
- For odd `x = 2k + 1`, one optimal allocation is `y1 = k`, `y2 = k`, `y5 = 1`.
- For even `x = 2k >= 2`, one optimal allocation is `y1 = k`, `y2 = k - 1`, `y5 = 1`.
- Other optimal allocations may exist. The code reports one optimizer found by exact enumeration.

### 3.5 My result for part (c)

The new directed link has base capacity `c8 = 1`. The extra-capacity budget `x` is separate from this base capacity. The code considered all candidate directed links between distinct nodes that were not already present.

The direct link `o -> d` was included because the assignment does not explicitly forbid source-to-sink links. Under that literal interpretation, `o -> d` is optimal because it crosses every `o-d` cut: `o` is always in `S`, and `d` is always outside `S`.

Assigning all extra capacity to `o -> d` gives:

`T_c(x) = 5 + 1 + x = 6 + x`.

The CSV chooses `b -> d` for `x = 0`, but `o -> d` also gives throughput `6` at `x = 0`; for `x >= 1`, the reported best link is `o -> d`.

Defensive answer if challenged:

> Under the literal statement, direct source-to-sink arcs are not excluded, so I included it in the candidate set. If the intended interpretation excluded direct source-to-sink links, the candidate set should be restricted and the optimization recomputed. My reported result is valid for the candidate set explicitly used.

### 3.6 Possible oral questions and model answers

1. Why are there 8 cuts?  
   There are three intermediate nodes, each either included or excluded from the source-side set, so `2^3 = 8`.

2. Why do we count only edges leaving `S`?  
   Because only those edges can carry flow from the source side to the sink side.

3. What is the max-flow/min-cut theorem?  
   The maximum feasible source-sink flow equals the minimum source-sink cut capacity.

4. Why is the min cut 5?  
   The cuts `{o,b,c}` and `{o,a,b,c}` both have capacity `5`, and all other cuts have capacity at least `6`.

5. Why is max flow also 5?  
   By max-flow/min-cut, and NetworkX also computes a maximum flow value of `5`.

6. What does feasible flow mean?  
   Nonnegative edge flows satisfying capacities and conservation at all non-source/non-sink nodes.

7. Why does part (b) grow like `ceil(x/2)`?  
   Because extra capacity must often be added to complementary bottlenecks before it raises the minimum cut.

8. Is the optimal allocation unique?  
   No. The exhaustive search returns one optimum; ties can exist.

9. Why is `x` integer?  
   The assignment treats the extra-capacity budget as integer units, and the code enumerates integer allocations.

10. Why is `o -> d` optimal in part (c)?  
    It crosses every source-sink cut and adds direct capacity from source to sink.

11. Is `o -> d` a loophole?  
    It may be considered too trivial, but it is valid under the literal candidate set unless direct source-sink arcs are forbidden.

12. What changes if `o -> d` is forbidden?  
    The candidate set changes, and the added-link optimization must be recomputed.

13. What is the difference between base capacity and extra capacity?  
    The new link gets base capacity `1` automatically; `x` is an additional budget allocated afterward.

14. How did the code enumerate allocations?  
    It generated all nonnegative integer compositions of each budget over the available edges.

15. How do you know the result is globally optimal?  
    For each budget and candidate link, the code exhaustively checks all integer allocations and uses min-cut values.

---

## 4. Exercise 2 oral-defense guide - Katz centrality and PageRank

### 4.1 Graph reconstruction

The graph is undirected and simple. The reconstructed edge set is:

```text
n1--n2, n1--n3, n1--n4, n1--n5, n1--n6
n2--n3, n2--n4, n2--n5, n2--n6
n3--n4, n3--n5, n3--n6
n4--n5, n4--n6
n5--n6
n6--n15, n6--n7
n7--n8
n8--n9
n9--n10, n9--n11, n9--n12, n9--n13, n9--n14
```

Key structure:

- Nodes `n1,...,n6` form a clique.
- `n6` connects the clique to the rest of the graph.
- `n9` is the center of a star-like part.
- The two parts are connected through `n6 - n7 - n8 - n9`.
- `n15` is attached to `n6`.

The graph reconstruction must match the assignment figure; all centrality values depend on this edge set.

### 4.2 Katz centrality

Katz centrality is:

`x = (I - beta W)^(-1) mu`.

Here `W` is the adjacency matrix and `mu` is intrinsic centrality. I used uniform intrinsic centrality, `mu_i = 1/n = 1/15`. If `mu_i = 1` were used instead, all Katz values would be multiplied by `15`, but the ranking would not change.

The spectral condition is:

`beta < 1 / lambda_max(W)`.

Numerically:

- `lambda_max(W) = 5.0710614673`
- `1/lambda_max(W) = 0.1971973731`
- `beta = 0.15`, so the condition is satisfied.

Key Katz results:

- `n6` ranked first with Katz `0.3177933766`.
- `n1,...,n5` tied for second with Katz `0.2858391829`.
- `n9` is below `n6`, with Katz `0.1498337716`.

Katz favors `n6` because Katz counts walks and is not degree-normalized. Node `n6` is inside or attached to the dense clique and also bridges to the rest of the graph, so many walks reach or pass through it.

### 4.3 Distributed PageRank

The distributed PageRank update is:

`x_i(t+1) = (1 - beta) mu_i + beta * sum_{j in N(i)} x_j(t) / deg(j)`.

Each node `j` sends `x_j(t)/deg(j)` to each neighbor. Node `i` only needs incoming messages from its neighbors and the common parameters `beta` and `mu_i`, so the algorithm is distributed.

The code uses tolerance `1e-12`. For `beta = 0.15`, it converged in `15` iterations with final `L1` difference `2.011e-13`.

Key PageRank result for `beta = 0.15`:

- `n9` ranked first with PageRank `0.1059574692`.
- `n6` ranked second with PageRank `0.0801149613`.
- Degree-one leaves send all their mass to `n9`.
- Clique nodes split their mass among many neighbors, so `n6` receives normalized contributions.

### 4.4 Comparison between Katz and PageRank

Katz ranks `n6` higher because it counts raw walk contributions. The dense clique creates many walks involving `n6`.

PageRank ranks `n9` higher because of degree normalization. The leaves `n10,...,n14` have degree one, so each sends all its neighbor mass to `n9`. This is the central conceptual point: Katz rewards dense walk structure, while PageRank rewards concentrated incoming normalized mass.

### 4.5 Sensitivity to beta

For `beta in {0, 1/4, 1/2, 3/4, 1}`:

| `beta` | `PR(n6)` | `PR(n9)` | `PR(n6)-PR(n9)` |
|---:|---:|---:|---:|
| 0 | 0.0666666667 | 0.0666666667 | 0.0000000000 |
| 1/4 | 0.0875066385 | 0.1269777272 | -0.0394710887 |
| 1/2 | 0.1014039855 | 0.1666925466 | -0.0652885611 |
| 3/4 | 0.1106044872 | 0.1911159378 | -0.0805114506 |
| 1 | 0.1458333333 | 0.1250000000 | 0.0208333333 |

At `beta = 0`, only intrinsic centrality matters, so all nodes are equal. At `beta = 1`, PageRank becomes the stationary distribution of the simple random walk. Because the graph is connected and non-bipartite, the stationary distribution is unique and the iteration converges. In an undirected graph:

`pi_i = deg(i) / (2|E|)`.

Here `|E| = 24`, so `2|E| = 48`. Since `deg(n6) = 7` and `deg(n9) = 6`, `PR(n6) = 7/48` and `PR(n9) = 6/48`, so the difference is `1/48 = 0.0208333333`.

The difference is not monotone: it starts at `0`, becomes negative for intermediate `beta`, then becomes positive at `beta = 1`.

### 4.6 Possible oral questions and model answers

1. What is the difference between Katz and PageRank?  
   Katz counts raw walks; PageRank uses degree-normalized random-walk contributions.

2. Why must Katz satisfy a spectral condition?  
   To make `(I - beta W)` invertible and the walk series convergent.

3. Why does PageRank not need the same condition for `beta < 1`?  
   The update is a contraction because of teleportation.

4. What happens at `beta = 0`?  
   The result equals `mu`, so all nodes are equal for uniform `mu`.

5. What happens at `beta = 1`?  
   PageRank becomes the stationary distribution of the simple random walk.

6. Why does `n9` beat `n6` for `beta = 0.15`?  
   Its degree-one leaves send all their PageRank mass to it.

7. Why does `n6` beat `n9` at `beta = 1`?  
   In an undirected graph, stationary mass is proportional to degree; `n6` has degree `7`, `n9` degree `6`.

8. Is PageRank iteration guaranteed to converge at `beta = 1`?  
   Not always; here it converges because the graph is connected and non-bipartite.

9. What role does non-bipartiteness play?  
   It removes period-two oscillations in the simple random walk.

10. Why is `mu_i = 1/15` used?  
    It represents uniform intrinsic centrality over 15 nodes.

11. Would Katz rankings change if `mu_i = 1`?  
    No, values scale by 15 but rankings stay the same.

12. What does "distributed algorithm" mean here?  
    Each node updates using only messages from neighbors.

13. What information does each node need?  
    Its own `mu_i`, `beta`, and incoming `x_j/deg(j)` messages.

14. Why do leaves matter so much for PageRank?  
    A degree-one leaf sends all its random-walk mass to its only neighbor.

15. What would happen if the graph were directed?  
    The transition matrix would use outgoing degrees and dangling nodes would need special handling.

16. Why is `n6` high in Katz?  
    It participates in many walks due to the clique and bridge position.

17. Why are `n1,...,n5` tied in Katz?  
    They are symmetric in the reconstructed clique.

18. Why is `n15` low in PageRank?  
    It is a degree-one leaf attached to `n6` and receives little normalized mass back.

19. What is `W`?  
    The adjacency matrix of the undirected graph.

20. What is the main lesson of Exercise 2?  
    Different centrality definitions capture different structural roles.

---

## 5. Exercise 3 oral-defense guide - Traffic assignment

### 5.1 Problem setup

Exercise 3 uses a directed traffic network with `17` nodes and `28` links. The data files provide:

- incidence matrix `B` from `traffic.mat`;
- capacity vector `C` from `capacities.mat`;
- free-flow travel-time vector `l` from `traveltime.mat`;
- measured flow vector `f` from `flow.mat`;
- source node `1`;
- sink node `17`.

### 5.2 Incidence matrix convention

Each column of `B` is a directed link. The convention in the code is:

- `+1` marks the tail;
- `-1` marks the head.

Flow conservation is:

`B f = nu`.

A positive `nu_i` means net exogenous inflow/source at node `i`. A negative `nu_i` means net exogenous outflow/sink.

### 5.3 Shortest path

The shortest path uses free-flow times `l_e` as edge weights. This is the fastest path in an empty network, before congestion effects.

Result:

- path: `1 -> 2 -> 3 -> 9 -> 13 -> 17`;
- link indices: `1, 2, 12, 9, 25`;
- total free-flow travel time: `0.559833`.

### 5.4 Maximum flow

Maximum flow uses capacities `C_e` and ignores congestion delays. The result is:

`22448`.

This is the largest feasible flow from node `1` to node `17` if links only have capacity constraints.

### 5.5 Vector nu = Bf

The original imbalance is computed as `nu = Bf` from the measured flow. The assignment then constructs a reduced single-origin/single-destination demand vector:

`nu_new_1 = 16806`  
`nu_new_17 = -16806`

All other entries are set to zero because the optimization models one origin-destination demand from node `1` to node `17`.

### 5.6 Delay function

The link delay is:

`tau_e(f_e) = l_e / (1 - f_e / C_e)`.

At zero flow, delay equals the free-flow time `l_e`. As flow approaches capacity, delay tends to infinity. The domain is:

`0 <= f_e < C_e`.

The numerical constraint:

`f_e <= 0.999 C_e`

keeps the solver inside the finite domain. It is a numerical safeguard; after solving, maximum utilization should be checked to make sure the artificial cap is not determining the result.

### 5.7 Social optimum

The social optimum minimizes total travel time:

`sum_e f_e tau_e(f_e)`

subject to:

`Bf = nu_new`, `0 <= f <= 0.999 C`.

The equivalent convex CVXPY objective is:

`sum_e (l_e C_e / (1 - f_e/C_e) - l_e C_e)`.

This is equivalent because:

`f_e l_e / (1 - f_e/C_e) = l_e C_e / (1 - f_e/C_e) - l_e C_e`.

Result:

- solver status: `optimal`;
- social optimum total travel time: `26142.669150`.

### 5.8 Wardrop equilibrium

Wardrop equilibrium means no driver can reduce their travel time by unilaterally changing route. For separable increasing link costs, it can be computed by minimizing the Beckmann potential:

`sum_e integral_0^{f_e} tau_e(s) ds`.

For this delay function:

`integral_0^{f_e} tau_e(s) ds = -l_e C_e log(1 - f_e/C_e)`.

Result:

- Wardrop total travel time: `26495.321759`;
- price of anarchy:

`PoA = Wardrop cost / Social optimum cost = 1.013490`.

Selfish routing is about `1.349%` worse than the social optimum.

### 5.9 Marginal-cost tolls

Individual users perceive `tau_e(f_e)`. The social planner cares about `f_e tau_e(f_e)`. The marginal social cost is:

`tau_e(f_e) + f_e tau'_e(f_e)`.

Therefore the toll is the externality:

`omega_e = f_e^* tau'_e(f_e^*)`.

The derivative is:

`tau'_e(f_e) = (l_e/C_e) / (1 - f_e/C_e)^2`.

Tolls make selfish users internalize the congestion externality, so the tolled Wardrop equilibrium recovers the social optimum.

Verification:

- absolute error: `||f_tolled - f^*||_2 = 2.993561e-1`;
- relative error: `1.325734e-5`.

Flows are thousands, so this is tiny.

### 5.10 Additional-delay objective

The additional-delay objective is total delay beyond free flow:

`psi_e(f_e) = f_e (tau_e(f_e) - l_e)`.

It ignores unavoidable free-flow travel time and penalizes only congestion delay.

Result:

- objective value: `15350.353857`.

For a general system objective `psi_e(f_e)`, the decentralizing toll is:

`omega_e = psi'_e(f_e^*) - tau_e(f_e^*)`.

For additional delay:

`omega_e^add = -l_e + f_e^* tau'_e(f_e^*)`.

Negative tolls are not necessarily errors. They can be interpreted as subsidies or offsets, and they appear because the formula subtracts the free-flow term `l_e`.

Verification:

- absolute error: `||f_tolled,add - f^*,add||_2 = 5.658613e-1`;
- relative error: `2.507184e-5`.

### 5.11 Possible oral questions and model answers

1. What does `Bf = nu` mean?  
   It is flow conservation with exogenous inflow/outflow at nodes.

2. Why is `+1` tail and `-1` head important?  
   It fixes the sign convention for interpreting `nu`.

3. What is positive `nu`?  
   Net exogenous inflow, so a source.

4. What is negative `nu`?  
   Net exogenous outflow, so a sink.

5. Why do we build `nu_new`?  
   To model a single OD demand from node 1 to node 17.

6. Why use free-flow times for shortest path?  
   The question asks for the fastest path in an empty network.

7. Difference between shortest path and Wardrop?  
   Shortest path ignores congestion; Wardrop accounts for flow-dependent delays.

8. Difference between max flow and traffic assignment?  
   Max flow only uses capacities; traffic assignment optimizes delay under demand.

9. Why does delay go to infinity near capacity?  
   The denominator `1 - f_e/C_e` tends to zero.

10. Why impose `0.999 C_e`?  
    To keep the numerical solver inside the finite domain.

11. What is the social optimum?  
    The feasible flow minimizing total travel time.

12. What is Wardrop equilibrium?  
    A selfish-routing state where no used route can be improved unilaterally.

13. Why use Beckmann potential?  
    Its first-order conditions match Wardrop conditions for separable costs.

14. Derive the Beckmann integral.  
    Integrating `l_e/(1-s/C_e)` from `0` to `f_e` gives `-l_e C_e log(1-f_e/C_e)`.

15. What is price of anarchy?  
    Wardrop total cost divided by social optimum total cost.

16. Why is Wardrop worse than social optimum?  
    Users ignore the congestion externality they impose on others.

17. What are marginal-cost tolls?  
    Charges equal to the external marginal delay `f_e tau'_e(f_e)`.

18. Why do tolls recover the social optimum?  
    They make perceived private cost equal marginal social cost.

19. Why can additional-delay tolls be negative?  
    The formula includes `-l_e`, so the toll may be a subsidy.

20. What does a negative toll mean?  
    A subsidy or offset, not a numerical mistake by itself.

21. How do you verify numerically that tolls worked?  
    Compare target optimal flow with tolled Wardrop flow using an `L2` norm.

22. What solver was used?  
    CVXPY, primarily CLARABEL, with ECOS fallback for tolled cases if needed.

23. Is the cap active?  
    It should be checked by capacity utilization; the cap is intended only as a domain safeguard.

24. What is the demand size?  
    `16806` units from node `1` to node `17`.

25. Why is the social objective convex?  
    The link cost expression is convex on `0 <= f_e < C_e`, and constraints are linear.

---

## 6. Formula cheat sheet

Exercise 1:

- Cut capacity: `C(S) = sum_{e: tail(e) in S, head(e) notin S} c_e`. This is the capacity available to leave the source side of a directed cut.
- Max-flow/min-cut theorem: `max feasible o-d flow = min_S C(S)`. The bottleneck cut determines the maximum throughput.
- Flow conservation: for every intermediate node, total incoming flow equals total outgoing flow. This is what prevents flow from being created or destroyed inside the network.
- Capacity constraint: `0 <= f_e <= c_e`. This is the local feasibility condition on each directed edge.
- Existing-link throughput: `T_b(x) = 5 + ceil(x/2)`. Integer capacity on existing links improves throughput roughly every two units after the base flow.
- Added direct-link throughput: `T_c(x) = 6 + x`. The direct base capacity gives one extra unit, and every extra unit on `o -> d` raises throughput by one.

Exercise 2:

- Katz centrality: `x = (I - beta W)^(-1) mu`. A node is central if many walks reach it.
- Katz walk expansion: `x = sum_{k>=0} beta^k W^k mu` when the spectral condition holds. This shows explicitly that Katz counts discounted walks of all lengths.
- Spectral condition: `beta < 1/lambda_max(W)`. This ensures the Katz inverse and walk expansion are well defined.
- PageRank update: `x_i(t+1) = (1-beta)mu_i + beta sum_{j in N(i)} x_j(t)/deg(j)`. Each node receives degree-normalized mass from neighbors.
- Matrix PageRank form: `x(t+1) = beta P x(t) + (1-beta)mu`. This is the random-walk-with-teleportation interpretation.
- PageRank at `beta = 1`: `pi_i = deg(i)/(2|E|)` for connected undirected graphs. Stationary probability is proportional to degree.

Exercise 3:

- Flow conservation: `Bf = nu`. Link flows must match node inflows and outflows.
- Delay: `tau_e(f_e) = l_e/(1 - f_e/C_e)`. Delay increases to infinity near capacity.
- Capacity-utilization ratio: `rho_e = f_e/C_e`. This tells how close a link is to its capacity boundary.
- Social optimum: `min sum_e f_e tau_e(f_e)`. This minimizes total travel time.
- CVXPY social objective: `sum_e (l_e C_e/(1-f_e/C_e) - l_e C_e)`. This is algebraically equal to `sum_e f_e tau_e(f_e)`.
- Beckmann potential: `sum_e -l_e C_e log(1-f_e/C_e)`. Its derivative is the private link delay.
- Price of anarchy: `PoA = cost_Wardrop / cost_social`. It measures inefficiency from selfish routing.
- Marginal social cost: `MSC_e(f_e) = tau_e(f_e) + f_e tau'_e(f_e)`. It is the system-level cost of adding one more unit of flow to a link.
- Marginal-cost toll: `omega_e = f_e^* tau'_e(f_e^*)`. This charges the congestion externality.
- Tolled perceived cost: `tau_e(f_e) + omega_e`. Wardrop users respond to this modified private cost.
- Additional-delay objective: `sum_e f_e(tau_e(f_e)-l_e)`. This counts only congestion delay above free flow.
- Additional-delay toll: `omega_e^add = -l_e + f_e^* tau'_e(f_e^*)`. This can be negative because free-flow time is subtracted.

---

## 7. Numerical results cheat sheet

| Quantity | Value |
|---|---|
| Exercise 1 min cut | `5` |
| Exercise 1 max flow | `5` |
| Exercise 1 part (b) throughput | `T_b(x) = 5 + ceil(x/2)` |
| Exercise 1 part (c) direct-link throughput | `T_c(x) = 6 + x` |
| Katz `lambda_max(W)` | `5.0710614673` |
| Katz threshold `1/lambda_max(W)` | `0.1971973731` |
| Katz top nodes | `n6` first; `n1,...,n5` tied second |
| PageRank top nodes, `beta=0.15` | `n9` first, `n6` second |
| `PR(n6)-PR(n9)`, beta `0` | `0.0000000000` |
| `PR(n6)-PR(n9)`, beta `1/4` | `-0.0394710887` |
| `PR(n6)-PR(n9)`, beta `1/2` | `-0.0652885611` |
| `PR(n6)-PR(n9)`, beta `3/4` | `-0.0805114506` |
| `PR(n6)-PR(n9)`, beta `1` | `0.0208333333` |
| Exercise 3 shortest path | `1 -> 2 -> 3 -> 9 -> 13 -> 17` |
| Shortest-path link indices | `1, 2, 12, 9, 25` |
| Shortest-path free-flow time | `0.559833` |
| Exercise 3 max flow | `22448` |
| `nu_new_1` | `16806` |
| `nu_new_17` | `-16806` |
| Social optimum objective | `26142.669150` |
| Wardrop objective | `26495.321759` |
| Price of anarchy | `1.013490` |
| Marginal-cost tolled flow error | `2.993561e-1`, relative `1.325734e-5` |
| Additional-delay objective | `15350.353857` |
| Additional-delay tolled flow error | `5.658613e-1`, relative `2.507184e-5` |

---

## 8. Weak points and how to defend them

### 8.1 Direct link `o -> d` in Exercise 1(c)

Risk: It may be considered too trivial.

Defense: The assignment does not forbid direct source-sink arcs. The code includes it under the literal candidate set. It is optimal because `o -> d` crosses every `o-d` cut. If the intended interpretation excluded direct source-sink arcs, the candidate set must be restricted and the optimization recomputed.

### 8.2 Graph reconstruction in Exercise 2

Risk: If the figure is ambiguous, edge reconstruction matters.

Defense: The report explicitly states the edge set. The centrality values are tied to that reconstructed graph. No visually ambiguous edges were included.

### 8.3 PageRank at beta = 1

Risk: PageRank iteration at `beta = 1` may not always converge.

Defense: For this undirected connected non-bipartite graph, the simple random walk has a unique stationary distribution and convergence is valid. The stationary distribution is proportional to degree.

### 8.4 Numerical optimization in Exercise 3

Risk: Solver approximation and the artificial cap `0.999 C`.

Defense: The delay is infinite at capacity, so the cap keeps the solver in the finite domain. Solver status was checked as `optimal`, and the toll-recovered flows are extremely close to their target flows.

### 8.5 Negative additional-delay tolls

Risk: Negative tolls may look wrong.

Defense: The toll formula contains `-l_e`. Negative tolls are subsidies or offsets and are mathematically valid for decentralizing the additional-delay optimum.

---

## 9. Five-minute oral defense script

In this homework I solved three network problems: max-flow and cuts, centrality measures, and traffic assignment.

In Exercise 1, the network is a directed graph with source `o`, sink `d`, and six capacitated links. I first enumerated all source-sink cuts. Since there are three intermediate nodes, there are `2^3 = 8` cuts. For each cut I used the directed convention: only edges leaving the source-side set count. The minimum cut capacity is `5`, and the maximum flow computed on the same network is also `5`, which verifies the max-flow/min-cut theorem. Then I studied capacity augmentation. If I can add `x` integer units only on existing links, the best throughput follows `T_b(x) = 5 + ceil(x/2)`. The reason is that extra capacity must often be added on complementary bottlenecks before it increases the full source-sink throughput. Finally, I allowed one new directed link with base capacity `1`. Under the literal interpretation of the assignment, I included `o -> d`, since direct source-sink links were not forbidden. That link crosses every source-sink cut, so with all extra capacity on it the throughput is `T_c(x) = 6 + x`.

In Exercise 2, I reconstructed the graph as an undirected simple graph with 15 nodes. The main structure is a clique on `n1,...,n6`, a path from `n6` to `n9` through `n7,n8`, a star centered at `n9`, and `n15` attached to `n6`. For Katz centrality I used `x = (I - beta W)^(-1) mu` with uniform `mu` and `beta = 0.15`. The largest eigenvalue is `5.0710614673`, so `1/lambda_max = 0.1971973731`, and the spectral condition is satisfied. Katz ranks `n6` first because Katz counts raw walks, and `n6` benefits from the dense clique and bridge position. For distributed PageRank, each node sends `x_j/deg(j)` to its neighbors. With `beta = 0.15`, PageRank ranks `n9` first and `n6` second. This happens because degree-one leaves send all their mass to `n9`, while the clique nodes split their mass among many neighbors. In the beta sensitivity test, `PR(n6)-PR(n9)` is not monotone: it is zero at `beta = 0`, negative for intermediate beta, and positive at `beta = 1`, where stationary mass is proportional to degree.

In Exercise 3, I used the traffic network data with 17 nodes and 28 links. The incidence matrix uses `+1` for the tail and `-1` for the head, so `Bf = nu` gives node imbalances. The shortest path from node 1 to node 17 using free-flow travel times is `1 -> 2 -> 3 -> 9 -> 13 -> 17`, with total free-flow time `0.559833`. The maximum flow using only capacities is `22448`. From the measured flow, I formed a reduced demand vector with `16806` entering at node 1 and `-16806` at node 17. The delay function is `tau_e(f_e) = l_e/(1-f_e/C_e)`, so it tends to infinity near capacity.

For routing, the social optimum minimizes total travel time and gives `26142.669150`. The Wardrop equilibrium is computed through the Beckmann potential and gives total travel time `26495.321759`, so the price of anarchy is `1.013490`, about `1.349%` worse than the social optimum. Then I computed marginal-cost tolls `omega_e = f_e^* tau'_e(f_e^*)`, and the tolled Wardrop flow matches the social optimum with `L2` error about `0.299`, tiny relative to flows of thousands. Finally, I solved the additional-delay objective, which minimizes only congestion delay above free-flow time. Its optimum is `15350.353857`, and its tolls can be negative because the formula includes `-l_e`.

---

## 10. Ten-minute oral defense script

This homework has three exercises, all centered on how network structure affects flows, rankings, and equilibria.

I start with Exercise 1. The graph has nodes `{o,a,b,c,d}`, source `o`, sink `d`, and directed edges `o -> a`, `a -> d`, `o -> b`, `b -> c`, `c -> d`, and `a -> b`, with capacities `3,3,3,3,2,1`. A feasible flow is nonnegative, respects these capacities, and satisfies conservation at `a,b,c`. An `o-d` cut is a set `S` containing `o` but not `d`. Since `a,b,c` can each be either inside or outside the source-side set, there are exactly `8` cuts. For directed cuts, I count only edges from `S` to `V \ S`.

The two minimum cuts are `{o,b,c}` with crossing edges `e1,e5`, and `{o,a,b,c}` with crossing edges `e2,e5`. Both have capacity `5`. The computed max flow is also `5`, matching the max-flow/min-cut theorem. So the minimum capacity that must be removed to disconnect the source from the sink is `5`.

For part (b), I add an integer extra-capacity budget `x` to existing links. The variables `y_i` represent extra units on edge `e_i`. The exact enumeration gives `T_b(x) = 5 + ceil(x/2)` for the tested budgets. This is not linear with slope one because an extra unit may relax only one side of a bottleneck. In part (c), I add one new directed link with base capacity `1`, then allocate the extra budget. The code considered all missing directed arcs. Since the assignment did not explicitly forbid `o -> d`, I included it. That link crosses every source-sink cut, so the direct-link formula is `T_c(x) = 6 + x`. At `x = 0`, the CSV reports `b -> d` as one best link, but `o -> d` also achieves throughput `6`; for `x >= 1`, the best reported link is `o -> d`.

In Exercise 2, I reconstructed an undirected graph with 15 nodes. The important structure is that `n1,...,n6` form a clique, `n6` connects to `n7` and `n15`, then the path `n6-n7-n8-n9` reaches a star centered at `n9` with leaves `n10,...,n14`.

For Katz centrality, the formula is `x = (I - beta W)^(-1) mu`, where `W` is the adjacency matrix and `mu` is uniform, `1/15`. The spectral condition is `beta < 1/lambda_max(W)`. In my graph, `lambda_max(W) = 5.0710614673`, hence the threshold is `0.1971973731`, and `beta = 0.15` is valid. Katz ranks `n6` first with value `0.3177933766`, while `n1,...,n5` are tied second. This makes sense because Katz counts raw walks. The clique creates many walks, and `n6` also bridges from the clique to the rest of the graph.

For PageRank, I implemented the distributed update `x_i(t+1) = (1-beta)mu_i + beta sum_{j in N(i)} x_j(t)/deg(j)`. Each node sends its current value divided by its degree to all neighbors. With `beta = 0.15`, the iteration converged in `15` iterations. PageRank ranks `n9` first with value `0.1059574692` and `n6` second with value `0.0801149613`. The reason is degree normalization: the degree-one leaves attached to `n9` send all their mass to `n9`, while nodes in the clique split their mass among several neighbors. This explains why Katz and PageRank disagree.

The beta sensitivity test also supports this interpretation. At `beta = 0`, all nodes have PageRank `1/15`. At beta `1/4`, `1/2`, and `3/4`, `PR(n6)-PR(n9)` is negative: `-0.0394710887`, `-0.0652885611`, and `-0.0805114506`. At `beta = 1`, PageRank becomes the stationary distribution of the simple random walk. Because the graph is connected and non-bipartite, convergence is valid; because the graph is undirected, stationary probability is proportional to degree. So `n6` gets `7/48`, `n9` gets `6/48`, and the difference is `1/48`.

In Exercise 3, I use the traffic network with 17 nodes and 28 links. The incidence matrix has one column per link, with `+1` at the tail and `-1` at the head, so `Bf = nu` gives net node inflows. The shortest path with free-flow times is `1 -> 2 -> 3 -> 9 -> 13 -> 17`, using links `1,2,12,9,25`, with time `0.559833`. The maximum flow from node 1 to node 17 using capacities is `22448`.

From the measured flow vector, I computed `nu = Bf`. For the assignment routing problems, I then used the reduced single-demand vector with `nu_new_1 = 16806`, `nu_new_17 = -16806`, and all intermediate entries zero. The link delay is `tau_e(f_e) = l_e/(1 - f_e/C_e)`, with domain `0 <= f_e < C_e`. Numerically I constrained `f_e <= 0.999 C_e` to avoid the infinite boundary.

The social optimum minimizes `sum_e f_e tau_e(f_e)` subject to conservation and capacity constraints. The equivalent convex expression is `sum_e (l_e C_e/(1-f_e/C_e) - l_e C_e)`. The solver status is `optimal`, and the total travel time is `26142.669150`.

The Wardrop equilibrium is the selfish-routing equilibrium where no user can improve by changing route. For separable link costs, it is computed by minimizing the Beckmann potential. Here the integral is `-l_e C_e log(1-f_e/C_e)`. The Wardrop total travel time is `26495.321759`, so the price of anarchy is `1.013490`. This means selfish routing is about `1.349%` worse than the social optimum.

For tolls, the social marginal cost is `tau_e(f_e) + f_e tau'_e(f_e)`, while users only perceive `tau_e(f_e)`. Therefore the marginal-cost toll is `omega_e = f_e^* tau'_e(f_e^*)`, with `tau'_e(f_e) = (l_e/C_e)/(1-f_e/C_e)^2`. Solving the tolled Wardrop problem gives an `L2` flow error of `2.993561e-1`, relative error `1.325734e-5`, which confirms the tolls recover the social optimum numerically.

Finally, I solved an additional-delay objective `psi_e(f_e) = f_e(tau_e(f_e)-l_e)`. The optimum value is `15350.353857`. The decentralizing toll is `omega_e^add = psi'_e(f_e^*) - tau_e(f_e^*) = -l_e + f_e^* tau'_e(f_e^*)`. Some tolls are negative; this is mathematically expected because of the `-l_e` term and can be interpreted as a subsidy or offset. The tolled flow matches the additional-delay optimum with `L2` error `5.658613e-1`, relative error `2.507184e-5`.

---

## 11. Rapid Q&A bank

Exercise 1:

1. What is the source? `o`.
2. What is the sink? `d`.
3. How many nodes are in Exercise 1? Five.
4. How many original directed edges? Six.
5. What is an `o-d` cut? A set containing `o` and excluding `d`.
6. Why are there 8 cuts? Three intermediate nodes give `2^3` choices.
7. What edges count in a directed cut? Only edges leaving `S`.
8. What is the min cut? `5`.
9. Which cuts attain capacity 5? `{o,b,c}` and `{o,a,b,c}`.
10. What is the max flow? `5`.
11. What theorem connects them? Max-flow/min-cut.
12. What is `T_b(x)`? `5 + ceil(x/2)`.
13. Why not `5 + x` in part (b)? Capacity must be added across bottleneck structures.
14. What is the best direct-link formula in part (c)? `6 + x`.
15. Why is `o -> d` defensible? It was not explicitly forbidden.
16. Is the part (b) optimizer unique? Not necessarily.
17. How was global optimality checked? Exhaustive integer allocation search.
18. What does base capacity mean in part (c)? The new link starts with capacity `1` before extra budget.

Exercise 2:

1. Is the graph directed? No, undirected.
2. Is the graph simple? Yes.
3. What nodes form a clique? `n1,...,n6`.
4. What is special about `n6`? It is in the clique and bridges to the rest.
5. What is special about `n9`? It is the center of a star-like part.
6. What path connects `n6` and `n9`? `n6-n7-n8-n9`.
7. What is Katz formula? `x = (I-beta W)^(-1)mu`.
8. What is `W`? The adjacency matrix.
9. What is `mu`? Uniform intrinsic centrality.
10. What is `lambda_max(W)`? `5.0710614673`.
11. What is the Katz beta threshold? `0.1971973731`.
12. Is `beta = 0.15` valid for Katz? Yes.
13. Who is top Katz node? `n6`.
14. Why is `n6` top in Katz? Many raw walks involve it.
15. What is PageRank update? `(1-beta)mu_i + beta sum x_j/deg(j)`.
16. Who is top PageRank node at `beta = 0.15`? `n9`.
17. Why is `n9` top in PageRank? Leaves send all mass to it.
18. What happens at `beta = 0`? All nodes equal `1/15`.
19. What happens at `beta = 1`? Stationary random-walk distribution.
20. Why is `PR(n6)` higher at `beta = 1`? Degree `7` versus degree `6`.
21. Is `PR(n6)-PR(n9)` monotone in tested beta values? No.
22. What does distributed mean? Nodes use local neighbor messages.
23. Would scaling `mu` change Katz ranking? No.
24. Why does PageRank degree-normalize? It models random-walk splitting among neighbors.

Exercise 3:

1. How many nodes? `17`.
2. How many links? `28`.
3. What is the source? Node `1`.
4. What is the sink? Node `17`.
5. What is `B`? Directed incidence matrix.
6. What does `+1` mean in `B`? Tail of a link.
7. What does `-1` mean in `B`? Head of a link.
8. What is `Bf = nu`? Flow conservation with external imbalance.
9. What is the shortest path? `1 -> 2 -> 3 -> 9 -> 13 -> 17`.
10. What are the shortest-path links? `1, 2, 12, 9, 25`.
11. What is shortest-path time? `0.559833`.
12. What is max flow? `22448`.
13. What is `nu_new_1`? `16806`.
14. What is `nu_new_17`? `-16806`.
15. Why set other `nu_new` entries to zero? Single OD model.
16. What is the delay function? `l_e/(1-f_e/C_e)`.
17. What happens near capacity? Delay tends to infinity.
18. Why use `0.999 C`? Numerical domain safeguard.
19. What is the social optimum objective? Total travel time.
20. What is social optimum value? `26142.669150`.
21. What is Wardrop value? `26495.321759`.
22. What is PoA? `1.013490`.
23. What does PoA mean? Selfish routing cost divided by optimal cost.
24. What is Beckmann potential? Integral of link delay functions.
25. What are marginal-cost tolls? `f_e^* tau'_e(f_e^*)`.
26. What is tolled social error? `2.993561e-1`.
27. What is additional-delay objective? `sum f_e(tau_e-l_e)`.
28. What is additional-delay optimum? `15350.353857`.
29. Why can additional-delay tolls be negative? They include `-l_e`.
30. What is additional-delay tolled error? `5.658613e-1`.

General/reproducibility:

1. Where are scripts? `src/exercise1.py`, `src/exercise2.py`, `src/exercise3.py`.
2. Where are CSV outputs? `results/`.
3. Where are plots? `figures/`.
4. Where are Exercise 3 data files? `data/*.mat`.
5. Where is the report source? `report/report.tex`.
6. What command runs Exercise 1? `python src/exercise1.py`.
7. What command runs Exercise 2? `python src/exercise2.py`.
8. What command runs Exercise 3? `python src/exercise3.py`.
9. Is randomness used? No intentional randomness.
10. What Python packages are needed? `numpy`, `scipy`, `networkx`, `cvxpy`, `matplotlib`, `pandas`.
11. What solver is used? CVXPY with CLARABEL, and ECOS fallback for some tolled cases.
12. How is the report compiled? `pdflatex -interaction=nonstopmode -output-directory=report report/report.tex`.

---

## 12. Reproducibility defense

If asked "How can I reproduce your results?", answer:

The repository is organized as follows:

```text
hw1/
  data/       Exercise 3 .mat input files
  src/        Python scripts for the three exercises
  results/    Generated CSV result tables
  figures/    Generated plots
  report/     LaTeX report source and compiled PDF if built
```

The data files are:

- `data/flow.mat`
- `data/capacities.mat`
- `data/traffic.mat`
- `data/traveltime.mat`

Install dependencies:

```bash
pip install numpy scipy networkx cvxpy matplotlib pandas
```

Run scripts from the project root:

```bash
python src/exercise1.py
python src/exercise2.py
python src/exercise3.py
```

These regenerate:

- Exercise 1 CSVs: `exercise1_cuts.csv`, `exercise1_b_results.csv`, `exercise1_c_results.csv`.
- Exercise 2 CSVs: `exercise2_katz.csv`, `exercise2_pagerank_beta015.csv`, `exercise2_pagerank_beta_sensitivity.csv`.
- Exercise 3 CSVs: `exercise3_shortest_path.csv`, `exercise3_maxflow.csv`, `exercise3_nu.csv`, `exercise3_social_optimum.csv`, `exercise3_wardrop.csv`, `exercise3_tolls_total_travel_time.csv`, `exercise3_additional_delay.csv`.
- Figures in `figures/`.

Compile the report:

```bash
pdflatex -interaction=nonstopmode -output-directory=report report/report.tex
```

The convex optimization is done through CVXPY. The code primarily uses CLARABEL and may fall back to ECOS for tolled Wardrop problems. There is no intentional randomness, so rerunning the scripts should reproduce the same results up to solver tolerances.

The repository already has a README, and this section is aligned with it.

---

## 13. Final pre-oral checklist

- [ ] Report compiles.
- [ ] Code runs from a clean environment.
- [ ] CSV files match report.
- [ ] Figures match report.
- [ ] `.mat` files are included in `data/`.
- [ ] Solver is installed and CVXPY can find CLARABEL/ECOS as needed.
- [ ] Direct link `o -> d` defense is prepared.
- [ ] Exercise 1 cut table is understood.
- [ ] Exercise 1 part (b) edge cases for odd/even `x` are understood.
- [ ] Exercise 2 graph edge set is memorized.
- [ ] Katz spectral condition is understood.
- [ ] PageRank beta extremes are understood.
- [ ] Difference between Katz and PageRank is clear.
- [ ] Incidence matrix sign convention is clear.
- [ ] Wardrop versus social optimum is understood.
- [ ] Beckmann potential derivation is understood.
- [ ] Price of anarchy interpretation is memorized.
- [ ] Marginal-cost toll derivation is understood.
- [ ] Additional-delay toll derivation is understood.
- [ ] Negative additional-delay toll defense is prepared.
