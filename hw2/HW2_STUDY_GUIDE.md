# HW2 Study Guide — Network Dynamics and Learning

## 1. Big picture

Homework 2 studies stochastic dynamics on directed weighted networks.

The first model is a continuous-time random walk on five nodes. A particle waits an exponential time at its current node, then jumps to a neighbor with probability proportional to the outgoing rate. This is a continuous-time Markov chain (CTMC), described by a transition-rate matrix `Lambda`, an exit-rate vector `omega`, and a generator `Q`.

The same graph is then interpreted as a French-DeGroot opinion network. Opinions evolve according to a linear differential equation. The graph structure, especially strongly connected components and sink components, determines whether all opinions converge to one consensus value.

Problem 2 extends the closed CTMC to many independent particles. In the particle perspective, one follows each particle separately. In the node perspective, the state is the count vector of particles at each node. The average final counts are compared with `N*pi`, where `pi` is the stationary distribution of the single-particle CTMC.

Problem 3 studies an open network. Particles arrive from outside at node `o`, move through an acyclic network, and exit at node `d`. Two service rules are compared. Under proportional service, total service rate grows with the number of particles at a node. Under fixed service, each node has a bounded service capacity, so bottlenecks can create blow-up.

## 2. Repository and codebase structure

The current HW2 folder contains:

- `HW2_CONTEXT.md`: project scope, assignment conventions, expected folder structure, model definitions, planned outputs, and final checklist.
- `RESULTS_LOG.md`: factual log of runs, seeds, assumptions, numerical results, generated figures, and remaining issues.
- `README.md`: how to set up Python, run scripts, and find outputs.
- `src/`: Python source code.
- `results/`: generated CSV files.
- `figures/`: generated PNG figures.
- `report/report.tex`: report content. Use it for mathematical/content consistency only, not compilation issues.

### `src/constants.py`

Defines shared constants:

- `NODES = ["o", "a", "b", "c", "d"]`.
- `NODE_TO_INDEX`.
- `LAMBDA`, the closed-network rate matrix for Problems 1 and 2.
- `LAMBDA_OPEN`, the open-network routing matrix for Problem 3.
- `DEFAULT_SEED = 20260604`.
- `PROJECT_ROOT`, `FIGURES_DIR`, and `RESULTS_DIR`.

Mathematically, this file fixes the node order and the matrix orientation used everywhere: `Lambda[i,j]` is the rate from source node `i` to destination node `j`.

### `src/problem1_single_ctmc.py`

Implements single-particle CTMC routines:

- `exit_rates(lambda_matrix)`: computes `omega_i = sum_j Lambda[i,j]`.
- `generator(lambda_matrix)`: computes `Q = Lambda - diag(omega)`.
- `jump_matrix(lambda_matrix)`: computes the embedded jump matrix `P_jump`.
- `stationary_distribution(lambda_matrix)`: solves `pi Q = 0`, `sum_i pi_i = 1`.
- `simulate_return_times(...)`: Monte Carlo return times to a start node, including the initial holding time.
- `theoretical_return_time(...)`: uses `1 / (pi_i omega_i)`.
- `simulate_hitting_times(...)`: Monte Carlo hitting times.
- `theoretical_hitting_times(...)`: solves the CTMC hitting-time linear system.

Generated CSVs through `run_problem1.py`:

- `results/problem1_return_time_a.csv`
- `results/problem1_hitting_o_to_d.csv`

### `src/problem1_opinion.py`

Implements continuous-time French-DeGroot routines:

- `laplacian(lambda_matrix)`: returns `L = diag(Lambda @ 1) - Lambda`.
- `left_consensus_vector(lambda_matrix, component=None)`: normalized left nullvector for consensus values.
- `matrix_exponential_action(matrix, vector, times)`: evaluates `exp(matrix*t) @ vector`.
- `simulate_fdg(lambda_matrix, x0, times)`: computes `exp(-Lt)x0`.
- `remove_edges(lambda_matrix, node_to_index, edges)`: sets directed edge weights to zero.
- `strongly_connected_components(lambda_matrix)`: Tarjan SCC algorithm.
- `sink_components(lambda_matrix)`: SCCs with no outgoing edge to another SCC.
- `graph_summary(lambda_matrix, node_names)`: SCC and consensus-condition summary.
- `asymptotic_state(lambda_matrix, x0)`: approximates the limit by evaluating at large time.

Generated CSVs and figures through `run_problem1.py`:

- `results/problem1_fdg_original_consensus.csv`
- `results/problem1_consensus_variance.csv`
- `results/problem1_removed_edges_g_summary.csv`
- `results/problem1_removed_edges_h_summary.csv`
- `figures/problem1_fdg_original_trajectory.png`
- `figures/problem1_consensus_value_histogram.png`
- `figures/problem1_fdg_removed_edges_g.png`
- `figures/problem1_fdg_removed_edges_h_case1.png`
- `figures/problem1_fdg_removed_edges_h_case2.png`

### `src/problem2_many_particles.py`

Implements many-particle closed-network simulations:

- `ctmc_components(lambda_matrix)`: returns `omega`, `Q`, and `P_jump`.
- `stationary_distribution(q_matrix)`: solves `pi Q = 0`.
- `simulate_return_time(...)`: one particle's return time to the start node.
- `simulate_particle_return_run(...)`: average return time across `N` particles.
- `simulate_particle_return_runs(...)`: repeated Monte Carlo batches.
- `simulate_node_count_run(...)`: event-based count-vector simulation.
- `simulate_node_count_runs(...)`: repeated node-count simulations, returning mean, standard deviation, and final counts.

Generated CSVs and figures through `run_problem2.py`:

- `results/problem2_particle_return_times.csv`
- `results/problem2_node_time_series.csv`
- `results/problem2_node_final_counts.csv`
- `figures/problem2_node_counts_over_time.png`
- `figures/problem2_stationary_distribution_comparison.png`

The current source, saved CSVs, `RESULTS_LOG.md`, and report are intended to agree on 5000 node-perspective Monte Carlo runs for Problem 2.

### `src/problem3_open_network.py`

Implements the open-network simulator:

- `SimulationResult`: dataclass storing sample times, count matrix, and event count.
- `open_network_components(lambda_open)`: computes `omega`, `P_jump`, and `P_jump_cdf`; overrides `omega_d = 7/4`.
- `service_rates(counts, omega, service_mode)`: proportional or fixed service rates.
- `simulate_open_network(...)`: event-based simulation with arrivals, services, routing, and removal at `d`.
- `simulate_replicates(...)`: repeated simulations.
- `total_population_slope(...)`: fitted slope of total population over the final part of a trajectory.
- `theoretical_loads(lambda_value)`: fixed-rate traffic loads.
- `fixed_rate_capacity(omega)`: bottleneck capacity.
- `mean_absorption_times(lambda_open, omega)`: expected remaining lifetime of one particle.

Generated CSVs and figures through `run_problem3.py`:

- `results/problem3_proportional_timeseries_lambda100.csv`
- `results/problem3_proportional_stability_scan.csv`
- `results/problem3_fixed_timeseries_lambda2.csv`
- `results/problem3_fixed_stability_scan.csv`
- `figures/problem3_proportional_lambda100_counts.png`
- `figures/problem3_proportional_stability_scan.png`
- `figures/problem3_fixed_lambda2_counts.png`
- `figures/problem3_fixed_stability_scan.png`

### Run scripts and helpers

- `src/run_problem1.py`: runs all Problem 1 computations, writes CSVs, figures, and Problem 1 log entries.
- `src/run_problem2.py`: runs Problem 2 simulations, writes CSVs, figures, and Problem 2 log entries.
- `src/run_problem3.py`: runs Problem 3 simulations and scans, writes CSVs, figures, and Problem 3 log entries.
- `src/utils.py`: CSV writers and normal 95 percent confidence intervals.
- `src/plotting.py`: `save_figure`, which saves Matplotlib figures with tight layout and 200 dpi.

## 3. Global notation and conventions

Node order is always:

```text
[o, a, b, c, d]
```

The matrix convention is:

```text
Lambda[i,j] = rate or weight from source node i to destination node j
```

Rows are source nodes. Columns are destination nodes. This matters for both CTMC transitions and edge removals.

The exit-rate vector is:

```text
omega = Lambda @ 1
```

The CTMC generator is:

```text
Q = Lambda - diag(omega)
```

So:

- `Q[i,j] = Lambda[i,j]` for `i != j`.
- `Q[i,i] = -omega_i`.
- each row of `Q` sums to zero.

The embedded jump matrix is:

```text
P_jump[i,j] = Lambda[i,j] / omega_i
```

for nodes with `omega_i > 0`. It is the transition matrix observed only at actual jump times.

Do not confuse `P_jump` with a uniformized transition matrix. A uniformized matrix would choose a global rate `Omega >= max_i omega_i` and define:

```text
P_uniformized = I + Q / Omega
```

This includes self-loops. `P_jump` usually has no self-loop unless `Lambda[i,i] > 0`. Mixing them changes holding-time behavior and gives wrong simulations.

## 4. Problem 1 theory

### 4.1 Continuous-time random walk

A CTMC random walk evolves as follows:

1. If the particle is at node `i`, it waits for a holding time with exponential distribution of rate `omega_i`.
2. After the holding time, it jumps to node `j` with probability `Lambda[i,j] / omega_i`.

The exponential holding time appears because independent Poisson clocks can be attached to outgoing edges. If edge `i -> j` has rate `Lambda[i,j]`, the first clock to ring determines the next jump. The minimum of independent exponential clocks is exponential with rate equal to the sum of rates, `omega_i`, and the probability that clock `j` rings first is `Lambda[i,j] / omega_i`.

The memoryless property means that after each jump the process restarts probabilistically from the new node. This is why the CTMC can be simulated event by event.

In the code:

- `exit_rates` computes `omega`.
- `jump_matrix` computes destination probabilities.
- `simulate_return_times` and `simulate_hitting_times` repeatedly draw exponential holding times and destinations.

### 4.2 Return time

`T_a^+` is the first positive return time to node `a`. In this codebase, the convention is:

- start at node `a`;
- include the initial holding time in `a`;
- require the particle to leave `a`;
- stop at the next entrance to `a`.

The saved result `results/problem1_return_time_a.csv` explicitly records `initial_holding_time_included = yes`.

For an irreducible CTMC, the mean time between entrances to state `i` is:

```text
E_i[T_i^+] = 1 / (pi_i omega_i)
```

Here `pi_i` is the stationary fraction of time spent in state `i`, and `1/omega_i` is the mean holding time during one visit to `i`. Since stationary time fraction equals mean time in `i` per cycle divided by mean cycle length, the mean cycle length is `1/(pi_i omega_i)`.

The code computes this in `theoretical_return_time`.

### 4.3 Hitting time

The hitting time of node `d` is:

```text
T_d = inf{t >= 0 : X(t) = d}
```

Let:

```text
h_i = E_i[T_d]
```

Then:

```text
h_d = 0
sum_j Q[i,j] h_j = -1,  for i != d
```

This comes from first-step analysis. From node `i`, the chain waits mean time `1/omega_i`, then jumps according to `P_jump`. Equivalently:

```text
h_i = 1/omega_i + sum_j P_jump[i,j] h_j
```

The code solves the generator linear system in `theoretical_hitting_times` and validates it with Monte Carlo in `simulate_hitting_times`.

### 4.4 French-DeGroot dynamics

The weighted graph interpretation uses the same `Lambda`. The Laplacian is:

```text
L = diag(Lambda @ 1) - Lambda
```

The continuous-time French-DeGroot dynamics are:

```text
dx/dt = -Lx
```

Node `i` updates toward a weighted average of its out-neighbors:

```text
dx_i/dt = sum_j Lambda[i,j] (x_j - x_i)
```

The solution is:

```text
x(t) = exp(-Lt) x(0)
```

Consensus means all components converge to the same value. For directed graphs under this convention, the asymptotic behavior is determined by sink strongly connected components. If there is a unique sink SCC reachable from all nodes, the system converges to global consensus. If there are multiple sink SCCs, global consensus is not guaranteed for every initial condition.

The consensus value is:

```text
x* = sum_i pi_i x_i(0)
```

where `pi` is a normalized nonnegative left nullvector:

```text
pi^T L = 0,  sum_i pi_i = 1
```

The code computes SCCs, sink SCCs, and consensus vectors in `problem1_opinion.py`.

### 4.5 Variance of consensus value

If the initial opinions are independent random variables, and:

```text
x* = sum_i pi_i x_i(0)
```

then:

```text
Var(x*) = sum_i pi_i^2 Var(x_i(0))
```

The cross-covariance terms vanish because the initial variables are independent. In this homework the variances are:

- node `o`: 1
- node `a`: 2
- node `b`: 2
- node `c`: 2
- node `d`: 1

The code validates the formula by sampling many random initial conditions, computing consensus values, and comparing the empirical variance with the theoretical variance.

### 4.6 Edge-removal cases

Removing directed edges changes the graph topology and therefore the limiting opinion behavior.

If edge removals create multiple sink SCCs, different parts of the graph may converge to different limiting values. Transient nodes can converge to weighted combinations determined by which sink components they can reach.

If after edge removals there is a unique sink SCC reachable from all nodes, global consensus still occurs, but the consensus value depends only on the invariant weights supported by that sink SCC.

The code checks:

- which directed edges are removed;
- the SCCs of the modified graph;
- the sink SCCs;
- whether there is one sink reachable from all nodes;
- representative asymptotic states from matrix exponential simulation.

## 5. Problem 2 theory

### 5.1 Many independent particles

Problem 2 uses `N = 100` particles on the same closed CTMC as Problem 1. The particles are independent copies of the same Markov process.

Because each particle follows the same CTMC, an average return-time estimate over many particles should match the single-particle return time from Problem 1, up to Monte Carlo error. More particles or more repetitions reduce the standard error.

### 5.2 Particle perspective

In the particle perspective, each particle is followed individually:

1. Start each particle at node `a`.
2. Simulate its CTMC trajectory.
3. Record its return time to `a`.
4. Average over 100 particles.
5. Repeat over Monte Carlo runs.

The current saved result is in `results/problem2_particle_return_times.csv`.

### 5.3 Node perspective

In the node perspective, the state is the count vector:

```text
N(t) = (N_o(t), N_a(t), N_b(t), N_c(t), N_d(t))
```

If node `i` contains `n_i(t)` particles, the total departure rate from node `i` is:

```text
n_i(t) omega_i
```

The total event rate is:

```text
sum_i n_i(t) omega_i
```

At an event:

1. Choose the source node with probability proportional to `n_i(t) omega_i`.
2. Choose the destination from row `i` of `P_jump`.
3. Decrease the source count by 1 and increase the destination count by 1.

The code simulates this event process and samples counts on a time grid until `T = 60`.

### 5.4 Stationary distribution comparison

The stationary distribution `pi` of the single-particle CTMC satisfies:

```text
pi Q = 0
sum_i pi_i = 1
```

For `N` independent particles at stationarity, the expected count at node `i` is:

```text
N pi_i
```

Final counts at finite time need not exactly equal `N*pi`, because:

- the simulation horizon is finite;
- the initial condition is not stationary;
- Monte Carlo estimates have sampling error.

The saved final-count comparison is in `results/problem2_node_final_counts.csv`.

## 6. Problem 3 theory

### 6.1 Open network

The open network has external arrivals at node `o`. The routing matrix is `LAMBDA_OPEN`, with row-source/column-destination convention. Particles move through nodes and eventually exit at node `d`.

The service-rate vector is:

```text
omega = Lambda_open @ 1
```

with the assignment override:

```text
omega_d = 7/4
```

At node `d`, a service event removes one particle from the system if `d` is nonempty.

### 6.2 Proportional-rate service

In proportional service:

```text
service rate at node i = omega_i N_i(t)
```

This means service capacity grows with the number of particles. If many particles accumulate, departures become faster.

The event simulator includes:

- external arrivals at rate `lambda`;
- service events at rates `omega_i N_i(t)`;
- routing for nodes `o,a,b,c`;
- removal at node `d`.

The theoretical intuition is that the network behaves like a linear infinite-server open network. Each arriving particle has a finite expected lifetime before exiting. Therefore any finite external arrival rate gives finite expected total population. The scan still provides useful numerical evidence for tested finite lambdas, but it is not the theoretical proof.

### 6.3 Fixed-rate service

In fixed service:

```text
service clock rate at node i = omega_i
```

This rate does not grow with the queue size. If a node is empty when its clock ticks, nothing happens. If arrivals and routed traffic exceed a node's fixed capacity, the queue at that node grows.

The theoretical loads per unit external arrival rate are:

```text
(1, 1/2, 5/8, 3/4, 1)
```

The fixed-rate capacity by node is `omega_i / load_i`. With `omega = (2,1,1,1,7/4)`, the capacities are:

```text
(2, 2, 1.6, 4/3, 1.75)
```

The bottleneck is node `c`, so the theoretical largest stable input rate is below:

```text
lambda = 4/3
```

The long horizon `T = 6000` is used because near a stability threshold, queue growth can be slow and noisy.

### 6.4 Blow-up criterion

In this homework, "blow-up" means sustained growth of total population or queue sizes rather than fluctuation around a stable finite scale.

Simulation-based indicators include:

- large final total population;
- positive fitted slope of total population over the second half of the run;
- repeated behavior over multiple replicates;
- concentration of growth at bottleneck nodes.

Theoretical indicators include:

- finite expected particle lifetime for proportional service;
- traffic intensity below or above service capacity for fixed service.

One single stochastic trajectory is not enough to prove stability. It can illustrate behavior, but conclusions should combine theory and replicated scans.

## 7. Generated results and figures

### Result files in `results/`

`problem1_return_time_a.csv`

- Supports Problem 1(a)-(b).
- Generated by `src/run_problem1.py`.
- Contains return-time convention, seed, number of runs, simulation mean, standard deviation, standard error, 95 percent CI, theoretical value, difference, relative error, and theory method.

`problem1_hitting_o_to_d.csv`

- Supports Problem 1(c)-(d).
- Generated by `src/run_problem1.py`.
- Contains start/target nodes, seed, number of runs, simulation summary, theoretical hitting time, difference, relative error, and theory method.

`problem1_fdg_original_consensus.csv`

- Supports Problem 1(e).
- Generated by `src/run_problem1.py`.
- Contains the French-DeGroot convention, initial condition, SCCs, sink components, consensus condition, consensus vector, consensus value, final error, and figure path.

`problem1_consensus_variance.csv`

- Supports Problem 1(f).
- Generated by `src/run_problem1.py`.
- Contains random initial-condition variances, consensus vector, theoretical variance, Monte Carlo variance, Monte Carlo mean, error, relative error, and figure path.

`problem1_removed_edges_g_summary.csv`

- Supports Problem 1(g).
- Generated by `src/run_problem1.py`.
- Contains removed edges, initial condition, SCCs, sink components, consensus classification, asymptotic state, interpretation, and figure path.

`problem1_removed_edges_h_summary.csv`

- Supports Problem 1(h).
- Generated by `src/run_problem1.py`.
- Contains two cases, removed edges, initial conditions, SCCs, sink components, consensus classification, asymptotic states, spread of asymptotic state, interpretation, and figure paths.

`problem2_particle_return_times.csv`

- Supports Problem 2(a).
- Generated by `src/run_problem2.py`.
- Contains particle count, Monte Carlo repetitions, total simulated return times, mean of run averages, standard error, CI, comparison with Problem 1 simulation, and comparison with Problem 1 theory.

`problem2_node_time_series.csv`

- Supports Problem 2(b).
- Generated by `src/run_problem2.py`.
- Long-format table with `time`, `node`, `mean_count`, `std_count`, `standard_error`, seed, particle count, and Monte Carlo run count.

`problem2_node_final_counts.csv`

- Supports Problem 2(b) stationary comparison.
- Generated by `src/run_problem2.py`.
- Contains final average counts at `T = 60`, standard deviations, standard errors, stationary probabilities, expected counts `N*pi`, differences, relative errors, seed, particle count, and run count.

`problem3_proportional_timeseries_lambda100.csv`

- Supports Problem 3(a) main proportional run.
- Generated by `src/run_problem3.py`.
- Long-format table with `time`, `node`, `count`, `total_count`, `lambda`, service mode, and seed.

`problem3_proportional_stability_scan.csv`

- Supports Problem 3(a) stability scan.
- Generated by `src/run_problem3.py`.
- Contains tested lambdas, horizon, replicate count, final totals, fitted total-population slopes, drift threshold, simulation classification, and theoretical classification.

`problem3_fixed_timeseries_lambda2.csv`

- Supports Problem 3(b) main fixed-rate run.
- Generated by `src/run_problem3.py`.
- Long-format table with `time`, `node`, `count`, `total_count`, `lambda`, service mode, and seed.

`problem3_fixed_stability_scan.csv`

- Supports Problem 3(b) stability scan and bottleneck threshold.
- Generated by `src/run_problem3.py`.
- Contains tested lambdas, horizon, replicate count, final totals, slopes, drift and final-total thresholds, simulation classification, theoretical classification, theoretical capacity, and bottleneck node.

### Figures in `figures/`

`problem1_fdg_original_trajectory.png`

- Shows opinion trajectories for the original graph.
- Supports the claim that the original graph converges to consensus.

`problem1_consensus_value_histogram.png`

- Shows the Monte Carlo distribution of the random consensus value.
- Supports the variance comparison in Problem 1(f).

`problem1_fdg_removed_edges_g.png`

- Shows French-DeGroot trajectories after removing `(d,a)`, `(d,c)`, `(a,c)`, `(b,c)`.
- Supports the non-global-consensus interpretation.

`problem1_fdg_removed_edges_h_case1.png`

- Shows the first representative trajectory after removing `(b,o)` and `(d,a)`.
- Supports convergence to consensus for one initial condition.

`problem1_fdg_removed_edges_h_case2.png`

- Shows the second representative trajectory after removing `(b,o)` and `(d,a)`.
- Supports dependence of the consensus value on the initial condition.

`problem2_node_counts_over_time.png`

- Shows average node counts over time for `N = 100`.
- Supports the node-perspective simulation in Problem 2(b).

`problem2_stationary_distribution_comparison.png`

- Compares simulated final counts at `T = 60` with `N*pi`.
- Supports the stationary-distribution interpretation.

`problem3_proportional_lambda100_counts.png`

- Shows one proportional-service trajectory for `lambda = 100`.
- Supports the finite-population behavior discussion.

`problem3_proportional_stability_scan.png`

- Shows proportional-service scan summaries over tested lambdas.
- Supports the scan evidence that tested finite lambdas did not blow up.

`problem3_fixed_lambda2_counts.png`

- Shows one fixed-service trajectory for `lambda = 2`.
- Supports the fixed-rate blow-up/bottleneck discussion.

`problem3_fixed_stability_scan.png`

- Shows fixed-service scan summaries over tested lambdas.
- Supports the threshold around `4/3`.

No expected result or figure file from the inspected HW2 context is missing. The remaining issue to remember is not a missing file, but the need to keep source constants, CSV metadata, log entries, and report text synchronized after any rerun.

## 8. Reproducibility

Run commands from the `hw2/` folder.

Problem 1:

```powershell
python src/run_problem1.py
```

Expected outputs:

- six Problem 1 CSV files in `results/`;
- five Problem 1 figures in `figures/`;
- updated Problem 1 section in `RESULTS_LOG.md`.

Problem 1 seeds and run counts are documented in `RESULTS_LOG.md`. Return and hitting simulations currently use 100000 runs, and the consensus variance Monte Carlo uses 50000 runs.

Problem 2:

```powershell
python src/run_problem2.py
```

Expected outputs:

- `results/problem2_particle_return_times.csv`
- `results/problem2_node_time_series.csv`
- `results/problem2_node_final_counts.csv`
- `figures/problem2_node_counts_over_time.png`
- `figures/problem2_stationary_distribution_comparison.png`

The node-perspective run count should be 5000, matching the saved CSV metadata and report text.

Problem 3:

```powershell
python src/run_problem3.py
```

Expected outputs:

- four Problem 3 CSV files in `results/`;
- four Problem 3 figures in `figures/`;
- updated Problem 3 section in `RESULTS_LOG.md`.

To verify consistency:

- compare CSV metadata with `RESULTS_LOG.md`;
- compare report numbers with CSV values;
- confirm figure files exist;
- confirm all random seeds are documented;
- confirm source constants match logged run parameters.

Do not focus on report compilation.

## 9. Oral exam preparation

**What is the difference between `Lambda`, `Q`, and `P_jump`?**

`Lambda` contains transition rates from source node `i` to destination node `j`. `Q` is the CTMC generator: off-diagonal entries are rates, diagonal entries are negative exit rates, and rows sum to zero. `P_jump` is the embedded discrete-time transition matrix used after a jump occurs: `P_jump[i,j] = Lambda[i,j] / omega_i`.

**Why are holding times exponential?**

In a CTMC, outgoing transitions are governed by Poisson clocks. The minimum of exponential clocks is exponential, with rate equal to the sum of outgoing rates. The exponential distribution is memoryless, which gives the Markov property in continuous time.

**How do you compute a hitting time in a CTMC?**

Let `h_i = E_i[T_d]`. Set `h_d = 0`. For every other node, solve `sum_j Q[i,j] h_j = -1`. This linear system comes from conditioning on the first holding time and first jump.

**What is the meaning of `T_a^+`?**

It is the first positive return time to node `a`. In this codebase, it includes the initial holding time at `a`, then stops at the next entrance to `a` after the particle has left.

**Why does French-DeGroot converge or not converge?**

Convergence is controlled by the directed graph structure. If there is a unique sink strongly connected component reachable from all nodes, all opinions converge to one consensus. If there are multiple sink components, different sink components can retain different limiting values, so global consensus is not guaranteed.

**What determines the consensus value?**

The consensus value is a weighted average of initial opinions: `x* = sum_i pi_i x_i(0)`. The weights `pi` are the normalized left nullvector of the Laplacian. In graphs with a unique sink component, only nodes in that sink component may have positive weight.

**Why does the consensus variance have that formula?**

Because the consensus value is a linear combination of independent initial random variables. Variance of a weighted sum of independent variables is the sum of squared weights times variances.

**What changes after removing edges?**

Removing edges changes reachability, SCCs, and sink components. That can create multiple sinks, destroy global consensus, or preserve consensus through a different sink component. The direction of each removed edge matters.

**Why does Problem 2 match Problem 1 in particle perspective?**

Each particle is an independent copy of the single-particle CTMC. Averaging over 100 particles and many repetitions estimates the same expected return time as Problem 1, with Monte Carlo error.

**Why compare node counts with `N*pi`?**

At stationarity, each independent particle has probability `pi_i` of being at node `i`. Therefore the expected number of particles at node `i` is `N*pi_i`.

**What is the difference between proportional-rate and fixed-rate open networks?**

In proportional service, the total service rate at a node grows with the number of particles there. In fixed service, the service clock rate is bounded and independent of queue size. Fixed service can have bottlenecks; proportional service can absorb larger finite populations because capacity scales with population.

**What does blow-up mean?**

It means sustained growth of total population or queues rather than fluctuation around a finite scale. In simulations, blow-up is indicated by positive long-run drift and large final totals.

**How did you determine the largest stable input rate?**

For proportional service, theory indicates no finite upper bound for finite input rates in this absorbing network. For fixed service, compute routing loads per unit input and compare them with service capacities. The bottleneck node is `c`, giving threshold `lambda = 4/3`; stable rates are below that threshold.

**What are the main possible sources of numerical error?**

Monte Carlo sampling error, finite simulation horizon, randomness near stability thresholds, matrix eigendecomposition roundoff in opinion dynamics, and inconsistencies between code constants and saved outputs.

## 10. Common mistakes and how to avoid them

- Confusing `Lambda` with `Q`: remember `Lambda` has nonnegative rates, while `Q` has negative diagonal entries and row sums zero.
- Using the wrong edge direction: rows are sources and columns are destinations.
- Mixing `P_jump` and a uniformized transition matrix: `P_jump` is used only after real jumps.
- Using one simulation run instead of averages: stochastic claims need many repetitions or theoretical support.
- Forgetting Monte Carlo error: compare differences with standard errors and confidence intervals.
- Claiming consensus without graph proof: check SCCs and sink components.
- Deriving stability thresholds only from plots: use traffic loads and service capacities.
- Reporting numbers not saved in CSVs: every numerical claim should trace to `results/` or a documented formula.
- Writing a report that depends too much on code: state models, methods, equations, values, and interpretations directly in the report.
- Ignoring source/result drift: source constants, CSV metadata, report text, and `RESULTS_LOG.md` should be checked together after every rerun.

## 11. Final checklist

Mathematical correctness:

- [ ] Node order `[o, a, b, c, d]` is stated and used everywhere.
- [ ] `Lambda[i,j]` row-source convention is stated and used.
- [ ] `omega`, `Q`, and `P_jump` are correctly distinguished.
- [ ] Return-time convention is explicit.
- [ ] Hitting-time linear system is correct.
- [ ] French-DeGroot consensus claims are supported by SCC/sink reasoning.
- [ ] Open-network fixed-rate threshold is supported by load/capacity analysis.

Code correctness:

- [ ] Run scripts match documented seeds and run counts.
- [ ] Problem 2 `NODE_COUNT_RUNS` matches CSV metadata, report text, and `RESULTS_LOG.md`.
- [ ] Simulations use embedded jump probabilities only after actual departures.
- [ ] Node `d` removal rule is implemented correctly in Problem 3.
- [ ] Results are generated by scripts, not manual invention.

Result consistency:

- [ ] Every report number appears in a CSV or documented formula.
- [ ] `RESULTS_LOG.md` agrees with CSVs.
- [ ] CSV metadata agrees with source constants.
- [ ] Monte Carlo standard errors and confidence intervals are documented.

Figure consistency:

- [ ] Every referenced figure exists in `figures/`.
- [ ] Figure captions explain what the figure supports.
- [ ] Figures have readable axes and legends.

Report-content readiness:

- [ ] Every subquestion is answered.
- [ ] Methods are explained independently of code.
- [ ] Simulation/theory comparisons are interpreted.
- [ ] Conclusions are not overclaimed.
- [ ] Compilation issues are ignored for this audit.

Oral-defense readiness:

- [ ] You can explain `Lambda`, `Q`, `omega`, and `P_jump`.
- [ ] You can derive return and hitting time equations.
- [ ] You can explain consensus and sink components.
- [ ] You can justify the variance formula.
- [ ] You can explain particle versus node perspectives.
- [ ] You can justify proportional versus fixed stability behavior.
- [ ] You can identify unresolved issues honestly, especially the unverified lecture-specific French-DeGroot convention.
