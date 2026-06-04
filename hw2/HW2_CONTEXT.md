# HW2_CONTEXT.md — Network Dynamics and Learning Homework 2

## 0. Scope

This folder contains the full solution workflow for Homework 2 of the Network Dynamics and Learning course.

Work only inside:

```text
network-dynamics-learning/hw2/
```

Do not modify:

```text
hw1/
hw3/
shared/
```

Workflow:

1. Create and maintain `HW2_CONTEXT.md`.
2. Create and maintain `RESULTS_LOG.md`.
3. Solve one problem at a time.
4. Generate figures and CSV/result files.
5. Update `RESULTS_LOG.md` after each problem.
6. Write the LaTeX report section by section.
7. Perform a strict final audit before submission.

The final submission must contain:

- PDF lab report;
- source code;
- generated figures;
- generated numerical result files;
- clear explanation of theoretical computations;
- readable plots with axis labels;
- enough detail for the report to be understandable independently of the code.

---

## 1. Official submission requirements

- Upload the solution on the course website under the name `Homework2`.
- Submit a PDF lab report together with the code.
- Numerical codes or hand computations are allowed.
- Python and MATLAB are supported by the TA.
- Code should be general enough that small changes in the question require only small code modifications.
- Code must be clearly commented.
- The PDF must read like a standard lab report.
- The report must contain descriptions of what is being done and proper presentation of results.
- Figures must be readable and have axis labels.
- The report must contain all important information and be readable independently of the code.
- LaTeX is strongly encouraged.
- Clarity and synthesis are part of the evaluation.
- If collaboration is involved, specify collaborators and the parts discussed/worked on.

---

## 2. Required folder structure

Expected structure:

```text
hw2/
├── HW2_CONTEXT.md
├── RESULTS_LOG.md
├── data/
├── src/
├── figures/
├── results/
├── report/
│   ├── report.tex
│   └── report.pdf
└── README.md
```

Suggested source files:

```text
hw2/src/
├── constants.py
├── utils.py
├── plotting.py
├── problem1_single_ctmc.py
├── problem1_opinion.py
├── problem2_many_particles.py
├── problem3_open_network.py
├── run_problem1.py
├── run_problem2.py
└── run_problem3.py
```

---

## 3. Global implementation conventions

### Node order

Always use:

```python
NODES = ["o", "a", "b", "c", "d"]
```

Index mapping:

```text
o -> 0
a -> 1
b -> 2
c -> 3
d -> 4
```

### Matrix orientation

Use:

```text
Lambda[i, j] = rate/weight from node i to node j
```

Rows are source nodes. Columns are destination nodes.

This must be stated clearly in the report.

---

## 4. Closed network matrix

Node order:

```text
[o, a, b, c, d]
```

Matrix:

```text
Lambda =
[
 [0,   2/5, 1/5, 0,   0  ],
 [0,   0,   3/4, 1/4, 0  ],
 [1/2, 0,   0,   1/3, 0  ],
 [0,   0,   1/3, 0,   2/3],
 [0,   1/3, 0,   1/3, 0  ]
]
```

Continuous-time random walk convention:

```text
omega = Lambda @ 1
Q = Lambda - diag(omega)
```

where:

- `omega_i` is the total rate of leaving node `i`;
- `Q` is the CTMC generator.

Simulation rule:

1. if particle is in node `i`, draw holding time:
   ```text
   Delta t ~ Exp(omega_i)
   ```
2. choose next node `j` with probability:
   ```text
   Lambda[i,j] / omega_i
   ```

---

# Problem 1 — Single particle and French-DeGroot dynamics

## P1(a)

Simulate a single particle starting from node `a`.

Estimate the average time to leave node `a` and then return to node `a`.

Important:

- explicitly state whether the initial holding time in `a` is included;
- use the same convention in simulation and theory;
- average over many Monte Carlo repetitions.

Output:

```text
results/problem1_return_time_a.csv
```

---

## P1(b)

Compute the theoretical return time:

```text
E_a[T_a^+]
```

Possible methods:

- stationary distribution formula, if consistent with lecture definition;
- first-step equations / linear system.

The report must include:

- definition of `Q`;
- definition of return time;
- theoretical computation;
- comparison with simulation.

---

## P1(c)

Simulate the hitting time from node `o` to node `d`.

Output:

```text
results/problem1_hitting_o_to_d.csv
```

---

## P1(d)

Compute theoretical hitting time:

```text
E_o[T_d]
```

Let:

```text
h_i = E_i[T_d]
```

Boundary condition:

```text
h_d = 0
```

For every `i != d`, solve:

```text
sum_j Q[i,j] h_j = -1
```

Equivalent first-step form:

```text
h_i = 1/omega_i + sum_j P_jump[i,j] h_j
```

where:

```text
P_jump[i,j] = Lambda[i,j] / omega_i
```

---

## P1(e)

Interpret `Lambda` as weighted adjacency matrix of graph:

```text
G = (V, E, Lambda)
```

Simulate French-DeGroot dynamics.

If lecture notes do not impose another convention, use continuous-time:

```text
dx/dt = -Lx
L = diag(Lambda @ 1) - Lambda
```

Tasks:

- simulate from arbitrary initial condition;
- decide whether the dynamics converges to consensus for every initial condition;
- motivate using graph structure / eigenvectors.

Suggested figure:

```text
figures/problem1_fdg_original_trajectory.png
```

---

## P1(f)

Initial states:

```text
x_i(0) = xi_i
```

Independent random variables with variances:

```text
Var(x_a(0)) = Var(x_b(0)) = Var(x_c(0)) = 2
Var(x_o(0)) = Var(x_d(0)) = 1
```

If the consensus value is:

```text
x* = sum_i pi_i x_i(0)
```

then:

```text
Var(x*) = sum_i pi_i^2 Var(x_i(0))
```

where `pi` is the normalized left eigenvector associated with the consensus value.

Output:

```text
results/problem1_consensus_variance.csv
figures/problem1_consensus_value_histogram.png
```

---

## P1(g)

Remove edges:

```text
(d,a), (d,c), (a,c), (b,c)
```

Tasks:

- analyze asymptotic behavior;
- decide whether the dynamics converges;
- explain relation between asymptotic state and initial condition.

Use row-source, column-destination convention.

Suggested outputs:

```text
results/problem1_removed_edges_g_summary.csv
figures/problem1_fdg_removed_edges_g.png
```

---

## P1(h)

Remove edges:

```text
(b,o), (d,a)
```

Tasks:

- analyze French-DeGroot dynamics;
- describe how asymptotic behavior varies with initial condition;
- motivate using graph structure.

Suggested outputs:

```text
results/problem1_removed_edges_h_summary.csv
figures/problem1_fdg_removed_edges_h_case1.png
figures/problem1_fdg_removed_edges_h_case2.png
```

---

# Problem 2 — Many particles on closed network

Same closed network and same `Lambda`.

Definitions:

```text
omega = Lambda @ 1
P_jump = diag(omega)^(-1) Lambda
```

Important distinction:

```text
P_jump = embedded jump-chain transition matrix
```

If using uniformization, define:

```text
omega_star >= max_i omega_i
P_uniformized = I - diag(omega / omega_star) + Lambda / omega_star
```

Do not mix `P_jump` and `P_uniformized`.

---

## P2(a) Particle perspective

Tasks:

- simulate `N = 100` particles;
- all start in node `a`;
- estimate average return time to node `a`;
- compare with Problem 1.

Output:

```text
results/problem2_particle_return_times.csv
```

Explanation expected:

Each particle follows the same single-particle CTMC, so the average return time should agree with Problem 1 up to Monte Carlo error.

---

## P2(b) Node perspective

Tasks:

- `N = 100` particles start in node `a`;
- simulate for `T = 60`;
- compute average number of particles in each node at final time;
- plot number of particles in each node over time;
- compare with stationary distribution of single-particle CTMC.

At time `t`, if node `i` has `n_i(t)` particles, node `i` forwards particles at rate:

```text
n_i(t) omega_i
```

Destination sampled with:

```text
P_jump[i,j] = Lambda[i,j] / omega_i
```

Outputs:

```text
results/problem2_node_final_counts.csv
results/problem2_node_time_series.csv
figures/problem2_node_counts_over_time.png
figures/problem2_stationary_distribution_comparison.png
```

Stationary comparison:

```text
average final counts should be compared with N * pi
```

where `pi` is the stationary distribution of the CTMC.

---

# Problem 3 — Open network

## Open-network matrix

Node order:

```text
[o, a, b, c, d]
```

Matrix:

```text
Lambda_open =
[
 [0, 1,   1,   0, 0],
 [0, 0,   1/4, 1/4, 2/4],
 [0, 0,   0,   1, 0],
 [0, 0,   0,   0, 1],
 [0, 0,   0,   0, 0]
]
```

Particles enter node `o` according to a Poisson process with rate `lambda`.

Let:

```text
omega = Lambda_open @ 1
```

But the assignment states:

```text
omega_d = 7/4
```

When node `d` clock ticks:

- if `N_d(t) > 0`, remove one particle from `d`;
- otherwise do nothing.

---

## P3(a) Proportional rate

Rule:

```text
node i forwards particles at rate omega_i N_i(t)
```

Tasks:

- simulate for `T = 60`;
- use `lambda = 100`;
- plot evolution of particles in each node;
- determine largest input rate before blow-up.

Outputs:

```text
results/problem3_proportional_timeseries_lambda100.csv
results/problem3_proportional_stability_scan.csv
figures/problem3_proportional_lambda100_counts.png
figures/problem3_proportional_stability_scan.png
```

---

## P3(b) Fixed rate

Rule:

```text
node i forwards particles at fixed rate omega_i
```

Tasks:

- simulate for `T = 6000`;
- use `lambda = 2`;
- plot evolution of particles in each node;
- determine largest input rate before blow-up;
- motivate analytically.

Outputs:

```text
results/problem3_fixed_timeseries_lambda2.csv
results/problem3_fixed_stability_scan.csv
figures/problem3_fixed_lambda2_counts.png
figures/problem3_fixed_stability_scan.png
```

---

# Report structure

Suggested LaTeX structure:

```text
report/report.tex
```

Sections:

1. Introduction
2. Problem 1 — Single-particle CT random walk and opinion dynamics
3. Problem 2 — Many-particle CT random walks
4. Problem 3 — Open network
5. Conclusion
6. Collaboration statement
7. Appendix / implementation details if needed

Each section must contain:

- model definition;
- method;
- result;
- interpretation;
- comparison between simulation and theory where required.

---

# Final checklist

General:

- [ ] Work is entirely inside `hw2/`.
- [ ] `hw1/`, `hw3/`, and `shared/` were not modified.
- [ ] `HW2_CONTEXT.md` exists.
- [ ] `RESULTS_LOG.md` exists and is complete.
- [ ] `README.md` exists.
- [ ] `report/report.tex` exists.
- [ ] `report/report.pdf` compiles.
- [ ] Code is commented.
- [ ] Stochastic simulations use fixed seeds.
- [ ] Numerical values in the report come from saved result files.
- [ ] No invented numerical result appears in the report.
- [ ] Figures have readable labels, legends, and captions.
- [ ] Report is readable independently of the code.

Ambiguities to resolve:

- [ ] Confirm French-DeGroot convention from lecture notes.
- [ ] Confirm return-time convention: include or exclude initial holding time.
- [ ] Confirm edge direction convention after edge removals.
- [ ] Keep embedded jump chain and uniformized chain separate.