# Codebase Context - Homework 2

## Purpose

This repository contains the complete working folder for Homework 2 of Network Dynamics and Learning. It includes Python implementations, fixed-seed numerical simulations, generated CSV results, generated figures, a LaTeX report, and support Markdown files used to track assumptions, fixes, and study notes.

The project studies:

- single-particle continuous-time Markov chains on a directed weighted network;
- continuous-time French-DeGroot opinion dynamics;
- many-particle continuous-time random walks on the same closed network;
- open particle networks under proportional and fixed service rules.

## Current Status

The homework workflow is essentially finalized as of 2026-06-04.

- Problems 1, 2, and 3 have Python implementations in `src/`.
- Result CSVs exist in `results/`.
- Plot PNGs exist in `figures/`.
- The report source is `report/report.tex`.
- `report/report.pdf` exists, but the final log says LaTeX tools were not available locally, so the PDF could not be regenerated during finalization.
- `README.md`, `HW2_CONTEXT.md`, `RESULTS_LOG.md`, `HW2_STUDY_GUIDE.md`, and `FIX_PROMPTS.md` document workflow, assumptions, fixes, and explanations.

Remaining known issues:

- The continuous-time French-DeGroot convention is used consistently, but the course lecture-note convention was not independently available in this workspace.
- Local LaTeX compilation was not possible because no LaTeX engine was installed.

## Repository Layout

```text
hw2/
|-- CODEBASE_CONTEXT.md
|-- FIX_PROMPTS.md
|-- HW2_CONTEXT.md
|-- HW2_STUDY_GUIDE.md
|-- README.md
|-- RESULTS_LOG.md
|-- data/
|   `-- hw2_en_2025-26.pdf
|-- figures/
|   |-- problem1_consensus_value_histogram.png
|   |-- problem1_fdg_original_trajectory.png
|   |-- problem1_fdg_removed_edges_g.png
|   |-- problem1_fdg_removed_edges_h_case1.png
|   |-- problem1_fdg_removed_edges_h_case2.png
|   |-- problem2_node_counts_over_time.png
|   |-- problem2_stationary_distribution_comparison.png
|   |-- problem3_fixed_lambda2_counts.png
|   |-- problem3_fixed_stability_scan.png
|   |-- problem3_proportional_lambda100_counts.png
|   `-- problem3_proportional_stability_scan.png
|-- report/
|   |-- report.pdf
|   `-- report.tex
|-- results/
|   |-- problem1_consensus_variance.csv
|   |-- problem1_fdg_original_consensus.csv
|   |-- problem1_hitting_o_to_d.csv
|   |-- problem1_removed_edges_g_summary.csv
|   |-- problem1_removed_edges_h_summary.csv
|   |-- problem1_return_time_a.csv
|   |-- problem2_node_final_counts.csv
|   |-- problem2_node_time_series.csv
|   |-- problem2_particle_return_times.csv
|   |-- problem3_fixed_stability_scan.csv
|   |-- problem3_fixed_timeseries_lambda2.csv
|   |-- problem3_proportional_stability_scan.csv
|   `-- problem3_proportional_timeseries_lambda100.csv
`-- src/
    |-- constants.py
    |-- plotting.py
    |-- problem1_opinion.py
    |-- problem1_single_ctmc.py
    |-- problem2_many_particles.py
    |-- problem3_open_network.py
    |-- run_problem1.py
    |-- run_problem2.py
    |-- run_problem3.py
    `-- utils.py
```

## Global Modeling Conventions

Node order is fixed everywhere:

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

Matrix orientation:

```text
Lambda[i, j] = rate or weight from source node i to destination node j
```

Rows are source nodes and columns are destination nodes. This convention is important for all CTMC, graph, edge-removal, and opinion-dynamics calculations.

For the closed network:

```text
omega = Lambda @ 1
Q = Lambda - diag(omega)
P_jump[i, j] = Lambda[i, j] / omega_i
```

`P_jump` is the embedded jump-chain transition matrix after a particle leaves a node. It is not a uniformized transition matrix.

## Source Files

### `src/constants.py`

Defines reusable project constants:

- `NODES` and `NODE_TO_INDEX`;
- closed-network matrix `LAMBDA`;
- open-network matrix `LAMBDA_OPEN`;
- `DEFAULT_SEED = 20260604`;
- project paths `PROJECT_ROOT`, `FIGURES_DIR`, and `RESULTS_DIR`.

The closed-network `LAMBDA` is used for Problems 1 and 2. The open-network `LAMBDA_OPEN` is used for Problem 3.

### `src/problem1_single_ctmc.py`

Implements single-particle CTMC routines:

- exit rates;
- CTMC generator;
- embedded jump matrix;
- stationary distribution;
- Monte Carlo return-time simulation;
- theoretical return time using `1 / (pi_i omega_i)`;
- Monte Carlo hitting-time simulation;
- theoretical hitting times by solving the linear system with target boundary condition.

Important convention: return time to `a` includes the initial holding time in node `a`.

### `src/problem1_opinion.py`

Implements continuous-time French-DeGroot tools:

- Laplacian `L = diag(Lambda @ 1) - Lambda`;
- left consensus vector;
- matrix exponential action by eigendecomposition;
- simulation of `dx/dt = -Lx`;
- directed edge removal;
- adjacency lists;
- strongly connected components by Tarjan's algorithm;
- sink-component detection;
- graph summary;
- large-time asymptotic-state approximation.

This module is used for original graph consensus and both edge-removal cases.

### `src/problem2_many_particles.py`

Implements many-particle closed-network simulations:

- computes `omega`, `Q`, and `P_jump`;
- solves stationary distribution;
- simulates one particle return time;
- simulates repeated averages across `N = 100` particles;
- simulates the node-count CTMC directly, where node `i` has departure rate `n_i(t) omega_i`;
- returns mean, standard deviation, and final counts across Monte Carlo runs.

Problem 2 uses `N = 100`, horizon `T = 60`, 2000 particle-perspective Monte Carlo batches, and 5000 node-perspective Monte Carlo runs.

### `src/problem3_open_network.py`

Implements the open-network event simulator:

- open-network exit rates and embedded destination probabilities;
- assignment override `omega_d = 7/4`;
- service modes `proportional` and `fixed`;
- event simulation with external arrivals into node `o`;
- replicated simulations;
- total-population slope estimation for stability scans;
- theoretical fixed-rate routing loads;
- fixed-rate bottleneck capacity;
- mean absorption times for proportional-rate stability reasoning.

In proportional mode, service rate at node `i` is `omega_i N_i(t)`.

In fixed mode, node `i` has clock rate `omega_i`; empty service ticks do nothing.

### `src/run_problem1.py`

Runs all Problem 1 computations and writes CSV/figure outputs.

Key parameters:

- `RETURN_RUNS = 100_000`;
- `HITTING_RUNS = 100_000`;
- `CONSENSUS_MC_RUNS = 50_000`;
- `FDG_TIME_MAX = 60.0`;
- `FDG_TIME_POINTS = 600`.

Outputs include return-time and hitting-time CSVs, French-DeGroot trajectory figures, consensus variance results, and edge-removal summaries.

### `src/run_problem2.py`

Runs all Problem 2 simulations and writes outputs.

Key parameters:

- `NUM_PARTICLES = 100`;
- `PARTICLE_RETURN_RUNS = 2000`;
- `NODE_COUNT_RUNS = 5000`;
- `HORIZON = 60.0`;
- `TIME_STEP = 0.1`;
- `PARTICLE_SEED = 20260624`;
- `NODE_SEED = 20260625`.

It also reads the saved Problem 1 return-time result to compare the many-particle particle-perspective estimate against Problem 1.

### `src/run_problem3.py`

Runs Problem 3 open-network simulations and stability scans.

Main proportional case:

- `lambda = 100`;
- horizon `T = 60`;
- output `problem3_proportional_timeseries_lambda100.csv`.

Proportional scan:

- lambdas `[10, 50, 100, 200, 500]`;
- horizon `T = 120`;
- 2 replicates per lambda.

Main fixed case:

- `lambda = 2`;
- horizon `T = 6000`;
- output `problem3_fixed_timeseries_lambda2.csv`.

Fixed scan:

- lambdas `[0.5, 1.0, 1.2, 1.3, 4/3, 1.35, 1.5, 2.0]`;
- horizon `T = 6000`;
- 3 replicates per lambda.

The fixed-rate theoretical capacity is `4/3`, with node `c` as bottleneck.

### `src/plotting.py`

Contains `save_figure(fig, path)`, which creates parent folders, applies `tight_layout()`, saves at `dpi=200`, and closes the figure.

### `src/utils.py`

Contains shared helpers:

- create parent folders;
- write one-row CSVs;
- write multi-row CSVs;
- compute mean, sample standard deviation, standard error, and normal 95 percent confidence interval.

## Results Summary

### Problem 1

Problem 1 covers a single-particle CTMC and French-DeGroot dynamics.

Generated result files:

- `results/problem1_return_time_a.csv`;
- `results/problem1_hitting_o_to_d.csv`;
- `results/problem1_fdg_original_consensus.csv`;
- `results/problem1_consensus_variance.csv`;
- `results/problem1_removed_edges_g_summary.csv`;
- `results/problem1_removed_edges_h_summary.csv`.

Generated figures:

- `figures/problem1_fdg_original_trajectory.png`;
- `figures/problem1_consensus_value_histogram.png`;
- `figures/problem1_fdg_removed_edges_g.png`;
- `figures/problem1_fdg_removed_edges_h_case1.png`;
- `figures/problem1_fdg_removed_edges_h_case2.png`.

Key logged values:

- return time to `a`: simulation mean `6.69375284456`, theory `6.70833333333`;
- hitting time from `o` to `d`: simulation mean `10.7608365852`, theory `10.7666666667`;
- original graph reaches global consensus for every initial condition;
- consensus vector is approximately `(0.217391, 0.149068, 0.260870, 0.186335, 0.186335)`;
- consensus-value theoretical variance is `0.331970217198`;
- edge-removal case P1(g) has two sink components and does not converge to global consensus for every initial condition;
- edge-removal case P1(h) has one reachable sink component and does converge to global consensus.

### Problem 2

Problem 2 studies `N = 100` particles on the closed network.

Generated result files:

- `results/problem2_particle_return_times.csv`;
- `results/problem2_node_final_counts.csv`;
- `results/problem2_node_time_series.csv`.

Generated figures:

- `figures/problem2_node_counts_over_time.png`;
- `figures/problem2_stationary_distribution_comparison.png`.

Key logged values:

- particle-perspective average return time: `6.70478013087`;
- Problem 1 theoretical return time comparison: `6.70833333333`;
- stationary distribution: `o=0.217391304348`, `a=0.149068322981`, `b=0.260869565217`, `c=0.186335403727`, `d=0.186335403727`;
- final simulated counts at `T = 60` are close to `100*pi`.

Final average counts versus stationary expected counts:

```text
o: 21.7254 versus 21.7391304348
a: 14.9592 versus 14.9068322981
b: 25.991  versus 26.0869565217
c: 18.6466 versus 18.6335403727
d: 18.6778 versus 18.6335403727
```

### Problem 3

Problem 3 studies the open network with external arrivals into `o`.

Generated result files:

- `results/problem3_proportional_timeseries_lambda100.csv`;
- `results/problem3_proportional_stability_scan.csv`;
- `results/problem3_fixed_timeseries_lambda2.csv`;
- `results/problem3_fixed_stability_scan.csv`.

Generated figures:

- `figures/problem3_proportional_lambda100_counts.png`;
- `figures/problem3_proportional_stability_scan.png`;
- `figures/problem3_fixed_lambda2_counts.png`;
- `figures/problem3_fixed_stability_scan.png`.

Open-network service rates:

```text
omega_o = 2
omega_a = 1
omega_b = 1
omega_c = 1
omega_d = 7/4
```

Proportional service:

- main case uses `lambda = 100`, `T = 60`;
- stability scan treats proportional service as stable for every finite input rate;
- final logged proportional counts at `lambda = 100`, `T = 60`: `(48, 39, 72, 62, 63)`, total `284`.

Fixed service:

- main case uses `lambda = 2`, `T = 6000`;
- theoretical loads per unit lambda are `[1, 1/2, 5/8, 3/4, 1]`;
- theoretical largest stable input rate is `4/3`;
- bottleneck node is `c`;
- final logged fixed counts at `lambda = 2`, `T = 6000`: `(233, 15, 1240, 1577, 6)`, total `3071`.

## Documentation Files

### `README.md`

Operational entry point. It explains project purpose, folder structure, Python dependencies, commands for each problem, output locations, report compilation instructions, and reproducibility notes.

### `HW2_CONTEXT.md`

Original assignment/workflow context. It documents submission requirements, expected folder structure, modeling conventions, problem statements, output expectations, and a final checklist.

Some box-drawing characters appear mojibake-encoded in this file, but the content is still understandable.

### `RESULTS_LOG.md`

Authoritative work log for simulations and report values. It records:

- methods and assumptions;
- seeds;
- result file names;
- figure file names;
- numerical summaries;
- report-writing status;
- finalization notes;
- content-fix pass notes;
- remaining unresolved issues.

Use this file when checking whether report values are traceable to CSV outputs.

### `HW2_STUDY_GUIDE.md`

Study and explanation guide for the homework. It appears to consolidate conceptual explanations and corrected notes for understanding the assignment and code.

### `FIX_PROMPTS.md`

Prompt/history-style file for fixes and audits. It likely records requested correction passes and is useful for understanding why certain consistency checks were made.

## How To Run

From the repository root:

```powershell
python src/run_problem1.py
python src/run_problem2.py
python src/run_problem3.py
```

Python dependencies:

```text
numpy
matplotlib
```

The scripts also use standard-library modules such as `csv`, `json`, `pathlib`, `platform`, `sys`, `math`, and `dataclasses`.

If using a virtual environment:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install numpy matplotlib
```

## Report

Report source:

```text
report/report.tex
```

Expected PDF:

```text
report/report.pdf
```

Compile from `report/` if LaTeX is installed:

```powershell
pdflatex report.tex
pdflatex report.tex
```

The logs say the local environment did not have `pdflatex`, `latexmk`, `xelatex`, or `lualatex`, so the final PDF was not regenerated during finalization.

## Reproducibility And Maintenance Notes

- Do not manually edit CSV outputs.
- If a simulation/model bug is found, update the source, rerun the relevant script, and update `RESULTS_LOG.md` and `report/report.tex` together.
- Keep the row-source, column-destination matrix convention explicit.
- Keep `P_jump` separate from any possible uniformized transition matrix.
- Keep Problem 2 `NODE_COUNT_RUNS = 5000`, because the saved CSVs and report values use 5000 node-perspective Monte Carlo runs.
- If report values change, trace them to the exact CSV files in `results/`.
- If figures are regenerated, check that all `\includegraphics` paths in `report/report.tex` still resolve from the `report/` directory.

