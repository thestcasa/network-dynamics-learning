# Homework 3 - Network Dynamics and Learning

## Purpose

This folder contains the code, saved numerical results, figures, and LaTeX
report for Homework 3. Part 1 simulates the H1N1 2009 Sweden pandemic with a
discrete-time SIR model on a `k`-regular graph and on preferential-attachment
(PA) random graphs, adds a vaccination campaign, and estimates the network and
disease parameters `(k, beta, rho)` that best fit the real scaled data. Part 2
studies a 3-player coordination / anti-coordination network game: its Nash
equilibria and the asynchronous best-response and noisy best-response (logit)
dynamics.

## Folder Structure

```text
hw3/
|-- README.md
|-- RESULTS_LOG.md
|-- data/
|   `-- hw-3-2025-26_en.pdf
|-- src/
|   |-- constants.py            parameters, schedules, seeds, paths
|   |-- utils.py                CSV helpers
|   |-- plotting.py             figure helpers (save_figure, plot_lines, plot_histogram)
|   |-- graphs.py               k_regular, preferential_attachment, degree_stats
|   |-- sir.py                  simulate_sir + run_many (Monte Carlo averaging)
|   |-- games.py                utilities, Nash, BR / logit transition matrices, limits
|   |-- run_problem1.py         1.1 k-regular epidemic + 1.2 PA generator validation
|   |-- run_problem2.py         PA epidemic, no vaccination
|   |-- run_problem3.py         PA epidemic, with vaccination
|   |-- run_problem4.py         parameter estimation (k, beta, rho)
|   `-- run_games.py            Part 2: Nash + best-response / logit dynamics
|-- figures/
|-- results/
`-- report/
    `-- report.tex
```

## Python Environment

Use Python 3 with:

- `numpy`
- `scipy` (sparse matrices for the epidemic)
- `matplotlib`
- `networkx` (only for drawing the Part 2 transition graphs)

The scripts also use Python standard-library modules such as `csv`, `pathlib`,
and `math`.

Example setup from the `hw3/` folder:

```bash
pip install numpy scipy matplotlib networkx
```

## Running the Scripts

Run all commands from the `hw3/` folder.

```bash
python src/run_problem1.py   # 1.1 k-regular epidemic + 1.2 PA graph validation
python src/run_problem2.py   # PA-graph epidemic, no vaccination
python src/run_problem3.py   # PA-graph epidemic, with vaccination
python src/run_problem4.py   # parameter estimation for Sweden 2009 (~20 s)
python src/run_games.py      # Part 2: Nash equilibria + dynamics
```

Each script prints a short summary, writes CSV files to `results/`, and saves
PNG figures to `figures/`.

## Outputs

Figures are saved in `figures/` and numerical results as CSV files in
`results/`. `RESULTS_LOG.md` records the assumptions, seeds, result files,
figure files, and key numbers. The numerical values in `report/report.tex` are
traceable to the CSV files in `results/`.

## Compiling the Report

The report source is `report/report.tex`. With a LaTeX distribution installed,
compile from the `hw3/report/` folder:

```bash
pdflatex report.tex
pdflatex report.tex
```

Before submission, fill in the author name and the collaboration statement in
`report/report.tex` (currently placeholders).

## Reproducibility Notes

- All run scripts set fixed random seeds (base seed `20260608`, with per-task
  offsets); results are reproducible.
- Reported numerical values match the CSV files in `results/`.
- Modeling conventions (documented in `RESULTS_LOG.md`): synchronous SIR update
  from the start-of-week state; weekly vectors indexed 0..15 with week 0 holding
  the seeds; vaccination applied at the start of each week among the not-yet
  vaccinated; Problem 4 starts from a single infected seed (`I0(0) = 1`).
- Problem 4 uses `N = 100` Monte Carlo runs per evaluation (instead of the
  suggested `N = 10`) to stabilize the noisy RMSE objective; see the report and
  `RESULTS_LOG.md`.
- Do not manually edit CSV outputs. If a model or simulation bug is found,
  document it in `RESULTS_LOG.md` and rerun the relevant script.
```
