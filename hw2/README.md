# Homework 2 - Network Dynamics and Learning

## Purpose

This folder contains the code, saved numerical results, figures, and LaTeX
report for Homework 2. The homework studies continuous-time random walks,
French-DeGroot opinion dynamics, many-particle continuous-time random walks,
and open particle networks with proportional and fixed service rates.

## Folder Structure

```text
hw2/
|-- HW2_CONTEXT.md
|-- RESULTS_LOG.md
|-- README.md
|-- data/
|-- src/
|   |-- constants.py
|   |-- problem1_single_ctmc.py
|   |-- problem1_opinion.py
|   |-- problem2_many_particles.py
|   |-- problem3_open_network.py
|   |-- run_problem1.py
|   |-- run_problem2.py
|   |-- run_problem3.py
|   |-- plotting.py
|   `-- utils.py
|-- figures/
|-- results/
`-- report/
    |-- report.tex
    `-- report.pdf
```

## Python Environment

Use Python 3 with:

- `numpy`
- `matplotlib`

The scripts also use Python standard-library modules such as `csv`,
`pathlib`, and `math`.

Example setup from the `hw2/` folder:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install numpy matplotlib
```

## Running the Scripts

Run all commands from the `hw2/` folder.

Problem 1:

```powershell
python src/run_problem1.py
```

This generates the single-particle return-time and hitting-time results, the
French-DeGroot consensus results, and the Problem 1 figures.

Problem 2:

```powershell
python src/run_problem2.py
```

This generates the many-particle return-time comparison, node-count time
series, final node-count comparison, and Problem 2 figures.

Problem 3:

```powershell
python src/run_problem3.py
```

This generates the open-network proportional-service and fixed-service time
series, stability scans, and Problem 3 figures.

## Outputs

Figures are saved in:

```text
figures/
```

Numerical results are saved as CSV files in:

```text
results/
```

`RESULTS_LOG.md` records the simulation assumptions, seeds, result files,
figure files, and final report status.

The numerical values in `report/report.tex` are intended to be traceable to
the CSV files in `results/` or to theoretical formulas documented in
`RESULTS_LOG.md`. If a script is rerun and changes a CSV or figure, update the
corresponding report text and `RESULTS_LOG.md` in the same pass.

## Compiling the Report

The report source is:

```text
report/report.tex
```

If LaTeX is installed, compile from the `hw2/report/` folder:

```powershell
pdflatex report.tex
pdflatex report.tex
```

The expected PDF output is:

```text
report/report.pdf
```

## Reproducibility Notes

- The run scripts set fixed random seeds, documented in `RESULTS_LOG.md`.
- Reported numerical values should match the CSV files in `results/`.
- Figures are generated from the same scripts and saved in `figures/`.
- Problem 2 uses `N = 100` particles, `T = 60`, 2000 particle-perspective
  Monte Carlo batches, and 5000 node-perspective Monte Carlo runs.
- Problem 3 saves both main trajectories and stability-scan summaries; the
  report interprets the scan together with the analytical stability criteria.
- Do not manually edit CSV outputs. If a model or simulation bug is found,
  document the issue in `RESULTS_LOG.md` and rerun the relevant script.
