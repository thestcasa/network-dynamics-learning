# Homework 1 -- Network Dynamics and Learning

This project contains Python solutions, generated numerical results, figures,
and a LaTeX report for Homework 1. The exercises cover max-flow and cuts, Katz
centrality and PageRank, and traffic assignment with social optimum, Wardrop
equilibrium, and tolling computations.

## Dependencies

The scripts use:

- Python 3
- numpy
- scipy
- networkx
- cvxpy
- matplotlib
- pandas

Install the Python dependencies from the `hw1/` root folder with:

```bash
pip install numpy scipy networkx cvxpy matplotlib pandas
```

To compile the report, a LaTeX distribution with `pdflatex` is also required.

## How to Run

Run all commands from the project root folder, `hw1/`.

Exercise 1:

```bash
python src/exercise1.py
```

Exercise 2:

```bash
python src/exercise2.py
```

Exercise 3:

```bash
python src/exercise3.py
```

Compile the report:

```bash
pdflatex -interaction=nonstopmode -output-directory=report report/report.tex
```

## Regenerating Results and Figures

Running the three exercise scripts regenerates the numerical CSV outputs in
`results/` and the plots in `figures/`.

Exercise 1 outputs:

- `results/exercise1_cuts.csv`
- `results/exercise1_b_results.csv`
- `results/exercise1_c_results.csv`
- `figures/exercise1_b_throughput.png`
- `figures/exercise1_c_added_link_throughput.png`

Exercise 2 outputs:

- `results/exercise2_katz.csv`
- `results/exercise2_pagerank_beta015.csv`
- `results/exercise2_pagerank_beta_sensitivity.csv`
- `figures/exercise2_katz.png`
- `figures/exercise2_pagerank_beta015.png`
- `figures/exercise2_pagerank_difference_beta.png`

Exercise 3 outputs:

- `results/exercise3_shortest_path.csv`
- `results/exercise3_maxflow.csv`
- `results/exercise3_nu.csv`
- `results/exercise3_social_optimum.csv`
- `results/exercise3_wardrop.csv`
- `results/exercise3_tolls_total_travel_time.csv`
- `results/exercise3_additional_delay.csv`
- `figures/exercise3_flow_comparison.png`

The final report source is `report/report.tex`. When LaTeX compilation
succeeds, the PDF is saved as `report/report.pdf`.

## Data

Exercise 3 reads the MATLAB data files from `data/`:

- `data/flow.mat`
- `data/capacities.mat`
- `data/traffic.mat`
- `data/traveltime.mat`
