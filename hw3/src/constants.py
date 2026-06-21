"""Reusable constants for Homework 3.

Part 1 is a discrete-time SIR epidemic on graphs (k-regular and preferential
attachment); Part 2 is a 3-player coordination / anti-coordination network game.
All vaccination schedules are cumulative percentages of the population, indexed by
week (entry t-1 is the cumulative fraction targeted by the start of week t).
"""

from pathlib import Path

# Disease-dynamics parameters for Problems 1.1, 2, and 3.
BETA = 0.3          # per-link infection probability per week
RHO = 0.7           # per-week recovery probability
WEEKS = 15          # number of simulated weeks (states recorded for weeks 0..15)
N_RUNS = 100        # Monte Carlo repetitions for Problems 1.1, 2, 3
N_INIT_INFECTED = 10  # initial infected seeds for Problems 1.1, 2, 3

# Graph sizes and degrees.
N_KREGULAR = 500    # nodes for the k-regular graph (Problem 1.1)
K_KREGULAR = 4      # degree of the k-regular graph (Problem 1.1)
N_PA = 500          # nodes for the preferential-attachment graph (Problems 2, 3)
K_PA = 6            # target average degree for Problems 2, 3

# Cumulative vaccination percentage by week 1..15 for Problem 3.
VACC_PROBLEM3 = [0, 5, 15, 25, 35, 45, 55, 60, 60, 60, 60, 60, 60, 60, 60]

# Problem 4 uses the H1N1 Sweden 2009 population scaled down by 1e4.
N_SWEDEN = 934      # scaled population
# Newly infected per week, weeks 42/2009 -> 5/2010. First entry = initial infected.
I0_SWEDEN = [1, 1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0]
# Cumulative vaccinated % by week (already 5% at the start, week 42).
VACC_SWEDEN = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60]

# Coordinate-search starting guess and step sizes for Problem 4.
K0_SEARCH = 10
BETA0_SEARCH = 0.3
RHO0_SEARCH = 0.6
DELTA_K = 1
DELTA_BETA = 0.1
DELTA_RHO = 0.1
# Monte Carlo repetitions per parameter evaluation (Problem 4). The assignment suggests
# N=10, but with a single initial seed the RMSE objective is noisy (std ~1.2 at N=10);
# since the search minimizes over 27 noisy evaluations, that selection bias can lock
# onto spurious minima (true RMSE far from the apparent one). We use N=100 for a stable,
# reproducible objective (std ~0.4), consistent with the assignment's invitation to tune.
N_RUNS_SEARCH = 100

# Reproducibility.
DEFAULT_SEED = 20260608

# Project paths.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"
