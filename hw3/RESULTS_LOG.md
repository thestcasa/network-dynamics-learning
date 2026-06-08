# Results Log - Homework 3

Base random seed: `20260608` (constant `DEFAULT_SEED`), with per-task offsets noted
below. All numbers come from the CSV files in `results/`; figures are in `figures/`.
Report status: `report/report.tex` written; author name and collaboration statement are
placeholders to fill before submission.

## Modeling conventions (apply to all of Part 1)

- **SIR update**: synchronous, computed from the start-of-week state `X(t)`. New
  infections (from the current infected set) and recoveries (of currently infected nodes)
  are evaluated together; a node infected this week cannot recover the same week.
- **Weekly vectors** have length 16, indexed by week `0..15`. Week 0 is the initial
  configuration; `newly_infected[0]` = number of seeds. Totals `S/I/R(/V)` are recorded at
  every week. Invariant `S+I+R+V = n` is asserted in code each week.
- **Vaccination** (Problems 3, 4): applied at the start of each week before infection
  spreads; cumulative target `round(Vacc[week-1]/100 * n)`; new vaccinees chosen uniformly
  among the not-yet-vaccinated (any S/I/R state) and moved to an absorbing `V`.
- **Infected-neighbour counts** via a single `scipy.sparse` matrix-vector product per week.

## Problem 1.1 - Epidemic on a known k-regular graph

- Graph: symmetric k-regular, `n=500`, `k=4` (realized degree exactly 4 for all nodes).
- Params: `beta=0.3`, `rho=0.7`, 15 weeks, 10 random seeds, `N=100`. Epidemic seed base
  `DEFAULT_SEED + 1000`.
- Result: newly infected peaks in week 1 (~11.3) then decays monotonically; final totals
  `S ~ 439.8`, `R ~ 60.0` (small, locally-saturating outbreak).
- Files: `results/problem1_kregular_newly_infected.csv`,
  `results/problem1_kregular_sir_totals.csv`;
  figures `figures/problem1_kregular_newly_infected.png`,
  `figures/problem1_kregular_sir_totals.png`.

## Problem 1.2 - Preferential-attachment graph generator

- Generator: start from complete `K_{k+1}`; each new node attaches to `c=k/2` existing
  nodes without replacement, probability proportional to current degree; odd `k` alternates
  floor/ceil. Parameterized by `(n,k)` only.
- Validation graph: `n=500`, `k=6`, rng seed `DEFAULT_SEED + 1100`. Realized mean degree
  exactly `6.0`, degrees range `3..74` (heavy tail with hubs).
- Files: `results/problem1_pa_degree_summary.csv`;
  figure `figures/problem1_pa_degree_distribution.png`.

## Problem 2 - Epidemic on a PA graph, no vaccination

- Graph: PA, `n=500`, `k=6`, rng seed `DEFAULT_SEED + 2100`.
- Params: `beta=0.3`, `rho=0.7`, 15 weeks, 10 seeds, `N=100`. Epidemic seed base
  `DEFAULT_SEED + 2000`.
- Result: newly infected peaks in week 3 (~106); final `R ~ 408.8` (~82% of population) -
  much larger and faster than the k-regular case (hubs as super-spreaders).
- Files: `results/problem2_newly_infected.csv`, `results/problem2_sir_totals.csv`;
  figures `figures/problem2_newly_infected.png`, `figures/problem2_sir_totals.png`.

## Problem 3 - Epidemic on a PA graph, with vaccination

- Graph: PA, `n=500`, `k=6`, rng seed `DEFAULT_SEED + 3100`.
- Vacc schedule: `[0,5,15,25,35,45,55,60,...,60]`. Epidemic seed base `DEFAULT_SEED + 3000`.
- Result: newly infected peaks in week 3 (~75, vs 106 without vaccination); final
  `R ~ 127.6` (~26%, vs 82%); `V` rises to exactly `300` (60% of 500) by week 8.
- Files: `results/problem3_newly_infected_vaccinated.csv`,
  `results/problem3_sirv_totals.csv`;
  figures `figures/problem3_newly_infected_vaccinated.png`,
  `figures/problem3_sirv_totals.png`.

## Problem 4 - Parameter estimation (Sweden 2009)

- Population `n=934` (scaled by 1e4); initial infected `I0(0)=1`; real
  `I0(t)=[1,1,3,5,9,17,32,32,17,5,2,1,0,0,0,0]`;
  `Vacc(t)=[5,9,16,24,32,40,47,54,59,60,...,60]`.
- Method: coordinate search over the 3x3x3 neighbourhood of `(k,beta,rho)` from
  `(10,0.3,0.6)`, steps `(1,0.1,0.1)`; halve `(d_beta, d_rho)` twice after convergence.
  PA graph cached per `k` (seed `DEFAULT_SEED + 4000 + k`); RMSE over weeks 1..15.
  Search epidemic seed base `DEFAULT_SEED + 4100`; final presentation run seed base
  `DEFAULT_SEED + 4200`.
- **Monte Carlo size**: `N=100` per evaluation (not the suggested `N=10`). With a single
  initial seed the RMSE objective is noisy (std ~1.2 at N=10), and minimizing over 27 noisy
  points can lock onto spurious minima (observed apparent RMSE 5.6 whose true value was
  13.6). At N=100 the objective std drops to ~0.4 and the result is reproducible.
- **Best parameters**: `k=8`, `beta=0.2`, `rho=0.625`; search RMSE `5.57`, final RMSE
  (N=100, different seed) `6.73`. Fitted curve peaks ~30 (real 32), leads the data by ~1-2
  weeks, matches the overall bell shape and total attack rate (~13%).
- Files: `results/problem4_best_params.csv`, `results/problem4_search_log.csv`,
  `results/problem4_fit_vs_real.csv`, `results/problem4_sirv_totals.csv`;
  figures `figures/problem4_fit_vs_real.png`, `figures/problem4_sirv_totals.png`.

## Part 2 - Network games and dynamics

- Setup: `n=3`, actions `{-1,+1}`, complete-graph interaction; players `1..n1` coordinate,
  the rest anti-coordinate. Computation is exact over all 8 profiles (no randomness).
- **Nash equilibria** (pure): `n1=3` -> {`---`,`+++`}; `n1=2` -> {`--+`,`++-`};
  `n1=1` -> {`--+`,`-+-`,`+-+`,`++-`}; `n1=0` -> the six non-consensus `2-1` splits.
  Files: `results/games_nash_n1_{0,1,2,3}.csv`, `results/games_nash_summary.csv`.
- **Dynamics** (discrete-time, asynchronous): best-response (strict-improvement tie-break,
  absorbing states = Nash) and logit `P(x_i=+1) ~ exp(u_i(+1)/eps)`.
  Transition-graph figures `figures/games_transition_n1_3.png`,
  `figures/games_transition_n1_0.png`.
- **Limit distributions** from `X(0)=(+1,-1,+1)`:
  - Noiseless BR (depends on X0): `n1=3` -> `+++` w.p. 1; `n1=0` -> `+-+` w.p. 1 (already NE).
  - Noisy BR, vanishing noise (independent of X0, exact Gibbs/potential limit):
    `n1=3` -> 1/2 on `+++`, 1/2 on `---`; `n1=0` -> 1/6 on each of the six splits.
    Confirmed numerically via the logit stationary distribution as eps -> 0
    (`figures/games_logit_convergence_n1_{0,3}.png`).
  File: `results/games_limit_distributions.csv`.

## Notes

- The logit stationary distribution is computed via a least-squares solve of
  `pi(P-I)=0` with `sum(pi)=1`, which stays symmetry-preserving at small `eps` where an
  eigenvector solve becomes unreliable (near-reducible chain).
- LaTeX was not available in the build environment; compile `report/report.tex` locally
  with `pdflatex` (run twice for references).
```
