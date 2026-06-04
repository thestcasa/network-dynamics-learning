# RESULTS_LOG.md - Homework 2

This file tracks all numerical results, generated figures, assumptions, and report-writing status for HW2.

Do not insert invented values. Every numerical result must come from a script output or a checked theoretical computation.

---

## Initialization - 2026-06-04

- Folder structure created for `data/`, `src/`, `figures/`, `results/`, and `report/`.
- Placeholder source files created in `src/`.
- Minimal LaTeX report skeleton created at `report/report.tex`.
- Placeholder `report/report.pdf` created because `pdflatex` is not available in this environment.
- No simulations were run during initialization.

---

## Problem 1 run - 2026-06-04

### Global methods and assumptions

- Work performed only inside `hw2/`.
- Node order: `['o', 'a', 'b', 'c', 'd']`.
- Matrix convention: `Lambda[i,j] = rate/weight from source node i to destination node j`.
- CTMC generator: `Q = Lambda - diag(Lambda @ 1)`.
- Return-time convention for P1(a)-P1(b): initial holding time in node `a` is included. The simulated time starts immediately at node `a`, samples the first holding time in `a`, and stops at the next entrance to `a` after leaving.
- French-DeGroot convention for P1(e)-P1(h): continuous-time `dx/dt = -Lx`, with `L = diag(Lambda @ 1) - Lambda`, because lecture-specific notes were not available in this workspace.
- Python version: `3.8.2`.
- OS: `Windows 10`.
- Main libraries: numpy, matplotlib, Python standard library.

### P1(a)-P1(b): return time to node a

- Method: Monte Carlo CTMC simulation for P1(a); stationary CTMC cycle formula `E_a[T_a^+] = 1 / (pi_a omega_a)` for P1(b).
- Seed: `20260615`.
- Number of Monte Carlo runs: `100000`.
- Result file: `results/problem1_return_time_a.csv`.
- Simulation mean: `6.69375284456`.
- Simulation standard deviation: `4.9131262948`.
- Simulation standard error: `0.0155366695236`.
- 95 percent CI: [`6.66330097229`, `6.72420471683`].
- Theoretical value: `6.70833333333`.
- Simulation minus theory: `-0.0145804887734`.
- Relative error: `-0.00217348900969`.

### P1(c)-P1(d): hitting time from node o to node d

- Method: Monte Carlo CTMC simulation for P1(c); linear system `sum_j Q[i,j] h_j = -1` for `i != d`, with `h_d = 0`, for P1(d).
- Seed: `20260616`.
- Number of Monte Carlo runs: `100000`.
- Result file: `results/problem1_hitting_o_to_d.csv`.
- Simulation mean: `10.7608365852`.
- Simulation standard deviation: `8.94426026682`.
- Simulation standard error: `0.0282842344285`.
- 95 percent CI: [`10.7053994857`, `10.8162736847`].
- Theoretical value: `10.7666666667`.
- Simulation minus theory: `-0.00583008148547`.
- Relative error: `-0.000541493636421`.

### P1(e): original French-DeGroot dynamics

- Method: exact matrix exponential via eigendecomposition for `exp(-Lt)x(0)`; SCC/sink-component graph check.
- Result file: `results/problem1_fdg_original_consensus.csv`.
- Figure: `figures/problem1_fdg_original_trajectory.png`.
- SCCs: `[["o", "a", "b", "c", "d"]]`.
- Sink components: `[["o", "a", "b", "c", "d"]]`.
- Converges to consensus for every initial condition: `yes`.
- Consensus vector pi: `{"o": 0.21739130434782614, "a": 0.14906832298136644, "b": 0.26086956521739135, "c": 0.18633540372670807, "d": 0.186335403726708}`.
- Consensus value for saved initial condition: `0.990683229814`.
- Max absolute final-minus-consensus on plotted horizon: `8.48210390814e-14`.

### P1(f): variance of consensus value

- Method: theoretical variance `sum_i pi_i^2 Var(x_i(0))`; Monte Carlo sampling of independent zero-mean Gaussian initial states with requested variances.
- Seed: `20260617`.
- Number of Monte Carlo runs: `50000`.
- Result file: `results/problem1_consensus_variance.csv`.
- Figure: `figures/problem1_consensus_value_histogram.png`.
- Initial variances by node: `{"o": 1.0, "a": 2.0, "b": 2.0, "c": 2.0, "d": 1.0}`.
- Theoretical variance: `0.331970217198`.
- Monte Carlo variance: `0.331277875538`.
- Monte Carlo mean: `-0.00270137580843`.
- Monte Carlo minus theory: `-0.000692341660048`.
- Relative error: `-0.00208555353517`.

### P1(g): removed edges `(d,a), (d,c), (a,c), (b,c)`

- Method: set requested directed entries of Lambda to zero; compute SCCs and sink SCCs; simulate a representative trajectory.
- Result file: `results/problem1_removed_edges_g_summary.csv`.
- Figure: `figures/problem1_fdg_removed_edges_g.png`.
- SCCs: `[["o", "a", "b"], ["c"], ["d"]]`.
- Sink components: `[["o", "a", "b"], ["d"]]`.
- Converges to global consensus for every initial condition: `no`.
- Asymptotic state for saved initial condition: `{"o": 0.14634146341463458, "a": 0.1463414634146345, "b": 0.1463414634146345, "c": 3.3821138211382116, "d": 5.0}`.
- Interpretation note: `two_sink_components_so_limit_depends_on_sink_component_initial_values_and_transient_routing`.

### P1(h): removed edges `(b,o), (d,a)`

- Method: set requested directed entries of Lambda to zero; compute SCCs and sink SCCs; simulate two representative trajectories.
- Result file: `results/problem1_removed_edges_h_summary.csv`.
- Figures: `figures/problem1_fdg_removed_edges_h_case1.png`, `figures/problem1_fdg_removed_edges_h_case2.png`.
- Case 1 asymptotic state: `{"o": 2.749999999999979, "a": 2.749999999999979, "b": 2.7499999999999787, "c": 2.749999999999979, "d": 2.7499999999999796}`.
- Case 2 asymptotic state: `{"o": 1.9428902930940094e-16, "a": 1.9428902930940094e-16, "b": 1.9428902930940092e-16, "c": 1.9428902930940094e-16, "d": 1.94289029309401e-16}`.
- SCCs: `[["o"], ["a"], ["b", "c", "d"]]`.
- Sink components: `[["b", "c", "d"]]`.
- Converges to global consensus for every initial condition: `yes`.
- Interpretation note: `single_sink_component_reachable_from_all_so_global_consensus_depends_on_sink_component_initial_values`.

### Report-writing notes and ambiguities

- Confirm whether course lecture notes use the same continuous-time French-DeGroot convention.
- State explicitly in the report that return times include the initial holding time in node `a`.
- These factual notes were used as source material for the final report prose.

### Problem 1 report section status - 2026-06-04

- Problem 1 report section written in `report/report.tex`.
- Figures inserted:
  - `figures/problem1_fdg_original_trajectory.png`
  - `figures/problem1_consensus_value_histogram.png`
  - `figures/problem1_fdg_removed_edges_g.png`
  - `figures/problem1_fdg_removed_edges_h_case1.png`
  - `figures/problem1_fdg_removed_edges_h_case2.png`
- Numerical values checked against CSV files in `results/`.
- Missing items: none among the requested Problem 1 result files and figures.
- Problem 2 and Problem 3 report sections were not written.

---

## Problem 2 - Many particles

Solved on 2026-06-04.

### Methods and assumptions

- Work performed only inside `hw2/`.
- Same closed network as Problem 1; node order `['o', 'a', 'b', 'c', 'd']`.
- Matrix convention: `Lambda[i,j]` is the rate from source node `i` to destination node `j`.
- CTMC definitions: `omega = Lambda @ 1`, `Q = Lambda - diag(omega)`, `P_jump = diag(omega)^(-1) Lambda`.
- `P_jump` was used only as the embedded jump-chain transition matrix after a departure event.
- No uniformized transition matrix was used.
- P2(a) uses the same return-time convention as Problem 1: the initial holding time in node `a` is included.
- P2(b) uses a direct event-based node-count CTMC: node `i` has departure rate `n_i(t) omega_i`.

### P2(a): particle perspective

- Seed: `20260624`.
- Number of particles per Monte Carlo run: `100`.
- Number of Monte Carlo runs: `2000`.
- Total particle return times simulated: `200000`.
- Result file: `results/problem2_particle_return_times.csv`.
- Mean of run particle averages: `6.70478013087`.
- Standard deviation of run particle averages: `0.480086936306`.
- Standard error of run particle averages: `0.0107350702469`.
- 95 percent CI: [`6.68373939319`, `6.72582086855`].
- Problem 1 simulation mean: `6.69375284456`.
- Problem 1 theoretical value: `6.70833333333`.
- Difference from Problem 1 simulation mean: `0.0110272863097`.
- Difference from Problem 1 theory: `-0.00355320246031`.

### P2(b): node perspective

- Seed: `20260625`.
- Number of particles: `100`.
- Initial condition: all particles in node `a`.
- Simulation horizon: `60`.
- Number of Monte Carlo runs: `5000`.
- Time-series sampling step: `0.1`.
- Result files:
  - `results/problem2_node_final_counts.csv`
  - `results/problem2_node_time_series.csv`
- Figures:
  - `figures/problem2_node_counts_over_time.png`
  - `figures/problem2_stationary_distribution_comparison.png`
- Stationary distribution pi: `o: 0.217391304348, a: 0.149068322981, b: 0.260869565217, c: 0.186335403727, d: 0.186335403727`.
- Final count comparison:
- Node `o`: simulated average final count `21.7254`, stationary expected count `21.7391304348`, difference `-0.0137304347826`.
- Node `a`: simulated average final count `14.9592`, stationary expected count `14.9068322981`, difference `0.0523677018634`.
- Node `b`: simulated average final count `25.991`, stationary expected count `26.0869565217`, difference `-0.0959565217391`.
- Node `c`: simulated average final count `18.6466`, stationary expected count `18.6335403727`, difference `0.0130596273292`.
- Node `d`: simulated average final count `18.6778`, stationary expected count `18.6335403727`, difference `0.0442596273292`.

### Generated artifacts

- `C:/Users/aless/OneDrive/Desktop/network dynamics and learning/hw2/results/problem2_particle_return_times.csv`
- `C:/Users/aless/OneDrive/Desktop/network dynamics and learning/hw2/results/problem2_node_final_counts.csv`
- `C:/Users/aless/OneDrive/Desktop/network dynamics and learning/hw2/results/problem2_node_time_series.csv`
- `C:/Users/aless/OneDrive/Desktop/network dynamics and learning/hw2/figures/problem2_node_counts_over_time.png`
- `C:/Users/aless/OneDrive/Desktop/network dynamics and learning/hw2/figures/problem2_stationary_distribution_comparison.png`

### Problem 2 report section status - 2026-06-04

- Problem 2 report section written in `report/report.tex`.
- Particle-perspective agreement with Problem 1 explained.
- Node-count final averages compared with `N*pi`, with finite-time and Monte Carlo deviations explained.
- Figures inserted:
  - `figures/problem2_node_counts_over_time.png`
  - `figures/problem2_stationary_distribution_comparison.png`
- Numerical values checked against CSV files in `results/`:
  - `results/problem2_particle_return_times.csv`
  - `results/problem1_return_time_a.csv`
  - `results/problem2_node_final_counts.csv`
- Missing items: none among the requested Problem 2 result files and figures.

---

# Problem 3 - Open network

Solved on 2026-06-04.

### Methods and assumptions

- Work performed only inside `hw2/`.
- Node order: `['o', 'a', 'b', 'c', 'd']`.
- Matrix convention: `Lambda_open[i,j]` is the rate from source node `i` to destination node `j`.
- External arrivals enter node `o` according to a Poisson process with rate `lambda`.
- Service clock rates use `omega = Lambda_open @ 1`, with the assignment override `omega_d = 7/4`.
- Omega values used: `o: 2, a: 1, b: 1, c: 1, d: 1.75`.
- Event simulator supports service modes `proportional` and `fixed`.
- In proportional mode, node `i` has service rate `omega_i N_i(t)`.
- In fixed mode, node `i` has clock rate `omega_i`; empty-node clock ticks do nothing.
- For nodes `o,a,b,c`, a nonempty service event forwards one particle using normalized outgoing rates from `Lambda_open`.
- For node `d`, a nonempty service event removes one particle from the system.

### P3(a): proportional-rate scenario

- Main simulation seed: `20260904`.
- Main simulation lambda: `100`.
- Main simulation horizon: `60`.
- Result file: `results/problem3_proportional_timeseries_lambda100.csv`.
- Figure: `figures/problem3_proportional_lambda100_counts.png`.
- Stability scan seed base: `20260906`.
- Stability scan lambdas tested: `10`, `50`, `100`, `200`, `500`.
- Stability scan horizon: `120`.
- Stability scan replicates per lambda: `2`.
- Blow-up criterion for scan: classify as blow-up if the mean fitted slope of total population over the second half of the horizon exceeds `0.05` times the tested lambda.
- Result file: `results/problem3_proportional_stability_scan.csv`.
- Figure: `figures/problem3_proportional_stability_scan.png`.
- Theoretical motivation note: proportional service makes the open network a linear infinite-server network; the expected single-particle remaining lifetimes are `o: 2.94642857143, a: 2.32142857143, b: 2.57142857143, c: 1.57142857143, d: 0.571428571429`, so finite lambda gives finite expected population.
- Proportional scan summary:
- Lambda `10`: mean final total `30`, mean second-half drift `0.0136088709677`, classification `stable_by_scan`.
- Lambda `50`: mean final total `165.5`, mean second-half drift `-0.00756048387097`, classification `stable_by_scan`.
- Lambda `100`: mean final total `308`, mean second-half drift `-0.0405241935484`, classification `stable_by_scan`.
- Lambda `200`: mean final total `607`, mean second-half drift `-0.407056451613`, classification `stable_by_scan`.
- Lambda `500`: mean final total `1445.5`, mean second-half drift `-0.571370967742`, classification `stable_by_scan`.

### P3(b): fixed-rate scenario

- Main simulation seed: `20260905`.
- Main simulation lambda: `2`.
- Main simulation horizon: `6000`.
- Result file: `results/problem3_fixed_timeseries_lambda2.csv`.
- Figure: `figures/problem3_fixed_lambda2_counts.png`.
- Stability scan seed base: `20260907`.
- Stability scan lambdas tested: `0.5`, `1`, `1.2`, `1.3`, `1.33333333333`, `1.35`, `1.5`, `2`.
- Stability scan horizon: `6000`.
- Stability scan replicates per lambda: `3`.
- Blow-up criterion for scan: classify as blow-up if the mean fitted slope of total population over the second half of the horizon exceeds `0.005` or if the mean final total exceeds `1000`.
- Result file: `results/problem3_fixed_stability_scan.csv`.
- Figure: `figures/problem3_fixed_stability_scan.png`.
- Theoretical routing loads per unit lambda are saved in the run code as `[1, 1/2, 5/8, 3/4, 1]`.
- Fixed-rate capacity by node: `o: 2, a: 2, b: 1.6, c: 1.33333333333, d: 1.75`.
- Theoretical fixed-rate bottleneck node: `c`.
- Theoretical largest stable input rate: `1.33333333333`.
- Fixed scan summary:
- Lambda `0.5`: mean final total `2.66666666667`, mean second-half drift `-5.39096084582e-05`, classification `stable_by_scan`, theoretical classification `stable_below_bottleneck`.
- Lambda `1`: mean final total `7.33333333333`, mean second-half drift `-0.000332810503079`, classification `stable_by_scan`, theoretical classification `stable_below_bottleneck`.
- Lambda `1.2`: mean final total `13`, mean second-half drift `-0.00193249680493`, classification `stable_by_scan`, theoretical classification `stable_below_bottleneck`.
- Lambda `1.3`: mean final total `27`, mean second-half drift `0.00343644707796`, classification `stable_by_scan`, theoretical classification `stable_below_bottleneck`.
- Lambda `1.33333333333`: mean final total `125.333333333`, mean second-half drift `0.017733066109`, classification `blow_up_by_scan`, theoretical classification `not_stable_at_or_above_bottleneck`.
- Lambda `1.35`: mean final total `115.666666667`, mean second-half drift `0.00748809108865`, classification `blow_up_by_scan`, theoretical classification `not_stable_at_or_above_bottleneck`.
- Lambda `1.5`: mean final total `677.333333333`, mean second-half drift `0.107035029627`, classification `blow_up_by_scan`, theoretical classification `not_stable_at_or_above_bottleneck`.
- Lambda `2`: mean final total `3144.66666667`, mean second-half drift `0.518004182642`, classification `blow_up_by_scan`, theoretical classification `not_stable_at_or_above_bottleneck`.

### Generated artifacts

- `results/problem3_proportional_timeseries_lambda100.csv`
- `results/problem3_proportional_stability_scan.csv`
- `results/problem3_fixed_timeseries_lambda2.csv`
- `results/problem3_fixed_stability_scan.csv`
- `figures/problem3_proportional_lambda100_counts.png`
- `figures/problem3_proportional_stability_scan.png`
- `figures/problem3_fixed_lambda2_counts.png`
- `figures/problem3_fixed_stability_scan.png`

### Problem 3 report section status - 2026-06-04

- Problem 3 report section written in `report/report.tex`.
- Figures inserted:
  - `figures/problem3_proportional_lambda100_counts.png`
  - `figures/problem3_proportional_stability_scan.png`
  - `figures/problem3_fixed_lambda2_counts.png`
  - `figures/problem3_fixed_stability_scan.png`
- Numerical values checked against CSV files in `results/`:
  - `results/problem3_proportional_timeseries_lambda100.csv`
  - `results/problem3_proportional_stability_scan.csv`
  - `results/problem3_fixed_timeseries_lambda2.csv`
  - `results/problem3_fixed_stability_scan.csv`
- Missing items: none among the requested Problem 3 result files and figures.

---

## Finalization - 2026-06-04

### Report status

- `report/report.tex` finalized with title, author field, Introduction, Problem 1, Problem 2, Problem 3, Conclusion, and Collaboration statement.
- Collaboration statement used: `This work was completed individually. No collaboration was used beyond general course material and standard references.`
- No references section was added because the report does not cite external references beyond standard course material and formulas derived in the text.
- Problem 2 node-count values were corrected in the report and this log to match `results/problem2_node_final_counts.csv`. No simulation outputs were changed.

### Compilation status

- LaTeX compilation was attempted by checking for `pdflatex`, `latexmk`, `xelatex`, and `lualatex`.
- None of those commands is available in this environment, so `report/report.pdf` was not regenerated during finalization.
- Existing `report/report.pdf` was left unchanged.

### README status

- `README.md` finalized with homework purpose, folder structure, Python requirements, commands for Problems 1--3, output locations, report compilation instructions, and reproducibility notes.

### Checked result files

- `results/problem1_return_time_a.csv`
- `results/problem1_hitting_o_to_d.csv`
- `results/problem1_fdg_original_consensus.csv`
- `results/problem1_consensus_variance.csv`
- `results/problem1_removed_edges_g_summary.csv`
- `results/problem1_removed_edges_h_summary.csv`
- `results/problem2_particle_return_times.csv`
- `results/problem2_node_final_counts.csv`
- `results/problem2_node_time_series.csv`
- `results/problem3_proportional_timeseries_lambda100.csv`
- `results/problem3_proportional_stability_scan.csv`
- `results/problem3_fixed_timeseries_lambda2.csv`
- `results/problem3_fixed_stability_scan.csv`

### Checked figure files

- `figures/problem1_fdg_original_trajectory.png`
- `figures/problem1_consensus_value_histogram.png`
- `figures/problem1_fdg_removed_edges_g.png`
- `figures/problem1_fdg_removed_edges_h_case1.png`
- `figures/problem1_fdg_removed_edges_h_case2.png`
- `figures/problem2_node_counts_over_time.png`
- `figures/problem2_stationary_distribution_comparison.png`
- `figures/problem3_proportional_lambda100_counts.png`
- `figures/problem3_proportional_stability_scan.png`
- `figures/problem3_fixed_lambda2_counts.png`
- `figures/problem3_fixed_stability_scan.png`

### Cleanup checks

- Searched for the requested cleanup keywords; no unresolved matches remained after finalization.
- Verified that every `\includegraphics` path in `report/report.tex` resolves to an existing file from the `report/` directory.
- Verified that stale Problem 2 final-count values from the older 500-run text no longer appear in `report/report.tex` or `RESULTS_LOG.md`.

### Remaining issues

- LaTeX is not installed in this environment, so the final PDF could not be regenerated here.
- The continuous-time French-DeGroot convention was documented and used consistently, but the course lecture-note convention was not independently available in this workspace.

---

## Final content-fix pass - 2026-06-04

### What was checked

- Checked `HW2_CONTEXT.md`, `RESULTS_LOG.md`, `README.md`, `src/`, `results/`, `figures/`, `report/report.tex`, `HW2_STUDY_GUIDE.md`, and `FIX_PROMPTS.md`.
- Checked that the report answers the requested Homework 2 content: CTMC definitions, return and hitting times, French-DeGroot consensus, edge-removal cases, many-particle particle/node perspectives, stationary count comparison, and open-network stability.
- Checked that all `\includegraphics` entries in `report/report.tex` point to existing files in `figures/`.
- Checked that the Problem 2 final counts in `report/report.tex` match `results/problem2_node_final_counts.csv`.
- Checked that the Problem 3 proportional final counts at `T=60` match `results/problem3_proportional_timeseries_lambda100.csv`.
- Checked that the Problem 3 fixed final counts at `T=6000` match `results/problem3_fixed_timeseries_lambda2.csv`.
- Checked the previously detected Problem 2 run-count mismatch between source code and saved outputs.

### What was fixed

- Updated `src/run_problem2.py` so `NODE_COUNT_RUNS = 5000`, matching:
  - `results/problem2_node_time_series.csv`;
  - `results/problem2_node_final_counts.csv`;
  - `report/report.tex`;
  - this log.
- Updated `report/report.tex` content to make the row-source convention, `P_jump` versus uniformization distinction, and proportional-rate stability interpretation clearer.
- Updated `README.md` to state how report values connect to CSV results and to document key reproducibility parameters.
- Updated `HW2_STUDY_GUIDE.md` to remove the now-resolved Problem 2 run-count inconsistency and replace it with a general source/result synchronization warning.

### Code/results/figures regenerated

- No simulation scripts were rerun.
- No CSV/result files were regenerated.
- No figure files were regenerated.
- Reason: the saved Problem 2 outputs already record 5000 node-perspective Monte Carlo runs; the required fix was to align the source constant with those saved outputs.

### Verification performed

- Ran `python -m py_compile src\run_problem2.py`; syntax check passed.
- Verified all report figure references exist:
  - `figures/problem1_fdg_original_trajectory.png`
  - `figures/problem1_consensus_value_histogram.png`
  - `figures/problem1_fdg_removed_edges_g.png`
  - `figures/problem1_fdg_removed_edges_h_case1.png`
  - `figures/problem1_fdg_removed_edges_h_case2.png`
  - `figures/problem2_node_counts_over_time.png`
  - `figures/problem2_stationary_distribution_comparison.png`
  - `figures/problem3_proportional_lambda100_counts.png`
  - `figures/problem3_proportional_stability_scan.png`
  - `figures/problem3_fixed_lambda2_counts.png`
  - `figures/problem3_fixed_stability_scan.png`
- Verified Problem 2 final average counts and stationary expected counts:
  - `o`: `21.7254` versus `21.7391304348`, with `5000` node-perspective Monte Carlo runs.
  - `a`: `14.9592` versus `14.9068322981`, with `5000` node-perspective Monte Carlo runs.
  - `b`: `25.991` versus `26.0869565217`, with `5000` node-perspective Monte Carlo runs.
  - `c`: `18.6466` versus `18.6335403727`, with `5000` node-perspective Monte Carlo runs.
  - `d`: `18.6778` versus `18.6335403727`, with `5000` node-perspective Monte Carlo runs.
- Verified Problem 3 proportional final counts at `lambda = 100`, `T = 60`: `(48, 39, 72, 62, 63)`, total `284`.
- Verified Problem 3 fixed final counts at `lambda = 2`, `T = 6000`: `(233, 15, 1240, 1577, 6)`, total `3071`.

### Remaining unresolved issues

- The continuous-time French-DeGroot convention is used consistently in code, results, log, study guide, and report, but the course lecture-note convention was not independently available in this workspace.
- LaTeX compilation issues were deliberately ignored in this pass, as requested. No missing packages, compiler warnings, or PDF build issues were treated as problems.

---

## Second content-fix pass - 2026-06-04

### Changes applied to `report/report.tex`

- Fixed author placeholder: replaced "Student Name" with "Alessandro Casadei".
- Added closed-network omega vector $\omega = (3/5, 1, 5/6, 1, 2/3)$ after the Lambda definition in Problem 1.
- Strengthened CTMC return-time explanation in Problem 1(b): added paragraph explaining why $\omega_a$ appears in $1/(\pi_a \omega_a)$.
- Specified random initialisation distribution for Problem 1(f): confirmed Gaussian from `rng.normal()` in `src/run_problem1.py` line 150; added explicit statement of $\mathcal{N}(0,1)$ and $\mathcal{N}(0,2)$ distributions.
- Improved Figure 2 caption to mention comparison of empirical consensus distribution with theoretical variance.
- Added explicit asymptotic formula for Problem 1(g): $x_c^\infty = (1/3) z_{oab} + (2/3) x_d(0) = 3.3821$, verified against saved asymptotic state.
- Added explicit left-eigenvector formula for Problem 1(h): $\pi^{(h)} = (0, 0, 1/4, 1/4, 1/2)$, with consensus value formulas for both representative initial conditions verified: $2.75$ and $0$.
- Added Little's-law stability reasoning for Problem 3(a) proportional service; clarified that the stability scan is implementation evidence, not the theoretical proof.
- Replaced fixed-rate threshold wording in Problem 3(b): "bottleneck node is $c$, giving critical threshold $\lambda_c = 4/3$" and "$4/3$ is the supremum of sustainable input rates".
- Clarified Problem 3 figure captions: both proportional ($\lambda=100$) and fixed ($\lambda=2$) trajectory figures are single sample paths, verified from `src/run_problem3.py` which calls `simulate_open_network` (not `simulate_replicates`) for the main cases.
- Reduced excessive decimal precision in prose values to 4--5 significant digits where safe; CSVs not modified.

### Simulations rerun

- None. All changes are text-only in `report/report.tex`.

### Compile status

- `pdflatex` is not available in this environment. PDF must be regenerated externally.
