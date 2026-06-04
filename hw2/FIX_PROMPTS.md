# FIX_PROMPTS.md - Homework 2

This file contains standalone prompts to give to a future coding/report agent. Each prompt is scoped to technical, mathematical, code, result, and report-content consistency only. Ignore LaTeX compilation problems.

## Prompt A - Fix Problem 1 technical/mathematical issues

### Title
Fix Problem 1 CTMC and French-DeGroot consistency.

### Objective
Audit and fix the technical, mathematical, code, result, figure, and report-content consistency for Problem 1.

### Files to inspect
- `hw2/HW2_CONTEXT.md`
- `hw2/RESULTS_LOG.md`
- `hw2/README.md`
- `hw2/src/constants.py`
- `hw2/src/problem1_single_ctmc.py`
- `hw2/src/problem1_opinion.py`
- `hw2/src/run_problem1.py`
- `hw2/src/utils.py`
- `hw2/src/plotting.py`
- `hw2/results/problem1_return_time_a.csv`
- `hw2/results/problem1_hitting_o_to_d.csv`
- `hw2/results/problem1_fdg_original_consensus.csv`
- `hw2/results/problem1_consensus_variance.csv`
- `hw2/results/problem1_removed_edges_g_summary.csv`
- `hw2/results/problem1_removed_edges_h_summary.csv`
- `hw2/figures/problem1_fdg_original_trajectory.png`
- `hw2/figures/problem1_consensus_value_histogram.png`
- `hw2/figures/problem1_fdg_removed_edges_g.png`
- `hw2/figures/problem1_fdg_removed_edges_h_case1.png`
- `hw2/figures/problem1_fdg_removed_edges_h_case2.png`
- `hw2/report/report.tex`, for content consistency only.

### Files allowed to modify
- `hw2/src/constants.py`
- `hw2/src/problem1_single_ctmc.py`
- `hw2/src/problem1_opinion.py`
- `hw2/src/run_problem1.py`
- `hw2/src/utils.py`
- `hw2/src/plotting.py`
- `hw2/results/*.csv`
- `hw2/figures/*.png`
- `hw2/RESULTS_LOG.md`
- `hw2/report/report.tex`
- `hw2/README.md`, only if commands or output descriptions change.

### Exact task
Work only inside `hw2/`. Do not modify `hw1/`, `hw3/`, or `shared/`. Check the CTMC generator convention, return-time convention, theoretical return time, theoretical hitting time, simulation/theory comparison, French-DeGroot convention, consensus condition, consensus variance formula, edge-removal cases, SCC/sink-component reasoning, and consistency between code, CSVs, figures, and report text.

In particular:
- Verify that the matrix orientation is consistently `Lambda[i,j] = rate/weight from source node i to destination node j`.
- Verify that `omega = Lambda @ 1`, `Q = Lambda - diag(omega)`, and `P_jump[i,j] = Lambda[i,j] / omega_i` are used consistently.
- Verify that the return time `T_a^+` includes the initial holding time in node `a` everywhere, or change all code/results/report text consistently if another convention is required by the course.
- Verify that the theoretical return-time formula `1 / (pi_a omega_a)` matches the chosen return-time convention.
- Verify that the hitting-time linear system uses `sum_j Q[i,j] h_j = -1` for `i != d` and `h_d = 0`.
- Verify that simulation means, standard errors, confidence intervals, theoretical values, and relative errors in the report match the saved CSV files.
- Verify that continuous-time French-DeGroot is intended: `dx/dt = -Lx`, `L = diag(Lambda @ 1) - Lambda`. If course notes require a different convention, update code/results/report consistently and document it.
- Verify that the original graph consensus claim is justified by SCC/sink-component structure and the invariant left eigenvector.
- Verify the consensus variance formula `Var(x*) = sum_i pi_i^2 Var(x_i(0))` and the Monte Carlo validation.
- Verify edge removals in parts (g) and (h) with row-source/column-destination convention.
- Recompute SCCs, sink SCCs, asymptotic states, and figures if needed.
- Save figures in `figures/`.
- Save CSV/tables in `results/`.
- Update `RESULTS_LOG.md`.
- Keep the report readable independently of code.
- Ignore report compilation problems.

### Expected outputs
- Corrected Problem 1 code, if inconsistencies are found.
- Regenerated Problem 1 CSVs and figures, if code or model conventions change.
- Updated `RESULTS_LOG.md` with factual methods, seeds, result paths, and values.
- Updated Problem 1 report prose in `report/report.tex` using only values from CSVs.
- A short final audit summary listing what changed and what was verified.

### Checks to perform
- Run `python src/run_problem1.py` from `hw2/` if code or numerical outputs change.
- Confirm every Problem 1 numerical value in `report/report.tex` appears in the corresponding CSV or is explicitly derived from it.
- Confirm every Problem 1 figure referenced by the report exists in `figures/`.
- Confirm no report text overclaims consensus without graph-structure proof.
- Confirm edge-removal descriptions match the actual removed directed edges.
- Confirm `RESULTS_LOG.md` agrees with source constants and generated outputs.

### What not to do
- Do not work outside `hw2/`.
- Do not modify `hw1/`, `hw3/`, or `shared/`.
- Do not invent numerical results.
- Do not manually edit CSV values unless they are produced by a checked computation and the computation is documented.
- Do not save figures outside `figures/`.
- Do not save CSV/tables outside `results/`.
- Do not remove documented seeds or reproducibility information.
- Do not make the report depend on reading the code.
- Do not fix, discuss, or prioritize LaTeX compilation problems.

## Prompt B - Fix Problem 2 technical/mathematical issues

### Title
Fix Problem 2 many-particle CTMC consistency.

### Objective
Audit and fix Problem 2 particle-perspective and node-perspective consistency, especially the number of node-count Monte Carlo runs.

### Files to inspect
- `hw2/HW2_CONTEXT.md`
- `hw2/RESULTS_LOG.md`
- `hw2/README.md`
- `hw2/src/constants.py`
- `hw2/src/problem1_single_ctmc.py`
- `hw2/src/problem2_many_particles.py`
- `hw2/src/run_problem2.py`
- `hw2/results/problem1_return_time_a.csv`
- `hw2/results/problem2_particle_return_times.csv`
- `hw2/results/problem2_node_time_series.csv`
- `hw2/results/problem2_node_final_counts.csv`
- `hw2/figures/problem2_node_counts_over_time.png`
- `hw2/figures/problem2_stationary_distribution_comparison.png`
- `hw2/report/report.tex`, for content consistency only.

### Files allowed to modify
- `hw2/src/problem2_many_particles.py`
- `hw2/src/run_problem2.py`
- `hw2/src/constants.py`, only if a shared convention is wrong.
- `hw2/results/problem2_particle_return_times.csv`
- `hw2/results/problem2_node_time_series.csv`
- `hw2/results/problem2_node_final_counts.csv`
- `hw2/figures/problem2_node_counts_over_time.png`
- `hw2/figures/problem2_stationary_distribution_comparison.png`
- `hw2/RESULTS_LOG.md`
- `hw2/report/report.tex`
- `hw2/README.md`, only if commands or output descriptions change.

### Exact task
Work only inside `hw2/`. Do not modify `hw1/`, `hw3/`, or `shared/`. Check the particle perspective, node perspective, `N = 100` particles, `T = 60`, embedded jump chain versus uniformized transition matrix, stationary distribution of the CTMC, final-count comparison with `N*pi`, and consistency between simulation, CSVs, figures, and report text.

Resolve the detected consistency issue: `src/run_problem2.py` currently states `NODE_COUNT_RUNS = 500`, while `results/problem2_node_time_series.csv`, `results/problem2_node_final_counts.csv`, `RESULTS_LOG.md`, and `report/report.tex` state/use `num_monte_carlo_runs = 5000`. Determine which is intended. If the saved outputs came from 5000 runs, update the source constant and rerun only if needed. If 500 runs is intended, regenerate outputs and update the report/log. Do not invent numerical results.

Also verify:
- `NUM_PARTICLES = 100` for both perspectives.
- Particle-perspective return times include the same initial holding-time convention as Problem 1.
- The particle-perspective comparison uses values from `results/problem1_return_time_a.csv`.
- `P_jump` is used only as the embedded jump-chain transition matrix after actual departures.
- No uniformized transition matrix is accidentally used or described as `P_jump`.
- Node-count simulation starts all 100 particles in node `a`.
- Node departure rate is `n_i(t) omega_i`.
- Final counts at `T = 60` are compared with `N*pi`, where `pi Q = 0` and `sum_i pi_i = 1`.
- Report text explains that finite-time counts may differ from stationarity because of finite horizon and Monte Carlo error.
- Save figures in `figures/`.
- Save CSV/tables in `results/`.
- Update `RESULTS_LOG.md`.
- Keep the report readable independently of code.
- Ignore report compilation problems.

### Expected outputs
- Corrected `src/run_problem2.py` and/or regenerated Problem 2 results, depending on the chosen run-count resolution.
- Updated Problem 2 CSVs and figures if rerun.
- Updated `RESULTS_LOG.md` documenting seeds, run counts, and file paths.
- Updated Problem 2 report prose that matches the CSVs exactly.
- A final summary of the resolved run-count inconsistency.

### Checks to perform
- Run `python src/run_problem2.py` from `hw2/` if source constants or outputs change.
- Confirm `NODE_COUNT_RUNS`, CSV `num_monte_carlo_runs`, `RESULTS_LOG.md`, and `report/report.tex` all agree.
- Confirm every Problem 2 numerical value in the report comes from `results/`.
- Confirm figures referenced in the report exist.
- Confirm final counts sum close to 100 for each reported average vector.
- Confirm `P_jump` rows sum to one for positive exit-rate states.

### What not to do
- Do not work outside `hw2/`.
- Do not modify `hw1/`, `hw3/`, or `shared/`.
- Do not invent numerical results.
- Do not silently change the Monte Carlo run count without documenting it.
- Do not confuse `P_jump` with a uniformized transition matrix.
- Do not save figures outside `figures/`.
- Do not save CSV/tables outside `results/`.
- Do not leave report values that do not match CSV files.
- Do not fix, discuss, or prioritize LaTeX compilation problems.

## Prompt C - Fix Problem 3 technical/mathematical issues

### Title
Fix Problem 3 open-network consistency.

### Objective
Audit and fix Problem 3 open-network matrix, simulation rules, stability reasoning, blow-up classification, and report/result consistency.

### Files to inspect
- `hw2/HW2_CONTEXT.md`
- `hw2/RESULTS_LOG.md`
- `hw2/README.md`
- `hw2/src/constants.py`
- `hw2/src/problem3_open_network.py`
- `hw2/src/run_problem3.py`
- `hw2/results/problem3_proportional_timeseries_lambda100.csv`
- `hw2/results/problem3_proportional_stability_scan.csv`
- `hw2/results/problem3_fixed_timeseries_lambda2.csv`
- `hw2/results/problem3_fixed_stability_scan.csv`
- `hw2/figures/problem3_proportional_lambda100_counts.png`
- `hw2/figures/problem3_proportional_stability_scan.png`
- `hw2/figures/problem3_fixed_lambda2_counts.png`
- `hw2/figures/problem3_fixed_stability_scan.png`
- `hw2/report/report.tex`, for content consistency only.

### Files allowed to modify
- `hw2/src/constants.py`
- `hw2/src/problem3_open_network.py`
- `hw2/src/run_problem3.py`
- `hw2/results/problem3_*.csv`
- `hw2/figures/problem3_*.png`
- `hw2/RESULTS_LOG.md`
- `hw2/report/report.tex`
- `hw2/README.md`, only if commands or output descriptions change.

### Exact task
Work only inside `hw2/`. Do not modify `hw1/`, `hw3/`, or `shared/`. Check the open-network matrix, external Poisson arrivals, proportional-rate simulation, fixed-rate simulation, `omega_d = 7/4`, removal rule at node `d`, stability/blow-up thresholds, analytical motivation of the largest sustainable input rate, simulation scan consistency, and consistency between CSVs, figures, and report text.

In particular:
- Verify `LAMBDA_OPEN` uses the row-source/column-destination convention.
- Verify external arrivals enter node `o` as a Poisson process with rate `lambda`.
- Verify `omega = Lambda_open @ 1` with the assignment override `omega_d = 7/4`.
- Verify node `d` removes one particle on a nonempty service event and does not route particles onward.
- Verify proportional mode uses service rate `omega_i N_i(t)`.
- Verify fixed mode uses service clock rate `omega_i`, with empty-node ticks doing nothing.
- Verify the proportional-rate stability claim: if the model is an infinite-server open network with finite absorption times, explain why every finite input rate has finite expected population. If the course expects a finite scan-based largest tested stable rate instead, state that separately and do not overclaim.
- Verify the fixed-rate bottleneck calculation from routing loads `[1, 1/2, 5/8, 3/4, 1]`, capacities `omega_i / load_i`, bottleneck node `c`, and threshold `4/3`.
- Verify the scan classification rule is documented and not presented as proof by itself.
- Verify the main time-series final counts quoted in the report match the last rows of their CSV files.
- Save figures in `figures/`.
- Save CSV/tables in `results/`.
- Update `RESULTS_LOG.md`.
- Keep the report readable independently of code.
- Ignore report compilation problems.

### Expected outputs
- Corrected Problem 3 code if inconsistencies are found.
- Regenerated Problem 3 CSVs and figures if code or scan rules change.
- Updated `RESULTS_LOG.md` with factual methods, seeds, thresholds, and results.
- Updated Problem 3 report prose using only saved result values and checked theoretical calculations.
- Clear distinction between theoretical stability threshold and simulation scan evidence.

### Checks to perform
- Run `python src/run_problem3.py` from `hw2/` if code or numerical outputs change.
- Confirm every Problem 3 numerical value in `report/report.tex` appears in CSVs or is a documented theoretical calculation.
- Confirm `omega_d = 7/4` is used in both service modes.
- Confirm all Problem 3 figures referenced in the report exist.
- Confirm scan CSV classifications match the stated blow-up criterion.
- Confirm `RESULTS_LOG.md` and report text agree on lambdas, horizons, seeds, and thresholds.

### What not to do
- Do not work outside `hw2/`.
- Do not modify `hw1/`, `hw3/`, or `shared/`.
- Do not invent numerical results.
- Do not derive the stability threshold only from plots.
- Do not present a single stochastic run as conclusive proof of stability or blow-up.
- Do not save figures outside `figures/`.
- Do not save CSV/tables outside `results/`.
- Do not make the report depend on reading the code.
- Do not fix, discuss, or prioritize LaTeX compilation problems.

## Prompt D - Fix report content only, excluding compilation

### Title
Fix report content consistency without touching compilation.

### Objective
Improve only the mathematical and explanatory content of `report/report.tex`, excluding all LaTeX compilation concerns.

### Files to inspect
- `hw2/HW2_CONTEXT.md`
- `hw2/RESULTS_LOG.md`
- `hw2/README.md`
- `hw2/report/report.tex`
- `hw2/src/*.py`
- `hw2/results/*.csv`
- `hw2/figures/*.png`

### Files allowed to modify
- `hw2/report/report.tex`
- `hw2/RESULTS_LOG.md`, only to document report-content changes or discovered content inconsistencies.
- `hw2/README.md`, only if report/output descriptions are factually wrong.

### Exact task
Work only inside `hw2/`. Do not modify `hw1/`, `hw3/`, or `shared/`. Fix report content only. Do not modify code, regenerate numerical outputs, or edit figures unless explicitly asked in a separate task. Check whether all subquestions are answered, methods are explained clearly, all numerical values come from result files, figures are referenced and interpreted, conclusions are not overclaimed, and the report is readable independently of code.

Required checks:
- Verify every problem and sub-question has a clear model, method, result, and interpretation.
- Verify all numerical values in the report can be traced to `results/*.csv` or to formulas documented in `RESULTS_LOG.md`.
- Verify figure captions and surrounding text explain what each figure supports.
- Verify report text clearly states the row-source/column-destination convention.
- Verify report text clearly distinguishes `Lambda`, `Q`, `P_jump`, and any uniformized matrix concept.
- Verify report text clearly states the return-time convention.
- Verify report text does not overclaim proportional-rate stability beyond the mathematical argument and finite scan evidence.
- Verify report text mentions the remaining ambiguity that the course lecture-note French-DeGroot convention was not independently available, if that remains true.
- Save figures in `figures/` if any report-only figure reference correction requires a new figure request, but do not generate figures unless explicitly needed.
- Save CSV/tables in `results/` if any table must be produced by a checked computation, but do not invent numerical results.
- Update `RESULTS_LOG.md`.
- Keep the report readable independently of code.
- Ignore report compilation problems.

### Expected outputs
- Updated `report/report.tex` content only.
- Updated `RESULTS_LOG.md` note describing content edits and any remaining factual ambiguity.
- A short checklist of report-content issues fixed.

### Checks to perform
- Search `report/report.tex` for every number and verify it against CSV/log/theoretical formulas.
- Confirm every `\includegraphics` path points to an existing generated figure, but do not address compilation warnings or package errors.
- Confirm all report conclusions are supported by code/results/math.
- Confirm no report prose says values were visually estimated if CSV values exist.

### What not to do
- Do not work outside `hw2/`.
- Do not modify `hw1/`, `hw3/`, or `shared/`.
- Do not invent numerical results.
- Do not alter Python code or rerun simulations in this report-only task unless explicitly asked.
- Do not save figures outside `figures/`.
- Do not save CSV/tables outside `results/`.
- Do not add LaTeX package fixes, compilation fixes, build-system fixes, or PDF build troubleshooting.
- Do not treat missing TeX packages, PDF build errors, or LaTeX warnings as problems to solve.

## Prompt E - Final consistency audit, excluding compilation

### Title
Perform final HW2 consistency audit excluding compilation.

### Objective
Audit the entire HW2 submission for code/result/report consistency, reproducibility, missing explanations, and oral-defense readiness.

### Files to inspect
- `hw2/HW2_CONTEXT.md`
- `hw2/RESULTS_LOG.md`
- `hw2/README.md`
- `hw2/src/`
- `hw2/results/`
- `hw2/figures/`
- `hw2/report/report.tex`, for content consistency only.

### Files allowed to modify
- `hw2/RESULTS_LOG.md`
- `hw2/README.md`
- `hw2/report/report.tex`
- `hw2/src/*.py`, only if a real code/result/report inconsistency must be fixed.
- `hw2/results/*.csv`, only by rerunning checked scripts or writing documented checked tables.
- `hw2/figures/*.png`, only by rerunning checked scripts.

### Exact task
Work only inside `hw2/`. Do not modify `hw1/`, `hw3/`, or `shared/`. Perform a final consistency audit excluding compilation. Focus on code/result/report consistency, numerical consistency, figure consistency, reproducibility, missing explanations, and oral-defense readiness.

Audit checklist:
- Compare source constants and run-script constants against CSV metadata, `RESULTS_LOG.md`, and report text.
- Specifically check the detected Problem 2 inconsistency: `src/run_problem2.py` says `NODE_COUNT_RUNS = 500`, while saved outputs/log/report say 5000. Resolve or document it.
- Verify every generated result file in `results/` is documented and has an explained purpose.
- Verify every generated figure in `figures/` is documented, referenced if relevant, and interpretable.
- Verify every numerical claim in `report/report.tex` comes from CSVs or checked formulas.
- Verify `README.md` commands reproduce the saved outputs.
- Verify random seeds and Monte Carlo repetitions are documented.
- Verify mathematical conventions are consistent across all files: node order, matrix orientation, `omega`, `Q`, `P_jump`, return time, hitting time, French-DeGroot, open-network service rules.
- Verify conclusions are supported and not overclaimed.
- Save figures in `figures/`.
- Save CSV/tables in `results/`.
- Update `RESULTS_LOG.md`.
- Keep the report readable independently of code.
- Ignore report compilation problems.

### Expected outputs
- A final audit entry in `RESULTS_LOG.md`.
- Any necessary targeted corrections to code/results/report/README.
- A short final audit summary with remaining unresolved issues.

### Checks to perform
- Run the relevant problem scripts only if a correction requires regenerated results.
- Confirm all CSV files named in `RESULTS_LOG.md` exist.
- Confirm all figure files named in `RESULTS_LOG.md` and `report/report.tex` exist.
- Confirm no stale numerical values remain in report/log.
- Confirm oral-defense explanations are present for the core theory: CTMC generator, return/hitting times, consensus, particle/node perspectives, stationary counts, open-network stability.

### What not to do
- Do not work outside `hw2/`.
- Do not modify `hw1/`, `hw3/`, or `shared/`.
- Do not invent numerical results.
- Do not use the final audit to perform broad refactors.
- Do not save figures outside `figures/`.
- Do not save CSV/tables outside `results/`.
- Do not treat LaTeX compilation problems as audit findings.
- Do not fix, discuss, or prioritize missing TeX packages, PDF build errors, or LaTeX warnings.
