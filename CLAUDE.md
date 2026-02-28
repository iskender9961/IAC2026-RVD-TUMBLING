# CLAUDE.md -- Project conventions for AI assistants

## Project Overview
MPC-based rendezvous and docking (RVD) simulation with a tumbling target spacecraft.
Body-fixed LOS tetrahedral corridor constraint, nonlinear ECI truth dynamics, LTV MPC prediction in target body frame, OSQP solver.

## Architecture
- **Plant (truth):** Nonlinear ECI propagation with 2-body + J2 (ode113, tight tolerances)
- **MPC prediction:** LTV discrete model in target body frame, re-linearized every step via finite differences
- **Solver:** OSQP (rebuilt each step due to sparsity pattern changes in Ad/Bd)
- **LOS constraint:** Polyhedral (tetrahedral) inner approximation of circular cone, body-fixed to tumbling target

## Key Conventions
- Quaternion: scalar-first `[qw; qx; qy; qz]`
- `R_eci_tb`: rotates from target body to ECI (`v_eci = R_eci_tb * v_tb`)
- Transport theorem: `v_tb = R_tb_eci * dv_eci - omega_body x r_tb`
- LVLH: x=radial (outward), y=along-track, z=orbit normal
- Docking axis: +yT (target body y-axis)
- Decision variable ordering: `z = [x_0..x_Np; u_0..u_{Np-1}]`
- All units SI (meters, seconds, radians)

## File Structure
```
params.m             -- central config (all tuning knobs here)
main_sim.m           -- entry point, simulation loop, logging
run_sim_headless.m   -- headless simulation (used by MC and comparisons)
run_tests.m          -- 5 sanity tests
run_rdu_comparison.m -- Rdu sweep (4 scenarios) with comparison plots/GIFs
run_monte_carlo.m    -- MC feasibility sweep (parfor, saves .mat + plots)
plot_mc_figures.m    -- regenerate MC plots from saved mc_sweep_all.mat

dynamics/            -- ECI equations of motion
  accel_2body.m, accel_J2.m, eom_eci.m, propagate_truth_step.m

frames/              -- coordinate transformations
  lvlh_from_rv.m, quat_propagate.m, rotm_from_quat.m,
  target_attitude_model.m, transform_rel_to_TB.m, control_TB_to_ECI.m

mpc/                 -- MPC formulation and solver
  build_qp_osqp.m, update_qp_osqp.m, linearize_discrete_model.m,
  los_tetra_constraints.m

viz/                 -- visualization
  plot_all_results.m, plot_time_histories.m, generate_gifs.m,
  generate_comparison_gifs.m, draw_los_tetra.m, draw_cube.m,
  draw_triad.m, animate_scene.m

utils/               -- small helpers
  skew.m, clamp.m, assert_units.m

docs/                -- LaTeX paper and README
results/             -- simulation outputs (.mat, .png, .gif)
```

## Known Issues
- OSQP problem must be rebuilt (not updated) each step because finite-difference
  linearization can produce different sparsity patterns in Ad/Bd across steps.
- MATLAB `-batch` mode does not support local functions in scripts; keep helpers
  as separate .m files or inside function-based scripts.
- `gh` CLI may need manual auth setup on Windows; use `curl` + git credential
  token as fallback for GitHub API calls.

## Tuning Guide
All tuning parameters are in `params.m`. Key knobs:
- `Q`: Stage state weight. Heavy on x_TB/z_TB drives docking axis tracking.
- `QN`: Terminal weight (multiplier of Q). Higher = stronger terminal pull.
- `Rdu`: Delta-u penalty. Controls smoothness vs responsiveness.
- `Ru`: Absolute input penalty. Keep small (regularization only).
- `Np`: Prediction horizon. Longer = better anticipation of tumble, slower solve.
- `u_max`: Max thrust. Must be large enough to track tumbling target.
- `omega_body`: Target tumble rate. Higher = harder problem.
- `y_hold_tau`: Exponential decay time constant for approach distance.

## Running
```matlab
>> main_sim          % full simulation + plots + GIFs
>> run_tests         % sanity checks (should print "All tests PASSED")
>> run_rdu_comparison  % Rdu sweep + comparison plots + GIFs
>> run_monte_carlo     % MC feasibility sweep (omega x a_max, parfor)
>> plot_mc_figures     % regenerate MC plots from saved data (no re-run)
```

## Performance
- Always use `parfor` for Monte Carlo / batch simulations (MATLAB Parallel Computing Toolbox)
- Use `parallel.pool.DataQueue` for progress reporting inside `parfor`
- For MC speed: reduce `Np`, `Tsim`, relax ODE tolerances, reduce `osqp_max_iter`

## MC Compatibility Contract
Any changes to the simulation code MUST keep MC scripts working. The chain is:
```
params.m → run_sim_headless.m → run_monte_carlo.m → plot_mc_figures.m
```
- `run_sim_headless(p)` is the MC workhorse. Its signature and output struct `lg` must not break.
- `lg` must always contain: `r_tb_hist`, `v_tb_hist`, `t_hist`, `u_hist`, `cost_hist`, `status_hist`
- `params.m` fields overridden by MC: `omega_body`, `u_max`, `Rdu`, `Np`, `Tsim`, `dr_lvlh0`, `dv_lvlh0`, `y_hold_start`, `y_hold_tau`, `cone_draw_L`, `ode_opts`, `osqp_max_iter`
- After any code change, verify: `run_tests`, then `run_monte_carlo` (or at minimum `run_sim_headless(params())`)
- `plot_mc_figures.m` reads from `results/mc_sweep_all.mat` -- re-run `run_monte_carlo` if `run_sim_headless` output format changes

## Do NOT
- Do not rename `los_tetra_constraints.m` or `draw_los_tetra.m` (old names were los_polyhedral_constraints / draw_cone_poly)
- Do not use `update('Ax', ...)` on OSQP -- always rebuild
- Do not add dependencies beyond MATLAB base + OSQP
- Do not change `run_sim_headless` return signature `[lg, p, sim_terminated, term_reason]`
- Do not remove fields from the `lg` struct -- only add new ones
