# MPC 90-Case Parameter Sweep

## Overview

This sweep evaluates a receding-horizon MPC controller for proximity approach
to a tumbling target using Clohessy-Wiltshire-Hill (CWH/HCW) linearised dynamics.

**Total: 90 cases x 10 initial conditions = 900 closed-loop simulations.**

## Sweep Grid

| Parameter          | Values                          | Count |
|--------------------|---------------------------------|-------|
| `a_max` (m/s^2)   | 0.10, 0.06, 0.02               | 3     |
| `y_limit` (m)      | 30, 100, 150                    | 3     |
| `omega` (deg/s)    | 1, 2, 3, 4, 5                  | 5     |
| Cost variant       | A (pos+du), B (pos+vel+du)      | 2     |

## Cost Function Variants

- **Cost A (baseline):** Penalises position error in body frame + delta-u
  smoothing + DARE terminal cost. No velocity penalty. No direct-u penalty.
- **Cost B (augmented):** Same as A, plus velocity error penalty in body frame.
  DARE terminal cost includes velocity weights.

Both cost functions use the same `r_du = 30.0` delta-u smoothing weight.

## Initial Conditions (10 per case)

All ICs place the chaser at the approach distance `y_limit` along the
target's +y_B axis (body-frame docking direction), with varying lateral
offsets in x_B:

| IC  | x_B (m) | y_B (m)     | z_B (m) |
|-----|---------|-------------|---------|
| 01  | 0       | y_limit     | 0       |
| 02  | -10     | y_limit     | 0       |
| 03  | +10     | y_limit     | 0       |
| 04  | -20     | y_limit     | 0       |
| 05  | +20     | y_limit     | 0       |
| 06  | -5      | y_limit     | 0       |
| 07  | +5      | y_limit     | 0       |
| 08  | -15     | y_limit     | 0       |
| 09  | +15     | y_limit     | 0       |
| 10  | 0       | y_limit + 5 | 0       |

Body-frame velocity is zero (co-rotation with target). The body-frame IC
is properly transformed to LVLH using `r_L = R(theta) r_B`,
`v_L = R(theta)(v_B + omega x r_B)`.

## Dynamics & Model

- **Truth model:** Discrete CWH at dt = 0.1s (truth sub-stepping, 10 sub-steps
  per MPC interval)
- **MPC model:** Discrete CWH at dt = 1.0s
- **Prediction horizon:** N = 30 steps (30s lookahead)
- **Simulation duration:** 600s (10 minutes)
- **Note:** Both truth and MPC models use linearised HCW. For nonlinear
  (J2+drag) truth propagation: **TODO**.

## Constraints

- **LOS cone:** Polyhedral (5 half-planes) in target body frame, rotating
  with the target at omega_t. Soft constraints via slack variables with
  scaled quadratic+linear penalty.
- **Thrust saturation:** Component-wise box `|u_i| <= a_max`.
- **LOS margin:** 0.05 m tightening for inter-sample rotation safety.
- **Docking:** Range <= 0.5 m and speed <= 0.05 m/s.

## Solver

OSQP (Operator Splitting QP Solver):
- eps_abs = 1e-4, eps_rel = 1e-4
- max_iter = 4000, polish = True
- adaptive_rho = True, scaling = 10
- Hessian normalisation applied for conditioning

Fallback braking controller for QP infeasibility (brake towards origin
+ velocity damping).

## Output Structure

```
results/sweep_mpc_90cases/
  summary.csv                           # Aggregated stats (1 row per case)
  README_sweep.md                       # This file
  a{a_max}_y{y_limit}_w{omega}_cost{A|B}/
    case_config.json                    # Full parameter set
    body_position_vs_time.pdf           # Body-frame [x_B, y_B, z_B] vs time
    control_and_dv_proxy.pdf            # Control + cumulative dV
    mpc_control_vs_time.pdf             # Control components vs time
    mpc_cost_breakdown.pdf              # Position/velocity/du/constraint costs
    mpc_range_vs_time.pdf               # Range to target vs time
    mpc_traj_body_3d.pdf                # 3D body-frame trajectory
    range_vs_time.pdf                   # Range vs time (duplicate for compatibility)
    approach_body_xy.gif                # Animated body x-y trajectory
    approach_lvlh_3d.gif                # Animated 3D LVLH trajectory
    approach_lvlh_xy_r0_50m.gif         # Animated LVLH x-y with 50m circle
```

## How to Rerun

```bash
# Full sweep (90 cases x 10 ICs)
python IAC/run_sweep_mpc_90cases.py

# Validation only (1 case, 2 ICs)
python IAC/run_sweep_mpc_90cases.py --validate
```

## Interpreting Results

### summary.csv columns

| Column                       | Description                                |
|------------------------------|--------------------------------------------|
| `success_count`              | Number of ICs (out of 10) that docked       |
| `min_range`                  | Minimum range achieved across all ICs       |
| `max_lateral_excursion`      | Maximum |x_B| across all ICs                |
| `total_dv_mean/std`          | Mean/std of total delta-V across 10 ICs     |
| `constraint_violations_total`| Total LOS violations across all ICs         |
| `qp_failures_total`          | Total QP solver failures across all ICs     |
| `r_final_mean/std`           | Mean/std of final range across 10 ICs       |

### Key physics

- **r_sync = 2*a_max / omega^2:** Synchronisation radius. If y_limit > r_sync,
  the chaser cannot sustain co-rotation and will likely fail to dock.
- **Cost A vs B:** Cost B penalises velocity errors, producing smoother
  trajectories but potentially more fuel. Cost A allows the velocity to
  evolve freely as long as position converges.
- **QP failures:** OSQP may hit max iterations for difficult cases (large
  range, low thrust). The fallback braking controller maintains safety.

## Key Assumptions

1. Planar analysis (z = 0 throughout)
2. Target tumbles about z-axis only (simple spin)
3. No sensor noise, no navigation errors
4. No thruster dynamics or minimum impulse bit
5. HCW validity assumed (circular reference orbit, small separations)
