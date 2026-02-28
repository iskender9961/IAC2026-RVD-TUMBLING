# MPC Rendezvous with Tumbling Target -- Body-Fixed LOS Constraints

## Overview

Model Predictive Control (MPC) simulation for spacecraft rendezvous and docking with a tumbling target in LEO. The chaser approaches along the target's **+y body axis** (docking axis) while staying inside a **body-fixed line-of-sight (LOS) cone** that rotates with the target.

### Architecture

| Component | Description |
|-----------|-------------|
| **Plant (truth)** | Nonlinear ECI dynamics: two-body + J2, integrated with `ode113` |
| **MPC prediction** | LTV discrete model in target-body frame, re-linearized each step |
| **Solver** | OSQP (sparse QP, warm-started) |
| **LOS constraint** | Polyhedral cone approximation, body-fixed to tumbling target |

## Prerequisites

- MATLAB R2020b or later
- [OSQP for MATLAB](https://osqp.org/docs/get_started/matlab.html) (install via `osqpinstall` or add to path)

## Quick Start

```matlab
% From the project root directory:
main_sim        % Run the full simulation
run_tests       % Run sanity checks
```

## File Structure

```
.
├── main_sim.m                  % Entry point
├── params.m                    % All tunable parameters
├── run_tests.m                 % Sanity checks
├── dynamics/
│   ├── accel_2body.m           % Two-body gravity
│   ├── accel_J2.m              % J2 perturbation
│   ├── eom_eci.m               % ECI equations of motion
│   └── propagate_truth_step.m  % One-step nonlinear propagation
├── frames/
│   ├── lvlh_from_rv.m          % LVLH rotation from orbit state
│   ├── quat_propagate.m        % Quaternion integration
│   ├── rotm_from_quat.m        % Quaternion to rotation matrix
│   ├── target_attitude_model.m % Target tumble model
│   ├── transform_rel_to_TB.m   % ECI to target-body relative state
│   └── control_TB_to_ECI.m     % Control transform TB -> ECI
├── mpc/
│   ├── build_qp_osqp.m         % Build initial OSQP problem
│   ├── update_qp_osqp.m        % Update and solve QP each step
│   ├── linearize_discrete_model.m  % Finite-difference Ad, Bd
│   └── los_polyhedral_constraints.m % Polyhedral LOS cone
├── viz/
│   ├── animate_scene.m         % 3D animation with cubes/triads/cone
│   ├── draw_cube.m             % Draw oriented cube
│   ├── draw_triad.m            % Draw coordinate frame triad
│   ├── draw_cone_poly.m        % Draw polyhedral cone wireframe
│   └── plot_time_histories.m   % Static time-history plots
├── utils/
│   ├── skew.m                  % Skew-symmetric matrix
│   ├── clamp.m                 % Element-wise clamping
│   └── assert_units.m          % Order-of-magnitude sanity check
└── docs/
    ├── README.md               % This file
    └── mpc_rvd_bodyfixed_los.tex % LaTeX documentation
```

## Key Parameters (in `params.m`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p.omega_body` | `[0; 0; 0.02]` | Target body tumble rate [rad/s] |
| `p.cone_half_angle_deg` | `30` | LOS cone half-angle [deg] |
| `p.u_max` | `0.1` | Max thrust acceleration per axis [m/s^2] |
| `p.Np` | `30` | MPC prediction horizon [steps] |
| `p.dt` | `1.0` | Control step size [s] |
| `p.Tsim` | `300` | Total simulation time [s] |
| `p.Q` | `diag([1,1,1,0.1,0.1,0.1])` | State tracking weight |
| `p.Rdu` | `1e4*I_3` | Delta-u (smoothness) weight |
| `p.y_hold_start` | `200` | Initial hold distance on +yT [m] |
| `p.y_hold_end` | `5` | Final approach distance [m] |
| `p.y_hold_tau` | `150` | Exponential decay time constant [s] |

## Tuning Tips

- **Increase `Rdu`** for smoother control (less jitter, slower response).
- **Decrease `cone_half_angle_deg`** for tighter approach corridor (may cause infeasibility if too small).
- **Increase `omega_body`** for faster target tumble (more challenging scenario).
- **Increase `Np`** for better look-ahead but slower solve times.
- **Set `p.save_mp4 = true`** to save the animation as an MP4 file.

## Coordinate Frames

- **ECI**: Earth-Centered Inertial (truth propagation frame).
- **LVLH**: Local Vertical Local Horizontal. x = radial, y = along-track, z = orbit normal.
- **Target Body (TB)**: Starts aligned with LVLH at t=0, then rotates with the target's tumble rate. The docking axis is **+y_TB**. All MPC states and inputs are expressed in this frame.

## Termination Conditions

The simulation stops immediately if:
1. The chaser violates the body-fixed LOS cone at any plant step.
2. OSQP returns an infeasible status.

Both conditions print diagnostic information (time, state, margin).
