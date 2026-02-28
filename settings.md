# Simulation Settings Reference

All parameters are defined in `params.m`. This file documents their meaning and default values.

## Earth Constants
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `mu` | 3.986004418e14 | m^3/s^2 | Gravitational parameter |
| `Re` | 6378137.0 | m | Mean equatorial radius |
| `J2` | 1.08263e-3 | - | J2 zonal harmonic |

## Target Orbit (Circular LEO)
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `alt` | 500e3 | m | Altitude |
| `a` | Re + alt | m | Semi-major axis |
| `n` | sqrt(mu/a^3) | rad/s | Mean motion |

## Initial Conditions
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `dr_lvlh0` | [0; 200; 0] | m | Chaser initial relative position (LVLH) |
| `dv_lvlh0` | [0; 0; 0] | m/s | Chaser initial relative velocity (LVLH) |
| `omega_body` | [0; 0; 0.015] | rad/s | Target body angular velocity (~0.86 deg/s about z) |

## Simulation Timing
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `Tsim` | 400 | s | Total simulation time |
| `dt` | 1.0 | s | Control time step |

## MPC Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `Np` | 40 | Prediction horizon (steps) |
| `nx` | 6 | State dimension [r_TB; v_TB] |
| `nu` | 3 | Input dimension (accel in TB) |
| `u_max` | 0.15 | Max accel per axis [m/s^2] |

## Cost Weights
| Parameter | Value | Description |
|-----------|-------|-------------|
| `Q` | diag([15, 1, 15, 1, 1, 1]) | Stage state cost (heavy on x_TB, z_TB) |
| `QN` | 30 * Q | Terminal state cost |
| `Rdu` | 300 * I_3 | Delta-u smoothness penalty |
| `Ru` | 0.01 * I_3 | Input magnitude regularization |

## Reference / Approach Strategy
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `y_hold_start` | 200 | m | Initial hold distance on +yT |
| `y_hold_end` | 5 | m | Final approach distance |
| `y_hold_tau` | 200 | s | Exponential decay time constant |

## LOS Tetrahedral Corridor
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `cone_half_angle_deg` | 30 | deg | Half-angle of LOS cone |
| `cone_k` | tan(30 deg) | - | Cone slope parameter |
| `cone_nfaces` | 8 | - | Number of polyhedral faces |
| `y_min` | 1.0 | m | Corridor floor (minimum y_T) |

## Integrator Tolerances
| Parameter | Value | Description |
|-----------|-------|-------------|
| `ode_opts.RelTol` | 1e-10 | Relative tolerance for ode113 |
| `ode_opts.AbsTol` | 1e-12 | Absolute tolerance for ode113 |

## Visualization
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| `cube_size` | 3 | m | Half-edge of drawn cubes |
| `triad_len` | 20 | m | Triad axis length |
| `cone_draw_L` | 0 (auto) | m | Cone draw length (0 = set from initial y_TB) |

## OSQP Solver Settings
| Parameter | Value | Description |
|-----------|-------|-------------|
| `osqp_max_iter` | 4000 | Maximum OSQP iterations |
| `osqp_eps_abs` | 1e-5 | Absolute tolerance |
| `osqp_eps_rel` | 1e-5 | Relative tolerance |
| `osqp_warm_start` | true | Enable warm-starting |
| `osqp_verbose` | false | Suppress solver output |
