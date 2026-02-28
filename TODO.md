# TODO

## High Priority
- [ ] Validate monotonic cost decay after sync phase (tune Q, Rdu, Ru)
- [ ] Add chaser-body-frame attitude dynamics (currently assumed aligned with TB)
- [ ] Add terminal constraint set for guaranteed docking safety
- [ ] Test multi-axis tumble scenarios (omega_x, omega_y nonzero)

## Medium Priority
- [ ] Add fuel-optimal or minimum-time fallback cost modes
- [ ] Implement variable prediction horizon (shrink Np near terminal)
- [ ] Add sensor noise / navigation filter (EKF or UKF)
- [ ] LVLH-based constraint fallback when TB frame loses observability
- [ ] Profile OSQP solve time vs horizon length; consider warm-start tuning
- [ ] Add approach waypoints (hold-points at intermediate y_T distances)

## Low Priority / Future Work
- [ ] GPU-accelerated batch Monte Carlo for robustness analysis
- [ ] ROS2 / Simulink integration for hardware-in-the-loop testing
- [ ] Add solar radiation pressure perturbation
- [ ] Extend to elliptical target orbits
- [ ] Add collision avoidance (keep-out zone around target body)
- [ ] Compare polyhedral vs SOCP cone formulation (fidelity vs speed)

## Completed
- [x] Replace CVX with OSQP (sparse QP with warm-start)
- [x] Replace HCW plant with nonlinear ECI 2-body + J2 dynamics
- [x] Body-fixed LOS tetrahedral corridor on tumbling target (+yT axis)
- [x] LTV MPC in target body frame with finite-difference linearization
- [x] Delta-u smoothness penalty as primary cost
- [x] Comprehensive plotting (18 PNGs) and GIF generation (5 GIFs)
- [x] Sanity test suite (5 tests, all passing)
- [x] LaTeX documentation
