"""Phase 6-7 diagnostics: rerun suspicious case after fixes.

Validates that the CWH truth model + no-blending dynamics produce
physically honest results.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from IAC.config import get_config
from IAC.constraints import LOS3D, body_to_lvlh, lvlh_to_body
from IAC.controller_mpc import SpiralMPCController
from IAC.dynamics import simulate_closed_loop
from IAC.reference_spiral import generate_spiral_reference


def run_case(omega_deg: float, a_max_val: float, label: str) -> dict:
    """Run one scenario and print diagnostics."""
    cfg = get_config("dev")
    cfg.tumble_omega_rad_s = np.radians(omega_deg)
    cfg.thrust_to_mass_ratio_mps2 = a_max_val
    np.random.seed(cfg.seed)
    cfg.ensure_output_dirs()

    omega_rad = cfg.tumble_omega_rad_s
    a_max = cfg.thrust_to_mass_ratio_mps2

    print(f"\n{'=' * 76}")
    print(f"CASE: {label}")
    print(f"  omega_t = {omega_deg:.1f} deg/s = {omega_rad:.6f} rad/s")
    print(f"  a_max   = {a_max:.4f} m/s^2")
    print(f"  n_orbit = {cfg.mean_motion_rad_s:.6e} rad/s")
    r_sync = 2.0 * a_max / omega_rad**2
    print(f"  r_sync  = 2*a_max/omega^2 = {r_sync:.3f} m")
    print(f"  r_init  = {cfg.initial_range_m:.1f} m")
    if cfg.initial_range_m > r_sync:
        print(f"  >>> EXPECTED: infeasible (r_init >> r_sync)")
    print(f"{'=' * 76}")

    t_grid = np.arange(0.0, cfg.sim_duration_s + cfg.ts, cfg.ts)
    reference = generate_spiral_reference(t_grid=t_grid, cfg=cfg)
    r0 = reference["r_lvlh"][0]
    x0 = np.hstack((r0, np.zeros(3)))

    los = LOS3D.from_config(cfg)
    controller = SpiralMPCController(cfg=cfg, los=los, ref=reference)

    sim = simulate_closed_loop(
        cfg=cfg, initial_state_lvlh=x0, controller=controller,
        los=los, reference=reference,
    )

    # ── Summary ──
    controls = sim["control"]
    states = sim["state_lvlh"]
    t_arr = sim["time"]
    body_states = sim["state_body"]

    ctrl_norms = np.linalg.norm(controls, axis=1)
    dv_total = float(np.sum(ctrl_norms) * cfg.ts)

    # Velocity change validation
    velocities = states[:, 3:]
    dv_actual = np.diff(velocities, axis=0)
    dv_actual_norms = np.linalg.norm(dv_actual, axis=1)
    dv_from_thrust = ctrl_norms * cfg.ts

    # For CWH, gravity gradient also changes velocity, so we allow
    # dv_actual to exceed dv_thrust by a CWH-consistent amount.
    # HCW max velocity change from gravity: ~3n^2*r_max*dt for radial
    n_orb = cfg.mean_motion_rad_s
    r_max = float(np.max(sim["range_m"]))
    dv_cwh_bound = 3.0 * n_orb**2 * r_max * cfg.ts  # max HCW gravity dv per step

    dv_excess = dv_actual_norms - dv_from_thrust - dv_cwh_bound
    n_phys_violations = int(np.sum(dv_excess > 0.001))

    print(f"\n  DYNAMICS INTEGRITY:")
    print(f"    max ||u||                  = {np.max(ctrl_norms):.6f} m/s^2")
    print(f"    mean ||u||                 = {np.mean(ctrl_norms):.6f} m/s^2")
    print(f"    total delta-v              = {dv_total:.3f} m/s")
    print(f"    max achievable delta-v     = {a_max * cfg.sim_duration_s:.3f} m/s")
    print(f"    CWH gravity dv bound/step  = {dv_cwh_bound:.6f} m/s")
    print(f"    physics violations (>1mm)  = {n_phys_violations} / {len(dv_excess)}")
    print(f"    dv_excess_total (reported) = {sim.get('dv_excess_total', 'N/A')}")

    print(f"\n  RESULT:")
    print(f"    final range                = {sim['range_m'][-1]:.4f} m")
    print(f"    final speed                = {sim['speed_mps'][-1]:.4f} m/s")
    print(f"    hold_success               = {sim['hold_success']}")
    print(f"    hold_mode_enter_time       = {sim['hold_mode_enter_time_s']}")
    print(f"    stays_feasible             = {sim['stays_feasible']}")
    print(f"    min LOS margin             = {sim['min_constraint_margin_m']:.4f} m")
    print(f"    keepout_violations         = {sim.get('keepout_violations', 'N/A')}")
    print(f"    los_violations             = {sim.get('los_violations', 'N/A')}")

    # LOS violations detail
    los_margins = sim["los_margin"]
    n_los_neg = int(np.sum(los_margins < -1e-3))
    if n_los_neg > 0:
        print(f"    LOS margin < -1e-3 count   = {n_los_neg}")
        print(f"    worst LOS margin           = {np.min(los_margins):.4f} m")

    # Range evolution
    ranges = sim["range_m"]
    print(f"\n  RANGE EVOLUTION:")
    print(f"    initial                    = {ranges[0]:.3f} m")
    print(f"    min                        = {np.min(ranges):.3f} m")
    print(f"    max                        = {np.max(ranges):.3f} m")
    print(f"    final                      = {ranges[-1]:.3f} m")
    for pct in [25, 50, 75]:
        idx = int(len(ranges) * pct / 100)
        print(f"    at {pct}% of sim (t={t_arr[idx]:.0f}s) = {ranges[idx]:.3f} m")

    return sim


def main() -> None:
    print("PHASE 7: POST-FIX VALIDATION")
    print("=" * 76)
    print("Changes applied:")
    print("  1. Truth dynamics: double integrator -> CWH (HCW)")
    print("  2. Reference blending: REMOVED (approach + hold)")
    print("  3. Safety projections: REMOVED (diagnostic only)")
    print("  4. Velocity/position overrides: REMOVED")
    print("=" * 76)

    # ── Case 1: Default scenario (should still work) ──
    run_case(1.72, 0.10, "DEFAULT: 1.72 deg/s, 0.10 m/s^2")

    # ── Case 2: Moderate challenge ──
    run_case(3.0, 0.05, "MODERATE: 3 deg/s, 0.05 m/s^2")

    # ── Case 3: The suspicious extreme case ──
    run_case(10.0, 0.02, "EXTREME: 10 deg/s, 0.02 m/s^2 (SHOULD FAIL)")

    # ── Case 4: Borderline ──
    run_case(5.0, 0.05, "BORDERLINE: 5 deg/s, 0.05 m/s^2")

    print(f"\n{'=' * 76}")
    print("DIAGNOSTICS COMPLETE")
    print("=" * 76)


if __name__ == "__main__":
    main()
