"""Run approach from multiple initial ranges with default tumble/thrust.

Usage:
    python IAC/run_custom.py
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
from IAC.constraints import LOS3D
from IAC.controller_mpc import SpiralMPCController
from IAC.dynamics import simulate_closed_loop
from IAC.plots import plot_body_position_time, plot_mpc_cost_breakdown
from IAC.animate import create_body_plane_animation, create_lvlh_plane_animation
from IAC.reference_spiral import generate_spiral_reference


def run_one_range(initial_range_m: float) -> None:
    """Run default scenario from a given initial range and produce 4 outputs."""
    cfg = get_config("dev")
    cfg.initial_range_m = initial_range_m
    # Default tumble and thrust (1.72 deg/s, 0.10 m/s^2)
    np.random.seed(cfg.seed)
    cfg.ensure_output_dirs()

    tag = f"r0={initial_range_m:.0f}m"
    print(f"\n{'=' * 70}")
    print(f"CASE: {tag}  ({cfg.scenario_tag})")
    print(f"{'=' * 70}")

    omega_rad = cfg.tumble_omega_rad_s
    a_max = cfg.thrust_to_mass_ratio_mps2
    r_sync = 2.0 * a_max / max(omega_rad**2, 1e-12)
    print(f"  r_sync = {r_sync:.1f} m")
    print(f"  initial_range = {initial_range_m:.0f} m")

    t_grid = np.arange(0.0, cfg.sim_duration_s + cfg.ts, cfg.ts)
    reference = generate_spiral_reference(t_grid=t_grid, cfg=cfg)
    r0 = reference["r_lvlh"][0]
    x0 = np.hstack((r0, np.zeros(3)))

    los = LOS3D.from_config(cfg)
    controller = SpiralMPCController(cfg=cfg, los=los, ref=reference)

    print("  Running closed-loop simulation...")
    sim = simulate_closed_loop(
        cfg=cfg,
        initial_state_lvlh=x0,
        controller=controller,
        los=los,
        reference=reference,
    )

    hold_str = "YES" if sim["hold_success"] else "NO"
    print(f"  hold_success     = {hold_str}")
    print(f"  final range      = {sim['range_m'][-1]:.3f} m")
    print(f"  min range        = {float(np.min(sim['range_m'])):.3f} m")
    print(f"  final speed      = {sim['speed_mps'][-1]:.4f} m/s")
    print(f"  delta-v proxy    = {sim['delta_v_proxy_mps']:.3f} m/s")
    print(f"  min LOS margin   = {sim['min_constraint_margin_m']:.4f} m")
    print(f"  LOS violations   = {sim.get('los_violations', 'N/A')}")

    suffix = f"_r0_{int(initial_range_m)}m"

    print("  Generating plots...")
    p1 = plot_body_position_time(cfg, sim, cfg.fig_dir / f"body_position_vs_time{suffix}.png")
    print(f"    {p1}")
    p2 = plot_mpc_cost_breakdown(cfg, sim, cfg.fig_dir / f"mpc_cost_breakdown{suffix}.png")
    print(f"    {p2}")

    print("  Generating animations...")
    a1 = create_body_plane_animation(cfg, sim, los, plane="xy",
                                     stem=f"approach_body_xy{suffix}")
    print(f"    approach_body_xy: {a1.get('gif', a1.get('frames_dir'))}")
    a2 = create_lvlh_plane_animation(cfg, sim, los, plane="xy",
                                     stem=f"approach_lvlh_xy{suffix}")
    print(f"    approach_lvlh_xy: {a2.get('gif', a2.get('frames_dir'))}")


def main() -> None:
    for r0 in [150.0, 100.0, 50.0]:
        run_one_range(r0)

    print(f"\n{'=' * 70}")
    print("ALL CASES COMPLETE")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
