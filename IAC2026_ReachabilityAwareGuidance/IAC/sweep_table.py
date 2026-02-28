"""Parameter sweep: omega_t vs a_max, with summary table + GIFs for default.

Sweeps omega_t = 1..5 deg/s and a_max = 0.10, 0.05, 0.02 m/s^2.
Produces:
  - Terminal + LaTeX table of results
  - approach_body_xy.gif, approach_lvlh_xy.gif for the default case
  - body_position_vs_time, mpc_cost_breakdown for the default case
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
from IAC.reference_spiral import generate_spiral_reference


def run_one(omega_deg: float, a_max_val: float) -> dict:
    """Run one case and return summary metrics."""
    cfg = get_config("dev")
    cfg.tumble_omega_rad_s = np.radians(omega_deg)
    cfg.thrust_to_mass_ratio_mps2 = a_max_val
    np.random.seed(cfg.seed)
    cfg.ensure_output_dirs()

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

    omega_rad = cfg.tumble_omega_rad_s
    r_sync = 2.0 * a_max_val / max(omega_rad**2, 1e-12)
    ctrl_norms = np.linalg.norm(sim["control"], axis=1)
    dv_total = float(np.sum(ctrl_norms) * cfg.ts)

    return {
        "omega_deg": omega_deg,
        "a_max": a_max_val,
        "r_sync": r_sync,
        "hold_success": sim["hold_success"],
        "final_range": sim["range_m"][-1],
        "min_range": float(np.min(sim["range_m"])),
        "final_speed": sim["speed_mps"][-1],
        "delta_v": dv_total,
        "min_los_margin": sim["min_constraint_margin_m"],
        "los_violations": sim.get("los_violations", -1),
        "stays_feasible": sim["stays_feasible"],
    }


def main() -> None:
    omega_values = [1.0, 2.0, 3.0, 4.0, 5.0]
    a_max_values = [0.10, 0.05, 0.02]

    print("=" * 90)
    print("PARAMETER SWEEP: omega_t x a_max")
    print("=" * 90)

    results = []
    for omega_deg in omega_values:
        for a_max_val in a_max_values:
            tag = f"omega={omega_deg:.0f}deg/s, a_max={a_max_val:.2f}"
            print(f"  Running {tag} ...", end=" ", flush=True)
            row = run_one(omega_deg, a_max_val)
            status = "HOLD" if row["hold_success"] else "FAIL"
            print(f"{status} | range=[{row['min_range']:.1f}, {row['final_range']:.1f}] m | dv={row['delta_v']:.1f} m/s")
            results.append(row)

    # ── Terminal table ──
    print("\n" + "=" * 90)
    print(f"{'omega':>6} {'a_max':>6} {'r_sync':>8} {'hold':>5} {'r_min':>7} {'r_final':>8} "
          f"{'v_final':>8} {'dv':>7} {'LOS_min':>8} {'LOS_viol':>9}")
    print("-" * 90)
    for r in results:
        hold_str = "YES" if r["hold_success"] else "NO"
        print(f"{r['omega_deg']:5.1f}° {r['a_max']:6.2f} {r['r_sync']:8.1f} {hold_str:>5} "
              f"{r['min_range']:7.1f} {r['final_range']:8.1f} "
              f"{r['final_speed']:8.3f} {r['delta_v']:7.1f} {r['min_los_margin']:8.2f} "
              f"{r['los_violations']:>9}")
    print("=" * 90)

    # ── LaTeX table ──
    latex_lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\caption{Parameter sweep: approach outcome vs.\ tumble rate and thrust authority.}",
        r"\label{tab:sweep}",
        r"\begin{tabular}{c c c c c c c c}",
        r"\hline",
        r"$\omega_t$ (°/s) & $a_{\max}$ (m/s$^2$) & $r_{\mathrm{sync}}$ (m) & Hold & "
        r"$r_{\min}$ (m) & $r_{\mathrm{final}}$ (m) & $\Delta v$ (m/s) & LOS viol. \\",
        r"\hline",
    ]
    for r in results:
        hold_str = "Yes" if r["hold_success"] else "No"
        latex_lines.append(
            f"{r['omega_deg']:.0f} & {r['a_max']:.2f} & {r['r_sync']:.1f} & {hold_str} & "
            f"{r['min_range']:.1f} & {r['final_range']:.1f} & {r['delta_v']:.1f} & "
            f"{r['los_violations']} \\\\"
        )
    latex_lines += [
        r"\hline",
        r"\end{tabular}",
        r"\end{table}",
    ]
    latex_str = "\n".join(latex_lines)

    cfg = get_config("dev")
    cfg.ensure_output_dirs()
    tex_path = cfg.data_dir / "sweep_table.tex"
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write(latex_str)
    print(f"\nLaTeX table saved to: {tex_path}")

    # ── Generate GIFs + plots for default case (1.72 deg/s, 0.10 m/s^2) ──
    print("\nGenerating GIFs + plots for default case ...")
    cfg = get_config("dev")
    np.random.seed(cfg.seed)
    cfg.ensure_output_dirs()

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

    from IAC.plots import plot_body_position_time, plot_mpc_cost_breakdown
    from IAC.animate import create_body_plane_animation, create_lvlh_plane_animation

    p1 = plot_body_position_time(cfg, sim, cfg.fig_dir / "body_position_vs_time.png")
    print(f"  {p1}")
    p2 = plot_mpc_cost_breakdown(cfg, sim, cfg.fig_dir / "mpc_cost_breakdown.png")
    print(f"  {p2}")

    a1 = create_body_plane_animation(cfg, sim, los, plane="xy", stem="approach_body_xy")
    print(f"  approach_body_xy: {a1.get('gif', a1.get('frames_dir'))}")
    a2 = create_lvlh_plane_animation(cfg, sim, los, plane="xy", stem="approach_lvlh_xy")
    print(f"  approach_lvlh_xy: {a2.get('gif', a2.get('frames_dir'))}")

    print("\nDone!")


if __name__ == "__main__":
    main()
