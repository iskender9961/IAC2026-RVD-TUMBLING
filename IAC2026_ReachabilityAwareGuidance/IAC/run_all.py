"""Run the upgraded IAC 2026 spiral-approach demo.

Usage:
    python IAC/run_all.py
    python IAC/run_all.py --mode full
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from IAC.animate import create_animation_suite
from IAC.config import get_config
from IAC.constraints import LOS3D
from IAC.controller_mpc import SpiralMPCController
from IAC.dynamics import simulate_closed_loop
from IAC.plots import create_all_plots
from IAC.reachability import (
    build_reachability_latex_table,
    format_reachability_table_terminal,
    run_reachability_study,
)
from IAC.reference_spiral import generate_spiral_reference


def _initial_state_from_reference(ref: Dict[str, np.ndarray]) -> np.ndarray:
    r0 = ref["r_lvlh"][0]
    v0 = np.zeros(3)
    return np.hstack((r0, v0))


def _format_performance_table(summary: Dict) -> str:
    rows = [
        ("Initial range", f"{summary['initial_range_m']:.3f}", "m"),
        ("Final range", f"{summary['final_range_m']:.3f}", "m"),
        ("Final speed", f"{summary['final_speed_mps']:.4f}", "m/s"),
        ("Min LOS margin", f"{summary['min_constraint_margin_m']:.5f}", "m"),
        ("Delta-v proxy", f"{summary['total_delta_v_proxy_mps']:.3f}", "m/s"),
        ("Hold entered at", "N/A" if summary["hold_mode_enter_time_s"] is None else f"{summary['hold_mode_enter_time_s']:.2f}", "s"),
        ("Hold success", "Yes" if summary["hold_success"] else "No", "-"),
        ("x_B overshoot", "Yes" if summary["x_body_overshoot"] else "No", "-"),
        ("max |x_B|", f"{summary['max_abs_x_body_m']:.3f}", "m"),
        ("max |z_B|", f"{summary['max_abs_z_body_m']:.3f}", "m"),
    ]
    w0 = max(len(r[0]) for r in rows) + 2
    w1 = max(len(r[1]) for r in rows) + 2
    header = f"{'Metric'.ljust(w0)}{'Value'.ljust(w1)}Unit"
    sep = "-" * (w0 + w1 + 4)
    body = [f"{k.ljust(w0)}{v.ljust(w1)}{u}" for k, v, u in rows]
    return "\n".join([header, sep, *body])


def _build_performance_latex_table(summary: Dict) -> str:
    hold_t = "N/A" if summary["hold_mode_enter_time_s"] is None else f"{summary['hold_mode_enter_time_s']:.2f} s"
    return "\n".join(
        [
            r"\begin{table}[H]",
            r"\centering",
            r"\caption{Single-run performance summary for the nominal case.}",
            r"\label{tab:performance_summary}",
            r"\begin{tabular}{l c}",
            r"\hline",
            r"Metric & Value \\",
            r"\hline",
            f"Initial range & {summary['initial_range_m']:.3f} m \\\\",
            f"Final range & {summary['final_range_m']:.3f} m \\\\",
            f"Final speed & {summary['final_speed_mps']:.4f} m/s \\\\",
            f"Minimum LOS margin & {summary['min_constraint_margin_m']:.5f} m \\\\",
            f"Total $\\Delta v$ proxy & {summary['total_delta_v_proxy_mps']:.3f} m/s \\\\",
            f"Hold switch time & {hold_t} \\\\",
            f"Hold success & {'Yes' if summary['hold_success'] else 'No'} \\\\",
            f"$x_B$ overshoot observed & {'Yes' if summary['x_body_overshoot'] else 'No'} \\\\",
            f"Max $|x_B|$ & {summary['max_abs_x_body_m']:.3f} m \\\\",
            f"Max $|z_B|$ & {summary['max_abs_z_body_m']:.3f} m \\\\",
            r"\hline",
            r"\end{tabular}",
            r"\end{table}",
        ]
    )


def run_pipeline(mode: str = "dev", verbose: bool = True, include_reachability: bool = True) -> Dict:
    cfg = get_config(mode)
    np.random.seed(cfg.seed)
    cfg.ensure_output_dirs()

    if verbose:
        print("=" * 76)
        print("IAC 2026 Reachability-Aware Guidance Demo (Spiral + Hold)")
        print("=" * 76)
        print(f"Mode: {mode}")
        print(f"Sampling time ts: {cfg.ts:.3f} s")
        print(f"Initial range target: <= 200 m (configured {cfg.initial_range_m:.1f} m)")
        print(f"Hold radius: {cfg.hold_radius_m:.2f} m")
        print(f"Thrust-to-mass bound: {cfg.thrust_to_mass_ratio_mps2:.3f} m/s^2")
        print()

    t_start = time.perf_counter()

    t_grid = np.arange(0.0, cfg.sim_duration_s + cfg.ts, cfg.ts)
    reference = generate_spiral_reference(t_grid=t_grid, cfg=cfg)
    x0 = _initial_state_from_reference(reference)

    if np.linalg.norm(x0[:3]) > 200.0 + 1e-9:
        raise RuntimeError("Initial range exceeds 200 m requirement.")

    los = LOS3D.from_config(cfg)
    controller = SpiralMPCController(cfg=cfg, los=los, ref=reference)

    sim = simulate_closed_loop(
        cfg=cfg,
        initial_state_lvlh=x0,
        controller=controller,
        los=los,
        reference=reference,
    )

    plot_paths = create_all_plots(cfg, sim, los)
    anim_paths = create_animation_suite(cfg, sim, los, stem_prefix="approach")

    reachability = None
    if include_reachability:
        if verbose:
            print()
            print("Running reachability sweep (body-frame x-y maps)...")
        reachability = run_reachability_study(cfg, verbose=verbose)

    runtime = time.perf_counter() - t_start

    monotonic_ok = bool(not np.any(sim["monotonic_violation"]))
    x_body = np.asarray(sim["state_body"][:, 0], dtype=float)
    z_body = np.asarray(sim["state_body"][:, 2], dtype=float)
    x_overshoot = bool(
        np.any(
            np.logical_and(
                np.abs(x_body[:-1]) <= cfg.lateral_brake_band_m,
                x_body[:-1] * x_body[1:] < -1e-6,
            )
        )
    )
    summary = {
        "config": cfg.to_dict(),
        "sampling_time_s": cfg.ts,
        "thrust_to_mass_used_mps2": cfg.thrust_to_mass_ratio_mps2,
        "initial_range_m": float(np.linalg.norm(x0[:3])),
        "feasibility": bool(sim["stays_feasible"]),
        "los_feasibility": bool(np.all(sim["los_margin"] >= -1e-6)),
        "monotonic_range_decrease_satisfied": monotonic_ok,
        "hold_success": bool(sim["hold_success"]),
        "hold_mode_enter_time_s": sim["hold_mode_enter_time_s"],
        "final_range_m": float(sim["range_m"][-1]),
        "final_speed_mps": float(sim["speed_mps"][-1]),
        "total_delta_v_proxy_mps": float(sim["delta_v_proxy_mps"]),
        "min_constraint_margin_m": float(sim["min_constraint_margin_m"]),
        "x_body_overshoot": x_overshoot,
        "max_abs_x_body_m": float(np.max(np.abs(x_body))),
        "max_abs_z_body_m": float(np.max(np.abs(z_body))),
        "runtime_s": float(runtime),
        "figures": plot_paths,
        "animations": anim_paths,
        "reachability": reachability,
    }

    perf_table_terminal = _format_performance_table(summary)
    perf_table_latex = _build_performance_latex_table(summary)
    perf_table_path = cfg.data_dir / "performance_summary_table.tex"
    perf_table_path.parent.mkdir(parents=True, exist_ok=True)
    with open(perf_table_path, "w", encoding="utf-8") as f:
        f.write(perf_table_latex)

    reach_table_terminal = None
    reach_table_latex = None
    reach_table_path = None
    if reachability is not None:
        reach_table_terminal = format_reachability_table_terminal(reachability)
        reach_table_latex = build_reachability_latex_table(reachability)
        reach_table_path = cfg.data_dir / "reachability_table.tex"
        with open(reach_table_path, "w", encoding="utf-8") as f:
            f.write(reach_table_latex)
        summary["reachability_table_tex"] = str(reach_table_path)
    summary["performance_table_tex"] = str(perf_table_path)

    summary_path = cfg.data_dir / "summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    np.savez(
        cfg.data_dir / "trajectory_demo.npz",
        time=sim["time"],
        state_lvlh=sim["state_lvlh"],
        state_body=sim["state_body"],
        control=sim["control"],
        mpc_cost_time=sim["mpc_cost_time"],
        mpc_cost=sim["mpc_cost"],
        cost_breakdown=sim["cost_breakdown"],
        los_margin=sim["los_margin"],
        range_m=sim["range_m"],
        speed_mps=sim["speed_mps"],
        phase_error_m=sim["phase_error_m"],
        reference_r_lvlh=reference["r_lvlh"],
        reference_v_lvlh=reference["v_lvlh"],
        reference_r_body=reference["r_body"],
        reference_v_body=reference["v_body"],
    )

    if verbose:
        print("Run completed.")
        print(f"Summary JSON: {summary_path}")
        for key, path in plot_paths.items():
            print(f"Figure [{key}]: {path}")

        for anim_name, anim_info in anim_paths.items():
            if anim_info.get("gif"):
                print(f"Animation [{anim_name}] (gif): {anim_info['gif']}")
            elif anim_info.get("mp4"):
                print(f"Animation [{anim_name}] (mp4): {anim_info['mp4']}")
            else:
                print(f"Animation [{anim_name}] frames: {anim_info.get('frames_dir')}")

        print()
        print("Performance Summary Table")
        print(perf_table_terminal)

        if reachability is not None and reach_table_terminal is not None:
            print()
            print("Reachability Feasible Fraction Table (inside cone)")
            print(reach_table_terminal)
            for key, path in reachability["figure_paths"].items():
                print(f"Reachability figure [{key}]: {path}")
            if reach_table_path is not None:
                print(f"LaTeX reachability table: {reach_table_path}")

        print(f"LaTeX performance table: {perf_table_path}")
        print(f"Data archive: {cfg.data_dir / 'trajectory_demo.npz'}")
        print("=" * 76)

    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run upgraded IAC spiral-approach demo.")
    parser.add_argument("--mode", type=str, default="dev", choices=["dev", "full"])
    parser.add_argument("--skip-reachability", action="store_true", help="Skip body-frame reachability sweep.")
    args = parser.parse_args()
    run_pipeline(mode=args.mode, verbose=True, include_reachability=not args.skip_reachability)


if __name__ == "__main__":
    main()
