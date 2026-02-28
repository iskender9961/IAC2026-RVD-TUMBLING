"""Configuration for the IAC 2026 reachability-aware guidance demo."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict

import numpy as np


@dataclass
class ScenarioConfig:
    """Global scenario configuration.

    Units are SI unless otherwise noted.
    """

    # ── Run mode ──────────────────────────────────────────────────────────
    # "dev" = fast iteration with reduced outputs
    # "full" = publication-ready sweeps, denser grids, longer horizon
    mode: str = "dev"

    # ── Timing ────────────────────────────────────────────────────────────
    ts: float = 0.1          # controller/output sample time (s)
    truth_substeps: int = 3
    sim_duration_s: float = 300.0  # 5 minutes

    # ── Quasi-inertial relative model in LVLH ─────────────────────────────
    mean_motion_rad_s: float = 1.1e-3

    # ── Initial / terminal geometry ───────────────────────────────────────
    initial_range_m: float = 195.0
    hold_radius_m: float = 0.50
    hold_band_min_m: float = 0.45
    hold_band_max_m: float = 0.55
    keepout_radius_m: float = 0.45

    # ── Safety shaping ────────────────────────────────────────────────────
    radial_decay_eps_m: float = 0.02
    max_rel_speed_mps: float = 20.0
    hold_switch_speed_mps: float = 0.08
    hold_confirm_time_s: float = 6.0

    # ── Target tumble (body frame wrt LVLH) ───────────────────────────────
    tumble_omega_rad_s: float = 0.03   # ~1.72 deg/s; feasible at 195m with a_max=0.1

    # ── LOS corridor parameters (body frame) ─────────────────────────────
    # Wider cone: cx=1.5 gives ~56° half-angle (rise/run = 1.5)
    los_x0_m: float = 2.5
    los_z0_m: float = 2.5
    los_cx: float = 1.5
    los_cz: float = 1.5

    # ── Control limits (N/kg == m/s^2) ────────────────────────────────────
    thrust_to_mass_ratio_mps2: float = 0.10

    # ── MPC-style tracker / safety filter ─────────────────────────────────
    mpc_horizon_steps: int = 6
    tracker_kp: float = 0.04
    tracker_kd: float = 0.30
    hold_kp: float = 0.50
    hold_kd: float = 1.80

    # Body-axis tracking emphasis (approach mode)
    # Drive x_B, z_B to centerline fast; keep y_B progression slower.
    approach_pos_weight_xz: float = 4.0
    approach_pos_weight_y: float = 1.0
    approach_vel_weight_xz: float = 4.0
    approach_vel_weight_y: float = 1.0

    # ── Input-rate penalty (mandatory for smooth control) ─────────────────
    # Penalises du_k = u_k - u_{k-1} in the QP objective.
    input_rate_weight: float = 2.0       # weight on ||du||^2
    input_magnitude_weight: float = 0.1  # weight on ||u||^2 in QP

    use_centerline_reference: bool = True
    align_lateral_tol_m: float = 0.15
    align_y_weight_scale: float = 0.0
    y_reference_blend_after_align: float = 0.50
    enforce_x_body_monotone: bool = True
    x_body_monotone_eps_m: float = 0.01
    lateral_brake_band_m: float = 2.5
    lateral_brake_gain: float = 5.0
    safety_iterations: int = 10
    safety_margin_m: float = 1e-3
    skip_los_projection: bool = False   # disable hard LOS projection (for reachability)

    # ── Spiral reference shaping ──────────────────────────────────────────
    spiral_rate_rad_s: float = 0.12
    spiral_approach_time_s: float = 270.0
    lateral_ratio: float = 0.28
    lateral_decay_power: float = 1.35
    lateral_ratio_z: float = 0.0
    straight_line_body_approach: bool = True
    lock_z_axis: bool = True

    # ── Plot / animation controls ─────────────────────────────────────────
    corridor_far_y_m: float = 220.0
    animation_stride: int = 20
    animation_fps: int = 15

    # ── Reachability sweep (body-frame x-y feasibility maps) ──────────────
    reachability_tumble_deg_s: tuple[float, ...] = (1.0, 2.0, 3.0, 4.0, 5.0)
    reachability_accel_mps2: tuple[float, ...] = (0.20, 0.10, 0.05, 0.01)
    reachability_y_min_m: float = 0.50
    reachability_y_max_m: float = 150.0
    reachability_nx: int = 200
    reachability_ny: int = 120

    # ── Reproducibility ───────────────────────────────────────────────────
    seed: int = 2026

    # ── Output paths ──────────────────────────────────────────────────────
    root_dir: Path = Path("IAC")
    data_dir: Path = Path("IAC/data")
    fig_dir: Path = Path("IAC/figures")
    anim_dir: Path = Path("IAC/animations")
    frames_dir: Path = Path("IAC/animations/frames")

    def ensure_output_dirs(self) -> None:
        for folder in (self.root_dir, self.data_dir, self.fig_dir,
                       self.anim_dir, self.frames_dir):
            folder.mkdir(parents=True, exist_ok=True)
        # Reachability data directory
        (self.data_dir / "reachability").mkdir(parents=True, exist_ok=True)

    def to_dict(self) -> Dict:
        out = asdict(self)
        for key in ("root_dir", "data_dir", "fig_dir", "anim_dir", "frames_dir"):
            out[key] = str(out[key])
        return out

    @property
    def n_steps(self) -> int:
        return int(np.floor(self.sim_duration_s / self.ts))

    @property
    def scenario_tag(self) -> str:
        """Compact string for plot titles: tumble rate + thrust authority."""
        omega_deg = np.degrees(self.tumble_omega_rad_s)
        return (
            f"$\\omega_t$={omega_deg:.2f}°/s, "
            f"$a_{{max}}$={self.thrust_to_mass_ratio_mps2:.2f} m/s²"
        )


def get_config(mode: str = "dev") -> ScenarioConfig:
    """Return runtime configuration.

    Both modes preserve ts=0.1 s; full mode extends runtime/horizon and
    densifies reachability grids.
    """
    cfg = ScenarioConfig()
    mode_key = mode.lower().strip()
    if mode_key not in {"dev", "full"}:
        raise ValueError(f"Unknown mode '{mode}'. Expected 'dev' or 'full'.")

    cfg.mode = mode_key

    if mode_key == "dev":
        cfg.sim_duration_s = 300.0
        cfg.mpc_horizon_steps = 6
    else:
        cfg.sim_duration_s = 360.0
        cfg.mpc_horizon_steps = 8
        cfg.spiral_approach_time_s = 320.0
        cfg.reachability_nx = 300
        cfg.reachability_ny = 200
        cfg.reachability_y_max_m = 200.0
    return cfg
