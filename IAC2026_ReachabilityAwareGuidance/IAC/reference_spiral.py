"""Spiral reference generation in LVLH and target body frames."""

from __future__ import annotations

from typing import Dict

import numpy as np

from IAC.config import ScenarioConfig
from IAC.constraints import body_to_lvlh


def _clamp01(x: np.ndarray) -> np.ndarray:
    return np.clip(x, 0.0, 1.0)


def generate_spiral_reference(
    t_grid: np.ndarray,
    cfg: ScenarioConfig,
    r_start: float | None = None,
    r_hold: float | None = None,
    omega_spiral: float | None = None,
) -> Dict[str, np.ndarray]:
    """Generate a tightening approach reference in the target body frame.

    By default this is planar (z_B locked to zero). 3D swirl can be enabled
    by setting cfg.lateral_ratio_z > 0 and cfg.lock_z_axis = False. If
    cfg.straight_line_body_approach is True, approach is centerline only:
    x_B=z_B=0 and only y_B decreases toward hold radius.
    """

    t = np.asarray(t_grid, dtype=float)
    if t.ndim != 1 or t.size < 3:
        raise ValueError("t_grid must be a 1D array with at least 3 samples.")

    r0 = float(cfg.initial_range_m if r_start is None else r_start)
    if r0 > 200.0 + 1e-9:
        raise ValueError("r_start must be <= 200 m for this scenario.")
    rh = float(cfg.hold_radius_m if r_hold is None else r_hold)
    if rh <= 0.0:
        raise ValueError("r_hold must be positive.")

    w_sp = float(cfg.spiral_rate_rad_s if omega_spiral is None else omega_spiral)
    t_app = float(max(cfg.spiral_approach_time_s, cfg.ts))

    s = _clamp01(t / t_app)
    rho = rh + (r0 - rh) * (1.0 - s) ** 2

    swirl_weight = (1.0 - s) ** cfg.lateral_decay_power
    lateral_mag_x = cfg.lateral_ratio * rho * swirl_weight
    lateral_mag_z = cfg.lateral_ratio_z * rho * swirl_weight

    psi = w_sp * t

    x_b = lateral_mag_x * np.cos(psi)
    z_b = lateral_mag_z * np.sin(psi)
    y_sq = np.maximum(rho**2 - x_b**2 - z_b**2, rh**2)
    y_b = np.sqrt(y_sq)

    if cfg.straight_line_body_approach:
        x_b[:] = 0.0
        z_b[:] = 0.0
        y_b = rho.copy()

    hold_mask = s >= 1.0
    x_b[hold_mask] = 0.0
    z_b[hold_mask] = 0.0
    y_b[hold_mask] = rh

    r_body = np.column_stack((x_b, y_b, z_b))
    v_body = np.gradient(r_body, t, axis=0)
    a_body = np.gradient(v_body, t, axis=0)

    if cfg.lock_z_axis:
        r_body[:, 2] = 0.0
        v_body[:, 2] = 0.0
        a_body[:, 2] = 0.0

    r_lvlh = np.zeros_like(r_body)
    v_lvlh = np.zeros_like(v_body)
    a_lvlh = np.zeros_like(a_body)

    for i, ti in enumerate(t):
        r_lvlh[i], v_lvlh[i] = body_to_lvlh(
            r_body[i],
            v_body[i],
            theta=cfg.tumble_omega_rad_s * float(ti),
            omega_z=cfg.tumble_omega_rad_s,
        )
        _, a_lvlh[i] = body_to_lvlh(
            np.zeros(3),
            a_body[i],
            theta=cfg.tumble_omega_rad_s * float(ti),
            omega_z=0.0,
        )

    if cfg.lock_z_axis:
        r_lvlh[:, 2] = 0.0
        v_lvlh[:, 2] = 0.0
        a_lvlh[:, 2] = 0.0

    speed = np.linalg.norm(v_lvlh, axis=1)
    vmax = max(cfg.max_rel_speed_mps, 1e-6)
    if np.max(speed) > vmax:
        scale = vmax / np.max(speed)
        v_lvlh *= scale
        v_body *= scale
        a_lvlh *= scale
        a_body *= scale

    return {
        "time": t,
        "r_body": r_body,
        "v_body": v_body,
        "a_body": a_body,
        "r_lvlh": r_lvlh,
        "v_lvlh": v_lvlh,
        "a_lvlh": a_lvlh,
        "rho": np.linalg.norm(r_body, axis=1),
        "psi": psi,
    }
