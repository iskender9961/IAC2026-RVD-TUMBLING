"""3D relative dynamics and closed-loop simulation.

Truth propagation uses linearised Clohessy-Wiltshire (HCW) equations in LVLH.
Control is applied as acceleration in the LVLH frame.  No reference blending
or position projection is performed — the state evolves purely from dynamics
plus the commanded acceleration, so feasibility claims are physically honest.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np

from IAC.config import ScenarioConfig
from IAC.constraints import LOS3D, body_to_lvlh, lvlh_to_body


# ─────────────────────────────────────────────────────────────────────────
# Discrete-time models
# ─────────────────────────────────────────────────────────────────────────

@dataclass
class DiscreteModel:
    a_d: np.ndarray
    b_d: np.ndarray


def discrete_cwh_3d(dt: float, n: float) -> DiscreteModel:
    """Discrete-time Clohessy-Wiltshire-Hill (CWH/HCW) state-transition.

    State: x = [x, y, z, vx, vy, vz]  in LVLH
        x  = radial (outward)
        y  = along-track
        z  = cross-track

    Continuous CWH:
        ẍ  = 3n²x  + 2nẏ  + a_x
        ÿ  =       - 2nẋ  + a_y
        z̈  = -n²z          + a_z

    Uses the exact matrix exponential (closed-form).
    """
    nt = n * dt
    c = float(np.cos(nt))
    s = float(np.sin(nt))

    # State transition Φ(dt)
    a_d = np.array([
        [4.0 - 3.0*c,    0.0, 0.0,  s/n,         2.0*(1.0 - c)/n,  0.0],
        [6.0*(s - nt),   1.0, 0.0, -2.0*(1.0-c)/n, (4.0*s - 3.0*nt)/n, 0.0],
        [0.0,            0.0, c,    0.0,           0.0,              s/n],
        [3.0*n*s,        0.0, 0.0,  c,             2.0*s,            0.0],
        [-6.0*n*(1.0-c), 0.0, 0.0, -2.0*s,         4.0*c - 3.0,      0.0],
        [0.0,            0.0, -n*s, 0.0,           0.0,              c],
    ], dtype=float)

    # Input matrix (zero-order hold on constant acceleration u over [0,dt]):
    # B = integral_0^dt Φ(τ) dτ * B_c  where B_c = [0;0;0;1;0;0;...] etc.
    # We integrate analytically for the exact ZOH input matrix.
    b_d = np.array([
        [(1.0 - c)/n**2,             2.0*(nt - s)/n**2,            0.0],
        [-2.0*(nt - s)/n**2,         (4.0*(1.0 - c) - 1.5*nt**2)/n**2, 0.0],
        [0.0,                        0.0,                          (1.0 - c)/n**2],
        [s/n,                        2.0*(1.0 - c)/n,              0.0],
        [-2.0*(1.0 - c)/n,           (4.0*s - 3.0*nt)/n,          0.0],
        [0.0,                        0.0,                          s/n],
    ], dtype=float)

    return DiscreteModel(a_d=a_d, b_d=b_d)


def discrete_double_integrator_3d(dt: float) -> DiscreteModel:
    """Discrete-time double integrator (no orbital coupling).

    Retained for comparison / ablation studies.
    """
    i3 = np.eye(3)
    z3 = np.zeros((3, 3))
    a_d = np.block([
        [i3, dt * i3],
        [z3, i3],
    ])
    b_d = np.block([
        [0.5 * dt * dt * i3],
        [dt * i3],
    ])
    return DiscreteModel(a_d=a_d, b_d=b_d)


# ─────────────────────────────────────────────────────────────────────────
# Hold reference
# ─────────────────────────────────────────────────────────────────────────

def _hold_reference_lvlh(cfg: ScenarioConfig, t: float) -> np.ndarray:
    """Station-keeping reference in LVLH for fixed body-frame docking direction."""
    r_b = np.array([0.0, cfg.hold_radius_m, 0.0])
    v_b = np.zeros(3)
    r_l, v_l = body_to_lvlh(r_b, v_b,
                             theta=cfg.tumble_omega_rad_s * t,
                             omega_z=cfg.tumble_omega_rad_s)
    return np.hstack((r_l, v_l))


# ─────────────────────────────────────────────────────────────────────────
# Closed-loop simulation
# ─────────────────────────────────────────────────────────────────────────

def simulate_closed_loop(
    cfg: ScenarioConfig,
    initial_state_lvlh: np.ndarray,
    controller,
    los: LOS3D,
    reference: Dict[str, np.ndarray],
) -> Dict:
    """Run closed-loop simulation with physically honest propagation.

    The state evolves via CWH dynamics + commanded acceleration only.
    No reference blending or state projection is applied.
    Constraint satisfaction is *checked*, not *enforced* post-hoc.
    """

    dt = cfg.ts
    n_steps = cfg.n_steps
    a_max = float(cfg.thrust_to_mass_ratio_mps2)
    n_orbit = float(cfg.mean_motion_rad_s)

    # CWH truth model (with substeps for accuracy)
    n_sub = max(1, cfg.truth_substeps)
    sub_dt = dt / n_sub
    model_sub = discrete_cwh_3d(sub_dt, n_orbit)

    x = np.asarray(initial_state_lvlh, dtype=float).copy()
    if cfg.lock_z_axis:
        x[2] = 0.0
        x[5] = 0.0

    # ── History storage ─────────────────────────────────────────────────
    time_hist = []
    state_hist = []
    body_state_hist = []
    ctrl_hist = []
    ctrl_info_hist = []
    los_margin_hist = []
    feasible_hist = []
    radial_hist = []
    speed_hist = []
    monotonic_violation_hist = []
    mode_hist = []
    phase_err_hist = []
    mpc_cost_hist = []
    cost_breakdown_hist = []

    # Integrity counters
    dv_excess_hist = []       # per-step excess dv (should all be ~0)
    keepout_violation_count = 0
    los_violation_count = 0

    hold_mode_enter_time = None
    hold_success = False
    hold_counter = 0.0

    prev_r = float(np.linalg.norm(x[:3]))

    for k in range(n_steps + 1):
        t = k * dt
        idx_ref = int(np.clip(np.round(t / dt), 0, len(reference["time"]) - 1))

        x_body = los.lvlh_to_body_state(x, t)
        los_margin = los.min_los_margin(x, t)
        r_norm = float(np.linalg.norm(x[:3]))
        v_norm = float(np.linalg.norm(x[3:]))
        phase_err = float(np.linalg.norm(x_body[:3] - reference["r_body"][idx_ref]))

        # ── Diagnostic checks (no state modification) ──────────────────
        keepout_ok = r_norm >= cfg.keepout_radius_m - 1e-9
        if not keepout_ok:
            keepout_violation_count += 1
        if los_margin < -1e-6:
            los_violation_count += 1

        monotonic_ok = True
        if r_norm <= 5.0 and r_norm > cfg.hold_band_max_m and k > 0:
            monotonic_ok = r_norm <= prev_r + 5e-3

        feasible = bool(
            los_margin >= -1e-3
            and keepout_ok
            and v_norm <= cfg.max_rel_speed_mps + 1e-6
        )

        # ── Record ──────────────────────────────────────────────────────
        time_hist.append(t)
        state_hist.append(x.copy())
        body_state_hist.append(x_body)
        los_margin_hist.append(los_margin)
        feasible_hist.append(feasible)
        radial_hist.append(r_norm)
        speed_hist.append(v_norm)
        monotonic_violation_hist.append(not monotonic_ok)
        phase_err_hist.append(phase_err)

        # ── Hold detection ──────────────────────────────────────────────
        is_hold = (cfg.hold_band_min_m <= r_norm <= cfg.hold_band_max_m
                   and v_norm <= cfg.hold_switch_speed_mps)
        if is_hold:
            if hold_mode_enter_time is None:
                hold_mode_enter_time = t
            hold_counter += dt
            if hold_counter >= cfg.hold_confirm_time_s:
                hold_success = True
        else:
            hold_counter = 0.0

        if k == n_steps:
            break

        # ── Control computation ─────────────────────────────────────────
        u, info = controller.compute_control(x, t, a_max)
        u = np.asarray(u, dtype=float)

        # Hard-enforce thrust magnitude limit
        u_norm = float(np.linalg.norm(u))
        if u_norm > a_max:
            u = u * (a_max / u_norm)
        if cfg.lock_z_axis:
            u[2] = 0.0

        mpc_cost_hist.append(float(info.get("mpc_cost", np.nan)))
        cbd = info.get("cost_breakdown", {})
        cost_breakdown_hist.append([
            float(cbd.get("position", 0.0)),
            float(cbd.get("velocity", 0.0)),
            float(cbd.get("input_magnitude", 0.0)),
            float(cbd.get("input_rate", 0.0)),
        ])

        # ── Truth propagation: CWH + control only ──────────────────────
        x_next = x.copy()
        for _ in range(n_sub):
            x_next = model_sub.a_d @ x_next + model_sub.b_d @ u

        if cfg.lock_z_axis:
            x_next[2] = 0.0
            x_next[5] = 0.0

        # ── Integrity check: verify dv is physical ─────────────────────
        dv_actual = np.linalg.norm(x_next[3:] - x[3:])
        dv_thrust = u_norm * dt
        dv_excess_hist.append(max(0.0, dv_actual - dv_thrust - 1e-9))

        # ── Record control ──────────────────────────────────────────────
        mode_now = info.get("mode", "approach")
        ctrl_hist.append(u)
        ctrl_info_hist.append(info)
        mode_hist.append(mode_now)

        prev_r = r_norm
        x = x_next

    # ── Package results ─────────────────────────────────────────────────
    states = np.asarray(state_hist)
    controls = np.asarray(ctrl_hist) if ctrl_hist else np.zeros((0, 3))
    body_states = np.asarray(body_state_hist)
    t_arr = np.asarray(time_hist)

    accel_mag = np.linalg.norm(controls, axis=1) if controls.size else np.zeros(0)
    dv_proxy = float(np.sum(accel_mag) * dt)

    out = {
        "time": t_arr,
        "state_lvlh": states,
        "state_body": body_states,
        "control": controls,
        "mpc_cost_time": t_arr[:-1].copy(),
        "mpc_cost": np.asarray(mpc_cost_hist, dtype=float),
        "cost_breakdown": np.asarray(cost_breakdown_hist, dtype=float)
            if cost_breakdown_hist else np.zeros((0, 4)),
        "control_mode": mode_hist,
        "control_info": ctrl_info_hist,
        "los_margin": np.asarray(los_margin_hist),
        "feasible": np.asarray(feasible_hist, dtype=bool),
        "range_m": np.asarray(radial_hist),
        "speed_mps": np.asarray(speed_hist),
        "monotonic_violation": np.asarray(monotonic_violation_hist, dtype=bool),
        "phase_error_m": np.asarray(phase_err_hist),
        "hold_mode_enter_time_s": (None if hold_mode_enter_time is None
                                   else float(hold_mode_enter_time)),
        "hold_success": bool(hold_success),
        "stays_feasible": bool(np.all(feasible_hist)),
        "min_constraint_margin_m": float(np.min(los_margin_hist)),
        "delta_v_proxy_mps": dv_proxy,
        "thrust_to_mass_mps2": a_max,
        # Integrity diagnostics
        "keepout_violations": keepout_violation_count,
        "los_violations": los_violation_count,
        "dv_excess_total": float(np.sum(dv_excess_hist)),
    }
    return out
