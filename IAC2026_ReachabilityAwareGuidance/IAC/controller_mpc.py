"""MPC-style tracking controller with LOS and safety filtering.

The controller solves a small QP at each step:
    min  w_track * ||e||^2  +  w_u * ||u||^2  +  w_du * ||u - u_prev||^2
    s.t. H u <= g   (LOS, radial, keepout, velocity constraints)
         ||u|| <= a_max  (handled via box approx + clip)

Horizon predictions use CWH (Clohessy-Wiltshire-Hill) free-drift and
input matrices so that the controller is consistent with the truth model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np

from IAC.config import ScenarioConfig
from IAC.constraints import LOS3D

# Try importing scipy for QP solver
try:
    from scipy.optimize import minimize as _scipy_minimize
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


@dataclass
class CostBreakdown:
    """Per-term cost breakdown for diagnostics."""
    position: float = 0.0
    velocity: float = 0.0
    input_magnitude: float = 0.0
    input_rate: float = 0.0
    total: float = 0.0


@dataclass
class ControlInfo:
    mode: str
    success: bool
    min_pred_margin: float
    iterations: int
    radial_bound_m: float
    mpc_cost: float
    cost_breakdown: dict

    def to_dict(self) -> Dict:
        return {
            "mode": self.mode,
            "success": self.success,
            "min_pred_margin": float(self.min_pred_margin),
            "iterations": int(self.iterations),
            "radial_bound_m": float(self.radial_bound_m),
            "mpc_cost": float(self.mpc_cost),
            "cost_breakdown": self.cost_breakdown,
        }


def _cwh_matrices(dt: float, n: float):
    """Return (Phi, B_d) for CWH constant-input over interval dt."""
    nt = n * dt
    c = float(np.cos(nt))
    s = float(np.sin(nt))

    phi = np.array([
        [4.0 - 3.0*c,    0.0, 0.0,  s/n,             2.0*(1.0 - c)/n,    0.0],
        [6.0*(s - nt),   1.0, 0.0, -2.0*(1.0 - c)/n,  (4.0*s - 3.0*nt)/n, 0.0],
        [0.0,            0.0, c,    0.0,               0.0,                s/n],
        [3.0*n*s,        0.0, 0.0,  c,                 2.0*s,              0.0],
        [-6.0*n*(1.0-c), 0.0, 0.0, -2.0*s,             4.0*c - 3.0,        0.0],
        [0.0,            0.0, -n*s, 0.0,               0.0,                c],
    ], dtype=float)

    b_d = np.array([
        [(1.0 - c)/n**2,           2.0*(nt - s)/n**2,               0.0],
        [-2.0*(nt - s)/n**2,       (4.0*(1.0 - c) - 1.5*nt**2)/n**2, 0.0],
        [0.0,                      0.0,                             (1.0 - c)/n**2],
        [s/n,                      2.0*(1.0 - c)/n,                 0.0],
        [-2.0*(1.0 - c)/n,        (4.0*s - 3.0*nt)/n,              0.0],
        [0.0,                      0.0,                             s/n],
    ], dtype=float)

    return phi, b_d


class SpiralMPCController:
    """Two-layer low-thrust controller with CWH-consistent predictions.

    Layer 1: nominal tracking acceleration (PD + CWH feedforward).
    Layer 2: QP-based safety filter with input-rate penalty.
    """

    def __init__(self, cfg: ScenarioConfig, los: LOS3D, ref: Dict[str, np.ndarray],
                 use_fast_solver: bool = False):
        self.cfg = cfg
        self.los = los
        self.ref = ref
        self.ts = cfg.ts
        self.horizon = int(cfg.mpc_horizon_steps)
        self._hold_mode = False
        self._u_prev = np.zeros(3)
        self._use_fast_solver = use_fast_solver

        # Pre-compute CWH matrices for each horizon step
        self._n = float(cfg.mean_motion_rad_s)
        self._horizon_phi = []
        self._horizon_bd = []
        for j in range(1, self.horizon + 1):
            tau = j * self.ts
            phi, bd = _cwh_matrices(tau, self._n)
            self._horizon_phi.append(phi)
            self._horizon_bd.append(bd)

        # Single-step matrices for constraints
        self._phi1, self._bd1 = _cwh_matrices(self.ts, self._n)

    @staticmethod
    def _clip_ball(u: np.ndarray, a_max: float) -> np.ndarray:
        nrm = float(np.linalg.norm(u))
        if nrm <= a_max:
            return u
        return u * (a_max / max(nrm, 1e-9))

    def _idx(self, t: float) -> int:
        idx = int(np.clip(np.round(t / self.ts), 0, len(self.ref["time"]) - 1))
        return idx

    def _approach_axis_weights(self) -> Tuple[np.ndarray, np.ndarray]:
        pos_w = np.array([
            self.cfg.approach_pos_weight_xz,
            self.cfg.approach_pos_weight_y,
            self.cfg.approach_pos_weight_xz,
        ], dtype=float)
        vel_w = np.array([
            self.cfg.approach_vel_weight_xz,
            self.cfg.approach_vel_weight_y,
            self.cfg.approach_vel_weight_xz,
        ], dtype=float)
        return pos_w, vel_w

    def _cwh_gravity_accel(self, x: np.ndarray) -> np.ndarray:
        """CWH gravity-gradient acceleration at state x (LVLH)."""
        n = self._n
        return np.array([
            3.0 * n**2 * x[0] + 2.0 * n * x[4],
            -2.0 * n * x[3],
            -n**2 * x[2],
        ])

    def _nominal_tracking(self, x: np.ndarray, t: float) -> np.ndarray:
        """Compute nominal (unconstrained) tracking acceleration.

        Three regimes:
        1. Hold mode: track co-rotating hold point in LVLH
        2. Far approach (r > r_sync): LVLH-frame approach toward target,
           aiming for a velocity that decelerates to zero at the hold point
        3. Close approach (r <= r_sync): body-frame tracking of LOS centerline

        Includes CWH feedforward so the controller compensates gravity gradient.
        """
        idx = self._idx(t)

        # CWH feedforward: compensate for gravity-gradient drift
        a_ff = self._cwh_gravity_accel(x)

        r = x[:3]
        v = x[3:]
        r_norm = float(np.linalg.norm(r))

        # Body-frame tracking threshold: the co-rotation velocity at range r
        # is omega*r.  We can only track the rotating frame if we can build
        # that velocity in a few timesteps.  Threshold: omega*r < a_max * T_ramp.
        omega = self.los.omega_rad_s
        a_max = self.cfg.thrust_to_mass_ratio_mps2
        t_ramp = 10.0 * self.ts  # 1 second ramp-up budget
        r_sync = a_max * t_ramp / max(omega, 1e-12)

        if self._hold_mode:
            # ── Hold mode: track co-rotating hold point ────────────────
            kp = self.cfg.hold_kp
            kd = self.cfg.hold_kd
            r_ref = self.ref["r_lvlh"][idx]
            v_ref = self.ref["v_lvlh"][idx]
            a_ref = self.ref.get("a_lvlh", np.zeros_like(self.ref["r_lvlh"]))[idx]
            u_nom = a_ref + kp * (r_ref - r) + kd * (v_ref - v) - a_ff

        elif r_norm > r_sync:
            # ── Far approach: LVLH PD tracking toward body-frame ref ───
            # Use LVLH reference (already body->LVLH transformed) so we
            # approach the correct spot but without demanding co-rotation.
            # Only track position; use damping to control velocity.
            kp = self.cfg.tracker_kp
            kd = self.cfg.tracker_kd
            r_ref = self.ref["r_lvlh"][idx]
            v_ref = self.ref["v_lvlh"][idx]

            # Limit reference velocity to what's physically achievable:
            # don't demand the full co-rotation velocity at long range.
            v_ref_mag = float(np.linalg.norm(v_ref))
            v_achievable = a_max * t * 0.5  # half the accumulated dv budget
            v_achievable = min(v_achievable, self.cfg.max_rel_speed_mps)
            if v_ref_mag > v_achievable and v_ref_mag > 1e-9:
                v_ref = v_ref * (v_achievable / v_ref_mag)

            u_nom = kp * (r_ref - r) + kd * (v_ref - v) - a_ff

        else:
            # ── Close approach: body-frame tracking ────────────────────
            kp = self.cfg.tracker_kp
            kd = self.cfg.tracker_kd
            pos_w, vel_w = self._approach_axis_weights()

            x_body = self.los.lvlh_to_body_state(x, t)
            r_body = x_body[:3]
            v_body = x_body[3:]
            lateral_norm = float(np.linalg.norm(r_body[[0, 2]]))
            align_mode = lateral_norm > self.cfg.align_lateral_tol_m

            r_ref_body = self.ref["r_body"][idx].copy()
            v_ref_body = self.ref["v_body"][idx].copy()
            a_ref_body = self.ref.get("a_body", np.zeros_like(self.ref["r_body"]))[idx].copy()
            if self.cfg.use_centerline_reference:
                r_ref_body[0] = 0.0
                r_ref_body[2] = 0.0
                v_ref_body[0] = 0.0
                v_ref_body[2] = 0.0
                a_ref_body[0] = 0.0
                a_ref_body[2] = 0.0

            if align_mode:
                r_ref_body[1] = r_body[1]
                v_ref_body[1] = 0.0
                a_ref_body[1] = 0.0
                pos_w[1] *= self.cfg.align_y_weight_scale
                vel_w[1] *= self.cfg.align_y_weight_scale
            else:
                y_blend = float(np.clip(self.cfg.y_reference_blend_after_align, 0.0, 1.0))
                r_ref_body[1] = (1.0 - y_blend) * r_body[1] + y_blend * r_ref_body[1]
                v_ref_body[1] = y_blend * v_ref_body[1]
                a_ref_body[1] = y_blend * a_ref_body[1]

            a_cmd_body = (a_ref_body
                          + kp * pos_w * (r_ref_body - r_body)
                          + kd * vel_w * (v_ref_body - v_body))

            # Local anti-overshoot damping near centerline for lateral axes
            for axis in (0, 2):
                p_ax = float(r_body[axis])
                v_ax = float(v_body[axis])
                if abs(p_ax) <= self.cfg.lateral_brake_band_m and p_ax * v_ax < 0.0:
                    gain = self.cfg.lateral_brake_gain
                    if abs(v_ax) * self.ts > abs(p_ax):
                        gain *= 1.6
                    a_cmd_body[axis] -= gain * v_ax

            theta = self.los.angle(t)
            r_mat = np.array([
                [np.cos(theta), -np.sin(theta), 0.0],
                [np.sin(theta),  np.cos(theta), 0.0],
                [0.0,            0.0,           1.0],
            ])
            u_nom = r_mat @ a_cmd_body - a_ff

        if self.cfg.lock_z_axis:
            u_nom[2] = 0.0
        return u_nom

    def maybe_switch_hold_mode(self, x: np.ndarray, t: float) -> bool:
        if self._hold_mode:
            return True
        idx = self._idx(t)
        r = float(np.linalg.norm(x[:3]))
        v = float(np.linalg.norm(x[3:]))
        rb = self.los.lvlh_to_body_state(x, t)
        phase_err = float(np.linalg.norm(rb[:3] - self.ref["r_body"][idx]))
        in_band = self.cfg.hold_band_min_m <= r <= self.cfg.hold_band_max_m
        if in_band and v <= self.cfg.hold_switch_speed_mps and phase_err <= 0.08:
            self._hold_mode = True
        return self._hold_mode

    def _build_horizon_constraints(
        self,
        x: np.ndarray,
        t: float,
        radial_upper: float,
        radial_lower: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Build H u <= g using CWH horizon predictions.

        For constant input u, x(t+tau) = Phi(tau) x + Bd(tau) u.
        Position: p(tau) = Phi_pos(tau) x + Bd_pos(tau) u
        """
        e_r = x[:3] / max(np.linalg.norm(x[:3]), 1e-9)

        h_rows = []
        g_vals = []

        for j in range(self.horizon):
            phi = self._horizon_phi[j]
            bd = self._horizon_bd[j]
            tau = (j + 1) * self.ts

            # Free-drift position (u=0): Phi[:3,:] @ x
            p_free = phi[:3, :] @ x
            # Input effect on position: Bd[:3,:] @ u
            bd_pos = bd[:3, :]

            theta = self.los.angle(t + tau)
            r_minus = np.array([
                [ np.cos(theta), np.sin(theta), 0.0],
                [-np.sin(theta), np.cos(theta), 0.0],
                [0.0,            0.0,           1.0],
            ])
            a_pos = self.los.A_c[:, :3] @ r_minus
            h_los = a_pos @ bd_pos
            g_los = self.los.b_c - a_pos @ p_free - self.cfg.safety_margin_m

            h_rows.append(h_los)
            g_vals.append(g_los)

        # Single-step predictions for radial/keepout/velocity
        p1_free = self._phi1[:3, :] @ x
        bd1_pos = self._bd1[:3, :]
        v1_free = self._phi1[3:, :] @ x
        bd1_vel = self._bd1[3:, :]

        if self.cfg.enforce_x_body_monotone and not self._hold_mode:
            x_body_now = self.los.lvlh_to_body_state(x, t)
            x_now = float(x_body_now[0])
            if abs(x_now) > 2e-3:
                theta1 = self.los.angle(t + self.ts)
                r_minus1 = np.array([
                    [ np.cos(theta1), np.sin(theta1), 0.0],
                    [-np.sin(theta1), np.cos(theta1), 0.0],
                    [0.0,             0.0,            1.0],
                ])
                p1_body_free = r_minus1 @ p1_free
                h_body = r_minus1 @ bd1_pos

                sgn = 1.0 if x_now >= 0.0 else -1.0
                x_target = max(0.0, sgn * x_now - self.cfg.x_body_monotone_eps_m)
                h_rows.append(np.array([sgn * h_body[0]]))
                g_vals.append(np.array([x_target - sgn * p1_body_free[0]]))

        # r_{k+1} <= radial_upper (linearised)
        h_rows.append(np.array([e_r @ bd1_pos]))
        g_vals.append(np.array([radial_upper - e_r @ p1_free]))

        # keepout: r_{k+1} >= radial_lower
        h_rows.append(np.array([-(e_r @ bd1_pos)]))
        g_vals.append(np.array([-radial_lower + e_r @ p1_free]))

        # One-step velocity box: |v_i(k+1)| <= v_max
        vmax = self.cfg.max_rel_speed_mps
        for i in range(3):
            # v1 = v1_free + bd1_vel @ u, so bd1_vel[i,:] @ u <= vmax - v1_free[i]
            h_rows.append(np.array([bd1_vel[i, :]]))
            g_vals.append(np.array([vmax - v1_free[i]]))

            h_rows.append(np.array([-bd1_vel[i, :]]))
            g_vals.append(np.array([vmax + v1_free[i]]))

        h = np.vstack(h_rows)
        g = np.hstack(g_vals)
        return h, g

    def _solve_qp(
        self,
        u_nom: np.ndarray,
        h: np.ndarray,
        g: np.ndarray,
        a_max: float,
    ) -> Tuple[np.ndarray, int]:
        """Solve the smoothed QP: min rate-penalised cost s.t. constraints."""
        w_u = max(self.cfg.input_magnitude_weight, 0.0)
        w_du = max(self.cfg.input_rate_weight, 0.0)

        # Add box constraints for thrust magnitude
        h_box = np.vstack([np.eye(3), -np.eye(3)])
        g_box = a_max * np.ones(6)
        H_all = np.vstack([h, h_box])
        g_all = np.hstack([g, g_box])

        if _HAS_SCIPY:
            def objective(u):
                du = u - self._u_prev
                return (float(np.sum((u - u_nom)**2))
                        + w_u * float(u @ u)
                        + w_du * float(du @ du))

            def gradient(u):
                du = u - self._u_prev
                return (2.0 * (u - u_nom)
                        + 2.0 * w_u * u
                        + 2.0 * w_du * du)

            u0 = self._clip_ball(u_nom.copy(), a_max)
            if self.cfg.lock_z_axis:
                u0[2] = 0.0

            bounds = [(-a_max, a_max), (-a_max, a_max),
                      (0.0, 0.0) if self.cfg.lock_z_axis else (-a_max, a_max)]

            constraints_list = []
            for i in range(H_all.shape[0]):
                hi = H_all[i].copy()
                gi = float(g_all[i])
                constraints_list.append({
                    'type': 'ineq',
                    'fun': lambda u, _h=hi, _g=gi: _g - float(_h @ u),
                    'jac': lambda u, _h=hi: -_h,
                })

            try:
                result = _scipy_minimize(
                    objective, u0, jac=gradient,
                    method='SLSQP',
                    bounds=bounds,
                    constraints=constraints_list,
                    options={'maxiter': 50, 'ftol': 1e-10},
                )
                u_sol = result.x.copy()
                n_iter = int(result.nit)
                u_sol = self._clip_ball(u_sol, a_max)
                if self.cfg.lock_z_axis:
                    u_sol[2] = 0.0
                return u_sol, n_iter
            except Exception:
                pass

        return self._project_with_rate(u_nom, h, g, a_max)

    def _project_with_rate(
        self,
        u_nom: np.ndarray,
        h: np.ndarray,
        g: np.ndarray,
        a_max: float,
    ) -> Tuple[np.ndarray, int]:
        """Projection fallback with rate-penalty blending."""
        w_du = max(self.cfg.input_rate_weight, 0.0)
        blend = w_du / (1.0 + w_du) if w_du > 0 else 0.0

        u_target = (1.0 - blend) * u_nom + blend * self._u_prev
        u = self._clip_ball(u_target.copy(), a_max)
        if self.cfg.lock_z_axis:
            u[2] = 0.0

        n_iter = 0
        for _ in range(self.cfg.safety_iterations):
            n_iter += 1
            changed = False
            for i in range(h.shape[0]):
                hi = h[i]
                lhs = float(hi @ u)
                if lhs > g[i]:
                    denom = float(hi @ hi) + 1e-12
                    u = u - ((lhs - g[i]) / denom) * hi
                    u = self._clip_ball(u, a_max)
                    if self.cfg.lock_z_axis:
                        u[2] = 0.0
                    changed = True
            if not changed:
                break
        return u, n_iter

    def _horizon_tracking_cost(
        self, x: np.ndarray, t: float, u: np.ndarray
    ) -> Tuple[float, CostBreakdown]:
        """Compute MPC-like horizon cost for diagnostics (CWH-based)."""
        idx0 = self._idx(t)
        dt = self.ts

        w_r = max(self.cfg.tracker_kp, 1e-3)
        w_v = max(self.cfg.tracker_kd, 1e-3)
        w_u = max(self.cfg.input_magnitude_weight, 0.0)
        w_du = max(self.cfg.input_rate_weight, 0.0)
        pos_w, vel_w = self._approach_axis_weights()

        j_pos = 0.0
        j_vel = 0.0
        j_input = 0.0
        n_ref = len(self.ref["time"]) - 1

        for h_idx in range(self.horizon):
            phi = self._horizon_phi[h_idx]
            bd = self._horizon_bd[h_idx]
            tau = (h_idx + 1) * dt

            x_pred = phi @ x + bd @ u
            r_pred = x_pred[:3]
            v_pred = x_pred[3:]
            idx = int(np.clip(idx0 + h_idx + 1, 0, n_ref))

            if self._hold_mode:
                r_ref = self.ref["r_lvlh"][idx]
                v_ref = self.ref["v_lvlh"][idx]
                e_r = r_pred - r_ref
                e_v = v_pred - v_ref
                j_pos += w_r * float(e_r @ e_r)
                j_vel += w_v * float(e_v @ e_v)
                j_input += w_u * float(u @ u)
            else:
                x_pred_body = self.los.lvlh_to_body_state(x_pred, t + tau)
                r_ref_body = self.ref["r_body"][idx].copy()
                v_ref_body = self.ref["v_body"][idx].copy()
                if self.cfg.use_centerline_reference:
                    r_ref_body[0] = 0.0
                    r_ref_body[2] = 0.0
                    v_ref_body[0] = 0.0
                    v_ref_body[2] = 0.0

                e_r_body = x_pred_body[:3] - r_ref_body
                e_v_body = x_pred_body[3:] - v_ref_body
                j_pos += w_r * float(np.sum(pos_w * (e_r_body ** 2)))
                j_vel += w_v * float(np.sum(vel_w * (e_v_body ** 2)))
                j_input += w_u * float(u @ u)

        du = u - self._u_prev
        j_rate = w_du * float(du @ du)

        total = j_pos + j_vel + j_input + j_rate
        breakdown = CostBreakdown(
            position=j_pos,
            velocity=j_vel,
            input_magnitude=j_input,
            input_rate=j_rate,
            total=total,
        )
        return total, breakdown

    def compute_control(
        self, x: np.ndarray, t: float, a_max: float
    ) -> Tuple[np.ndarray, Dict]:
        x = np.asarray(x, dtype=float)
        self.maybe_switch_hold_mode(x, t)

        r_now = float(np.linalg.norm(x[:3]))
        if self._hold_mode:
            radial_upper = self.cfg.hold_band_max_m
            radial_lower = self.cfg.hold_band_min_m
            mode = "hold"
        else:
            if r_now > 5.0:
                radial_upper = r_now
            else:
                radial_upper = max(self.cfg.hold_radius_m,
                                   r_now - self.cfg.radial_decay_eps_m)
            radial_lower = self.cfg.keepout_radius_m
            mode = "approach"

        u_nom = self._nominal_tracking(x, t)
        h, g = self._build_horizon_constraints(
            x=x, t=t,
            radial_upper=radial_upper,
            radial_lower=radial_lower,
        )

        if self._use_fast_solver:
            u, n_iter = self._project_with_rate(u_nom, h, g, a_max)
        else:
            u, n_iter = self._solve_qp(u_nom, h, g, a_max)
        if self.cfg.lock_z_axis:
            u[2] = 0.0

        residual = g - h @ u
        min_margin = float(np.min(residual)) if residual.size else 0.0
        mpc_cost, breakdown = self._horizon_tracking_cost(x, t, u)

        self._u_prev = u.copy()

        info = ControlInfo(
            mode=mode,
            success=min_margin >= -1e-6,
            min_pred_margin=min_margin,
            iterations=n_iter,
            radial_bound_m=radial_upper,
            mpc_cost=mpc_cost,
            cost_breakdown={
                "position": breakdown.position,
                "velocity": breakdown.velocity,
                "input_magnitude": breakdown.input_magnitude,
                "input_rate": breakdown.input_rate,
                "total": breakdown.total,
            },
        )
        return u, info.to_dict()
