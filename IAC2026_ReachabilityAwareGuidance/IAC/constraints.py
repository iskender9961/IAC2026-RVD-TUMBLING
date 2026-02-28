"""LOS constraints and frame transforms for the IAC 2026 demo."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from IAC.config import ScenarioConfig


def rotz(theta: float) -> np.ndarray:
    """Rotation matrix Rz(theta)."""

    c = float(np.cos(theta))
    s = float(np.sin(theta))
    return np.array(
        [
            [c, -s, 0.0],
            [s, c, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )


def omega_cross_matrix(omega_z: float) -> np.ndarray:
    """Skew matrix such that Omega @ r == omega x r for omega=[0,0,omega_z]."""

    return np.array(
        [
            [0.0, -omega_z, 0.0],
            [omega_z, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )


def body_to_lvlh(
    r_body: np.ndarray,
    v_body: np.ndarray,
    theta: float,
    omega_z: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Map body-frame relative state to LVLH.

    r_L = R(theta) r_B
    v_L = R(theta) (v_B + omega x r_B)
    """

    r_body = np.asarray(r_body, dtype=float)
    v_body = np.asarray(v_body, dtype=float)
    r_mat = rotz(theta)
    omega_cross_r = np.array([-omega_z * r_body[1], omega_z * r_body[0], 0.0])
    r_lvlh = r_mat @ r_body
    v_lvlh = r_mat @ (v_body + omega_cross_r)
    return r_lvlh, v_lvlh


def lvlh_to_body(
    r_lvlh: np.ndarray,
    v_lvlh: np.ndarray,
    theta: float,
    omega_z: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Map LVLH relative state to body frame.

    r_B = R(-theta) r_L
    v_B = R(-theta) v_L - omega x r_B
    """

    r_lvlh = np.asarray(r_lvlh, dtype=float)
    v_lvlh = np.asarray(v_lvlh, dtype=float)
    r_minus = rotz(-theta)
    r_body = r_minus @ r_lvlh
    v_body = r_minus @ v_lvlh - np.array([-omega_z * r_body[1], omega_z * r_body[0], 0.0])
    return r_body, v_body


@dataclass
class LOS3D:
    """Polyhedral LOS corridor in target body frame.

    State convention is x=[x,y,z,vx,vy,vz].

    Constraint template requested by scenario:
        A_c x <= b_c

    where
        A_c = [
          [  0,  -1,   0, 0,0,0],
          [  1, -cx,   0, 0,0,0],
          [ -1, -cx,   0, 0,0,0],
          [  0, -cz,   1, 0,0,0],
          [  0, -cz,  -1, 0,0,0],
        ]

    and b_c is defined here as:
        b1 = -r_hold
        b2 = x0 - cx * r_hold
        b3 = x0 - cx * r_hold
        b4 = z0 - cz * r_hold
        b5 = z0 - cz * r_hold

    This yields a wedge/tetrahedron-like corridor around +y_B:
        y_B >= r_hold                        (anti-overshoot)
        |x_B| <= x0 + cx * (y_B - r_hold)   (lateral cone)
        |z_B| <= z0 + cz * (y_B - r_hold)   (out-of-plane cone)
    """

    x0_m: float
    z0_m: float
    cx: float
    cz: float
    hold_radius_m: float
    keepout_radius_m: float
    omega_rad_s: float

    @classmethod
    def from_config(cls, cfg: ScenarioConfig) -> "LOS3D":
        return cls(
            x0_m=cfg.los_x0_m,
            z0_m=cfg.los_z0_m,
            cx=cfg.los_cx,
            cz=cfg.los_cz,
            hold_radius_m=cfg.hold_radius_m,
            keepout_radius_m=cfg.keepout_radius_m,
            omega_rad_s=cfg.tumble_omega_rad_s,
        )

    def angle(self, t: float) -> float:
        return self.omega_rad_s * t

    @property
    def A_c(self) -> np.ndarray:
        return np.array(
            [
                [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
                [1.0, -self.cx, 0.0, 0.0, 0.0, 0.0],
                [-1.0, -self.cx, 0.0, 0.0, 0.0, 0.0],
                [0.0, -self.cz, 1.0, 0.0, 0.0, 0.0],
                [0.0, -self.cz, -1.0, 0.0, 0.0, 0.0],
            ]
        )

    @property
    def b_c(self) -> np.ndarray:
        rh = self.hold_radius_m
        return np.array(
            [
                -rh,
                self.x0_m - self.cx * rh,
                self.x0_m - self.cx * rh,
                self.z0_m - self.cz * rh,
                self.z0_m - self.cz * rh,
            ]
        )

    def lvlh_to_body_state(self, x_lvlh: np.ndarray, t: float) -> np.ndarray:
        x_lvlh = np.asarray(x_lvlh, dtype=float)
        r_lvlh = x_lvlh[:3]
        v_lvlh = x_lvlh[3:]
        r_body, v_body = lvlh_to_body(r_lvlh, v_lvlh, theta=self.angle(t), omega_z=self.omega_rad_s)
        return np.hstack((r_body, v_body))

    def los_slack(self, x_lvlh: np.ndarray, t: float) -> np.ndarray:
        x_body = self.lvlh_to_body_state(x_lvlh, t)
        return self.b_c - self.A_c @ x_body

    def min_los_margin(self, x_lvlh: np.ndarray, t: float) -> float:
        return float(np.min(self.los_slack(x_lvlh, t)))

    def is_los_feasible(self, x_lvlh: np.ndarray, t: float) -> bool:
        return bool(np.all(self.los_slack(x_lvlh, t) >= -1e-9))

    def radial_distance(self, x_lvlh: np.ndarray) -> float:
        return float(np.linalg.norm(np.asarray(x_lvlh[:3], dtype=float)))

    def is_keepout_respected(self, x_lvlh: np.ndarray) -> bool:
        return self.radial_distance(x_lvlh) >= self.keepout_radius_m - 1e-9

    def corridor_vertices_body(self, y_far: float) -> np.ndarray:
        """Return 8 body-frame vertices for wedge visualization."""

        y_near = self.hold_radius_m
        y_far = max(y_far, y_near + 1.0)

        x_near = self.x0_m
        z_near = self.z0_m
        x_far = self.cx * (y_far - y_near) + self.x0_m
        z_far = self.cz * (y_far - y_near) + self.z0_m

        return np.array(
            [
                [-x_near, y_near, -z_near],
                [x_near, y_near, -z_near],
                [x_near, y_near, z_near],
                [-x_near, y_near, z_near],
                [-x_far, y_far, -z_far],
                [x_far, y_far, -z_far],
                [x_far, y_far, z_far],
                [-x_far, y_far, z_far],
            ],
            dtype=float,
        )

    def corridor_faces(self) -> list[list[int]]:
        """Face indices for wedge/tetrahedron-like frustum."""

        return [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        ]

    def corridor_vertices_lvlh(self, t: float, y_far: float) -> np.ndarray:
        verts_b = self.corridor_vertices_body(y_far)
        r_mat = rotz(self.angle(t))
        return (r_mat @ verts_b.T).T

    def los_margin_history(self, states_lvlh: np.ndarray, t_grid: np.ndarray) -> np.ndarray:
        margins = np.zeros(len(t_grid), dtype=float)
        for i, (x, t) in enumerate(zip(states_lvlh, t_grid)):
            margins[i] = self.min_los_margin(x, float(t))
        return margins

    def body_frame_history(self, states_lvlh: np.ndarray, t_grid: np.ndarray) -> np.ndarray:
        out = np.zeros_like(states_lvlh)
        for i, (x, t) in enumerate(zip(states_lvlh, t_grid)):
            out[i] = self.lvlh_to_body_state(x, float(t))
        return out
