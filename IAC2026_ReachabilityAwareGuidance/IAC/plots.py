"""Plot generation for the upgraded IAC scenario."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from IAC.config import ScenarioConfig
from IAC.constraints import LOS3D


def _style() -> None:
    plt.rcParams.update(
        {
            "figure.figsize": (9.0, 6.2),
            "font.size": 10,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "savefig.dpi": 220,
        }
    )


def _save_figure(fig: plt.Figure, out_path: Path) -> str:
    """Save figure as PNG and also as PDF (always)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    try:
        fig.savefig(out_path.with_suffix(".pdf"))
    except PermissionError:
        print(f"  [warn] PDF locked, skipped: {out_path.with_suffix('.pdf')}")
    plt.close(fig)
    return str(out_path)


def _try_pdf(fig: plt.Figure, out_path: Path) -> None:
    """Attempt PDF save, silently skip if file is locked."""
    try:
        fig.savefig(out_path.with_suffix(".pdf"))
    except PermissionError:
        print(f"  [warn] PDF locked, skipped: {out_path.with_suffix('.pdf')}")


def _normalize(vec: np.ndarray, fallback: np.ndarray) -> np.ndarray:
    nrm = float(np.linalg.norm(vec))
    if nrm < 1e-9:
        return fallback.copy()
    return vec / nrm


def _chaser_body_rotation(state_lvlh: np.ndarray) -> np.ndarray:
    vel = np.asarray(state_lvlh[3:6], dtype=float)
    x_c = _normalize(vel, np.array([1.0, 0.0, 0.0]))

    z_hint = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(x_c, z_hint))) > 0.95:
        z_hint = np.array([0.0, 1.0, 0.0])

    y_c = _normalize(np.cross(z_hint, x_c), np.array([0.0, 1.0, 0.0]))
    z_c = _normalize(np.cross(x_c, y_c), np.array([0.0, 0.0, 1.0]))
    return np.column_stack((x_c, y_c, z_c))


def _target_in_chaser_history(states_lvlh: np.ndarray) -> np.ndarray:
    out = np.zeros((len(states_lvlh), 3), dtype=float)
    for i in range(len(states_lvlh)):
        r_chaser_to_lvlh = _chaser_body_rotation(states_lvlh[i])
        out[i] = r_chaser_to_lvlh.T @ (-states_lvlh[i, :3])
    return out


def _plot_single_series(
    t: np.ndarray,
    y: np.ndarray,
    title: str,
    y_label: str,
    out_path: Path,
    color: str = "tab:blue",
    zero_line: bool = False,
) -> str:
    _style()
    n = min(len(t), len(y))
    t_s = np.asarray(t[:n], dtype=float)
    y_s = np.asarray(y[:n], dtype=float)

    fig, ax = plt.subplots(figsize=(9, 4.2))
    if n > 0:
        ax.plot(t_s, y_s, color=color, lw=1.9)
    if zero_line:
        ax.axhline(0.0, color="k", lw=1.0, ls="--")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    fig.tight_layout()
    return _save_figure(fig, out_path)


def _plot_component_triplet(
    t: np.ndarray,
    data_xyz: np.ndarray,
    y_labels: tuple[str, str, str],
    title: str,
    out_path: Path,
    colors: tuple[str, str, str] = ("tab:red", "tab:green", "tab:blue"),
) -> str:
    _style()
    n = min(len(t), len(data_xyz))
    t_s = np.asarray(t[:n], dtype=float)
    xyz = np.asarray(data_xyz[:n], dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(9, 7.6), sharex=True)
    for i in range(3):
        if n > 0:
            axes[i].plot(t_s, xyz[:, i], color=colors[i], lw=1.8)
        axes[i].set_ylabel(y_labels[i])
        axes[i].grid(True, alpha=0.28)
    axes[0].set_title(title)
    axes[-1].set_xlabel("Time (s)")
    fig.tight_layout()
    return _save_figure(fig, out_path)


def _add_corridor(ax, vertices: np.ndarray, faces: list[list[int]], color: str = "tab:orange", alpha: float = 0.12) -> None:
    poly_faces = [vertices[idx] for idx in faces]
    coll = Poly3DCollection(poly_faces, facecolor=color, edgecolor=color, linewidth=0.8, alpha=alpha)
    ax.add_collection3d(coll)


def _docking_box_vertices(half_x: float, half_y: float, half_z: float) -> np.ndarray:
    return np.array(
        [
            [-half_x, -half_y, -half_z],
            [half_x, -half_y, -half_z],
            [half_x, half_y, -half_z],
            [-half_x, half_y, -half_z],
            [-half_x, -half_y, half_z],
            [half_x, -half_y, half_z],
            [half_x, half_y, half_z],
            [-half_x, half_y, half_z],
        ],
        dtype=float,
    )


def _edges_from_faces(faces_idx: list[list[int]]) -> list[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for face in faces_idx:
        for i in range(len(face)):
            a = int(face[i])
            b = int(face[(i + 1) % len(face)])
            edges.add((a, b) if a < b else (b, a))
    return sorted(edges)


def _draw_wire_box(ax, vertices: np.ndarray, edges: list[tuple[int, int]], color: str = "k", lw: float = 1.4, alpha: float = 1.0) -> None:
    for a, b in edges:
        ax.plot(
            [vertices[a, 0], vertices[b, 0]],
            [vertices[a, 1], vertices[b, 1]],
            [vertices[a, 2], vertices[b, 2]],
            color=color,
            lw=lw,
            alpha=alpha,
        )


def plot_traj_lvlh_3d(cfg: ScenarioConfig, sim: Dict, los: LOS3D, out_path: Path) -> str:
    _style()
    states = sim["state_lvlh"]
    t = sim["time"]

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(states[:, 0], states[:, 1], states[:, 2], color="tab:blue", lw=2.0, label="Chaser path")

    box_half_x = max(0.55, 1.0 * cfg.hold_radius_m)
    box_half_y = max(0.45, 0.9 * cfg.hold_radius_m)
    box_half_z = max(0.35, 0.7 * cfg.hold_radius_m)
    target_box = _docking_box_vertices(box_half_x, box_half_y, box_half_z)
    box_faces = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ]
    _draw_wire_box(ax, target_box, _edges_from_faces(box_faces), color="k", lw=1.5, alpha=0.95)

    faces = los.corridor_faces()
    for ts in np.linspace(t[0], t[-1], 5):
        verts_l = los.corridor_vertices_lvlh(float(ts), y_far=cfg.corridor_far_y_m)
        _add_corridor(ax, verts_l, faces, color="tab:orange", alpha=0.08)

    lim = max(10.0, 1.05 * np.max(np.abs(states[:, :3])))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-0.45 * lim, 0.45 * lim)
    ax.set_xlabel("x_L (m)")
    ax.set_ylabel("y_L (m)")
    ax.set_zlabel("z_L (m)")
    ax.set_title(f"LVLH 3D Trajectory with Rotating LOS Corridor\n({cfg.scenario_tag})")
    ax.view_init(elev=24, azim=-46)
    ax.legend(loc="upper right")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_traj_body_3d(cfg: ScenarioConfig, sim: Dict, los: LOS3D, out_path: Path) -> str:
    _style()
    states_b = sim["state_body"]

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(states_b[:, 0], states_b[:, 1], states_b[:, 2], color="tab:green", lw=2.0, label="Chaser in body frame")

    box_half_x = max(0.55, 1.0 * cfg.hold_radius_m)
    box_half_y = max(0.45, 0.9 * cfg.hold_radius_m)
    box_half_z = max(0.35, 0.7 * cfg.hold_radius_m)
    target_box = _docking_box_vertices(box_half_x, box_half_y, box_half_z)
    box_faces = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ]
    _draw_wire_box(ax, target_box, _edges_from_faces(box_faces), color="k", lw=1.5, alpha=0.95)

    verts_b = los.corridor_vertices_body(y_far=cfg.corridor_far_y_m)
    _add_corridor(ax, verts_b, los.corridor_faces(), color="tab:red", alpha=0.10)

    # Fixed docking direction (+y_B)
    ax.quiver(0.0, 0.0, 0.0, 0.0, cfg.hold_radius_m * 2.8, 0.0, color="tab:red", linewidth=2.2)
    ax.text(0.0, cfg.hold_radius_m * 3.0, 0.0, "Docking direction (+y_B)")

    lim = max(10.0, 1.05 * np.max(np.abs(states_b[:, :3])))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(0.0, lim)
    ax.set_zlim(-0.45 * lim, 0.45 * lim)
    ax.set_xlabel("x_B (m)")
    ax.set_ylabel("y_B (m)")
    ax.set_zlabel("z_B (m)")
    ax.set_title(f"Target Body Frame Trajectory and Docking Corridor\n({cfg.scenario_tag})")
    ax.view_init(elev=21, azim=-56)
    ax.legend(loc="upper right")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_range_time(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t = sim["time"]
    r = sim["range_m"]
    hold_t = sim["hold_mode_enter_time_s"]

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.plot(t, r, color="tab:blue", lw=2.0, label="||r||")

    r_prev = r[:-1]
    r_next = r[1:]
    eps = cfg.radial_decay_eps_m
    mono_ok = np.ones_like(r, dtype=bool)
    mono_ok[1:] = r_next <= (r_prev - eps + 1e-6)
    ax.scatter(t[mono_ok], r[mono_ok], s=7, color="tab:green", label="Monotonic step ok")
    ax.scatter(t[~mono_ok], r[~mono_ok], s=12, color="tab:red", label="Monotonic violation")

    ax.axhline(cfg.hold_radius_m, color="k", ls="--", lw=1.0, label="Hold radius 0.50 m")
    ax.axhspan(cfg.hold_band_min_m, cfg.hold_band_max_m, color="tab:gray", alpha=0.15, label="Hold band")
    if hold_t is not None:
        ax.axvline(hold_t, color="tab:purple", lw=1.3, ls="-.", label="Hold mode switch")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Range ||r|| (m)")
    ax.set_title(f"Range vs Time (Monotonic Decrease + Hold Transition)\n({cfg.scenario_tag})")
    ax.legend(loc="upper right", ncols=2)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_constraint_margin(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t = sim["time"]
    m = sim["los_margin"]

    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.plot(t, m, color="tab:orange", lw=1.9)
    ax.axhline(0.0, color="k", lw=1.0, ls="--")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("min(b - A x)")
    ax.set_title(f"LOS Constraint Margin vs Time\n({cfg.scenario_tag})")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_speed_time(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t = sim["time"]
    v = sim["speed_mps"]

    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.plot(t, v, color="tab:cyan", lw=1.9)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("||v|| (m/s)")
    ax.set_title(f"Relative Speed vs Time\n({cfg.scenario_tag})")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_control_and_dv(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    u = sim["control"]
    if u.size == 0:
        u_norm = np.zeros(0)
        t_u = np.zeros(0)
    else:
        u_norm = np.linalg.norm(u, axis=1)
        t_u = sim["time"][:-1]

    dv = np.cumsum(u_norm) * cfg.ts if u_norm.size else np.zeros(0)

    fig, ax1 = plt.subplots(figsize=(9, 4.8))
    ax1.plot(t_u, u_norm, color="tab:red", lw=1.7, label="||a_cmd||")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Acceleration magnitude (m/s^2)", color="tab:red")
    ax1.tick_params(axis="y", labelcolor="tab:red")

    ax2 = ax1.twinx()
    ax2.plot(t_u, dv, color="tab:blue", lw=1.7, label="Cumulative Delta-v proxy")
    ax2.set_ylabel("Cumulative Delta-v proxy (m/s)", color="tab:blue")
    ax2.tick_params(axis="y", labelcolor="tab:blue")

    ax1.set_title(f"Control Effort and Cumulative Delta-v Proxy\n({cfg.scenario_tag})")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_phase_sync(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t = sim["time"]
    phase_err = sim["phase_error_m"]

    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.plot(t, phase_err, color="tab:purple", lw=1.9)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Body-frame phase error (m)")
    ax.set_title(f"Body-Frame Alignment / Phase Synchronization Metric\n({cfg.scenario_tag})")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_mpc_cost(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t_cost = np.asarray(sim.get("mpc_cost_time", sim["time"][:-1]), dtype=float)
    j_cost = np.asarray(sim.get("mpc_cost", np.zeros_like(t_cost)), dtype=float)
    n = min(len(t_cost), len(j_cost))
    t_cost = t_cost[:n]
    j_cost = j_cost[:n]

    fig, axes = plt.subplots(2, 1, figsize=(9, 7.0), sharex=True)

    # Top: log scale (handles large initial transient)
    if n > 0:
        j_pos = np.maximum(j_cost, 1e-6)
        axes[0].semilogy(t_cost, j_pos, color="tab:purple", lw=1.9)
    axes[0].set_ylabel("MPC cost J (log)")
    axes[0].set_title(f"MPC Cost vs Time (Log Scale) — {cfg.scenario_tag}")
    axes[0].grid(True, alpha=0.28)

    # Bottom: linear scale, clipped to exclude initial transient
    if n > 0:
        # Find where cost drops below 2x median (skip initial spike)
        median_cost = float(np.median(j_cost[min(20, n):]))
        clip_val = max(median_cost * 10.0, 1.0)
        j_clipped = np.minimum(j_cost, clip_val)
        axes[1].plot(t_cost, j_clipped, color="tab:purple", lw=1.9)
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("MPC cost J (linear, clipped)")
    axes[1].set_title("MPC Cost vs Time (Linear, Initial Transient Clipped)")
    axes[1].grid(True, alpha=0.28)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_mpc_cost_breakdown(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    """Per-term cost breakdown: position, velocity, input magnitude, input rate."""
    _style()
    t_cost = np.asarray(sim.get("mpc_cost_time", sim["time"][:-1]), dtype=float)
    cb = np.asarray(sim.get("cost_breakdown", np.zeros((0, 4))), dtype=float)
    n = min(len(t_cost), len(cb))

    fig, axes = plt.subplots(4, 1, figsize=(9, 9.5), sharex=True)
    labels = ["Position cost", "Velocity cost", "Input magnitude cost", "Input rate cost"]
    colors = ["tab:blue", "tab:orange", "tab:red", "tab:green"]

    for i in range(4):
        if n > 0:
            vals = cb[:n, i]
            axes[i].semilogy(t_cost[:n], np.maximum(vals, 1e-8), color=colors[i], lw=1.7)
        axes[i].set_ylabel(labels[i])
        axes[i].grid(True, alpha=0.28)

    axes[0].set_title(f"MPC Cost Per-Term Breakdown (Log Scale) — {cfg.scenario_tag}")
    axes[-1].set_xlabel("Time (s)")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_state_components_lvlh(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t = sim["time"]
    x = sim["state_lvlh"]

    fig, axes = plt.subplots(2, 1, figsize=(10, 7.2), sharex=True)
    for i, lbl in enumerate(("x_L", "y_L", "z_L")):
        axes[0].plot(t, x[:, i], lw=1.7, label=lbl)
    axes[0].set_ylabel("Position (m)")
    axes[0].set_title(f"LVLH Position Components vs Time — {cfg.scenario_tag}")
    axes[0].legend(loc="upper right", ncols=3)

    for i, lbl in enumerate(("vx_L", "vy_L", "vz_L")):
        axes[1].plot(t, x[:, i + 3], lw=1.7, label=lbl)
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Velocity (m/s)")
    axes[1].set_title("LVLH Velocity Components vs Time")
    axes[1].legend(loc="upper right", ncols=3)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_state_components_body(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t = sim["time"]
    x_b = sim["state_body"]

    fig, axes = plt.subplots(2, 1, figsize=(10, 7.2), sharex=True)
    for i, lbl in enumerate(("x_B", "y_B", "z_B")):
        axes[0].plot(t, x_b[:, i], lw=1.7, label=lbl)
    axes[0].set_ylabel("Position (m)")
    axes[0].set_title(f"Target-Body Position Components vs Time — {cfg.scenario_tag}")
    axes[0].legend(loc="upper right", ncols=3)

    for i, lbl in enumerate(("vx_B", "vy_B", "vz_B")):
        axes[1].plot(t, x_b[:, i + 3], lw=1.7, label=lbl)
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Velocity (m/s)")
    axes[1].set_title("Target-Body Velocity Components vs Time")
    axes[1].legend(loc="upper right", ncols=3)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_body_position_time(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t = sim["time"]
    x_b = sim["state_body"]

    fig, axes = plt.subplots(3, 1, figsize=(9, 7.4), sharex=True)
    labels = ("x_B (m)", "y_B (m)", "z_B (m)")
    colors = ("tab:red", "tab:green", "tab:blue")
    for i in range(3):
        axes[i].plot(t, x_b[:, i], color=colors[i], lw=1.8)
        axes[i].set_ylabel(labels[i])
        axes[i].grid(True, alpha=0.28)
    axes[0].set_title(f"Body-Frame Position Components vs Time — {cfg.scenario_tag}")
    axes[-1].set_xlabel("Time (s)")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_control_components(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t_u = sim["time"][:-1]
    u = sim["control"]

    fig, ax = plt.subplots(figsize=(9, 4.8))
    if u.size:
        ax.plot(t_u, u[:, 0], lw=1.6, label="a_x")
        ax.plot(t_u, u[:, 1], lw=1.6, label="a_y")
        ax.plot(t_u, u[:, 2], lw=1.6, label="a_z")
        ax.plot(t_u, np.linalg.norm(u, axis=1), lw=1.8, ls="--", color="k", label="||a||")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Acceleration (m/s^2)")
    ax.set_title(f"Control Components vs Time — {cfg.scenario_tag}")
    ax.legend(loc="upper right", ncols=4)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_fuel_use_time(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t_u = sim["time"][:-1]
    u = sim["control"]
    dv = np.cumsum(np.linalg.norm(u, axis=1)) * cfg.ts if u.size else np.zeros(0)

    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.plot(t_u, dv, color="tab:red", lw=1.9)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Cumulative Delta-v proxy (m/s)")
    ax.set_title(f"Fuel-Use Proxy vs Time — {cfg.scenario_tag}")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    _try_pdf(fig, out_path)
    plt.close(fig)
    return str(out_path)


def plot_lvlh_position_time(sim: Dict, out_path: Path) -> str:
    t = sim["time"]
    x = sim["state_lvlh"][:, :3]
    return _plot_component_triplet(
        t=t,
        data_xyz=x,
        y_labels=("x_L (m)", "y_L (m)", "z_L (m)"),
        title="LVLH Position Components vs Time (Separate Axes)",
        out_path=out_path,
    )


def plot_lvlh_velocity_time(sim: Dict, out_path: Path) -> str:
    t = sim["time"]
    v = sim["state_lvlh"][:, 3:6]
    return _plot_component_triplet(
        t=t,
        data_xyz=v,
        y_labels=("vx_L (m/s)", "vy_L (m/s)", "vz_L (m/s)"),
        title="LVLH Velocity Components vs Time (Separate Axes)",
        out_path=out_path,
        colors=("tab:cyan", "tab:blue", "tab:purple"),
    )


def plot_target_body_position_time(sim: Dict, out_path: Path) -> str:
    t = sim["time"]
    x_b = sim["state_body"][:, :3]
    return _plot_component_triplet(
        t=t,
        data_xyz=x_b,
        y_labels=("X_T_body (m)", "Y_T_body (m)", "Z_T_body (m)"),
        title="Target-Body Position Components vs Time",
        out_path=out_path,
    )


def plot_target_body_velocity_time(sim: Dict, out_path: Path) -> str:
    t = sim["time"]
    v_b = sim["state_body"][:, 3:6]
    return _plot_component_triplet(
        t=t,
        data_xyz=v_b,
        y_labels=("VX_T_body (m/s)", "VY_T_body (m/s)", "VZ_T_body (m/s)"),
        title="Target-Body Velocity Components vs Time",
        out_path=out_path,
        colors=("tab:olive", "tab:green", "tab:cyan"),
    )


def plot_target_in_chaser_position_time(sim: Dict, out_path: Path) -> str:
    t = sim["time"]
    x_c = _target_in_chaser_history(sim["state_lvlh"])
    return _plot_component_triplet(
        t=t,
        data_xyz=x_c,
        y_labels=("X_T_chaser (m)", "Y_T_chaser (m)", "Z_T_chaser (m)"),
        title="Target Position in Chaser Body Frame vs Time",
        out_path=out_path,
        colors=("tab:orange", "tab:brown", "tab:red"),
    )


def plot_control_norm_time(sim: Dict, out_path: Path) -> str:
    t_u = sim["time"][:-1]
    u = sim["control"]
    u_norm = np.linalg.norm(u, axis=1) if u.size else np.zeros(0)
    return _plot_single_series(
        t=t_u,
        y=u_norm,
        title="Control Norm vs Time",
        y_label="||a_cmd|| (m/s^2)",
        out_path=out_path,
        color="tab:brown",
    )


def plot_thrust_usage_time(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    t_u = sim["time"][:-1]
    u = sim["control"]
    u_norm = np.linalg.norm(u, axis=1) if u.size else np.zeros(0)
    usage = u_norm / max(float(cfg.thrust_to_mass_ratio_mps2), 1e-9) if u_norm.size else np.zeros(0)
    return _plot_single_series(
        t=t_u,
        y=usage,
        title="Thrust Usage Ratio vs Time",
        y_label="||a_cmd|| / a_max (-)",
        out_path=out_path,
        color="tab:gray",
    )


def plot_keepout_margin_time(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    t = sim["time"]
    margin = np.asarray(sim["range_m"], dtype=float) - float(cfg.keepout_radius_m)
    return _plot_single_series(
        t=t,
        y=margin,
        title="Keepout Margin vs Time",
        y_label="||r|| - r_keepout (m)",
        out_path=out_path,
        color="tab:orange",
        zero_line=True,
    )


def plot_hold_error_time(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    t = sim["time"]
    err = np.asarray(sim["range_m"], dtype=float) - float(cfg.hold_radius_m)
    return _plot_single_series(
        t=t,
        y=err,
        title="Hold-Radius Error vs Time",
        y_label="||r|| - r_hold (m)",
        out_path=out_path,
        color="tab:purple",
        zero_line=True,
    )


def plot_range_rate_time(sim: Dict, out_path: Path) -> str:
    t = np.asarray(sim["time"], dtype=float)
    r = np.asarray(sim["range_m"], dtype=float)
    if t.size > 1:
        r_dot = np.gradient(r, t)
    else:
        r_dot = np.zeros_like(r)
    return _plot_single_series(
        t=t,
        y=r_dot,
        title="Range Rate vs Time",
        y_label="d||r||/dt (m/s)",
        out_path=out_path,
        color="tab:pink",
        zero_line=True,
    )


def plot_body_speed_time(sim: Dict, out_path: Path) -> str:
    t = sim["time"]
    v_b = sim["state_body"][:, 3:6]
    v_b_norm = np.linalg.norm(v_b, axis=1)
    return _plot_single_series(
        t=t,
        y=v_b_norm,
        title="Target-Body Relative Speed vs Time",
        y_label="||v_B|| (m/s)",
        out_path=out_path,
        color="tab:olive",
    )


def plot_control_mode_time(sim: Dict, out_path: Path) -> str:
    _style()
    t_u = np.asarray(sim["time"][:-1], dtype=float)
    mode_hold = np.asarray(
        [1.0 if str(mode).lower() == "hold" else 0.0 for mode in sim.get("control_mode", [])],
        dtype=float,
    )
    n = min(len(t_u), len(mode_hold))

    fig, ax = plt.subplots(figsize=(9, 3.8))
    if n > 0:
        ax.step(t_u[:n], mode_hold[:n], where="post", color="tab:red", lw=1.8)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Hold mode (1/0)")
    ax.set_title("Control Mode Timeline")
    fig.tight_layout()
    return _save_figure(fig, out_path)


def plot_status_flags_time(sim: Dict, out_path: Path) -> str:
    _style()
    t = np.asarray(sim["time"], dtype=float)
    feasible = np.asarray(sim.get("feasible", np.ones(len(t), dtype=bool)), dtype=float)
    monotonic_ok = 1.0 - np.asarray(sim.get("monotonic_violation", np.zeros(len(t), dtype=bool)), dtype=float)
    los_margin = np.asarray(sim["los_margin"], dtype=float)

    n_f = min(len(t), len(feasible))
    n_m = min(len(t), len(monotonic_ok))
    n_l = min(len(t), len(los_margin))

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 6.2), sharex=True)
    if n_f > 0:
        axes[0].step(t[:n_f], feasible[:n_f], where="post", lw=1.8, color="tab:blue", label="Feasible")
    if n_m > 0:
        axes[0].step(t[:n_m], monotonic_ok[:n_m], where="post", lw=1.8, color="tab:green", label="Monotonic")
    axes[0].set_ylim(-0.05, 1.05)
    axes[0].set_ylabel("Flag (1/0)")
    axes[0].set_title("Feasibility and Monotonic Flags")
    axes[0].legend(loc="lower right")

    if n_l > 0:
        axes[1].plot(t[:n_l], los_margin[:n_l], lw=1.8, color="tab:orange")
    axes[1].axhline(0.0, color="k", lw=1.0, ls="--")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("min(b - A x)")
    axes[1].set_title("LOS Margin (for flag interpretation)")
    fig.tight_layout()
    return _save_figure(fig, out_path)


def plot_time_domain_overview(cfg: ScenarioConfig, sim: Dict, out_path: Path) -> str:
    _style()
    t = np.asarray(sim["time"], dtype=float)
    t_u = np.asarray(sim["time"][:-1], dtype=float)
    range_hist = np.asarray(sim["range_m"], dtype=float)
    speed_hist = np.asarray(sim["speed_mps"], dtype=float)
    los_margin = np.asarray(sim["los_margin"], dtype=float)
    phase_err = np.asarray(sim["phase_error_m"], dtype=float)
    mpc_t = np.asarray(sim.get("mpc_cost_time", t_u), dtype=float)
    mpc_j = np.asarray(sim.get("mpc_cost", np.zeros_like(mpc_t)), dtype=float)
    control = np.asarray(sim["control"], dtype=float)
    control_norm = np.linalg.norm(control, axis=1) if control.size else np.zeros(0)
    fuel_use = np.cumsum(control_norm) * cfg.ts if control_norm.size else np.zeros(0)
    hold_err = range_hist - float(cfg.hold_radius_m)

    fig, axes = plt.subplots(4, 2, figsize=(12.4, 11.2), sharex="col")
    ax = axes.ravel()

    n0 = min(len(t), len(range_hist))
    if n0 > 0:
        ax[0].plot(t[:n0], range_hist[:n0], color="tab:blue", lw=1.7)
    ax[0].set_title("Range")
    ax[0].set_ylabel("||r|| (m)")

    n1 = min(len(t), len(speed_hist))
    if n1 > 0:
        ax[1].plot(t[:n1], speed_hist[:n1], color="tab:cyan", lw=1.7)
    ax[1].set_title("Relative Speed")
    ax[1].set_ylabel("||v|| (m/s)")

    n2 = min(len(t), len(los_margin))
    if n2 > 0:
        ax[2].plot(t[:n2], los_margin[:n2], color="tab:orange", lw=1.7)
    ax[2].axhline(0.0, color="k", ls="--", lw=1.0)
    ax[2].set_title("LOS Margin")
    ax[2].set_ylabel("min(b - A x)")

    n3 = min(len(t), len(phase_err))
    if n3 > 0:
        ax[3].plot(t[:n3], phase_err[:n3], color="tab:green", lw=1.7)
    ax[3].set_title("Phase Error")
    ax[3].set_ylabel("(m)")

    n4 = min(len(mpc_t), len(mpc_j))
    if n4 > 0:
        ax[4].plot(mpc_t[:n4], mpc_j[:n4], color="tab:purple", lw=1.7)
    ax[4].set_title("MPC Cost")
    ax[4].set_ylabel("J")

    n5 = min(len(t_u), len(control_norm))
    if n5 > 0:
        ax[5].plot(t_u[:n5], control_norm[:n5], color="tab:brown", lw=1.7)
    ax[5].set_title("Control Norm")
    ax[5].set_ylabel("||a||")

    n6 = min(len(t_u), len(fuel_use))
    if n6 > 0:
        ax[6].plot(t_u[:n6], fuel_use[:n6], color="tab:red", lw=1.7)
    ax[6].set_title("Fuel-Use Proxy")
    ax[6].set_ylabel("Cum. Delta-v")
    ax[6].set_xlabel("Time (s)")

    n7 = min(len(t), len(hold_err))
    if n7 > 0:
        ax[7].plot(t[:n7], hold_err[:n7], color="tab:gray", lw=1.7)
    ax[7].axhline(0.0, color="k", ls="--", lw=1.0)
    ax[7].set_title("Hold Error")
    ax[7].set_ylabel("||r|| - r_hold")
    ax[7].set_xlabel("Time (s)")

    fig.suptitle("Time-Domain Overview Dashboard", y=0.995)
    fig.tight_layout()
    return _save_figure(fig, out_path)


def create_all_plots(cfg: ScenarioConfig, sim: Dict, los: LOS3D) -> Dict[str, str]:
    """Create all required scenario figures.

    All figures are saved as both PNG and PDF.
    """
    paths = {
        "trajectory_body_3d": plot_traj_body_3d(cfg, sim, los, cfg.fig_dir / "trajectory_body_3d.png"),
        "body_position_time": plot_body_position_time(cfg, sim, cfg.fig_dir / "body_position_vs_time.png"),
        "constraint_margin": plot_constraint_margin(cfg, sim, cfg.fig_dir / "constraint_margin_vs_time.png"),
        "control_dv": plot_control_and_dv(cfg, sim, cfg.fig_dir / "control_and_dv_proxy.png"),
        "control_components": plot_control_components(cfg, sim, cfg.fig_dir / "control_components_vs_time.png"),
        "mpc_cost": plot_mpc_cost(cfg, sim, cfg.fig_dir / "mpc_cost_vs_time.png"),
        "mpc_cost_breakdown": plot_mpc_cost_breakdown(cfg, sim, cfg.fig_dir / "mpc_cost_breakdown.png"),
        "range_time": plot_range_time(cfg, sim, cfg.fig_dir / "range_vs_time.png"),
        "state_components_body": plot_state_components_body(cfg, sim, cfg.fig_dir / "state_components_body.png"),
        "trajectory_lvlh_3d": plot_traj_lvlh_3d(cfg, sim, los, cfg.fig_dir / "trajectory_lvlh_3d.png"),
        "speed_time": plot_speed_time(cfg, sim, cfg.fig_dir / "relative_speed_vs_time.png"),
        "fuel_use_time": plot_fuel_use_time(cfg, sim, cfg.fig_dir / "fuel_use_vs_time.png"),
        "phase_sync": plot_phase_sync(cfg, sim, cfg.fig_dir / "phase_sync_metric.png"),
        "state_components_lvlh": plot_state_components_lvlh(cfg, sim, cfg.fig_dir / "state_components_lvlh.png"),
    }

    # Additional plots (full mode adds more)
    if cfg.mode == "full":
        paths.update({
            "range_rate_time": plot_range_rate_time(sim, cfg.fig_dir / "range_rate_vs_time.png"),
            "body_speed_time": plot_body_speed_time(sim, cfg.fig_dir / "body_speed_vs_time.png"),
            "control_norm_time": plot_control_norm_time(sim, cfg.fig_dir / "control_norm_vs_time.png"),
            "thrust_usage_time": plot_thrust_usage_time(cfg, sim, cfg.fig_dir / "thrust_usage_ratio_vs_time.png"),
            "control_mode_time": plot_control_mode_time(sim, cfg.fig_dir / "control_mode_vs_time.png"),
            "lvlh_position_time": plot_lvlh_position_time(sim, cfg.fig_dir / "lvlh_position_components_vs_time.png"),
            "lvlh_velocity_time": plot_lvlh_velocity_time(sim, cfg.fig_dir / "lvlh_velocity_components_vs_time.png"),
            "target_body_position_time": plot_target_body_position_time(sim, cfg.fig_dir / "target_body_position_vs_time.png"),
            "target_body_velocity_time": plot_target_body_velocity_time(sim, cfg.fig_dir / "target_body_velocity_vs_time.png"),
            "target_in_chaser_position_time": plot_target_in_chaser_position_time(
                sim, cfg.fig_dir / "target_in_chaser_position_vs_time.png",
            ),
            "keepout_margin_time": plot_keepout_margin_time(cfg, sim, cfg.fig_dir / "keepout_margin_vs_time.png"),
            "hold_error_time": plot_hold_error_time(cfg, sim, cfg.fig_dir / "hold_error_vs_time.png"),
            "status_flags_time": plot_status_flags_time(sim, cfg.fig_dir / "status_flags_vs_time.png"),
            "time_domain_overview": plot_time_domain_overview(cfg, sim, cfg.fig_dir / "time_domain_overview.png"),
        })

    return paths
