"""Animation utilities for the upgraded 3D approach scenario."""

from __future__ import annotations

from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from IAC.config import ScenarioConfig
from IAC.constraints import LOS3D, rotz


def _draw_docking_axis_3d(ax, theta: float, length: float, lw: float = 1.5) -> None:
    """Draw instantaneous docking axis (+y_B) as dashed black line in LVLH 3D."""
    r_mat = rotz(theta)
    y_b = r_mat @ np.array([0.0, 1.0, 0.0])
    end = length * y_b
    ax.plot([0.0, end[0]], [0.0, end[1]], [0.0, end[2]],
            color="k", ls="--", lw=lw, alpha=0.7)


def _draw_docking_axis_2d(ax, theta: float, length: float,
                          idx_x: int, idx_y: int, lw: float = 1.5) -> None:
    """Draw instantaneous docking axis (+y_B) as dashed black line in 2D projection."""
    r_mat = rotz(theta)
    y_b = r_mat @ np.array([0.0, 1.0, 0.0])
    end = length * y_b
    ax.plot([0.0, end[idx_x]], [0.0, end[idx_y]],
            color="k", ls="--", lw=lw, alpha=0.7)


def _draw_docking_axis_body_2d(ax, length: float,
                               idx_x: int, idx_y: int, lw: float = 1.5) -> None:
    """Draw docking axis in body frame (always +y_B direction)."""
    end = np.array([0.0, length, 0.0])
    ax.plot([0.0, end[idx_x]], [0.0, end[idx_y]],
            color="k", ls="--", lw=lw, alpha=0.7)


def _frame_indices(n: int, stride: int) -> np.ndarray:
    idx = np.arange(0, n, max(1, stride), dtype=int)
    if idx[-1] != n - 1:
        idx = np.hstack((idx, [n - 1]))
    return idx


def _save_animation(ani: animation.FuncAnimation, fig: plt.Figure, cfg: ScenarioConfig, stem: str) -> Dict[str, str | None]:
    cfg.anim_dir.mkdir(parents=True, exist_ok=True)
    cfg.frames_dir.mkdir(parents=True, exist_ok=True)

    mp4_path = cfg.anim_dir / f"{stem}.mp4"
    gif_path = cfg.anim_dir / f"{stem}.gif"
    frames_dir = cfg.frames_dir / stem
    frames_dir.mkdir(parents=True, exist_ok=True)

    out: Dict[str, str | None] = {"mp4": None, "gif": None, "frames_dir": None}

    # User-facing default: prefer GIF to make artifacts easy to view/share.
    try:
        writer = animation.PillowWriter(fps=max(1, cfg.animation_fps // 2))
        ani.save(gif_path, writer=writer)
        out["gif"] = str(gif_path)
        plt.close(fig)
        return out
    except Exception:
        pass

    try:
        writer = animation.FFMpegWriter(fps=cfg.animation_fps, bitrate=2200)
        ani.save(mp4_path, writer=writer)
        out["mp4"] = str(mp4_path)
        plt.close(fig)
        return out
    except Exception:
        pass

    for i, frame in enumerate(ani.new_frame_seq()):
        ani._func(frame)  # noqa: SLF001
        fig.savefig(frames_dir / f"frame_{i:04d}.png", dpi=180)
    out["frames_dir"] = str(frames_dir)
    plt.close(fig)
    return out


def _normalize(vec: np.ndarray, fallback: np.ndarray) -> np.ndarray:
    nrm = float(np.linalg.norm(vec))
    if nrm < 1e-9:
        return fallback.copy()
    return vec / nrm


def _chaser_body_rotation(state_lvlh: np.ndarray) -> np.ndarray:
    """Return a right-handed chaser body frame from LVLH state.

    x_C is aligned with the velocity direction when available, then a
    near-global z direction is used to build y_C, z_C.
    """

    vel = np.asarray(state_lvlh[3:6], dtype=float)
    x_c = _normalize(vel, np.array([1.0, 0.0, 0.0]))

    z_hint = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(x_c, z_hint))) > 0.95:
        z_hint = np.array([0.0, 1.0, 0.0])

    y_c = _normalize(np.cross(z_hint, x_c), np.array([0.0, 1.0, 0.0]))
    z_c = _normalize(np.cross(x_c, y_c), np.array([0.0, 0.0, 1.0]))
    return np.column_stack((x_c, y_c, z_c))


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


def _draw_chaser_box_3d(
    ax,
    center: np.ndarray,
    rot: np.ndarray,
    half_x: float = 0.25,
    half_y: float = 0.20,
    half_z: float = 0.15,
    color: str = "tab:blue",
) -> None:
    """Draw a small cubesat box for the chaser in 3D."""
    verts_local = _docking_box_vertices(half_x, half_y, half_z)
    verts_world = (rot @ verts_local.T).T + center
    _draw_edges_3d(ax, verts_world, _BOX_EDGES, color=color, lw=1.2, alpha=0.9)


def _draw_chaser_box_2d(
    ax,
    center: np.ndarray,
    rot: np.ndarray,
    idx_x: int,
    idx_y: int,
    half_x: float = 0.25,
    half_y: float = 0.20,
    half_z: float = 0.15,
    color: str = "tab:blue",
) -> None:
    """Draw a small cubesat box for the chaser in 2D projection."""
    verts_local = _docking_box_vertices(half_x, half_y, half_z)
    verts_world = (rot @ verts_local.T).T + center
    _draw_edges_2d(ax, verts_world, _BOX_EDGES, idx_x, idx_y, color=color, lw=1.2, alpha=0.9)


def _edges_from_faces(faces_idx: list[list[int]]) -> list[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for face in faces_idx:
        n_face = len(face)
        for i in range(n_face):
            a = int(face[i])
            b = int(face[(i + 1) % n_face])
            edges.add((a, b) if a < b else (b, a))
    return sorted(edges)


_BOX_FACES = [
    [0, 1, 2, 3],
    [4, 5, 6, 7],
    [0, 1, 5, 4],
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 0, 4, 7],
]
_BOX_EDGES = _edges_from_faces(_BOX_FACES)
_DOCKING_FACE_IDX = [3, 2, 6, 7]


def _draw_edges_3d(
    ax,
    verts: np.ndarray,
    edges: list[tuple[int, int]],
    color: str,
    lw: float = 1.2,
    alpha: float = 1.0,
) -> None:
    for a, b in edges:
        ax.plot(
            [verts[a, 0], verts[b, 0]],
            [verts[a, 1], verts[b, 1]],
            [verts[a, 2], verts[b, 2]],
            color=color,
            lw=lw,
            alpha=alpha,
        )


def _draw_edges_2d(
    ax,
    verts: np.ndarray,
    edges: list[tuple[int, int]],
    idx_x: int,
    idx_y: int,
    color: str,
    lw: float = 1.2,
    alpha: float = 1.0,
) -> None:
    for a, b in edges:
        ax.plot(
            [verts[a, idx_x], verts[b, idx_x]],
            [verts[a, idx_y], verts[b, idx_y]],
            color=color,
            lw=lw,
            alpha=alpha,
        )


def _draw_triad(
    ax,
    origin: np.ndarray,
    rot: np.ndarray,
    scale: float,
    lw: float = 3.0,
    alpha: float = 1.0,
    label_prefix: str | None = None,
) -> None:
    colors = ["tab:red", "tab:green", "tab:blue"]
    labels = ["x", "y", "z"]
    for i, unit in enumerate(np.eye(3)):
        vec = rot @ unit
        p2 = origin + scale * vec
        ax.plot(
            [origin[0], p2[0]],
            [origin[1], p2[1]],
            [origin[2], p2[2]],
            color=colors[i],
            lw=lw,
            alpha=alpha,
        )
        if label_prefix is not None:
            ax.text(p2[0], p2[1], p2[2], f"{label_prefix}{labels[i]}", color=colors[i], fontsize=13, fontweight="bold", alpha=alpha)


def _draw_triad_projection(
    ax,
    origin: np.ndarray,
    rot: np.ndarray,
    scale: float,
    idx_x: int,
    idx_y: int,
    lw: float = 2.8,
    alpha: float = 1.0,
    label_prefix: str | None = None,
) -> None:
    colors = ["tab:red", "tab:green", "tab:blue"]
    labels = ["x", "y", "z"]
    for i, unit in enumerate(np.eye(3)):
        vec = rot @ unit
        p2 = origin + scale * vec
        ax.plot(
            [origin[idx_x], p2[idx_x]],
            [origin[idx_y], p2[idx_y]],
            color=colors[i],
            lw=lw,
            alpha=alpha,
        )
        if label_prefix is not None:
            ax.text(p2[idx_x], p2[idx_y], f"{label_prefix}{labels[i]}", color=colors[i], fontsize=13, fontweight="bold", alpha=alpha)


def _draw_corridor_3d(
    ax,
    verts: np.ndarray,
    faces_idx: list[list[int]],
    face_color: str,
    edge_color: str,
    alpha: float,
) -> None:
    faces = [verts[idx] for idx in faces_idx]
    coll = Poly3DCollection(faces, facecolor=face_color, edgecolor=edge_color, alpha=alpha, linewidth=0.8)
    ax.add_collection3d(coll)


def _draw_docking_face(ax, verts: np.ndarray, face_idx: list[int], color: str = "dimgray", alpha: float = 0.35) -> None:
    face = [verts[face_idx]]
    coll = Poly3DCollection(face, facecolor=color, edgecolor=color, alpha=alpha, linewidth=1.0)
    ax.add_collection3d(coll)


def create_approach_animation(cfg: ScenarioConfig, sim: Dict, los: LOS3D, stem: str = "approach") -> Dict[str, str | None]:
    t = sim["time"]
    states = sim["state_lvlh"]
    controls = sim["control"]
    mpc_cost = np.asarray(sim.get("mpc_cost", np.zeros(max(0, len(t) - 1))), dtype=float)
    mpc_t = np.asarray(sim.get("mpc_cost_time", t[:-1]), dtype=float)
    if mpc_cost.size != mpc_t.size:
        n_cost = min(mpc_cost.size, mpc_t.size)
        mpc_cost = mpc_cost[:n_cost]
        mpc_t = mpc_t[:n_cost]
    if mpc_cost.size == 0 and len(t) > 1:
        mpc_cost = np.zeros(len(t) - 1, dtype=float)
        mpc_t = t[:-1]
    frame_ids = _frame_indices(len(t), cfg.animation_stride)

    fig = plt.figure(figsize=(14.5, 9.8))
    grid = fig.add_gridspec(2, 2, wspace=0.08, hspace=0.18)
    ax_l3d = fig.add_subplot(grid[0, 0], projection="3d")
    ax_xy = fig.add_subplot(grid[0, 1])
    ax_yz = fig.add_subplot(grid[1, 0])
    ax_cost = fig.add_subplot(grid[1, 1])

    lim_l = max(8.0, 1.06 * float(np.max(np.abs(states[:, :3]))))
    z_lim_l = 0.45 * lim_l

    faces_idx = los.corridor_faces()
    corridor_edges = _edges_from_faces(faces_idx)

    box_half_x = max(0.55, 1.0 * cfg.hold_radius_m)
    box_half_y = max(0.45, 0.9 * cfg.hold_radius_m)
    box_half_z = max(0.35, 0.7 * cfg.hold_radius_m)
    target_box_b = _docking_box_vertices(box_half_x, box_half_y, box_half_z)

    chaser_box_hx = 0.25 * max(1.0, cfg.hold_radius_m)
    chaser_box_hy = 0.20 * max(1.0, cfg.hold_radius_m)
    chaser_box_hz = 0.15 * max(1.0, cfg.hold_radius_m)

    triad_scale_target_l = max(6.0, 0.08 * lim_l)
    triad_scale_chaser_l = max(5.0, 0.07 * lim_l)
    triad_scale_lvlh_l = max(8.0, 0.10 * lim_l)
    finite_cost = mpc_cost[np.isfinite(mpc_cost)]
    if finite_cost.size > 0:
        cost_max = max(1e-6, 1.08 * float(np.max(finite_cost)))
    else:
        cost_max = 1.0

    def draw_scene(i: int) -> None:
        ax_l3d.cla()
        ax_xy.cla()
        ax_yz.cla()
        ax_cost.cla()

        ax_l3d.set_title("3D LVLH View")
        ax_l3d.set_xlim(-lim_l, lim_l)
        ax_l3d.set_ylim(-lim_l, lim_l)
        ax_l3d.set_zlim(-z_lim_l, z_lim_l)
        ax_l3d.set_xlabel("x_L (m)")
        ax_l3d.set_ylabel("y_L (m)")
        ax_l3d.set_zlabel("z_L (m)")
        ax_l3d.view_init(elev=24, azim=-52)

        ax_xy.set_title("Top View (x-y)")
        ax_xy.set_xlim(-lim_l, lim_l)
        ax_xy.set_ylim(-lim_l, lim_l)
        ax_xy.set_xlabel("x_L (m)")
        ax_xy.set_ylabel("y_L (m)")
        ax_xy.set_aspect("equal", adjustable="box")
        ax_xy.grid(True, alpha=0.28)

        # Requested side view: y is horizontal axis, z is vertical axis.
        ax_yz.set_title("Side View (y-z)")
        ax_yz.set_xlim(-lim_l, lim_l)
        ax_yz.set_ylim(-z_lim_l, z_lim_l)
        ax_yz.set_xlabel("y_L (m)")
        ax_yz.set_ylabel("z_L (m)")
        ax_yz.set_aspect("equal", adjustable="box")
        ax_yz.grid(True, alpha=0.28)

        ax_cost.set_title("MPC Cost")
        ax_cost.set_xlim(float(t[0]), float(t[-1]))
        ax_cost.set_ylim(0.0, cost_max)
        ax_cost.set_xlabel("Time (s)")
        ax_cost.set_ylabel("Cost J")
        ax_cost.grid(True, alpha=0.30)

        xyz_l = states[: i + 1, :3]
        chaser_l = states[i, :3]

        theta = los.angle(float(t[i]))
        r_body_to_lvlh = rotz(theta)
        r_chaser_to_lvlh = _chaser_body_rotation(states[i])

        verts_l_corridor = los.corridor_vertices_lvlh(float(t[i]), y_far=cfg.corridor_far_y_m)
        verts_l_box = (r_body_to_lvlh @ target_box_b.T).T

        _draw_corridor_3d(
            ax_l3d,
            verts_l_corridor,
            faces_idx=faces_idx,
            face_color="tab:orange",
            edge_color="tab:orange",
            alpha=0.11,
        )

        _draw_edges_3d(ax_l3d, verts_l_box, _BOX_EDGES, color="k", lw=1.4, alpha=0.95)
        _draw_docking_face(ax_l3d, verts_l_box, _DOCKING_FACE_IDX)

        _draw_edges_2d(ax_xy, verts_l_corridor, corridor_edges, 0, 1, color="tab:orange", lw=1.0, alpha=0.75)
        _draw_edges_2d(ax_yz, verts_l_corridor, corridor_edges, 1, 2, color="tab:orange", lw=1.0, alpha=0.75)
        _draw_edges_2d(ax_xy, verts_l_box, _BOX_EDGES, 0, 1, color="k", lw=1.3, alpha=0.95)
        _draw_edges_2d(ax_yz, verts_l_box, _BOX_EDGES, 1, 2, color="k", lw=1.3, alpha=0.95)

        ax_l3d.plot(xyz_l[:, 0], xyz_l[:, 1], xyz_l[:, 2], color="tab:blue", lw=2.0)
        ax_xy.plot(xyz_l[:, 0], xyz_l[:, 1], color="tab:blue", lw=2.0)
        ax_yz.plot(xyz_l[:, 1], xyz_l[:, 2], color="tab:blue", lw=2.0)

        # Chaser cubesat box
        _draw_chaser_box_3d(ax_l3d, chaser_l, r_chaser_to_lvlh,
                            half_x=chaser_box_hx, half_y=chaser_box_hy, half_z=chaser_box_hz)
        _draw_chaser_box_2d(ax_xy, chaser_l, r_chaser_to_lvlh, 0, 1,
                            half_x=chaser_box_hx, half_y=chaser_box_hy, half_z=chaser_box_hz)
        _draw_chaser_box_2d(ax_yz, chaser_l, r_chaser_to_lvlh, 1, 2,
                            half_x=chaser_box_hx, half_y=chaser_box_hy, half_z=chaser_box_hz)

        # LVLH inertial triad (fixed at origin)
        _draw_triad(
            ax_l3d, np.zeros(3), np.eye(3),
            scale=triad_scale_lvlh_l, lw=2.2, alpha=0.85, label_prefix="L",
        )
        _draw_triad_projection(
            ax_xy, np.zeros(3), np.eye(3),
            scale=triad_scale_lvlh_l, idx_x=0, idx_y=1,
            lw=2.0, alpha=0.85, label_prefix="L",
        )
        _draw_triad_projection(
            ax_yz, np.zeros(3), np.eye(3),
            scale=triad_scale_lvlh_l, idx_x=1, idx_y=2,
            lw=2.0, alpha=0.85, label_prefix="L",
        )

        # Target body triad (rotating)
        _draw_triad(
            ax_l3d,
            np.zeros(3),
            r_body_to_lvlh,
            scale=triad_scale_target_l,
            lw=1.9,
            alpha=0.95,
            label_prefix="B",
        )
        _draw_triad(
            ax_l3d,
            chaser_l,
            r_chaser_to_lvlh,
            scale=triad_scale_chaser_l,
            lw=1.3,
            alpha=0.95,
            label_prefix="C",
        )

        _draw_triad_projection(
            ax_xy,
            np.zeros(3),
            r_body_to_lvlh,
            scale=triad_scale_target_l,
            idx_x=0,
            idx_y=1,
            lw=1.6,
            alpha=0.95,
            label_prefix="B",
        )
        _draw_triad_projection(
            ax_xy,
            chaser_l,
            r_chaser_to_lvlh,
            scale=triad_scale_chaser_l,
            idx_x=0,
            idx_y=1,
            lw=1.3,
            alpha=0.95,
            label_prefix="C",
        )
        _draw_triad_projection(
            ax_yz,
            np.zeros(3),
            r_body_to_lvlh,
            scale=triad_scale_target_l,
            idx_x=1,
            idx_y=2,
            lw=1.6,
            alpha=0.95,
            label_prefix="B",
        )
        _draw_triad_projection(
            ax_yz,
            chaser_l,
            r_chaser_to_lvlh,
            scale=triad_scale_chaser_l,
            idx_x=1,
            idx_y=2,
            lw=1.3,
            alpha=0.95,
            label_prefix="C",
        )

        # Instantaneous docking axis (+y_B)
        _draw_docking_axis_3d(ax_l3d, theta, length=lim_l * 0.9)
        _draw_docking_axis_2d(ax_xy, theta, length=lim_l * 0.9, idx_x=0, idx_y=1)
        _draw_docking_axis_2d(ax_yz, theta, length=lim_l * 0.9, idx_x=1, idx_y=2)

        if mpc_cost.size > 0:
            k_cost = min(i, mpc_cost.size)
            if k_cost > 0:
                ax_cost.plot(mpc_t[:k_cost], mpc_cost[:k_cost], color="tab:purple", lw=2.0)
                ax_cost.plot([mpc_t[k_cost - 1]], [mpc_cost[k_cost - 1]], marker="o", color="tab:purple", ms=5)

        a_mag = 0.0 if i == 0 or controls.size == 0 else float(np.linalg.norm(controls[min(i - 1, len(controls) - 1)]))
        txt = f"t={t[i]:6.1f}s  ||r||={np.linalg.norm(chaser_l):6.2f}m  ||a||={a_mag:5.3f}m/s^2"
        ax_l3d.text2D(0.02, 0.96, txt, transform=ax_l3d.transAxes)
        ax_xy.text(0.02, 0.96, txt, transform=ax_xy.transAxes, va="top")
        ax_yz.text(0.02, 0.96, txt, transform=ax_yz.transAxes, va="top")
        ax_cost.text(0.02, 0.96, txt, transform=ax_cost.transAxes, va="top")

    def update(frame_idx: int):
        i = int(frame_ids[frame_idx])
        draw_scene(i)
        return []

    draw_scene(int(frame_ids[0]))
    fig.suptitle(f"Spiral Approach, Multi-View Geometry, and MPC Cost\n({cfg.scenario_tag})", fontsize=12)

    preview_path = cfg.fig_dir / f"{stem}_layout.png"
    preview_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(preview_path, dpi=220)

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(frame_ids),
        interval=1000.0 / max(1, cfg.animation_fps),
        blit=False,
    )
    out = _save_animation(ani, fig, cfg, stem=stem)
    out["layout_png"] = str(preview_path)
    return out


def _finite_limits(y: np.ndarray, pad_frac: float = 0.08, min_span: float = 1e-3) -> tuple[float, float]:
    vals = np.asarray(y, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return -1.0, 1.0
    y_min = float(np.min(vals))
    y_max = float(np.max(vals))
    span = max(y_max - y_min, min_span)
    pad = pad_frac * span
    return y_min - pad, y_max + pad


def _create_time_series_animation(
    cfg: ScenarioConfig,
    t_anim: np.ndarray,
    t_series: np.ndarray,
    y_series: np.ndarray,
    stem: str,
    title: str,
    y_label: str,
    color: str = "tab:blue",
) -> Dict[str, str | None]:
    n = min(len(t_series), len(y_series))
    if n == 0:
        return {"mp4": None, "gif": None, "frames_dir": None}

    t_s = np.asarray(t_series[:n], dtype=float)
    y_s = np.asarray(y_series[:n], dtype=float)
    frame_ids = _frame_indices(len(t_anim), cfg.animation_stride)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    y_lo, y_hi = _finite_limits(y_s)

    def draw_scene(i: int) -> None:
        ax.cla()
        t_now = float(t_anim[i])
        k = int(np.searchsorted(t_s, t_now, side="right"))

        ax.set_title(title)
        ax.set_xlim(float(t_anim[0]), float(t_anim[-1]))
        ax.set_ylim(y_lo, y_hi)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(y_label)
        ax.grid(True, alpha=0.3)

        if k > 0:
            ax.plot(t_s[:k], y_s[:k], color=color, lw=2.0)
            ax.plot([t_s[k - 1]], [y_s[k - 1]], marker="o", color=color, ms=5)
            txt = f"t={t_now:6.1f}s  value={y_s[k - 1]:.5g}"
        else:
            txt = f"t={t_now:6.1f}s"
        ax.text(0.02, 0.96, txt, transform=ax.transAxes, va="top")

    def update(frame_idx: int):
        i = int(frame_ids[frame_idx])
        draw_scene(i)
        return []

    draw_scene(int(frame_ids[0]))
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(frame_ids),
        interval=1000.0 / max(1, cfg.animation_fps),
        blit=False,
    )
    return _save_animation(ani, fig, cfg, stem=stem)


def create_lvlh_3d_animation(cfg: ScenarioConfig, sim: Dict, los: LOS3D, stem: str = "approach_lvlh_3d") -> Dict[str, str | None]:
    t = sim["time"]
    states = sim["state_lvlh"]
    controls = sim["control"]
    frame_ids = _frame_indices(len(t), cfg.animation_stride)

    fig = plt.figure(figsize=(8.8, 7.0))
    ax = fig.add_subplot(111, projection="3d")

    lim_l = max(8.0, 1.06 * float(np.max(np.abs(states[:, :3]))))
    z_lim_l = 0.45 * lim_l
    faces_idx = los.corridor_faces()

    box_half_x = max(0.55, 1.0 * cfg.hold_radius_m)
    box_half_y = max(0.45, 0.9 * cfg.hold_radius_m)
    box_half_z = max(0.35, 0.7 * cfg.hold_radius_m)
    target_box_b = _docking_box_vertices(box_half_x, box_half_y, box_half_z)

    chaser_box_hx = 0.25 * max(1.0, cfg.hold_radius_m)
    chaser_box_hy = 0.20 * max(1.0, cfg.hold_radius_m)
    chaser_box_hz = 0.15 * max(1.0, cfg.hold_radius_m)

    triad_scale_target_l = max(6.0, 0.08 * lim_l)
    triad_scale_chaser_l = max(5.0, 0.07 * lim_l)
    triad_scale_lvlh_l = max(8.0, 0.10 * lim_l)

    def draw_scene(i: int) -> None:
        ax.cla()
        ax.set_title(f"LVLH 3D View — {cfg.scenario_tag}")
        ax.set_xlim(-lim_l, lim_l)
        ax.set_ylim(-lim_l, lim_l)
        ax.set_zlim(-z_lim_l, z_lim_l)
        ax.set_xlabel("x_L (m)")
        ax.set_ylabel("y_L (m)")
        ax.set_zlabel("z_L (m)")
        ax.view_init(elev=24, azim=-52)

        xyz_l = states[: i + 1, :3]
        chaser_l = states[i, :3]
        theta = los.angle(float(t[i]))
        r_body_to_lvlh = rotz(theta)
        r_chaser_to_lvlh = _chaser_body_rotation(states[i])

        verts_l_corridor = los.corridor_vertices_lvlh(float(t[i]), y_far=cfg.corridor_far_y_m)
        verts_l_box = (r_body_to_lvlh @ target_box_b.T).T
        _draw_corridor_3d(
            ax,
            verts_l_corridor,
            faces_idx=faces_idx,
            face_color="tab:orange",
            edge_color="tab:orange",
            alpha=0.11,
        )
        _draw_edges_3d(ax, verts_l_box, _BOX_EDGES, color="k", lw=1.4, alpha=0.95)
        _draw_docking_face(ax, verts_l_box, _DOCKING_FACE_IDX)

        ax.plot(xyz_l[:, 0], xyz_l[:, 1], xyz_l[:, 2], color="tab:blue", lw=2.0)

        # Chaser cubesat box
        _draw_chaser_box_3d(ax, chaser_l, r_chaser_to_lvlh,
                            half_x=chaser_box_hx, half_y=chaser_box_hy, half_z=chaser_box_hz)

        # LVLH inertial triad (fixed at origin)
        _draw_triad(
            ax, np.zeros(3), np.eye(3),
            scale=triad_scale_lvlh_l, lw=2.2, alpha=0.85, label_prefix="L",
        )
        # Target body triad (rotating)
        _draw_triad(
            ax,
            np.zeros(3),
            r_body_to_lvlh,
            scale=triad_scale_target_l,
            lw=1.9,
            alpha=0.95,
            label_prefix="B",
        )
        # Chaser body triad
        _draw_triad(
            ax,
            chaser_l,
            r_chaser_to_lvlh,
            scale=triad_scale_chaser_l,
            lw=1.3,
            alpha=0.95,
            label_prefix="C",
        )

        # Instantaneous docking axis (+y_B)
        _draw_docking_axis_3d(ax, theta, length=lim_l * 0.9)

        a_mag = 0.0 if i == 0 or controls.size == 0 else float(np.linalg.norm(controls[min(i - 1, len(controls) - 1)]))
        txt = f"t={t[i]:6.1f}s  ||r||={np.linalg.norm(chaser_l):6.2f}m  ||a||={a_mag:5.3f}m/s^2"
        ax.text2D(0.02, 0.96, txt, transform=ax.transAxes)

    def update(frame_idx: int):
        i = int(frame_ids[frame_idx])
        draw_scene(i)
        return []

    draw_scene(int(frame_ids[0]))
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(frame_ids),
        interval=1000.0 / max(1, cfg.animation_fps),
        blit=False,
    )
    return _save_animation(ani, fig, cfg, stem=stem)


def create_lvlh_plane_animation(
    cfg: ScenarioConfig,
    sim: Dict,
    los: LOS3D,
    plane: str,
    stem: str,
) -> Dict[str, str | None]:
    t = sim["time"]
    states = sim["state_lvlh"]
    controls = sim["control"]
    frame_ids = _frame_indices(len(t), cfg.animation_stride)

    plane_key = plane.lower().strip()
    tag = cfg.scenario_tag
    if plane_key == "xy":
        idx_x, idx_y = 0, 1
        title = f"LVLH Top View (x-y) — {tag}"
        x_label, y_label = "x_L (m)", "y_L (m)"
    elif plane_key == "yz":
        idx_x, idx_y = 1, 2
        title = f"LVLH Side View (y-z) — {tag}"
        x_label, y_label = "y_L (m)", "z_L (m)"
    elif plane_key == "xz":
        idx_x, idx_y = 0, 2
        title = f"LVLH Side View (x-z) — {tag}"
        x_label, y_label = "x_L (m)", "z_L (m)"
    else:
        raise ValueError(f"Unknown plane '{plane}'. Expected one of: xy, yz, xz.")

    fig, ax = plt.subplots(figsize=(7.4, 6.2))

    lim_l = max(8.0, 1.06 * float(np.max(np.abs(states[:, :3]))))
    z_lim_l = 0.45 * lim_l
    faces_idx = los.corridor_faces()
    corridor_edges = _edges_from_faces(faces_idx)

    box_half_x = max(0.55, 1.0 * cfg.hold_radius_m)
    box_half_y = max(0.45, 0.9 * cfg.hold_radius_m)
    box_half_z = max(0.35, 0.7 * cfg.hold_radius_m)
    target_box_b = _docking_box_vertices(box_half_x, box_half_y, box_half_z)

    chaser_box_hx = 0.25 * max(1.0, cfg.hold_radius_m)
    chaser_box_hy = 0.20 * max(1.0, cfg.hold_radius_m)
    chaser_box_hz = 0.15 * max(1.0, cfg.hold_radius_m)

    triad_scale_target_l = max(6.0, 0.08 * lim_l)
    triad_scale_chaser_l = max(5.0, 0.07 * lim_l)
    triad_scale_lvlh_l = max(8.0, 0.10 * lim_l)

    def draw_scene(i: int) -> None:
        ax.cla()
        ax.set_title(title)
        if idx_x == 2 or idx_y == 2:
            x_lim = (-z_lim_l, z_lim_l) if idx_x == 2 else (-lim_l, lim_l)
            y_lim = (-z_lim_l, z_lim_l) if idx_y == 2 else (-lim_l, lim_l)
        else:
            x_lim = (-lim_l, lim_l)
            y_lim = (-lim_l, lim_l)
        ax.set_xlim(*x_lim)
        ax.set_ylim(*y_lim)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.28)

        xyz_l = states[: i + 1, :3]
        chaser_l = states[i, :3]
        theta = los.angle(float(t[i]))
        r_body_to_lvlh = rotz(theta)
        r_chaser_to_lvlh = _chaser_body_rotation(states[i])

        verts_l_corridor = los.corridor_vertices_lvlh(float(t[i]), y_far=cfg.corridor_far_y_m)
        verts_l_box = (r_body_to_lvlh @ target_box_b.T).T

        _draw_edges_2d(ax, verts_l_corridor, corridor_edges, idx_x, idx_y, color="tab:orange", lw=1.0, alpha=0.75)
        _draw_edges_2d(ax, verts_l_box, _BOX_EDGES, idx_x, idx_y, color="k", lw=1.3, alpha=0.95)

        ax.plot(xyz_l[:, idx_x], xyz_l[:, idx_y], color="tab:blue", lw=2.0)

        # Chaser cubesat box
        _draw_chaser_box_2d(ax, chaser_l, r_chaser_to_lvlh, idx_x, idx_y,
                            half_x=chaser_box_hx, half_y=chaser_box_hy, half_z=chaser_box_hz)

        # LVLH inertial triad (fixed at origin)
        _draw_triad_projection(
            ax, np.zeros(3), np.eye(3),
            scale=triad_scale_lvlh_l, idx_x=idx_x, idx_y=idx_y,
            lw=2.0, alpha=0.85, label_prefix="L",
        )
        # Target body triad (rotating)
        _draw_triad_projection(
            ax,
            np.zeros(3),
            r_body_to_lvlh,
            scale=triad_scale_target_l,
            idx_x=idx_x,
            idx_y=idx_y,
            lw=1.6,
            alpha=0.95,
            label_prefix="B",
        )
        # Chaser body triad
        _draw_triad_projection(
            ax,
            chaser_l,
            r_chaser_to_lvlh,
            scale=triad_scale_chaser_l,
            idx_x=idx_x,
            idx_y=idx_y,
            lw=1.3,
            alpha=0.95,
            label_prefix="C",
        )

        # Instantaneous docking axis (+y_B)
        _draw_docking_axis_2d(ax, theta, length=lim_l * 0.9, idx_x=idx_x, idx_y=idx_y)

        a_mag = 0.0 if i == 0 or controls.size == 0 else float(np.linalg.norm(controls[min(i - 1, len(controls) - 1)]))
        txt = f"t={t[i]:6.1f}s  ||r||={np.linalg.norm(chaser_l):6.2f}m  ||a||={a_mag:5.3f}m/s^2"
        ax.text(0.02, 0.96, txt, transform=ax.transAxes, va="top")

    def update(frame_idx: int):
        i = int(frame_ids[frame_idx])
        draw_scene(i)
        return []

    draw_scene(int(frame_ids[0]))
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(frame_ids),
        interval=1000.0 / max(1, cfg.animation_fps),
        blit=False,
    )
    return _save_animation(ani, fig, cfg, stem=stem)


def create_target_body_frame_animation(
    cfg: ScenarioConfig,
    sim: Dict,
    los: LOS3D,
    stem: str = "approach_target_body_frame",
) -> Dict[str, str | None]:
    t = sim["time"]
    states = sim["state_lvlh"]
    states_body = sim.get("state_body")
    if states_body is None:
        states_body = los.body_frame_history(states, t)
    controls = sim["control"]
    frame_ids = _frame_indices(len(t), cfg.animation_stride)

    fig = plt.figure(figsize=(8.8, 7.0))
    ax = fig.add_subplot(111, projection="3d")

    lim_b = max(8.0, 1.06 * float(np.max(np.abs(states_body[:, :3]))))
    z_lim_b = 0.45 * lim_b
    y_min_b = min(-0.2 * lim_b, 1.06 * float(np.min(states_body[:, 1])))
    y_max_b = max(0.8 * lim_b, 1.06 * float(np.max(states_body[:, 1])))

    faces_idx = los.corridor_faces()
    verts_b_corridor = los.corridor_vertices_body(y_far=cfg.corridor_far_y_m)

    box_half_x = max(0.55, 1.0 * cfg.hold_radius_m)
    box_half_y = max(0.45, 0.9 * cfg.hold_radius_m)
    box_half_z = max(0.35, 0.7 * cfg.hold_radius_m)
    target_box_b = _docking_box_vertices(box_half_x, box_half_y, box_half_z)

    chaser_box_hx = 0.25 * max(1.0, cfg.hold_radius_m)
    chaser_box_hy = 0.20 * max(1.0, cfg.hold_radius_m)
    chaser_box_hz = 0.15 * max(1.0, cfg.hold_radius_m)

    triad_scale_target_b = max(6.0, 0.08 * lim_b)
    triad_scale_chaser_b = max(5.0, 0.07 * lim_b)
    triad_scale_lvlh_b = max(8.0, 0.10 * lim_b)

    def draw_scene(i: int) -> None:
        ax.cla()
        ax.set_title(f"Target Body Frame — {cfg.scenario_tag}")
        ax.set_xlim(-lim_b, lim_b)
        ax.set_ylim(y_min_b, y_max_b)
        ax.set_zlim(-z_lim_b, z_lim_b)
        ax.set_xlabel("x_B (m)")
        ax.set_ylabel("y_B (m)")
        ax.set_zlabel("z_B (m)")
        ax.view_init(elev=22, azim=-54)

        xyz_b = states_body[: i + 1, :3]
        chaser_b = states_body[i, :3]
        theta = los.angle(float(t[i]))
        r_chaser_to_lvlh = _chaser_body_rotation(states[i])
        r_chaser_to_body = rotz(-theta) @ r_chaser_to_lvlh
        r_lvlh_to_body = rotz(-theta)

        _draw_corridor_3d(
            ax,
            verts_b_corridor,
            faces_idx=faces_idx,
            face_color="tab:red",
            edge_color="tab:red",
            alpha=0.10,
        )
        _draw_edges_3d(ax, target_box_b, _BOX_EDGES, color="k", lw=1.4, alpha=0.95)
        _draw_docking_face(ax, target_box_b, _DOCKING_FACE_IDX)

        ax.plot(xyz_b[:, 0], xyz_b[:, 1], xyz_b[:, 2], color="tab:green", lw=2.0)

        # Chaser cubesat box
        _draw_chaser_box_3d(ax, chaser_b, r_chaser_to_body,
                            half_x=chaser_box_hx, half_y=chaser_box_hy, half_z=chaser_box_hz)

        # LVLH inertial triad (rotates in body view as rotz(-theta))
        _draw_triad(
            ax, np.zeros(3), r_lvlh_to_body,
            scale=triad_scale_lvlh_b, lw=2.2, alpha=0.85, label_prefix="L",
        )
        # Target body triad (fixed in body frame = identity)
        _draw_triad(
            ax,
            np.zeros(3),
            np.eye(3),
            scale=triad_scale_target_b,
            lw=1.9,
            alpha=0.95,
            label_prefix="B",
        )
        # Chaser body triad
        _draw_triad(
            ax,
            chaser_b,
            r_chaser_to_body,
            scale=triad_scale_chaser_b,
            lw=1.3,
            alpha=0.95,
            label_prefix="C",
        )

        # Docking axis in body frame is always +y_B
        ax.plot([0, 0], [0, lim_b * 0.9], [0, 0],
                color="k", ls="--", lw=1.5, alpha=0.7)

        a_mag = 0.0 if i == 0 or controls.size == 0 else float(np.linalg.norm(controls[min(i - 1, len(controls) - 1)]))
        txt = f"t={t[i]:6.1f}s  ||r||={np.linalg.norm(states[i, :3]):6.2f}m  ||a||={a_mag:5.3f}m/s^2"
        ax.text2D(0.02, 0.96, txt, transform=ax.transAxes)

    def update(frame_idx: int):
        i = int(frame_ids[frame_idx])
        draw_scene(i)
        return []

    draw_scene(int(frame_ids[0]))
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(frame_ids),
        interval=1000.0 / max(1, cfg.animation_fps),
        blit=False,
    )
    return _save_animation(ani, fig, cfg, stem=stem)


def create_body_plane_animation(
    cfg: ScenarioConfig,
    sim: Dict,
    los: LOS3D,
    plane: str,
    stem: str,
) -> Dict[str, str | None]:
    t = sim["time"]
    states = sim["state_lvlh"]
    states_body = sim.get("state_body")
    if states_body is None:
        states_body = los.body_frame_history(states, t)
    controls = sim["control"]
    frame_ids = _frame_indices(len(t), cfg.animation_stride)

    plane_key = plane.lower().strip()
    tag = cfg.scenario_tag
    if plane_key == "xy":
        idx_x, idx_y = 0, 1
        title = f"Target Body Top View (x-y) — {tag}"
        x_label, y_label = "x_B (m)", "y_B (m)"
    elif plane_key == "yz":
        idx_x, idx_y = 1, 2
        title = f"Target Body Side View (y-z) — {tag}"
        x_label, y_label = "y_B (m)", "z_B (m)"
    elif plane_key == "xz":
        idx_x, idx_y = 0, 2
        title = f"Target Body Side View (x-z) — {tag}"
        x_label, y_label = "x_B (m)", "z_B (m)"
    else:
        raise ValueError(f"Unknown plane '{plane}'. Expected one of: xy, yz, xz.")

    fig, ax = plt.subplots(figsize=(7.4, 6.2))

    lim_b = max(8.0, 1.06 * float(np.max(np.abs(states_body[:, :3]))))
    z_lim_b = 0.45 * lim_b
    y_min_b = min(-0.2 * lim_b, 1.06 * float(np.min(states_body[:, 1])))
    y_max_b = max(0.8 * lim_b, 1.06 * float(np.max(states_body[:, 1])))
    faces_idx = los.corridor_faces()
    corridor_edges = _edges_from_faces(faces_idx)
    verts_b_corridor = los.corridor_vertices_body(y_far=cfg.corridor_far_y_m)

    box_half_x = max(0.55, 1.0 * cfg.hold_radius_m)
    box_half_y = max(0.45, 0.9 * cfg.hold_radius_m)
    box_half_z = max(0.35, 0.7 * cfg.hold_radius_m)
    target_box_b = _docking_box_vertices(box_half_x, box_half_y, box_half_z)

    chaser_box_hx = 0.25 * max(1.0, cfg.hold_radius_m)
    chaser_box_hy = 0.20 * max(1.0, cfg.hold_radius_m)
    chaser_box_hz = 0.15 * max(1.0, cfg.hold_radius_m)

    triad_scale_target_b = max(6.0, 0.08 * lim_b)
    triad_scale_chaser_b = max(5.0, 0.07 * lim_b)
    triad_scale_lvlh_b = max(8.0, 0.10 * lim_b)

    def _axis_limit(axis_idx: int) -> tuple[float, float]:
        if axis_idx == 2:
            return (-z_lim_b, z_lim_b)
        if axis_idx == 1:
            return (y_min_b, y_max_b)
        return (-lim_b, lim_b)

    def draw_scene(i: int) -> None:
        ax.cla()
        ax.set_title(title)
        ax.set_xlim(*_axis_limit(idx_x))
        ax.set_ylim(*_axis_limit(idx_y))
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.28)

        xyz_b = states_body[: i + 1, :3]
        chaser_b = states_body[i, :3]
        theta = los.angle(float(t[i]))
        r_chaser_to_lvlh = _chaser_body_rotation(states[i])
        r_chaser_to_body = rotz(-theta) @ r_chaser_to_lvlh
        r_lvlh_to_body = rotz(-theta)

        _draw_edges_2d(ax, verts_b_corridor, corridor_edges, idx_x, idx_y, color="tab:red", lw=1.0, alpha=0.75)
        _draw_edges_2d(ax, target_box_b, _BOX_EDGES, idx_x, idx_y, color="k", lw=1.3, alpha=0.95)

        ax.plot(xyz_b[:, idx_x], xyz_b[:, idx_y], color="tab:green", lw=2.0)

        # Chaser cubesat box
        _draw_chaser_box_2d(ax, chaser_b, r_chaser_to_body, idx_x, idx_y,
                            half_x=chaser_box_hx, half_y=chaser_box_hy, half_z=chaser_box_hz)

        # LVLH inertial triad (rotates in body view)
        _draw_triad_projection(
            ax, np.zeros(3), r_lvlh_to_body,
            scale=triad_scale_lvlh_b, idx_x=idx_x, idx_y=idx_y,
            lw=2.0, alpha=0.85, label_prefix="L",
        )
        # Target body triad (fixed = identity in body frame)
        _draw_triad_projection(
            ax,
            np.zeros(3),
            np.eye(3),
            scale=triad_scale_target_b,
            idx_x=idx_x,
            idx_y=idx_y,
            lw=1.6,
            alpha=0.95,
            label_prefix="B",
        )
        # Chaser body triad
        _draw_triad_projection(
            ax,
            chaser_b,
            r_chaser_to_body,
            scale=triad_scale_chaser_b,
            idx_x=idx_x,
            idx_y=idx_y,
            lw=1.3,
            alpha=0.95,
            label_prefix="C",
        )

        # Docking axis in body frame is always +y_B
        _draw_docking_axis_body_2d(ax, length=lim_b * 0.9, idx_x=idx_x, idx_y=idx_y)

        a_mag = 0.0 if i == 0 or controls.size == 0 else float(np.linalg.norm(controls[min(i - 1, len(controls) - 1)]))
        txt = f"t={t[i]:6.1f}s  ||r||={np.linalg.norm(states[i, :3]):6.2f}m  ||a||={a_mag:5.3f}m/s^2"
        ax.text(0.02, 0.96, txt, transform=ax.transAxes, va="top")

    def update(frame_idx: int):
        i = int(frame_ids[frame_idx])
        draw_scene(i)
        return []

    draw_scene(int(frame_ids[0]))
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(frame_ids),
        interval=1000.0 / max(1, cfg.animation_fps),
        blit=False,
    )
    return _save_animation(ani, fig, cfg, stem=stem)


def create_chaser_body_frame_animation(
    cfg: ScenarioConfig,
    sim: Dict,
    los: LOS3D,
    stem: str = "approach_chaser_body_frame",
) -> Dict[str, str | None]:
    t = sim["time"]
    states = sim["state_lvlh"]
    controls = sim["control"]
    frame_ids = _frame_indices(len(t), cfg.animation_stride)

    target_in_chaser = np.zeros((len(t), 3), dtype=float)
    for i in range(len(t)):
        r_chaser_to_lvlh = _chaser_body_rotation(states[i])
        r_lvlh_to_chaser = r_chaser_to_lvlh.T
        target_in_chaser[i] = r_lvlh_to_chaser @ (-states[i, :3])

    fig = plt.figure(figsize=(8.8, 7.0))
    ax = fig.add_subplot(111, projection="3d")

    lim_c = max(8.0, 1.12 * float(np.max(np.abs(target_in_chaser))))
    z_lim_c = 0.55 * lim_c

    faces_idx = los.corridor_faces()
    box_half_x = max(0.55, 1.0 * cfg.hold_radius_m)
    box_half_y = max(0.45, 0.9 * cfg.hold_radius_m)
    box_half_z = max(0.35, 0.7 * cfg.hold_radius_m)
    target_box_b = _docking_box_vertices(box_half_x, box_half_y, box_half_z)

    chaser_box_hx = 0.25 * max(1.0, cfg.hold_radius_m)
    chaser_box_hy = 0.20 * max(1.0, cfg.hold_radius_m)
    chaser_box_hz = 0.15 * max(1.0, cfg.hold_radius_m)

    triad_scale = max(6.0, 0.08 * lim_c)
    triad_scale_lvlh = max(8.0, 0.10 * lim_c)

    def draw_scene(i: int) -> None:
        ax.cla()
        ax.set_title(f"Chaser Body Frame — {cfg.scenario_tag}")
        ax.set_xlim(-lim_c, lim_c)
        ax.set_ylim(-lim_c, lim_c)
        ax.set_zlim(-z_lim_c, z_lim_c)
        ax.set_xlabel("x_C (m)")
        ax.set_ylabel("y_C (m)")
        ax.set_zlabel("z_C (m)")
        ax.view_init(elev=24, azim=-50)

        chaser_l = states[i, :3]
        theta = los.angle(float(t[i]))
        r_body_to_lvlh = rotz(theta)
        r_chaser_to_lvlh = _chaser_body_rotation(states[i])
        r_lvlh_to_chaser = r_chaser_to_lvlh.T

        verts_l_corridor = los.corridor_vertices_lvlh(float(t[i]), y_far=cfg.corridor_far_y_m)
        verts_l_box = (r_body_to_lvlh @ target_box_b.T).T
        verts_c_corridor = (r_lvlh_to_chaser @ (verts_l_corridor - chaser_l).T).T
        verts_c_box = (r_lvlh_to_chaser @ (verts_l_box - chaser_l).T).T

        _draw_corridor_3d(
            ax,
            verts_c_corridor,
            faces_idx=faces_idx,
            face_color="tab:orange",
            edge_color="tab:orange",
            alpha=0.10,
        )
        _draw_edges_3d(ax, verts_c_box, _BOX_EDGES, color="k", lw=1.4, alpha=0.95)
        _draw_docking_face(ax, verts_c_box, _DOCKING_FACE_IDX)

        xyz_c = target_in_chaser[: i + 1]
        tgt_c = target_in_chaser[i]
        ax.plot(xyz_c[:, 0], xyz_c[:, 1], xyz_c[:, 2], color="tab:blue", lw=2.0)

        # Chaser cubesat box at origin (chaser is at origin in its own frame)
        _draw_chaser_box_3d(ax, np.zeros(3), np.eye(3),
                            half_x=chaser_box_hx, half_y=chaser_box_hy, half_z=chaser_box_hz)
        # Target cubesat box at target position in chaser frame (already drawn as wire box above)

        # LVLH inertial triad in chaser frame
        _draw_triad(ax, np.zeros(3), r_lvlh_to_chaser,
                    scale=triad_scale_lvlh, lw=2.2, alpha=0.85, label_prefix="L")

        r_target_to_chaser = r_lvlh_to_chaser @ r_body_to_lvlh
        _draw_triad(ax, np.zeros(3), np.eye(3), scale=triad_scale, lw=1.9, alpha=0.95, label_prefix="C")
        _draw_triad(ax, tgt_c, r_target_to_chaser, scale=triad_scale, lw=1.3, alpha=0.95, label_prefix="B")

        # Docking axis (+y_B) in chaser frame
        y_b_chaser = r_lvlh_to_chaser @ r_body_to_lvlh @ np.array([0.0, 1.0, 0.0])
        dock_end = tgt_c + lim_c * 0.9 * y_b_chaser
        ax.plot([tgt_c[0], dock_end[0]], [tgt_c[1], dock_end[1]], [tgt_c[2], dock_end[2]],
                color="k", ls="--", lw=1.5, alpha=0.7)

        a_mag = 0.0 if i == 0 or controls.size == 0 else float(np.linalg.norm(controls[min(i - 1, len(controls) - 1)]))
        txt = f"t={t[i]:6.1f}s  ||r||={np.linalg.norm(states[i, :3]):6.2f}m  ||a||={a_mag:5.3f}m/s^2"
        ax.text2D(0.02, 0.96, txt, transform=ax.transAxes)

    def update(frame_idx: int):
        i = int(frame_ids[frame_idx])
        draw_scene(i)
        return []

    draw_scene(int(frame_ids[0]))
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(frame_ids),
        interval=1000.0 / max(1, cfg.animation_fps),
        blit=False,
    )
    return _save_animation(ani, fig, cfg, stem=stem)


def create_animation_suite(cfg: ScenarioConfig, sim: Dict, los: LOS3D, stem_prefix: str = "approach") -> Dict[str, Dict[str, str | None]]:
    """Create animations.

    Dev mode: only 5 GIFs (body_xy, body_yz, lvlh_3d, lvlh_xy, lvlh_xz).
    Full mode: broad set of separate animations.
    """
    out: Dict[str, Dict[str, str | None]] = {}

    # ── Core 6 GIFs (always generated) ───────────────────────────────────
    out["lvlh_3d"] = create_lvlh_3d_animation(cfg, sim, los, stem=f"{stem_prefix}_lvlh_3d")
    out["body_xy"] = create_body_plane_animation(cfg, sim, los, plane="xy", stem=f"{stem_prefix}_body_xy")
    out["body_yz"] = create_body_plane_animation(cfg, sim, los, plane="yz", stem=f"{stem_prefix}_body_yz")
    out["lvlh_xy"] = create_lvlh_plane_animation(cfg, sim, los, plane="xy", stem=f"{stem_prefix}_lvlh_xy")
    out["target_body_frame"] = create_target_body_frame_animation(cfg, sim, los, stem=f"{stem_prefix}_target_body_frame")

    # MPC cost time series
    t = sim["time"]
    t_u = t[:-1]
    out["mpc_cost_time"] = _create_time_series_animation(
        cfg, t_anim=t,
        t_series=sim.get("mpc_cost_time", t_u),
        y_series=sim.get("mpc_cost", np.zeros_like(t_u)),
        stem=f"{stem_prefix}_mpc_cost_time",
        title="MPC Cost vs Time", y_label="Cost J", color="tab:purple",
    )

    if cfg.mode != "dev":
        # Full mode adds all remaining animations
        states = sim["state_lvlh"]
        states_body = sim.get("state_body")
        if states_body is None:
            states_body = los.body_frame_history(states, t)
        controls = sim["control"]
        accel_norm = np.linalg.norm(controls, axis=1) if controls.size else np.zeros(0)
        fuel_use = np.cumsum(accel_norm) * cfg.ts if accel_norm.size else np.zeros(0)

        out["multiview"] = create_approach_animation(cfg, sim, los, stem=stem_prefix)
        out["lvlh_yz"] = create_lvlh_plane_animation(cfg, sim, los, plane="yz", stem=f"{stem_prefix}_lvlh_yz")
        out["lvlh_xz"] = create_lvlh_plane_animation(cfg, sim, los, plane="xz", stem=f"{stem_prefix}_lvlh_xz")
        out["body_xz"] = create_body_plane_animation(cfg, sim, los, plane="xz", stem=f"{stem_prefix}_body_xz")
        out["chaser_body_frame"] = create_chaser_body_frame_animation(cfg, sim, los, stem=f"{stem_prefix}_chaser_body_frame")
        out["range_time"] = _create_time_series_animation(
            cfg, t_anim=t, t_series=t, y_series=sim["range_m"],
            stem=f"{stem_prefix}_range_time",
            title="Range vs Time", y_label="||r|| (m)", color="tab:blue",
        )
        out["los_margin_time"] = _create_time_series_animation(
            cfg, t_anim=t, t_series=t, y_series=sim["los_margin"],
            stem=f"{stem_prefix}_los_margin_time",
            title="LOS Margin vs Time", y_label="min(b - A x)", color="tab:orange",
        )
        out["fuel_use_time"] = _create_time_series_animation(
            cfg, t_anim=t, t_series=t_u, y_series=fuel_use,
            stem=f"{stem_prefix}_fuel_use_time",
            title="Cumulative Fuel-Use Proxy vs Time",
            y_label="Cumulative Delta-v proxy (m/s)", color="tab:red",
        )

    return out
