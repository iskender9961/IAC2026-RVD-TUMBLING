"""Body-frame x-y reachability/feasibility sweeps for rotating LOS constraints.

Implements a reachability-inspired safe start-region analysis that goes beyond
brute-force Monte Carlo:

1. **Geometric erosion analysis** (primary, instant): For each grid point
   (x_B, y_B), compute the per-constraint LOS margin and the directional
   erosion due to the rotating cone.  Points where every constraint's margin
   exceeds its erosion are marked "safe".  This is a conservative inner
   approximation.

2. **Overlay plots**: For each tumble rate, produce a single overlay plot
   showing the LOS cone boundary and nested safe regions for all thrust
   levels (0.2, 0.1, 0.05, 0.01), layered largest-to-smallest.

Erosion model with synchronization range limit
------------------------------------------------
When the body frame rotates at omega_z about z, an inertially-stationary
chaser at (x_B, y_B, 0) sees apparent body-frame velocity:

    v_rot = [omega_z * y_B,  -omega_z * x_B,  0]

(This is -omega x r_B, the velocity from the inverse rotation transport.)

For each constraint face i with position normal a_i = A_c[i, :3], the rate
at which the constraint slack changes is:

    d(slack_i)/dt = -a_i @ v_rot

If this rate is negative (margin shrinking), the controller must apply
thrust a_max to arrest the drift.  The settling-time analysis gives:

    settling time:      t_s = |slack_rate| / a_max
    margin consumed:    erosion_i = 0.5 * slack_rate^2 / a_max

If the rate is non-negative (margin growing or constant), no erosion occurs.

This makes the safe region **asymmetric in x_B**: for positive omega_z
(CCW target rotation), the right-hand cone boundary erodes faster, tilting
the safe region toward -x_B.

In addition to the per-constraint erosion, a **synchronization range limit**
is enforced.  At range r the apparent rotational speed is v = omega * r.
To synchronize, the chaser must arrest this velocity, requiring braking
distance d_brake = 0.5 * v^2 / a_max = 0.5 * omega^2 * r^2 / a_max.
If d_brake exceeds the range r itself, the chaser cannot close to the
hold point without overshooting.  This gives:

    r_max = 2 * a_max / omega^2

Points beyond r_max are marked infeasible regardless of constraint margins.
This bounds the safe region in range even on the downstream side where
directional erosion is zero.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from IAC.config import ScenarioConfig
from IAC.constraints import LOS3D


# ── Geometric erosion analysis ──────────────────────────────────────────────

def _per_constraint_erosion(
    x_b: float,
    y_b: float,
    omega_rad_s: float,
    a_max: float,
    A_c: np.ndarray,
) -> np.ndarray:
    """Directional per-constraint erosion at (x_B, y_B).

    The rotation of the body frame creates an apparent body-frame velocity
    for a stationary chaser:
        v_rot = [omega * y_B, -omega * x_B, 0]

    For each constraint face i, the slack rate is:
        d(slack_i)/dt = -A_c[i, :3] @ v_rot

    If negative (margin shrinking), settling-time drift gives:
        erosion_i = 0.5 * slack_rate_i^2 / a_max

    If non-negative (margin growing), erosion_i = 0.
    """
    v_rot = np.array([omega_rad_s * y_b, -omega_rad_s * x_b, 0.0])
    n_c = A_c.shape[0]
    erosions = np.zeros(n_c)
    for i in range(n_c):
        slack_rate = -float(A_c[i, :3] @ v_rot)
        if slack_rate < 0.0:
            erosions[i] = 0.5 * slack_rate * slack_rate / max(a_max, 1e-12)
    return erosions


def _compute_safe_mask(
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    cone_mask: np.ndarray,
    los: LOS3D,
    omega_rad_s: float,
    a_max: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute safe mask and margin field for one (omega, a_max) pair.

    Two conditions must hold for a point to be marked safe:

    1. **Directional constraint erosion**: each constraint face is eroded
       based on the negative slack rate under the rotation-induced apparent
       velocity.  All per-constraint net margins must be positive.

    2. **Synchronization range limit**: the apparent rotational speed at
       range r is v_rot = omega * r.  To synchronize, the chaser must
       arrest this velocity.  The braking distance is 0.5 * v_rot^2 / a_max.
       If this exceeds the range itself, the chaser overshoots during
       synchronization and cannot close to the hold point.  The condition
       r < 2 * a_max / omega^2 bounds the safe region in range.

    Returns (safe_mask, net_margin) where:
        safe_mask[iy, ix] = True if both conditions are met
        net_margin[iy, ix] = min over constraints of (slack_i - erosion_i)
    """
    ny, nx = cone_mask.shape
    net_margin = np.full((ny, nx), np.nan, dtype=float)
    safe_mask = np.zeros((ny, nx), dtype=bool)
    A_c = los.A_c
    b_c = los.b_c

    # Synchronization range limit: r_max = 2 * a_max / omega^2
    w2 = omega_rad_s * omega_rad_s
    r_sync_max = 2.0 * a_max / max(w2, 1e-30)

    for iy in range(ny):
        for ix in range(nx):
            if not cone_mask[iy, ix]:
                continue
            xb = float(x_grid[ix])
            yb = float(y_grid[iy])
            rng = float(np.sqrt(xb * xb + yb * yb))

            p = np.array([xb, yb, 0.0, 0.0, 0.0, 0.0], dtype=float)
            per_slack = b_c - A_c @ p
            erosions = _per_constraint_erosion(xb, yb, omega_rad_s, a_max, A_c)
            net_per = per_slack - erosions
            net = float(np.min(net_per))
            net_margin[iy, ix] = net
            safe_mask[iy, ix] = (net > 0.0) and (rng < r_sync_max)

    return safe_mask, net_margin


def _grid_and_mask(
    cfg: ScenarioConfig, los: LOS3D
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build evaluation grid and inside-cone mask.

    The x-range is derived from the cone geometry so the grid always
    covers the full cone width at y_max:
        |x_B| <= x0 + cx * (y_max - r_hold)
    A 5 % padding is added so the forbidden boundary is visible.
    """
    y_grid = np.linspace(
        cfg.reachability_y_min_m, cfg.reachability_y_max_m, int(cfg.reachability_ny)
    )
    x_cone_max = los.x0_m + los.cx * (cfg.reachability_y_max_m - los.hold_radius_m)
    x_limit = x_cone_max * 1.05  # 5 % padding outside cone edge
    x_grid = np.linspace(-x_limit, x_limit, int(cfg.reachability_nx))
    mask = np.zeros((len(y_grid), len(x_grid)), dtype=bool)
    for iy, yb in enumerate(y_grid):
        for ix, xb in enumerate(x_grid):
            p = np.array([xb, yb, 0.0, 0.0, 0.0, 0.0], dtype=float)
            slacks = los.b_c - los.A_c @ p
            mask[iy, ix] = bool(np.all(slacks >= -1e-9))
    return x_grid, y_grid, mask


# ── Overlay plots ───────────────────────────────────────────────────────────

def _plot_overlay_for_tumble(
    cfg: ScenarioConfig,
    los_geom: LOS3D,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    cone_mask: np.ndarray,
    safe_masks: dict,
    tumble_deg_s: float,
    accel_levels: np.ndarray,
) -> str:
    """Create overlay plot for one tumble rate, all thrust levels."""
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    extent = (float(x_grid[0]), float(x_grid[-1]),
              float(y_grid[0]), float(y_grid[-1]))

    y_curve = np.linspace(
        max(float(y_grid[0]), los_geom.hold_radius_m), float(y_grid[-1]), 300
    )
    x_curve = los_geom.cx * (y_curve - los_geom.hold_radius_m) + los_geom.x0_m

    # Background: red for forbidden, strong blue for initial LOS cone
    bg = np.zeros((*cone_mask.shape, 4), dtype=float)
    bg[~cone_mask] = [0.85, 0.25, 0.22, 0.55]   # red forbidden zone
    bg[cone_mask] = [0.40, 0.60, 1.0, 0.50]      # visible blue for feasible
    ax.imshow(bg, origin="lower", extent=extent, aspect="auto",
              interpolation="nearest")

    colors = {
        0.20: "#a5d6a7",
        0.10: "#66bb6a",
        0.05: "#2e7d32",
        0.01: "#1b5e20",
    }
    inside_count = max(int(np.sum(cone_mask)), 1)

    legend_items = [
        Patch(facecolor=(0.85, 0.25, 0.22, 0.55), edgecolor="k", lw=0.5,
              label="Forbidden (outside LOS cone)"),
        Patch(facecolor=(0.40, 0.60, 1.0, 0.50), edgecolor="k", lw=0.5,
              label="LOS cone (initially feasible)"),
    ]

    for amax in sorted(accel_levels, reverse=True):
        mask_key = f"{amax:.2f}"
        if mask_key not in safe_masks:
            continue
        safe = safe_masks[mask_key]
        color = colors.get(round(amax, 2), "#388e3c")

        overlay = np.zeros((*cone_mask.shape, 4), dtype=float)
        r, g, b = (int(color[1:3], 16) / 255,
                    int(color[3:5], 16) / 255,
                    int(color[5:7], 16) / 255)
        overlay[safe] = [r, g, b, 0.7]
        ax.imshow(overlay, origin="lower", extent=extent, aspect="auto",
                  interpolation="nearest")

        frac = int(np.sum(safe)) / inside_count
        legend_items.append(
            Patch(facecolor=color, edgecolor="none",
                  label=f"Safe, $a_{{max}}$={amax:.2f} m/s$^2$ ({100*frac:.0f}%)")
        )

    ax.plot(x_curve, y_curve, "k--", lw=1.5)
    ax.plot(-x_curve, y_curve, "k--", lw=1.5)
    ax.axhline(los_geom.hold_radius_m, color="k", ls=":", lw=1.0)

    ax.set_xlabel("$x_B$ (m)", fontsize=11)
    ax.set_ylabel("$y_B$ (m)", fontsize=11)
    ax.set_title(
        f"Safe Start Region ($v_0$=0), Target Tumble = {tumble_deg_s:.0f} deg/s",
        fontsize=12,
    )
    ax.legend(handles=legend_items, loc="upper right", fontsize=8,
              framealpha=0.9)
    ax.grid(True, alpha=0.18)
    fig.tight_layout()

    stem = f"safe_maps_v0_zero_omega_{int(round(tumble_deg_s))}deg"
    png_path = cfg.fig_dir / f"{stem}.png"
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=220)
    try:
        fig.savefig(png_path.with_suffix(".pdf"))
    except PermissionError:
        print(f"  [warn] PDF locked, skipped: {png_path.with_suffix('.pdf')}")
    plt.close(fig)
    return str(png_path)


def _plot_combined_overview(
    cfg: ScenarioConfig,
    los_geom: LOS3D,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    cone_mask: np.ndarray,
    all_safe_masks: dict,
    tumble_rates: np.ndarray,
    accel_levels: np.ndarray,
) -> str:
    """Combined overview with one panel per tumble rate."""
    n_w = len(tumble_rates)
    fig, axes = plt.subplots(1, n_w, figsize=(3.8 * n_w, 5.2),
                             sharex=True, sharey=True)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])

    extent = (float(x_grid[0]), float(x_grid[-1]),
              float(y_grid[0]), float(y_grid[-1]))
    y_curve = np.linspace(
        max(float(y_grid[0]), los_geom.hold_radius_m), float(y_grid[-1]), 300
    )
    x_curve = los_geom.cx * (y_curve - los_geom.hold_radius_m) + los_geom.x0_m
    inside_count = max(int(np.sum(cone_mask)), 1)

    colors_ordered = ["#a5d6a7", "#66bb6a", "#2e7d32", "#1b5e20"]

    for iw, w_deg in enumerate(tumble_rates):
        ax = axes[iw]
        bg = np.zeros((*cone_mask.shape, 4), dtype=float)
        bg[~cone_mask] = [0.85, 0.25, 0.22, 0.55]   # red forbidden
        bg[cone_mask] = [0.40, 0.60, 1.0, 0.50]      # visible blue
        ax.imshow(bg, origin="lower", extent=extent, aspect="auto",
                  interpolation="nearest")

        w_key = f"{w_deg:.1f}"
        if w_key in all_safe_masks:
            for ja, amax in enumerate(sorted(accel_levels, reverse=True)):
                mask_key = f"{amax:.2f}"
                if mask_key in all_safe_masks[w_key]:
                    safe = all_safe_masks[w_key][mask_key]
                    color = colors_ordered[ja % len(colors_ordered)]
                    overlay = np.zeros((*cone_mask.shape, 4), dtype=float)
                    r, g, b = (int(color[1:3], 16) / 255,
                               int(color[3:5], 16) / 255,
                               int(color[5:7], 16) / 255)
                    overlay[safe] = [r, g, b, 0.7]
                    ax.imshow(overlay, origin="lower", extent=extent,
                              aspect="auto", interpolation="nearest")

        ax.plot(x_curve, y_curve, "k--", lw=1.2)
        ax.plot(-x_curve, y_curve, "k--", lw=1.2)
        ax.axhline(los_geom.hold_radius_m, color="k", ls=":", lw=0.8)
        ax.set_title(f"$\\omega$={w_deg:.0f} deg/s", fontsize=10)
        ax.set_xlabel("$x_B$ (m)")
        if iw == 0:
            ax.set_ylabel("$y_B$ (m)")
        ax.grid(True, alpha=0.18)

    legend_items = [
        Patch(facecolor=(0.85, 0.25, 0.22, 0.55), label="Forbidden"),
        Patch(facecolor=(0.40, 0.60, 1.0, 0.50), label="LOS cone"),
    ]
    for ja, amax in enumerate(sorted(accel_levels, reverse=True)):
        legend_items.append(
            Patch(facecolor=colors_ordered[ja % len(colors_ordered)],
                  label=f"$a_{{max}}$={amax:.2f}")
        )
    fig.legend(handles=legend_items, loc="lower center",
               ncol=len(legend_items), frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=8)
    fig.suptitle(
        "Safe Start Regions ($v_0$=0) — Erosion-Based Inner Approximation",
        y=1.02, fontsize=12,
    )
    fig.tight_layout()

    png_path = cfg.fig_dir / "safe_maps_v0_zero.png"
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, bbox_inches="tight", dpi=220)
    try:
        fig.savefig(png_path.with_suffix(".pdf"), bbox_inches="tight")
    except PermissionError:
        print(f"  [warn] PDF locked, skipped: {png_path.with_suffix('.pdf')}")
    plt.close(fig)
    return str(png_path)


# ── Table formatting ────────────────────────────────────────────────────────

def format_reachability_table_terminal(reach: Dict) -> str:
    tumble = np.asarray(reach["tumble_deg_s"], dtype=float)
    accel = np.asarray(reach["accel_mps2"], dtype=float)
    frac = np.asarray(reach["feasible_fraction"], dtype=float)

    cell_w = 10
    header = "omega/amax".ljust(cell_w) + "".join(
        f"{a:>{cell_w}.2f}" for a in accel
    )
    sep = "-" * len(header)
    lines = [header, sep]
    for i, w in enumerate(tumble):
        row = f"{w:>{cell_w}.0f} deg/s"
        row += "".join(
            f"{100.0 * frac[i, j]:>{cell_w}.1f}%" for j in range(len(accel))
        )
        lines.append(row)
    return "\n".join(lines)


def build_reachability_latex_table(reach: Dict) -> str:
    tumble = np.asarray(reach["tumble_deg_s"], dtype=float)
    accel = np.asarray(reach["accel_mps2"], dtype=float)
    frac = np.asarray(reach["feasible_fraction"], dtype=float)

    cols = "c|" + "c" * len(accel)
    header = " & ".join(
        [r"$\omega_t$ (deg/s)"] + [f"{a:.2f}" for a in accel]
    ) + r" \\"
    rows = []
    for i, w in enumerate(tumble):
        row = [f"{w:.0f}"] + [
            f"{100.0 * frac[i, j]:.1f}\\%" for j in range(len(accel))
        ]
        rows.append(" & ".join(row) + r" \\")

    return "\n".join([
        r"\begin{table}[H]",
        r"\centering",
        r"\caption{Safe fraction of the body-frame LOS cone for various tumble-rate "
        r"and thrust-authority combinations.  Criteria: directional per-constraint "
        r"erosion $\delta_i = \frac{1}{2}\dot{s}_i^2 / a_{\max}$ and synchronization "
        r"range bound $r < 2a_{\max}/\omega_t^2$.}",
        r"\label{tab:reachability_sweep}",
        rf"\begin{{tabular}}{{{cols}}}",
        r"\hline",
        header,
        r"\hline",
        *rows,
        r"\hline",
        r"\end{tabular}",
        r"\end{table}",
    ])


# ── Main entry point ───────────────────────────────────────────────────────

def run_reachability_study(
    cfg: ScenarioConfig,
    verbose: bool = True,
) -> Dict:
    """Compute erosion-based safe start regions for all (omega, a_max) pairs.

    Uses a purely geometric analysis (no closed-loop simulation required).
    Grid density is controlled by cfg.reachability_nx / ny.
    """
    cfg.ensure_output_dirs()
    tumble_deg = np.asarray(cfg.reachability_tumble_deg_s, dtype=float)
    accel = np.asarray(cfg.reachability_accel_mps2, dtype=float)
    los_geom = LOS3D.from_config(cfg)
    x_grid, y_grid, cone_mask = _grid_and_mask(cfg, los_geom)

    n_w = len(tumble_deg)
    n_a = len(accel)
    ny, nx = cone_mask.shape
    safe_eroded = np.zeros((n_w, n_a, ny, nx), dtype=bool)
    net_margin_all = np.full((n_w, n_a, ny, nx), np.nan, dtype=float)
    feasible_fraction = np.zeros((n_w, n_a), dtype=float)
    inside_count = max(int(np.sum(cone_mask)), 1)

    for iw, w_deg in enumerate(tumble_deg):
        omega_rad = float(np.deg2rad(w_deg))
        for ia, amax in enumerate(accel):
            if verbose:
                print(f"[Reachability] omega={w_deg:.0f} deg/s, "
                      f"a_max={amax:.2f} m/s^2")

            safe, net_m = _compute_safe_mask(
                x_grid, y_grid, cone_mask, los_geom,
                omega_rad, float(amax),
            )
            safe_eroded[iw, ia] = safe
            net_margin_all[iw, ia] = net_m
            safe_count = int(np.sum(np.logical_and(cone_mask, safe)))
            feasible_fraction[iw, ia] = safe_count / inside_count

            if verbose:
                print(f"  safe fraction: {100.0 * feasible_fraction[iw, ia]:.1f}%")

    # Build overlay plots
    figure_paths: Dict[str, str] = {}
    all_safe_masks: Dict[str, Dict[str, np.ndarray]] = {}

    for iw, w_deg in enumerate(tumble_deg):
        w_key = f"{w_deg:.1f}"
        safe_masks: Dict[str, np.ndarray] = {}
        for ia, amax in enumerate(accel):
            mask_key = f"{amax:.2f}"
            safe_masks[mask_key] = np.logical_and(
                cone_mask, safe_eroded[iw, ia]
            )
        all_safe_masks[w_key] = safe_masks

        key = f"tumble_{int(round(w_deg))}deg"
        figure_paths[key] = _plot_overlay_for_tumble(
            cfg=cfg, los_geom=los_geom,
            x_grid=x_grid, y_grid=y_grid, cone_mask=cone_mask,
            safe_masks=safe_masks, tumble_deg_s=float(w_deg),
            accel_levels=accel,
        )

    # Combined overview figure
    figure_paths["safe_maps_v0_zero"] = _plot_combined_overview(
        cfg=cfg, los_geom=los_geom,
        x_grid=x_grid, y_grid=y_grid, cone_mask=cone_mask,
        all_safe_masks=all_safe_masks, tumble_rates=tumble_deg,
        accel_levels=accel,
    )

    # Save numeric data
    npz_path = cfg.data_dir / "reachability" / "reachability_xy_sweep.npz"
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        npz_path,
        x_grid=x_grid, y_grid=y_grid, cone_mask=cone_mask,
        tumble_deg_s=tumble_deg, accel_mps2=accel,
        safe_eroded=safe_eroded, net_margin=net_margin_all,
        feasible_fraction=feasible_fraction,
        seed=cfg.seed,
    )

    # Save metadata JSON
    meta_path = cfg.data_dir / "reachability" / "metadata.json"
    meta = {
        "method": "geometric_erosion_continuous",
        "tumble_deg_s": tumble_deg.tolist(),
        "accel_mps2": accel.tolist(),
        "grid_nx": int(cfg.reachability_nx),
        "grid_ny": int(cfg.reachability_ny),
        "x_range_m": [float(x_grid[0]), float(x_grid[-1])],
        "y_range_m": [float(cfg.reachability_y_min_m),
                       float(cfg.reachability_y_max_m)],
        "seed": int(cfg.seed),
        "erosion_formula": "per-constraint directional erosion_i = 0.5 * slack_rate_i^2 / a_max (if slack_rate<0) + sync range r < 2*a_max/omega^2",
    }
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    out = {
        "x_grid": x_grid.tolist(),
        "y_grid": y_grid.tolist(),
        "inside_count": inside_count,
        "tumble_deg_s": tumble_deg.tolist(),
        "accel_mps2": accel.tolist(),
        "feasible_fraction": feasible_fraction.tolist(),
        "figure_paths": figure_paths,
        "npz_path": str(npz_path),
    }
    return out
