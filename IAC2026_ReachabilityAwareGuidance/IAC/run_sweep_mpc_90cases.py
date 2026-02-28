"""MPC-HCW 90-case parameter sweep with 10 initial conditions per case.

Sweep grid:
    a_max     = {0.10, 0.06, 0.02}  m/s²   (3 values)
    y_limit   = {30, 100, 150}      m       (3 values)
    omega     = {1,2,3,4,5}         deg/s   (5 values)
    cost      = {A, B}                      (2 variants)
    Total     = 3 × 3 × 5 × 2 = 90 cases

Cost variants:
    Cost A (baseline):  position error + delta-U  (no velocity penalty, no direct U)
    Cost B (augmented): position error + velocity error + delta-U  (no direct U)

Per-case:  10 initial conditions in target body frame, all at y_B = y_limit.

Outputs per case (aggregated 10-run plots):
    body_position_vs_time.pdf
    control_and_dv_proxy.pdf
    mpc_control_vs_time.pdf
    mpc_cost_breakdown.pdf
    mpc_range_vs_time.pdf
    mpc_traj_body_3d.pdf
    range_vs_time.pdf
    approach_body_xy.gif
    approach_lvlh_3d.gif
    approach_lvlh_xy_r0_50m.gif
    case_config.json

Plus: results/sweep_mpc_90cases/summary.csv

Usage
-----
    python IAC/run_sweep_mpc_90cases.py                 # full 90-case sweep
    python IAC/run_sweep_mpc_90cases.py --validate      # 1-case validation (2 ICs)
    python IAC/run_sweep_mpc_90cases.py --replot        # replot from saved .npz (no re-sim)
"""
from __future__ import annotations
import sys, os, json, csv, time as _timer, argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
THIS_DIR   = Path(__file__).resolve().parent
RESULTS    = THIS_DIR / "results" / "sweep_mpc_90cases"

# ---------------------------------------------------------------------------
# Physics constants
# ---------------------------------------------------------------------------
n_orbit   = 1.1e-3          # mean motion [rad/s]
ts_truth  = 0.1             # truth propagation step [s]
ts_mpc    = 1.0             # MPC control interval [s]
sub_steps = int(round(ts_mpc / ts_truth))  # 10
sim_dur   = 600.0           # 10 min per case
N_hor     = 30              # MPC horizon
nx, nu    = 6, 3
n_los     = 5

# LOS cone geometry (body frame)
y_min     = 0.50
x_0_los   = 2.50
z_0_los   = 2.50
c_x       = 1.5
c_z       = 1.5
A_los     = np.array([
    [ 0., -1.,  0.],
    [ 1., -c_x, 0.],
    [-1., -c_x, 0.],
    [ 0., -c_z, 1.],
    [ 0., -c_z,-1.],
])
b_los     = np.array([
    -y_min,
    x_0_los - c_x*y_min,
    x_0_los - c_x*y_min,
    z_0_los - c_z*y_min,
    z_0_los - c_z*y_min,
])
los_margin = 0.05
dock_range = 0.5
dock_speed = 0.05

# ---------------------------------------------------------------------------
# Weight matrices
# ---------------------------------------------------------------------------
# Body-frame position weights (y_B penalised 10×)
q_pos  = np.array([5.0, 50.0, 5.0])
q_vel  = np.array([2.0, 10.0, 2.0])
w_slack = 1e4

# Cost A: position + delta-u (no velocity, no direct u)
Q_body_A = np.diag(np.concatenate([q_pos, np.zeros(3)]))       # vel weight = 0
r_du_A   = 30.0

# Cost B: position + velocity + delta-u (no direct u)
Q_body_B = np.diag(np.concatenate([q_pos, q_vel]))              # vel weight on
r_du_B   = 30.0

# Sweep grid
a_maxes     = [0.10, 0.06, 0.02]
y_limits    = [30, 100, 150]
omegas_deg  = [1, 2, 3, 4, 5]

# 10 initial conditions (body frame offsets)
def build_ics(y0):
    """Return list of 10 (label, x_B, y_B, z_B) tuples."""
    return [
        ("IC01: x=0,y=Y0",      0.,  y0,    0.),
        ("IC02: x=-10,y=Y0",  -10.,  y0,    0.),
        ("IC03: x=+10,y=Y0",   10.,  y0,    0.),
        ("IC04: x=-20,y=Y0",  -20.,  y0,    0.),
        ("IC05: x=+20,y=Y0",   20.,  y0,    0.),
        ("IC06: x=-5,y=Y0",    -5.,  y0,    0.),
        ("IC07: x=+5,y=Y0",     5.,  y0,    0.),
        ("IC08: x=-15,y=Y0",  -15.,  y0,    0.),
        ("IC09: x=+15,y=Y0",   15.,  y0,    0.),
        ("IC10: x=0,y=Y0+5",    0.,  y0+5., 0.),
    ]

# ---------------------------------------------------------------------------
# CWH discretisation
# ---------------------------------------------------------------------------
def discretize_cwh(n, dt):
    Ac = np.zeros((6,6))
    Ac[0,3]=1; Ac[1,4]=1; Ac[2,5]=1
    Ac[3,0]=3*n**2; Ac[3,4]=2*n; Ac[4,3]=-2*n; Ac[5,2]=-n**2
    Bc = np.zeros((6,3)); Bc[3,0]=1; Bc[4,1]=1; Bc[5,2]=1
    M = np.zeros((9,9)); M[:6,:6]=Ac*dt; M[:6,6:]=Bc*dt
    eM = la.expm(M)
    return eM[:6,:6].copy(), eM[:6,6:].copy()

Ad_mpc,   Bd_mpc   = discretize_cwh(n_orbit, ts_mpc)
Ad_truth, Bd_truth  = discretize_cwh(n_orbit, ts_truth)

def rotz(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c,-s,0],[s,c,0],[0,0,1.]])

def T6(th):
    Rt = rotz(th).T
    T = np.zeros((6,6)); T[:3,:3]=Rt; T[3:,3:]=Rt
    return T

# ---------------------------------------------------------------------------
# DARE terminal cost (one per cost variant, computed once at import)
# ---------------------------------------------------------------------------
from scipy.linalg import solve_discrete_are

def _compute_dare(Q, R_du):
    """Compute DARE terminal cost: P = dare(A, B, Q, R)."""
    R = R_du * np.eye(3)
    try:
        P = solve_discrete_are(Ad_mpc, Bd_mpc, Q, R)
        eigs = la.eigvalsh(P)
        return P, eigs.min(), eigs.max()
    except Exception:
        return 10.0 * Q, 0., 0.

P_dare_A, eA_lo, eA_hi = _compute_dare(Q_body_A, r_du_A)
P_dare_B, eB_lo, eB_hi = _compute_dare(Q_body_B, r_du_B)

print(f"DARE-A (pos+du): eig [{eA_lo:.1f}, {eA_hi:.1f}]")
print(f"DARE-B (pos+vel+du): eig [{eB_lo:.1f}, {eB_hi:.1f}]")

# ---------------------------------------------------------------------------
# Body-frame IC → LVLH state
# ---------------------------------------------------------------------------
def body_ic_to_lvlh(x_B, y_B, z_B, omega_t, t0=0.0):
    """Convert body-frame position (zero body velocity) to LVLH state.

    Body-frame velocity is zero → relative velocity is purely from co-rotation.
    r_L = R(theta) @ r_B
    v_L = R(theta) @ (v_B + omega×r_B)  with v_B = 0
    """
    theta = omega_t * t0
    r_B = np.array([x_B, y_B, z_B])
    R = rotz(theta)
    r_L = R @ r_B
    omega_cross_rB = np.array([-omega_t * y_B, omega_t * x_B, 0.0])
    v_L = R @ omega_cross_rB
    return np.concatenate([r_L, v_L])

# ---------------------------------------------------------------------------
# MPC QP solver
# ---------------------------------------------------------------------------
def solve_mpc(x0, u_prev, t_now, omega_t, a_max_val, cost_id):
    """Solve one MPC step.

    cost_id='A': position + delta-u, no velocity cost, no direct-u
    cost_id='B': position + velocity + delta-u, no direct-u
    """
    import osqp
    N  = N_hor; ns = n_los
    nz = N*nx + N*nu + N*ns

    def xi(j): return j*nx
    def ui(j): return N*nx + j*nu
    def si(j): return N*nx + N*nu + j*ns

    Q_stage = Q_body_A if cost_id == 'A' else Q_body_B
    P_term  = P_dare_A if cost_id == 'A' else P_dare_B
    r_du    = r_du_A   if cost_id == 'A' else r_du_B

    P_mat = np.zeros((nz, nz)); q_vec = np.zeros(nz)

    # -- State cost (body frame) + DARE terminal --
    for j in range(N):
        th_j = omega_t * (t_now + (j+1)*ts_mpc)
        Tj = T6(th_j)
        Wj = P_term if j == N-1 else Q_stage
        P_mat[xi(j):xi(j)+nx, xi(j):xi(j)+nx] += Tj.T @ Wj @ Tj

    # -- Delta-u cost (no direct u penalty) --
    for j in range(N):
        s_ = ui(j)
        P_mat[s_:s_+nu, s_:s_+nu] += r_du * np.eye(nu)
        if j == 0:
            q_vec[s_:s_+nu] += -2 * r_du * u_prev
        else:
            sp_ = ui(j-1)
            P_mat[sp_:sp_+nu, sp_:sp_+nu] += r_du * np.eye(nu)
            P_mat[s_:s_+nu,  sp_:sp_+nu]  += -r_du * np.eye(nu)
            P_mat[sp_:sp_+nu, s_:s_+nu]   += -r_du * np.eye(nu)

    # -- Slack penalty (scaled by a_max for conditioning) --
    ws = w_slack * (a_max_val / 0.10)
    for j in range(N):
        s_ = si(j)
        P_mat[s_:s_+ns, s_:s_+ns] += ws * np.eye(ns)
        q_vec[s_:s_+ns]           += ws * np.ones(ns)

    P_mat = 0.5 * (P_mat + P_mat.T)

    # -- Equality + inequality constraints --
    n_con = N*nx + N*nu + N*ns + N*ns
    A_mat = np.zeros((n_con, nz))
    l_con = np.zeros(n_con); u_con = np.zeros(n_con)

    # Dynamics: x_{j} = A x_{j-1} + B u_j
    off = 0
    for j in range(N):
        r0 = off + j*nx
        A_mat[r0:r0+nx, xi(j):xi(j)+nx] = np.eye(nx)
        A_mat[r0:r0+nx, ui(j):ui(j)+nu] = -Bd_mpc
        if j == 0:
            rhs = Ad_mpc @ x0
        else:
            A_mat[r0:r0+nx, xi(j-1):xi(j-1)+nx] = -Ad_mpc
            rhs = np.zeros(nx)
        l_con[r0:r0+nx] = rhs; u_con[r0:r0+nx] = rhs

    # Thrust bounds (box)
    off = N*nx
    for j in range(N):
        r0 = off + j*nu
        A_mat[r0:r0+nu, ui(j):ui(j)+nu] = np.eye(nu)
        l_con[r0:r0+nu] = -a_max_val; u_con[r0:r0+nu] = a_max_val

    # LOS constraints (body-frame cone, rotated)
    off = N*nx + N*nu
    for j in range(N):
        th_j = omega_t * (t_now + (j+1)*ts_mpc)
        Hlj = A_los @ rotz(th_j).T          # project LVLH pos → body
        r0 = off + j*ns
        A_mat[r0:r0+ns, xi(j):xi(j)+3] = Hlj
        A_mat[r0:r0+ns, si(j):si(j)+ns] = -np.eye(ns)
        l_con[r0:r0+ns] = -1e30
        u_con[r0:r0+ns] = b_los - los_margin

    # Slack ≥ 0
    off = N*nx + N*nu + N*ns
    for j in range(N):
        r0 = off + j*ns
        A_mat[r0:r0+ns, si(j):si(j)+ns] = np.eye(ns)
        l_con[r0:r0+ns] = 0.0; u_con[r0:r0+ns] = 1e30

    # Normalise Hessian for OSQP conditioning
    p_scale = max(1.0, float(np.max(np.abs(np.diag(P_mat)))))
    P_mat /= p_scale
    q_vec /= p_scale

    m = osqp.OSQP()
    m.setup(sp.csc_matrix(P_mat), q_vec, sp.csc_matrix(A_mat), l_con, u_con,
            verbose=False, eps_abs=1e-4, eps_rel=1e-4,
            max_iter=4000, polish=True, warm_start=True,
            adaptive_rho=True, scaling=10)
    res = m.solve()

    if res.info.status not in ("solved", "solved_inaccurate"):
        return np.zeros(nu), np.inf, {"feasible": False, "status": res.info.status}

    u_opt = np.asarray(res.x[ui(0):ui(0)+nu])
    total_slack = sum(
        float(np.sum(np.maximum(res.x[si(j):si(j)+ns], 0))) for j in range(N))

    # Extract cost breakdown from solution
    x_pred = np.zeros((N, nx))
    u_seq  = np.zeros((N, nu))
    for j in range(N):
        x_pred[j] = res.x[xi(j):xi(j)+nx]
        u_seq[j]  = res.x[ui(j):ui(j)+nu]

    # Compute individual cost terms
    pos_cost = 0.; vel_cost = 0.; du_cost = 0.; constraint_pen = 0.
    for j in range(N):
        th_j = omega_t * (t_now + (j+1)*ts_mpc)
        xb_j = T6(th_j) @ x_pred[j]
        Wj = P_term if j == N-1 else Q_stage
        w_diag = np.diag(Wj)
        pos_cost += float(np.sum(w_diag[:3] * xb_j[:3]**2))
        vel_cost += float(np.sum(w_diag[3:] * xb_j[3:]**2))
        if j == 0:
            du_cost += float(r_du * np.sum((u_seq[j] - u_prev)**2))
        else:
            du_cost += float(r_du * np.sum((u_seq[j] - u_seq[j-1])**2))
    constraint_pen = float(total_slack * ws)

    info = {
        "feasible": True,
        "total_slack": total_slack,
        "pos_cost": pos_cost,
        "vel_cost": vel_cost,
        "du_cost": du_cost,
        "constraint_penalty": constraint_pen,
    }
    return u_opt, float(res.info.obj_val), info

# ---------------------------------------------------------------------------
# Run one simulation (one case, one IC)
# ---------------------------------------------------------------------------
def run_one_sim(omega_deg, a_max_val, y_limit, cost_id, ic_label, x_B0, y_B0, z_B0):
    """Run closed-loop MPC for a single IC.  Returns result dict."""
    omega_t = np.radians(omega_deg)
    r_sync  = 2 * a_max_val / omega_t**2
    n_sim   = int(sim_dur / ts_mpc)

    # Body IC → LVLH state  (zero body velocity → co-rotation)
    x = body_ic_to_lvlh(x_B0, y_B0, z_B0, omega_t)
    u_prev = np.zeros(nu)

    x_hist  = np.zeros((n_sim+1, nx)); x_hist[0] = x
    xb_hist = np.zeros((n_sim+1, nx)); xb_hist[0] = T6(0.) @ x
    u_hist  = np.zeros((n_sim, nu))
    slk_coarse = np.zeros(n_sim+1)
    slk_coarse[0] = float(np.min(b_los - A_los @ (rotz(0.).T @ x[:3])))

    n_fine   = n_sim * sub_steps
    slk_fine = np.zeros(n_fine + 1)
    slk_fine[0] = slk_coarse[0]

    # Cost breakdown per step
    pos_costs = np.zeros(n_sim)
    vel_costs = np.zeros(n_sim)
    du_costs  = np.zeros(n_sim)
    con_pens  = np.zeros(n_sim)

    dock_step = None; n_infeasible = 0; consec_grow = 0
    r0_m = float(np.linalg.norm(x[:3]))
    prev_rng = r0_m

    for k in range(n_sim):
        t_now = k * ts_mpc
        u_k, _, info = solve_mpc(x, u_prev, t_now, omega_t, a_max_val, cost_id)

        if not info.get("feasible", True):
            n_infeasible += 1
            # Fallback: brake towards origin + velocity damping
            r_vec = x[:3]; r_norm = np.linalg.norm(r_vec)
            v_vec = x[3:]
            if r_norm > 0.1:
                u_k = -a_max_val * r_vec / r_norm * 0.5
                v_norm = np.linalg.norm(v_vec)
                if v_norm > 1e-4:
                    u_k -= a_max_val * v_vec / v_norm * 0.5
        else:
            pos_costs[k] = info.get("pos_cost", 0.)
            vel_costs[k] = info.get("vel_cost", 0.)
            du_costs[k]  = info.get("du_cost", 0.)
            con_pens[k]  = info.get("constraint_penalty", 0.)

        u_k = np.clip(u_k, -a_max_val, a_max_val)
        u_hist[k] = u_k

        # Truth sub-stepping
        for ss in range(sub_steps):
            x = Ad_truth @ x + Bd_truth @ u_k
            idx = k * sub_steps + (ss + 1)
            th_f = omega_t * (t_now + (ss+1)*ts_truth)
            slk_fine[idx] = float(np.min(b_los - A_los @ (rotz(th_f).T @ x[:3])))

        th = omega_t * ((k+1) * ts_mpc)
        x_hist[k+1]  = x
        xb_hist[k+1] = T6(th) @ x
        slk_coarse[k+1] = float(np.min(b_los - A_los @ (rotz(th).T @ x[:3])))
        u_prev = u_k.copy()

        rng = float(np.linalg.norm(x[:3]))
        spd = float(np.linalg.norm(x[3:]))
        if rng <= dock_range + los_margin + 1e-6 and spd <= dock_speed and dock_step is None:
            dock_step = k + 1

        if rng > prev_rng + 0.5:
            consec_grow += 1
        else:
            consec_grow = 0
        prev_rng = rng

        # Early stop
        if (rng > 2.0 * max(r0_m, y_limit) and k > 30) or (consec_grow > 30 and k > 50):
            for kk in range(k+2, n_sim+1):
                x_hist[kk] = x; xb_hist[kk] = xb_hist[k+1]
                slk_coarse[kk] = slk_coarse[k+1]
            u_hist[k+1:] = 0.0
            for kk in range(k*sub_steps+sub_steps+1, n_fine+1):
                slk_fine[kk] = slk_fine[min(idx, n_fine)]
            break

    ranges = np.linalg.norm(x_hist[:,:3], axis=1)
    speeds = np.linalg.norm(x_hist[:,3:], axis=1)
    dv_cumulative = np.cumsum(np.linalg.norm(u_hist, axis=1) * ts_mpc)
    dv_total = float(dv_cumulative[-1]) if len(dv_cumulative) > 0 else 0.

    # Max lateral excursion in body frame
    max_lateral = float(np.max(np.abs(xb_hist[:, 0])))  # x_B

    return {
        "ic_label": ic_label,
        "x_B0": x_B0, "y_B0": y_B0, "z_B0": z_B0,
        "r_sync": r_sync,
        "r_final": float(ranges[-1]),
        "v_final": float(speeds[-1]),
        "r_min":   float(np.min(ranges)),
        "dv_total": dv_total,
        "min_slk": float(np.min(slk_fine)),
        "n_viol":  int(np.sum(slk_fine < -1e-6)),
        "dock":    dock_step is not None,
        "dock_step": dock_step,
        "n_infeasible": n_infeasible,
        "max_lateral": max_lateral,
        # Trajectories
        "xb": xb_hist, "x_lvlh": x_hist, "u": u_hist,
        "slk": slk_coarse, "ranges": ranges, "speeds": speeds,
        "dv_cumulative": np.concatenate([[0.], dv_cumulative]),
        "time": np.arange(n_sim+1) * ts_mpc,
        "pos_costs": pos_costs, "vel_costs": vel_costs,
        "du_costs": du_costs, "con_pens": con_pens,
    }

# ---------------------------------------------------------------------------
# Case tag and directory
# ---------------------------------------------------------------------------
def case_dir_name(a_max_val, y_limit, omega_deg, cost_id):
    return f"a{a_max_val:.2f}_y{y_limit}_w{omega_deg}_cost{cost_id}"

def case_dir(a_max_val, y_limit, omega_deg, cost_id):
    d = RESULTS / case_dir_name(a_max_val, y_limit, omega_deg, cost_id)
    d.mkdir(parents=True, exist_ok=True)
    return d

# ---------------------------------------------------------------------------
# Color palette for 10 ICs
# ---------------------------------------------------------------------------
IC_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
]

# ---------------------------------------------------------------------------
# Plotting: aggregated 10-run plots
# ---------------------------------------------------------------------------
def _case_title(a_max_val, y_limit, omega_deg, cost_id):
    return (f"$a_{{max}}$={a_max_val:.2f} m/s²   "
            f"$y_{{lim}}$={y_limit} m   "
            f"$\\omega$={omega_deg}°/s   "
            f"Cost {cost_id}")

def plot_body_position_vs_time(runs, outdir, a_max_val, y_limit, omega_deg, cost_id):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    labels = ['$x_B$ (m)', '$y_B$ (m)', '$z_B$ (m)']
    for ax, idx, lab in zip(axes, range(3), labels):
        for i, r in enumerate(runs):
            ax.plot(r["time"], r["xb"][:, idx], color=IC_COLORS[i % 10],
                    lw=0.8, alpha=0.8, label=r["ic_label"] if idx == 0 else None)
        ax.set_ylabel(lab); ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("Time (s)")
    axes[0].legend(fontsize=6, ncol=3, loc='upper right')
    fig.suptitle(_case_title(a_max_val, y_limit, omega_deg, cost_id), fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outdir / "body_position_vs_time.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_control_and_dv_proxy(runs, outdir, a_max_val, y_limit, omega_deg, cost_id):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    # Top: control components
    ax = axes[0]
    comp_labels = ['$a_x$', '$a_y$', '$a_z$']
    for i, r in enumerate(runs):
        t_u = r["time"][:-1]
        for c in range(3):
            lbl = f"{r['ic_label']} {comp_labels[c]}" if i == 0 else None
            ax.plot(t_u, r["u"][:, c], color=IC_COLORS[i % 10],
                    lw=0.6, alpha=0.6)
    ax.axhline(a_max_val, color='r', ls='--', lw=0.7, alpha=0.5)
    ax.axhline(-a_max_val, color='r', ls='--', lw=0.7, alpha=0.5)
    ax.set_ylabel("Control (m/s²)"); ax.grid(True, alpha=0.3)
    ax.set_title("Control components (all ICs overlaid)")

    # Bottom: cumulative ΔV
    ax2 = axes[1]
    for i, r in enumerate(runs):
        ax2.plot(r["time"], r["dv_cumulative"], color=IC_COLORS[i % 10],
                 lw=0.8, alpha=0.8, label=r["ic_label"])
        # Annotate total ΔV at end
        dv_t = r["dv_total"]
        ax2.annotate(f"{dv_t:.1f}", xy=(r["time"][-1], dv_t),
                     fontsize=5, color=IC_COLORS[i % 10], alpha=0.7)
    ax2.set_ylabel("Cumulative $\\Delta V$ (m/s)"); ax2.set_xlabel("Time (s)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=5, ncol=3, loc='upper left')

    fig.suptitle(_case_title(a_max_val, y_limit, omega_deg, cost_id), fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outdir / "control_and_dv_proxy.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_mpc_control_vs_time(runs, outdir, a_max_val, y_limit, omega_deg, cost_id):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    comp_labels = ['$a_x$ (m/s²)', '$a_y$ (m/s²)', '$a_z$ (m/s²)']
    for ax, c, lab in zip(axes, range(3), comp_labels):
        for i, r in enumerate(runs):
            t_u = r["time"][:-1]
            ax.plot(t_u, r["u"][:, c], color=IC_COLORS[i % 10],
                    lw=0.7, alpha=0.7, label=r["ic_label"] if c == 0 else None)
        ax.set_ylabel(lab); ax.grid(True, alpha=0.3)
        ax.axhline(a_max_val, color='r', ls='--', lw=0.5, alpha=0.4)
        ax.axhline(-a_max_val, color='r', ls='--', lw=0.5, alpha=0.4)
    axes[-1].set_xlabel("Time (s)")
    axes[0].legend(fontsize=5, ncol=3, loc='upper right')
    fig.suptitle(_case_title(a_max_val, y_limit, omega_deg, cost_id), fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outdir / "mpc_control_vs_time.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_mpc_cost_breakdown(runs, outdir, a_max_val, y_limit, omega_deg, cost_id):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 1, figsize=(10, 9), sharex=True)
    titles = ['Position cost', 'Velocity cost', 'Delta-u cost', 'Constraint penalty']
    keys   = ['pos_costs', 'vel_costs', 'du_costs', 'con_pens']
    for ax, title, key in zip(axes, titles, keys):
        for i, r in enumerate(runs):
            t_u = r["time"][:-1]
            ax.plot(t_u, r[key], color=IC_COLORS[i % 10],
                    lw=0.6, alpha=0.7, label=r["ic_label"] if key == 'pos_costs' else None)
        ax.set_ylabel(title); ax.grid(True, alpha=0.3)
        ax.set_yscale('symlog', linthresh=1.0)
    axes[-1].set_xlabel("Time (s)")
    axes[0].legend(fontsize=5, ncol=3, loc='upper right')
    fig.suptitle(_case_title(a_max_val, y_limit, omega_deg, cost_id), fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outdir / "mpc_cost_breakdown.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_range_vs_time(runs, outdir, a_max_val, y_limit, omega_deg, cost_id, fname="mpc_range_vs_time.pdf"):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, r in enumerate(runs):
        ax.plot(r["time"], r["ranges"], color=IC_COLORS[i % 10],
                lw=0.8, alpha=0.8, label=r["ic_label"])
    ax.axhline(dock_range, color='r', ls='--', lw=0.8, label="Dock range")
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Range (m)")
    ax.set_title(_case_title(a_max_val, y_limit, omega_deg, cost_id))
    ax.legend(fontsize=6, ncol=3); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / fname, bbox_inches="tight")
    plt.close(fig)

def plot_traj_body_3d(runs, outdir, a_max_val, y_limit, omega_deg, cost_id):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection='3d')
    for i, r in enumerate(runs):
        xb = r["xb"]
        ax.plot(xb[:,0], xb[:,1], xb[:,2], color=IC_COLORS[i % 10],
                lw=0.7, alpha=0.7, label=r["ic_label"])
        ax.scatter([xb[0,0]], [xb[0,1]], [xb[0,2]], c=IC_COLORS[i % 10],
                   marker='o', s=20, alpha=0.8)
    # Draw LOS cone edges
    yf = max(10., 1.2 * max(r["xb"][:, 1].max() for r in runs))
    for sx in [1, -1]:
        for sz in [1, -1]:
            xc = sx * (x_0_los + c_x * (yf - y_min))
            zc = sz * (z_0_los + c_z * (yf - y_min))
            ax.plot([sx*x_0_los, xc], [y_min, yf], [sz*z_0_los, zc],
                    'r--', lw=0.6, alpha=0.4)
    ax.scatter([0], [0], [0], c='k', marker='s', s=50)
    ax.set_xlabel("$x_B$ (m)"); ax.set_ylabel("$y_B$ (m)"); ax.set_zlabel("$z_B$ (m)")
    ax.set_title(_case_title(a_max_val, y_limit, omega_deg, cost_id), fontsize=9)
    ax.legend(fontsize=5, ncol=2, loc='upper left')
    fig.tight_layout()
    fig.savefig(outdir / "mpc_traj_body_3d.pdf", bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------------------
# GIF generators
# ---------------------------------------------------------------------------
def make_approach_body_xy_gif(runs, outdir, a_max_val, y_limit, omega_deg, cost_id):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim

    fig, ax = plt.subplots(figsize=(6, 5), dpi=72)
    yf = max(10., 1.2 * max(float(np.max(np.abs(r["xb"][:, :2]))) for r in runs))
    for sx in [1, -1]:
        xc = sx*(x_0_los + c_x*(yf - y_min))
        ax.plot([sx*x_0_los, xc], [y_min, yf], "r--", lw=1, alpha=0.5)
    ax.plot(0, 0, "ks", ms=8)
    ax.set_xlim(-yf*0.8, yf*0.8); ax.set_ylim(-yf*0.05, yf*1.05)
    ax.set_xlabel("$x_B$ (m)"); ax.set_ylabel("$y_B$ (m)")
    ax.set_title(f"Body x-y  {case_dir_name(a_max_val, y_limit, omega_deg, cost_id)}", fontsize=8)
    ax.set_aspect("equal"); ax.grid(True, alpha=0.2)

    trails = []; dots = []
    for i in range(len(runs)):
        tr, = ax.plot([], [], color=IC_COLORS[i % 10], lw=0.8, alpha=0.7)
        dt, = ax.plot([], [], 'o', color=IC_COLORS[i % 10], ms=4)
        trails.append(tr); dots.append(dt)
    txt = ax.text(0.02, 0.98, "", transform=ax.transAxes, va="top",
                  fontsize=7, family="monospace",
                  bbox=dict(fc="white", alpha=0.7, ec="none"))

    n_pts = len(runs[0]["time"])
    skip = max(1, n_pts // 120)
    idxs = list(range(0, n_pts, skip))

    def update(fi):
        i = idxs[fi]
        for j, r in enumerate(runs):
            trails[j].set_data(r["xb"][:i+1, 0], r["xb"][:i+1, 1])
            dots[j].set_data([r["xb"][i, 0]], [r["xb"][i, 1]])
        txt.set_text(f"t={runs[0]['time'][i]:5.0f}s")
        return trails + dots + [txt]

    ani = anim.FuncAnimation(fig, update, frames=len(idxs), blit=True, interval=67)
    ani.save(str(outdir / "approach_body_xy.gif"), writer="pillow", fps=15)
    plt.close(fig)

def _cone_edges_body(yf_cone):
    """Return 4 cone edge lines in body frame: [(x1,y1,z1),(x2,y2,z2)] each."""
    edges = []
    for sx in [1, -1]:
        for sz in [1, -1]:
            x1 = sx * x_0_los; z1 = sz * z_0_los
            x2 = sx * (x_0_los + c_x * (yf_cone - y_min))
            z2 = sz * (z_0_los + c_z * (yf_cone - y_min))
            edges.append(np.array([[x1, y_min, z1], [x2, yf_cone, z2]]))
    # base ring (4 edges connecting near-end corners)
    corners_near = [
        np.array([sx * x_0_los, y_min, sz * z_0_los])
        for sx in [1, -1] for sz in [1, -1]
    ]
    # order: (+,+), (+,-), (-,+), (-,-) -> connect as rectangle
    order = [0, 1, 3, 2, 0]
    for k in range(4):
        edges.append(np.array([corners_near[order[k]], corners_near[order[k+1]]]))
    # far ring
    corners_far = [
        np.array([sx * (x_0_los + c_x*(yf_cone-y_min)), yf_cone,
                   sz * (z_0_los + c_z*(yf_cone-y_min))])
        for sx in [1, -1] for sz in [1, -1]
    ]
    for k in range(4):
        edges.append(np.array([corners_far[order[k]], corners_far[order[k+1]]]))
    return edges

def _rotate_edges_to_lvlh(edges_body, theta):
    """Rotate body-frame edges to LVLH at angle theta."""
    R = rotz(theta)
    return [R @ e.T for e in edges_body]  # each is 3×2

def make_approach_lvlh_3d_gif(runs, outdir, a_max_val, y_limit, omega_deg, cost_id):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim

    omega_t = np.radians(omega_deg)
    fig = plt.figure(figsize=(7, 6), dpi=72)
    ax = fig.add_subplot(111, projection='3d')
    lim = max(10., 1.2 * max(float(np.max(np.abs(r["x_lvlh"][:, :3]))) for r in runs))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim*0.3, lim*0.3)
    ax.scatter([0], [0], [0], c='k', marker='s', s=40)
    ax.set_xlabel("$x_L$ (m)"); ax.set_ylabel("$y_L$ (m)"); ax.set_zlabel("$z_L$ (m)")
    ax.set_title(f"LVLH 3D  {case_dir_name(a_max_val, y_limit, omega_deg, cost_id)}", fontsize=8)

    # Cone edges in body frame
    yf_cone = min(lim * 0.8, max(y_limit * 1.2, 30.))
    edges_body = _cone_edges_body(yf_cone)

    # Pre-draw cone lines (will update each frame)
    cone_lines = []
    for _ in edges_body:
        ln, = ax.plot([], [], [], 'r-', lw=0.7, alpha=0.5)
        cone_lines.append(ln)

    trails = []; dots = []
    for i in range(len(runs)):
        tr, = ax.plot([], [], [], color=IC_COLORS[i % 10], lw=0.7, alpha=0.7)
        dt  = ax.plot([], [], [], 'o', color=IC_COLORS[i % 10], ms=4)[0]
        trails.append(tr); dots.append(dt)

    n_pts = len(runs[0]["time"])
    skip = max(1, n_pts // 120)
    idxs = list(range(0, n_pts, skip))

    def update(fi):
        i = idxs[fi]
        t = runs[0]["time"][i]
        theta = omega_t * t
        # Update rotating cone
        rotated = _rotate_edges_to_lvlh(edges_body, theta)
        for k, edge_lvlh in enumerate(rotated):
            cone_lines[k].set_data(edge_lvlh[0, :], edge_lvlh[1, :])
            cone_lines[k].set_3d_properties(edge_lvlh[2, :])
        # Update trajectories
        for j, r in enumerate(runs):
            xl = r["x_lvlh"]
            trails[j].set_data(xl[:i+1, 0], xl[:i+1, 1])
            trails[j].set_3d_properties(xl[:i+1, 2])
            dots[j].set_data([xl[i, 0]], [xl[i, 1]])
            dots[j].set_3d_properties([xl[i, 2]])
        return cone_lines + trails + dots

    ani = anim.FuncAnimation(fig, update, frames=len(idxs), blit=False, interval=67)
    ani.save(str(outdir / "approach_lvlh_3d.gif"), writer="pillow", fps=15)
    plt.close(fig)

def make_approach_lvlh_xy_r50_gif(runs, outdir, a_max_val, y_limit, omega_deg, cost_id):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim
    from matplotlib.patches import Polygon as MplPolygon

    omega_t = np.radians(omega_deg)
    fig, ax = plt.subplots(figsize=(6, 5), dpi=72)
    lim = max(60., 1.2 * max(float(np.max(np.abs(r["x_lvlh"][:, :2]))) for r in runs))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim*0.1, lim*1.1)
    ax.plot(0, 0, "ks", ms=8)
    # Draw 50m circle
    th_c = np.linspace(0, 2*np.pi, 100)
    ax.plot(50*np.cos(th_c), 50*np.sin(th_c), 'k--', lw=0.8, alpha=0.4, label="r=50m")
    ax.set_xlabel("$x_L$ (m)"); ax.set_ylabel("$y_L$ (m)")
    ax.set_title(f"LVLH x-y  {case_dir_name(a_max_val, y_limit, omega_deg, cost_id)}", fontsize=7)
    ax.set_aspect("equal"); ax.grid(True, alpha=0.2)

    # Cone x-y cross-section in body frame (z=0 slice)
    yf_cone = min(lim * 0.85, max(y_limit * 1.2, 30.))
    cone_body_xy = np.array([
        [-x_0_los, y_min],
        [x_0_los, y_min],
        [x_0_los + c_x*(yf_cone - y_min), yf_cone],
        [-(x_0_los + c_x*(yf_cone - y_min)), yf_cone],
    ])

    # Rotating cone polygon (filled, semi-transparent)
    cone_patch = MplPolygon(cone_body_xy, closed=True, fc='red', alpha=0.08,
                             ec='red', lw=1.0, ls='--')
    ax.add_patch(cone_patch)
    # Rotating cone edge lines (left & right)
    cone_l, = ax.plot([], [], 'r--', lw=1.0, alpha=0.6)
    cone_r, = ax.plot([], [], 'r--', lw=1.0, alpha=0.6)
    # Docking axis line (+y_B direction)
    dock_line, = ax.plot([], [], 'r-', lw=0.5, alpha=0.3)

    trails = []; dots = []
    for i in range(len(runs)):
        tr, = ax.plot([], [], color=IC_COLORS[i % 10], lw=0.8, alpha=0.7)
        dt, = ax.plot([], [], 'o', color=IC_COLORS[i % 10], ms=4)
        trails.append(tr); dots.append(dt)
    txt = ax.text(0.02, 0.98, "", transform=ax.transAxes, va="top",
                  fontsize=7, family="monospace",
                  bbox=dict(fc="white", alpha=0.7, ec="none"))

    n_pts = len(runs[0]["time"])
    skip = max(1, n_pts // 120)
    idxs = list(range(0, n_pts, skip))

    def update(fi):
        i = idxs[fi]
        t = runs[0]["time"][i]
        theta = omega_t * t
        R2 = np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta),  np.cos(theta)]])
        # Rotate cone vertices to LVLH
        cone_lvlh = (R2 @ cone_body_xy.T).T
        cone_patch.set_xy(cone_lvlh)
        # Cone edge lines
        cone_l.set_data([cone_lvlh[0,0], cone_lvlh[3,0]],
                        [cone_lvlh[0,1], cone_lvlh[3,1]])
        cone_r.set_data([cone_lvlh[1,0], cone_lvlh[2,0]],
                        [cone_lvlh[1,1], cone_lvlh[2,1]])
        # Docking axis (+y_B in LVLH)
        dock_dir = R2 @ np.array([0, yf_cone])
        dock_line.set_data([0, dock_dir[0]], [0, dock_dir[1]])
        # Update trajectories
        for j, r in enumerate(runs):
            xl = r["x_lvlh"]
            trails[j].set_data(xl[:i+1, 0], xl[:i+1, 1])
            dots[j].set_data([xl[i, 0]], [xl[i, 1]])
        txt.set_text(f"t={t:5.0f}s")
        return [cone_patch, cone_l, cone_r, dock_line] + trails + dots + [txt]

    ani = anim.FuncAnimation(fig, update, frames=len(idxs), blit=False, interval=67)
    ani.save(str(outdir / "approach_lvlh_xy_r0_50m.gif"), writer="pillow", fps=15)
    plt.close(fig)

# ---------------------------------------------------------------------------
# Save case config
# ---------------------------------------------------------------------------
def save_case_config(outdir, a_max_val, y_limit, omega_deg, cost_id, r_sync):
    cfg = {
        "a_max": a_max_val,
        "y_limit": y_limit,
        "omega_deg_s": omega_deg,
        "cost_variant": cost_id,
        "dt_truth_s": ts_truth,
        "dt_mpc_s": ts_mpc,
        "horizon": N_hor,
        "sim_duration_s": sim_dur,
        "n_orbit_rad_s": n_orbit,
        "r_sync_m": r_sync,
        "dock_range_m": dock_range,
        "dock_speed_mps": dock_speed,
        "los_margin_m": los_margin,
        "solver": "OSQP",
        "solver_settings": {
            "eps_abs": 1e-4, "eps_rel": 1e-4,
            "max_iter": 4000, "polish": True,
            "adaptive_rho": True, "scaling": 10,
        },
        "constraints": {
            "LOS_cone": True, "thrust_saturation": True,
            "x0_los": x_0_los, "z0_los": z_0_los,
            "cx": c_x, "cz": c_z, "y_min": y_min,
        },
        "cost_weights": {
            "q_pos": q_pos.tolist(),
            "q_vel": q_vel.tolist() if cost_id == 'B' else [0,0,0],
            "r_du": r_du_A if cost_id == 'A' else r_du_B,
            "r_u": 0.0,
            "w_slack": w_slack,
        },
        "model": "HCW (linearised CWH)",
        "truth_model": "HCW (same as prediction; TODO: add nonlinear for J2+drag)",
        "n_initial_conditions": 10,
    }
    with open(outdir / "case_config.json", "w") as f:
        json.dump(cfg, f, indent=2)

# ---------------------------------------------------------------------------
# Run one case (all 10 ICs) and generate all outputs
# ---------------------------------------------------------------------------
def _save_runs_npz(runs, outdir):
    """Save trajectory data as .npz for later replotting."""
    arrays = {}
    for i, r in enumerate(runs):
        arrays[f"xb_{i}"]    = r["xb"]
        arrays[f"xlvlh_{i}"] = r["x_lvlh"]
        arrays[f"u_{i}"]     = r["u"]
        arrays[f"slk_{i}"]   = r["slk"]
        arrays[f"ranges_{i}"]= r["ranges"]
        arrays[f"speeds_{i}"]= r["speeds"]
        arrays[f"dvcum_{i}"] = r["dv_cumulative"]
        arrays[f"time_{i}"]  = r["time"]
        arrays[f"poscost_{i}"] = r["pos_costs"]
        arrays[f"velcost_{i}"] = r["vel_costs"]
        arrays[f"ducost_{i}"]  = r["du_costs"]
        arrays[f"conpen_{i}"]  = r["con_pens"]
    # Scalars stored as small arrays
    ic_labels = [r["ic_label"] for r in runs]
    meta = {
        "n_runs": len(runs),
        "dv_totals": [r["dv_total"] for r in runs],
        "r_finals": [r["r_final"] for r in runs],
        "v_finals": [r["v_final"] for r in runs],
        "r_mins": [r["r_min"] for r in runs],
        "n_viols": [r["n_viol"] for r in runs],
        "docks": [r["dock"] for r in runs],
        "n_infeasibles": [r["n_infeasible"] for r in runs],
        "max_laterals": [r["max_lateral"] for r in runs],
        "ic_labels": ic_labels,
    }
    arrays["meta_json"] = np.array([json.dumps(meta)])
    np.savez_compressed(outdir / "trajectories.npz", **arrays)

def _load_runs_npz(outdir):
    """Load trajectory data from .npz."""
    data = np.load(outdir / "trajectories.npz", allow_pickle=True)
    meta = json.loads(str(data["meta_json"][0]))
    n_runs = meta["n_runs"]
    runs = []
    for i in range(n_runs):
        r = {
            "xb": data[f"xb_{i}"],
            "x_lvlh": data[f"xlvlh_{i}"],
            "u": data[f"u_{i}"],
            "slk": data[f"slk_{i}"],
            "ranges": data[f"ranges_{i}"],
            "speeds": data[f"speeds_{i}"],
            "dv_cumulative": data[f"dvcum_{i}"],
            "time": data[f"time_{i}"],
            "pos_costs": data[f"poscost_{i}"],
            "vel_costs": data[f"velcost_{i}"],
            "du_costs": data[f"ducost_{i}"],
            "con_pens": data[f"conpen_{i}"],
            "ic_label": meta["ic_labels"][i],
            "dv_total": meta["dv_totals"][i],
            "r_final": meta["r_finals"][i],
            "v_final": meta["v_finals"][i],
            "r_min": meta["r_mins"][i],
            "n_viol": meta["n_viols"][i],
            "dock": meta["docks"][i],
            "n_infeasible": meta["n_infeasibles"][i],
            "max_lateral": meta["max_laterals"][i],
        }
        runs.append(r)
    return runs

def _generate_all_plots_gifs(runs, outdir, a_max_val, y_limit, omega_deg, cost_id):
    """Generate all PDF plots and GIF animations for one case."""
    plot_body_position_vs_time(runs, outdir, a_max_val, y_limit, omega_deg, cost_id)
    plot_control_and_dv_proxy(runs, outdir, a_max_val, y_limit, omega_deg, cost_id)
    plot_mpc_control_vs_time(runs, outdir, a_max_val, y_limit, omega_deg, cost_id)
    plot_mpc_cost_breakdown(runs, outdir, a_max_val, y_limit, omega_deg, cost_id)
    plot_range_vs_time(runs, outdir, a_max_val, y_limit, omega_deg, cost_id,
                       fname="mpc_range_vs_time.pdf")
    plot_traj_body_3d(runs, outdir, a_max_val, y_limit, omega_deg, cost_id)
    plot_range_vs_time(runs, outdir, a_max_val, y_limit, omega_deg, cost_id,
                       fname="range_vs_time.pdf")
    make_approach_body_xy_gif(runs, outdir, a_max_val, y_limit, omega_deg, cost_id)
    make_approach_lvlh_3d_gif(runs, outdir, a_max_val, y_limit, omega_deg, cost_id)
    make_approach_lvlh_xy_r50_gif(runs, outdir, a_max_val, y_limit, omega_deg, cost_id)

def _aggregate_stats(runs, a_max_val, y_limit, omega_deg, cost_id, r_sync):
    """Compute aggregate statistics from a list of run dicts."""
    dock_count = sum(1 for r in runs if r["dock"])
    dv_vals = [r["dv_total"] for r in runs]
    return {
        "a_max": a_max_val, "y_limit": y_limit,
        "omega_deg": omega_deg, "cost_variant": cost_id,
        "r_sync": r_sync,
        "success_count": dock_count,
        "min_range": min(r["r_min"] for r in runs),
        "max_lateral_excursion": max(r["max_lateral"] for r in runs),
        "total_dv_mean": float(np.mean(dv_vals)),
        "total_dv_std": float(np.std(dv_vals)),
        "constraint_violations_total": sum(r["n_viol"] for r in runs),
        "qp_failures_total": sum(r["n_infeasible"] for r in runs),
        "r_final_mean": float(np.mean([r["r_final"] for r in runs])),
        "r_final_std": float(np.std([r["r_final"] for r in runs])),
    }

def run_one_case_full(a_max_val, y_limit, omega_deg, cost_id, n_ics=10):
    """Run all ICs for one parameter combination, generate plots and GIFs."""
    omega_t = np.radians(omega_deg)
    r_sync  = 2 * a_max_val / omega_t**2
    ics = build_ics(float(y_limit))[:n_ics]
    outdir = case_dir(a_max_val, y_limit, omega_deg, cost_id)

    runs = []
    for ic_label, x_B0, y_B0, z_B0 in ics:
        r = run_one_sim(omega_deg, a_max_val, y_limit, cost_id,
                        ic_label, x_B0, y_B0, z_B0)
        runs.append(r)

    # Save config + trajectory data
    save_case_config(outdir, a_max_val, y_limit, omega_deg, cost_id, r_sync)
    _save_runs_npz(runs, outdir)

    # Generate all plots and GIFs
    _generate_all_plots_gifs(runs, outdir, a_max_val, y_limit, omega_deg, cost_id)

    return _aggregate_stats(runs, a_max_val, y_limit, omega_deg, cost_id, r_sync)

# ---------------------------------------------------------------------------
# Summary CSV
# ---------------------------------------------------------------------------
def write_summary_csv(all_stats):
    csv_path = RESULTS / "summary.csv"
    fieldnames = [
        "a_max", "y_limit", "omega_deg", "cost_variant", "r_sync",
        "success_count", "min_range", "max_lateral_excursion",
        "total_dv_mean", "total_dv_std",
        "constraint_violations_total", "qp_failures_total",
        "r_final_mean", "r_final_std",
    ]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for s in all_stats:
            writer.writerow(s)
    print(f"  Summary CSV: {csv_path}")


def write_summary_latex(all_stats):
    """Generate LaTeX table from sweep stats and save to results dir."""
    tex_path = RESULTS / "summary_table.tex"
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{MPC 90-case parameter sweep results "
                 r"(10 ICs per case, 600\,s simulation).}")
    lines.append(r"\label{tab:sweep_results}")
    lines.append(r"\footnotesize")
    lines.append(r"\begin{tabular}{cc cc c rr rr rr}")
    lines.append(r"\toprule")
    lines.append(r"$a_{\max}$ & $y_{\mathrm{lim}}$ & $\omega$ & Cost & "
                 r"$r_{\mathrm{sync}}$ & Dock & $\Delta V_{\mathrm{mean}}$ & "
                 r"$\Delta V_{\mathrm{std}}$ & $r_{f,\mathrm{mean}}$ & "
                 r"$r_{f,\mathrm{std}}$ & Viol. & QP \\")
    lines.append(r"(m/s$^2$) & (m) & ($^\circ$/s) & & (m) & /10 & "
                 r"(m/s) & (m/s) & (m) & (m) & & fail \\")
    lines.append(r"\midrule")
    prev_cost = None
    for s in all_stats:
        if prev_cost is not None and s["cost_variant"] != prev_cost:
            lines.append(r"\midrule")
        prev_cost = s["cost_variant"]
        dock_str = f"{s['success_count']}"
        # Highlight full-success rows
        if s["success_count"] == 10:
            dock_str = r"\textbf{10}"
        line = (f"{s['a_max']:.2f} & {s['y_limit']} & "
                f"{s['omega_deg']} & {s['cost_variant']} & "
                f"{s['r_sync']:.0f} & {dock_str} & "
                f"{s['total_dv_mean']:.1f} & {s['total_dv_std']:.1f} & "
                f"{s['r_final_mean']:.1f} & {s['r_final_std']:.1f} & "
                f"{s['constraint_violations_total']} & "
                f"{s['qp_failures_total']} \\\\")
        lines.append(line)
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    tex = "\n".join(lines) + "\n"
    with open(tex_path, "w") as f:
        f.write(tex)
    print(f"  LaTeX table: {tex_path}")
    return tex_path

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate", action="store_true",
                        help="Run validation: 1 case, 2 ICs only")
    parser.add_argument("--replot", action="store_true",
                        help="Replot from saved .npz data (no re-simulation)")
    args = parser.parse_args()

    RESULTS.mkdir(parents=True, exist_ok=True)

    if args.validate:
        print("=" * 60)
        print("VALIDATION RUN: a=0.10, y=30, w=1, costA, 2 ICs")
        print("=" * 60)
        t0 = _timer.perf_counter()
        stats = run_one_case_full(0.10, 30, 1, 'A', n_ics=2)
        dt = _timer.perf_counter() - t0
        print(f"\nValidation complete in {dt:.1f}s")
        print(f"  Dock success: {stats['success_count']}/2")
        print(f"  dV mean: {stats['total_dv_mean']:.2f} +/- {stats['total_dv_std']:.2f} m/s")
        print(f"  Min range: {stats['min_range']:.2f} m")
        print(f"  Max lateral: {stats['max_lateral_excursion']:.2f} m")
        print(f"  LOS violations: {stats['constraint_violations_total']}")
        print(f"  QP failures: {stats['qp_failures_total']}")
        outdir = case_dir(0.10, 30, 1, 'A')
        print(f"\n  Output dir: {outdir}")
        for f in sorted(outdir.iterdir()):
            print(f"    {f.name}  ({f.stat().st_size/1024:.0f} KB)")
        return

    if args.replot:
        print("=" * 72)
        print("REPLOT MODE: regenerating plots/GIFs from saved .npz data")
        print("=" * 72)
        wall0 = _timer.perf_counter()
        all_stats = []
        case_num = 0
        total_cases = len(a_maxes) * len(y_limits) * len(omegas_deg) * 2
        n_replotted = 0; n_skipped = 0
        for cost_id in ['A', 'B']:
            for omega in omegas_deg:
                for am in a_maxes:
                    for yl in y_limits:
                        case_num += 1
                        tag = case_dir_name(am, yl, omega, cost_id)
                        outdir = case_dir(am, yl, omega, cost_id)
                        npz_file = outdir / "trajectories.npz"
                        if not npz_file.exists():
                            print(f"  [{case_num:2d}/{total_cases}] {tag}  "
                                  f"SKIP (no .npz)")
                            n_skipped += 1
                            continue
                        t0 = _timer.perf_counter()
                        runs = _load_runs_npz(outdir)
                        _generate_all_plots_gifs(runs, outdir, am, yl, omega, cost_id)
                        omega_t = np.radians(omega)
                        r_sync = 2*am/omega_t**2
                        stats = _aggregate_stats(runs, am, yl, omega, cost_id, r_sync)
                        all_stats.append(stats)
                        dt = _timer.perf_counter() - t0
                        print(f"  [{case_num:2d}/{total_cases}] {tag}  "
                              f"dock={stats['success_count']}/10  [{dt:.0f}s]")
                        n_replotted += 1
        if all_stats:
            write_summary_csv(all_stats)
            write_summary_latex(all_stats)
        wall_total = _timer.perf_counter() - wall0
        print(f"\nReplot done: {n_replotted} cases replotted, "
              f"{n_skipped} skipped (no .npz)")
        print(f"Total wall time: {wall_total:.0f}s ({wall_total/60:.1f} min)")
        return

    # Full sweep
    total_cases = len(a_maxes) * len(y_limits) * len(omegas_deg) * 2
    print("=" * 72)
    print(f"MPC SWEEP  ({total_cases} cases × 10 ICs = {total_cases*10} sims)")
    print(f"  a_max    = {a_maxes} m/s²")
    print(f"  y_limit  = {y_limits} m")
    print(f"  omega    = {omegas_deg} deg/s")
    print(f"  cost     = [A: pos+du, B: pos+vel+du]")
    print(f"  sim      = {sim_dur}s ({int(sim_dur/ts_mpc)} MPC steps)")
    print(f"  ICs/case = 10")
    print("=" * 72)

    wall0 = _timer.perf_counter()
    all_stats = []
    case_num = 0

    for cost_id in ['A', 'B']:
        for omega in omegas_deg:
            for am in a_maxes:
                for yl in y_limits:
                    case_num += 1
                    tag = case_dir_name(am, yl, omega, cost_id)
                    omega_t = np.radians(omega)
                    r_sync = 2*am/omega_t**2
                    print(f"\n  [{case_num:2d}/{total_cases}] {tag}  "
                          f"r_sync={r_sync:.1f}m", end="  ", flush=True)

                    t0 = _timer.perf_counter()
                    stats = run_one_case_full(am, yl, omega, cost_id, n_ics=10)
                    dt = _timer.perf_counter() - t0

                    print(f"dock={stats['success_count']}/10  "
                          f"dv={stats['total_dv_mean']:.1f}±{stats['total_dv_std']:.1f}  "
                          f"viol={stats['constraint_violations_total']}  "
                          f"[{dt:.0f}s]")
                    all_stats.append(stats)

    # Write summary
    print("\n" + "=" * 72)
    write_summary_csv(all_stats)
    write_summary_latex(all_stats)

    wall_total = _timer.perf_counter() - wall0
    print(f"\nTotal wall time: {wall_total:.0f}s ({wall_total/60:.1f} min)")

if __name__ == "__main__":
    main()
