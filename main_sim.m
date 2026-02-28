%MAIN_SIM  MPC-based rendezvous with tumbling target.
%
%  Architecture:
%    - Plant (truth): nonlinear ECI 2-body + J2, integrated with ode113.
%    - MPC prediction: LTV discrete model in target-body frame.
%    - Solver: OSQP (QP with polyhedral LOS cone).
%    - LOS cone is body-fixed to the tumbling target (+yT docking axis).
%
%  Run:  >> main_sim
%
%  Iskender O.B. -- IAC 2026 -- Refactored from HCW/CVX monolith.

clear; close all; clc;

%% ===== Add subfolders to path =====
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'dynamics'));
addpath(fullfile(this_dir, 'frames'));
addpath(fullfile(this_dir, 'mpc'));
addpath(fullfile(this_dir, 'viz'));
addpath(fullfile(this_dir, 'utils'));

%% ===== Load parameters =====
p = params();

%% ===== Initialization =====
Nsteps = floor(p.Tsim / p.dt);

% Target initial ECI state
r_tgt_eci = p.r_tgt0_eci;
v_tgt_eci = p.v_tgt0_eci;

% LVLH frame at t=0
R_eci_lvlh_0 = lvlh_from_rv(r_tgt_eci, v_tgt_eci);

% Target body frame starts aligned with LVLH
R_eci_tb = R_eci_lvlh_0;
q_tb = quat_from_rotm_main(R_eci_tb);

% Chaser initial ECI state (from LVLH relative state)
r_chs_eci = r_tgt_eci + R_eci_lvlh_0 * p.dr_lvlh0;
v_chs_eci = v_tgt_eci + R_eci_lvlh_0 * p.dv_lvlh0;

% Initial relative state in target-body frame
[r_tb, v_tb] = transform_rel_to_TB(r_chs_eci, v_chs_eci, ...
                                    r_tgt_eci, v_tgt_eci, ...
                                    R_eci_tb, p.omega_body);
x_tb = [r_tb; v_tb];

% Previous input (for delta-u penalty)
u_prev = zeros(3, 1);

%% ===== Logging arrays =====
% Target body frame
lg.r_tb_hist       = zeros(3, Nsteps+1);
lg.v_tb_hist       = zeros(3, Nsteps+1);
% LVLH frame
lg.r_lvlh_hist     = zeros(3, Nsteps+1);
lg.v_lvlh_hist     = zeros(3, Nsteps+1);
% Chaser body frame (assumed = TB for now; extensible)
lg.r_cb_hist       = zeros(3, Nsteps+1);
lg.v_cb_hist       = zeros(3, Nsteps+1);
% ECI (absolute)
lg.r_tgt_eci_hist  = zeros(3, Nsteps+1);
lg.v_tgt_eci_hist  = zeros(3, Nsteps+1);
lg.r_chs_eci_hist  = zeros(3, Nsteps+1);
lg.v_chs_eci_hist  = zeros(3, Nsteps+1);
% Inputs
lg.u_hist          = zeros(3, Nsteps);
lg.u_lvlh_hist     = zeros(3, Nsteps);
% Rotations
lg.R_eci_tb_hist   = zeros(3, 3, Nsteps+1);
lg.R_eci_lvlh_hist = zeros(3, 3, Nsteps+1);
% Docking axis in LVLH (for plotting the rotating axis)
lg.dock_axis_lvlh_hist = zeros(3, Nsteps+1);
% Time, solver, cost
lg.t_hist          = zeros(1, Nsteps+1);
lg.status_hist     = cell(1, Nsteps);
lg.solve_time_hist = zeros(1, Nsteps);
lg.cost_hist       = zeros(1, Nsteps);
lg.ref_hist        = zeros(6, Nsteps+1);

% Compute initial LVLH relative state
R_eci_lvlh = R_eci_lvlh_0;
dr_eci = r_chs_eci - r_tgt_eci;
dv_eci = v_chs_eci - v_tgt_eci;
r_lvlh = R_eci_lvlh' * dr_eci;
v_lvlh = R_eci_lvlh' * dv_eci;

% Docking axis (+yT) expressed in LVLH
dock_axis_lvlh = R_eci_lvlh' * R_eci_tb(:,2);

% Store initial states
lg.r_tb_hist(:,1)         = r_tb;
lg.v_tb_hist(:,1)         = v_tb;
lg.r_lvlh_hist(:,1)       = r_lvlh;
lg.v_lvlh_hist(:,1)       = v_lvlh;
lg.r_cb_hist(:,1)         = r_tb;   % chaser body = TB initially
lg.v_cb_hist(:,1)         = v_tb;
lg.r_tgt_eci_hist(:,1)    = r_tgt_eci;
lg.v_tgt_eci_hist(:,1)    = v_tgt_eci;
lg.r_chs_eci_hist(:,1)    = r_chs_eci;
lg.v_chs_eci_hist(:,1)    = v_chs_eci;
lg.R_eci_tb_hist(:,:,1)   = R_eci_tb;
lg.R_eci_lvlh_hist(:,:,1) = R_eci_lvlh;
lg.dock_axis_lvlh_hist(:,1) = dock_axis_lvlh;
lg.t_hist(1)              = 0;
lg.ref_hist(:,1)          = compute_reference(0, p);

%% ===== Build initial OSQP problem =====
fprintf('Linearizing initial model...\n');
[Ad, Bd] = linearize_discrete_model(x_tb, R_eci_tb, ...
                r_tgt_eci, v_tgt_eci, p.omega_body, ...
                zeros(3,1), p.dt, p);

x_ref = compute_reference(0, p);
fprintf('Building OSQP problem...\n');
[prob, qp_data] = build_qp_osqp(Ad, Bd, x_tb, x_ref, u_prev, p);

%% ===== Main simulation loop =====
fprintf('Starting simulation: %d steps, dt=%.1f s, Np=%d\n', Nsteps, p.dt, p.Np);
fprintf('==========================================================\n');

sim_terminated = false;
final_step = Nsteps;

for k = 1:Nsteps
    t_now = (k-1) * p.dt;

    % ---- Check LOS violation ----
    [los_ok, los_margin] = check_los(r_tb, p.cone_k, p.y_min);
    if ~los_ok
        fprintf('[TERMINATED] LOS violation at t=%.1f s\n', t_now);
        fprintf('  r_TB = [%.3f, %.3f, %.3f] m\n', r_tb(1), r_tb(2), r_tb(3));
        fprintf('  LOS margin = %.3f m\n', los_margin);
        final_step = k - 1;
        sim_terminated = true;
        break;
    end

    % ---- Reference ----
    x_ref = compute_reference(t_now, p);

    % ---- Linearize ----
    [Ad, Bd] = linearize_discrete_model(x_tb, R_eci_tb, ...
                    r_tgt_eci, v_tgt_eci, p.omega_body, ...
                    u_prev, p.dt, p);

    % ---- Solve MPC QP ----
    tic;
    [u_opt, status, x_pred] = update_qp_osqp(prob, qp_data, Ad, Bd, ...
                            x_tb, x_ref, u_prev, p);
    solve_time = toc;

    lg.status_hist{k}     = status;
    lg.solve_time_hist(k) = solve_time;

    if ~strcmp(status, 'solved') && ~strcmp(status, 'solved_inaccurate')
        fprintf('[TERMINATED] OSQP infeasible at t=%.1f s  (status: %s)\n', t_now, status);
        fprintf('  r_TB = [%.3f, %.3f, %.3f] m\n', r_tb(1), r_tb(2), r_tb(3));
        final_step = k - 1;
        sim_terminated = true;
        break;
    end

    % Compute stage cost for logging
    du = u_opt - u_prev;
    J_k = (x_tb - x_ref)' * p.Q * (x_tb - x_ref) + ...
          u_opt' * p.Ru * u_opt + du' * p.Rdu * du;
    lg.cost_hist(k) = J_k;

    % Clamp input
    u_applied = clamp(u_opt, -p.u_max, p.u_max);
    lg.u_hist(:, k) = u_applied;

    % Input in LVLH frame
    u_eci = control_TB_to_ECI(u_applied, R_eci_tb);
    lg.u_lvlh_hist(:, k) = R_eci_lvlh' * u_eci;

    % ---- Propagate plant (nonlinear, ECI) ----
    a_ctrl_eci = control_TB_to_ECI(u_applied, R_eci_tb);

    x_tgt_next = propagate_truth_step([r_tgt_eci; v_tgt_eci], ...
                    p.dt, p.mu, p.Re, p.J2, [0;0;0], p.ode_opts);
    r_tgt_eci = x_tgt_next(1:3);
    v_tgt_eci = x_tgt_next(4:6);

    x_chs_next = propagate_truth_step([r_chs_eci; v_chs_eci], ...
                    p.dt, p.mu, p.Re, p.J2, a_ctrl_eci, p.ode_opts);
    r_chs_eci = x_chs_next(1:3);
    v_chs_eci = x_chs_next(4:6);

    % ---- Propagate target attitude ----
    [q_tb, R_eci_tb] = target_attitude_model(q_tb, p.omega_body, p.dt);

    % ---- Compute LVLH ----
    R_eci_lvlh = lvlh_from_rv(r_tgt_eci, v_tgt_eci);

    % ---- Transform to target body frame ----
    [r_tb, v_tb] = transform_rel_to_TB(r_chs_eci, v_chs_eci, ...
                                        r_tgt_eci, v_tgt_eci, ...
                                        R_eci_tb, p.omega_body);
    x_tb = [r_tb; v_tb];

    % ---- Transform to LVLH frame ----
    dr_eci = r_chs_eci - r_tgt_eci;
    dv_eci = v_chs_eci - v_tgt_eci;
    r_lvlh = R_eci_lvlh' * dr_eci;
    v_lvlh = R_eci_lvlh' * dv_eci;

    % ---- Docking axis in LVLH ----
    dock_axis_lvlh = R_eci_lvlh' * R_eci_tb(:,2);

    % ---- Store for delta-u ----
    u_prev = u_applied;

    % ---- Log ----
    lg.r_tb_hist(:, k+1)         = r_tb;
    lg.v_tb_hist(:, k+1)         = v_tb;
    lg.r_lvlh_hist(:, k+1)       = r_lvlh;
    lg.v_lvlh_hist(:, k+1)       = v_lvlh;
    lg.r_cb_hist(:, k+1)         = r_tb;   % chaser body ~ TB
    lg.v_cb_hist(:, k+1)         = v_tb;
    lg.r_tgt_eci_hist(:, k+1)    = r_tgt_eci;
    lg.v_tgt_eci_hist(:, k+1)    = v_tgt_eci;
    lg.r_chs_eci_hist(:, k+1)    = r_chs_eci;
    lg.v_chs_eci_hist(:, k+1)    = v_chs_eci;
    lg.R_eci_tb_hist(:,:,k+1)    = R_eci_tb;
    lg.R_eci_lvlh_hist(:,:,k+1)  = R_eci_lvlh;
    lg.dock_axis_lvlh_hist(:,k+1) = dock_axis_lvlh;
    lg.t_hist(k+1)               = k * p.dt;
    lg.ref_hist(:, k+1)          = compute_reference(k*p.dt, p);

    % ---- Print progress ----
    if mod(k, 10) == 0 || k == 1
        fprintf('  t=%6.1f s | y_T=%8.2f m | rad_dev=%7.2f m | ||u||=%.4f | J=%.2f | %s (%.3f s)\n', ...
            t_now, r_tb(2), sqrt(r_tb(1)^2+r_tb(3)^2), norm(u_applied), J_k, status, solve_time);
    end
end

%% ===== Trim logs =====
N = final_step;
if sim_terminated
    lg.r_tb_hist   = lg.r_tb_hist(:, 1:N+1);
    lg.v_tb_hist   = lg.v_tb_hist(:, 1:N+1);
    lg.r_lvlh_hist = lg.r_lvlh_hist(:, 1:N+1);
    lg.v_lvlh_hist = lg.v_lvlh_hist(:, 1:N+1);
    lg.r_cb_hist   = lg.r_cb_hist(:, 1:N+1);
    lg.v_cb_hist   = lg.v_cb_hist(:, 1:N+1);
    lg.r_tgt_eci_hist = lg.r_tgt_eci_hist(:, 1:N+1);
    lg.v_tgt_eci_hist = lg.v_tgt_eci_hist(:, 1:N+1);
    lg.r_chs_eci_hist = lg.r_chs_eci_hist(:, 1:N+1);
    lg.v_chs_eci_hist = lg.v_chs_eci_hist(:, 1:N+1);
    lg.u_hist      = lg.u_hist(:, 1:N);
    lg.u_lvlh_hist = lg.u_lvlh_hist(:, 1:N);
    lg.t_hist      = lg.t_hist(1:N+1);
    lg.R_eci_tb_hist   = lg.R_eci_tb_hist(:,:,1:N+1);
    lg.R_eci_lvlh_hist = lg.R_eci_lvlh_hist(:,:,1:N+1);
    lg.dock_axis_lvlh_hist = lg.dock_axis_lvlh_hist(:,1:N+1);
    lg.status_hist = lg.status_hist(1:N);
    lg.solve_time_hist = lg.solve_time_hist(1:N);
    lg.cost_hist   = lg.cost_hist(1:N);
    lg.ref_hist    = lg.ref_hist(:, 1:N+1);
else
    fprintf('==========================================================\n');
    fprintf('Simulation completed successfully (t=%.1f s).\n', p.Tsim);
end

%% ===== Save results =====
results_dir = fullfile(this_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
save(fullfile(results_dir, 'sim_results.mat'), 'lg', 'p', 'sim_terminated');
fprintf('Results saved to results/sim_results.mat\n');

%% ===== Generate all plots and save =====
fprintf('Generating plots...\n');
plot_all_results(lg, p, results_dir);

%% ===== Generate GIFs =====
fprintf('Generating GIFs...\n');
generate_gifs(lg, p, results_dir);

fprintf('All done.\n');


%% ===== Local helper functions =====

function x_ref = compute_reference(t, p)
    y_hold = p.y_hold_end + (p.y_hold_start - p.y_hold_end) * exp(-t / p.y_hold_tau);
    x_ref = [0; y_hold; 0; 0; 0; 0];
end

function [ok, margin] = check_los(r_tb, cone_k, y_min)
    xT = r_tb(1); yT = r_tb(2); zT = r_tb(3);
    rad = sqrt(xT^2 + zT^2);
    if yT < y_min
        ok = false;
        margin = yT - y_min;
    elseif rad > cone_k * yT
        ok = false;
        margin = cone_k * yT - rad;
    else
        ok = true;
        margin = cone_k * yT - rad;
    end
end

function q = quat_from_rotm_main(R)
    tr = trace(R);
    if tr > 0
        s = 0.5 / sqrt(tr + 1);
        qw = 0.25 / s;
        qx = (R(3,2) - R(2,3)) * s;
        qy = (R(1,3) - R(3,1)) * s;
        qz = (R(2,1) - R(1,2)) * s;
    elseif R(1,1) > R(2,2) && R(1,1) > R(3,3)
        s = 2 * sqrt(1 + R(1,1) - R(2,2) - R(3,3));
        qw = (R(3,2) - R(2,3)) / s;
        qx = 0.25 * s;
        qy = (R(1,2) + R(2,1)) / s;
        qz = (R(1,3) + R(3,1)) / s;
    elseif R(2,2) > R(3,3)
        s = 2 * sqrt(1 + R(2,2) - R(1,1) - R(3,3));
        qw = (R(1,3) - R(3,1)) / s;
        qx = (R(1,2) + R(2,1)) / s;
        qy = 0.25 * s;
        qz = (R(2,3) + R(3,2)) / s;
    else
        s = 2 * sqrt(1 + R(3,3) - R(1,1) - R(2,2));
        qw = (R(2,1) - R(1,2)) / s;
        qx = (R(1,3) + R(3,1)) / s;
        qy = (R(2,3) + R(3,2)) / s;
        qz = 0.25 * s;
    end
    q = [qw; qx; qy; qz];
    q = q / norm(q);
end
