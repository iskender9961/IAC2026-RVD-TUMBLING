%% SIM_3D_TUMBLE_MPC  —  MPC approach to a 3D tumbling target
%
%  Self-contained simulation of a chaser approaching a torque-free
%  tumbling target that rotates about all three principal axes.
%  Extends the planar (z-axis only) MPC from main_sim.m to a full 3D
%  tumble scenario using Euler's equations for attitude propagation.
%
%  Outputs:
%    1) Static 3D trajectory plot (body frame) — saved as PNG
%    2) GIF animation in target body frame
%    3) GIF animation in LVLH frame showing tumbling target
%
%  Author : O. B. Iskender, NTU Singapore
%  Date   : Feb 2026
%  Ref    : Schaub & Junkins, "Analytical Mechanics of Space Systems", 4th ed.
%
%  Requirements: MATLAB R2020b+ (no external toolboxes)
%  -----------------------------------------------------------------------

clear; clc; close all;
fprintf('=== 3D Tumble MPC Simulation ===\n');

% Check for Optimization Toolbox (quadprog)
if ~exist('quadprog', 'file')
    error('sim_3d_tumble_mpc:noToolbox', ...
        'This script requires MATLAB Optimization Toolbox (quadprog).');
end

%% ====================================================================
%  SECTION 1 — PARAMETERS
%  ====================================================================

% --- Output directory ---
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', '3d_tumble');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% --- Orbit parameters (LEO, 500 km circular) ---
mu   = 3.986004418e14;          % Earth GM  [m^3/s^2]
Re   = 6378137.0;               % Earth equatorial radius [m]
J2   = 1.08263e-3;              % J2 zonal harmonic
alt  = 500e3;                   % Orbital altitude [m]
a_orb = Re + alt;               % Semi-major axis [m]
n_orb = sqrt(mu / a_orb^3);    % Mean motion [rad/s]

% --- Target inertia (asymmetric body, kg·m²) ---
I_body = diag([500, 800, 600]);

% --- Target initial angular velocity in body frame [rad/s] ---
omega0_body = [0.025; 0.015; 0.020];  % ~2 deg/s total, off all axes
omega_mag   = norm(omega0_body);      % ≈ 0.035 rad/s ≈ 2.0 deg/s

% --- Chaser parameters ---
a_max = 0.20;                   % Max thrust acceleration [m/s^2]

% --- Simulation timing ---
T_sim = 400;                    % Total simulation time [s]
dt    = 1.0;                    % Timestep [s]
N_steps = round(T_sim / dt);

% --- MPC parameters (matched to existing code) ---
Np    = 40;                     % Prediction horizon [steps]
nx    = 6;                      % State dimension [r; v] in body frame
nu    = 3;                      % Input dimension (body-frame accel)

% Cost matrices
Q     = diag([15, 30, 15, 1, 1, 1]);    % Stage cost (y-axis emphasised)
QN    = 30 * Q;                          % Terminal cost
Ru    = 1e-2 * eye(nu);                 % Input regularisation
Rdu   = 1e4 * eye(nu);                  % Input-rate penalty

% Reference trajectory (exponential approach along +y_B)
y_start = 100;                  % Initial hold distance [m]
y_end   = 5;                    % Final hold distance [m]
tau_ref = 100;                  % Exponential time constant [s]

% LOS cone
cone_half_angle = 30;                       % [deg]
cone_k = tan(deg2rad(cone_half_angle));     % slope parameter
y_min  = 1.0;                               % Minimum y_B [m]
n_faces = 8;                                % Polyhedral faces for cone

% Synchronisation range for 3D tumble
r_sync = 2 * a_max / omega_mag^2;

% --- Initial relative state in body frame ---
r0_tb = [10; 100; 5];           % Position [m] — offset in all axes
v0_tb = [0; 0; 0];              % Velocity [m/s]

% --- Target initial ECI state (circular orbit, x-direction at t=0) ---
r_tgt0_eci = [a_orb; 0; 0];
v_tgt0_eci = [0; sqrt(mu/a_orb); 0];

% --- Target initial quaternion (scalar-first, identity = body aligned with ECI) ---
q_tb = [1; 0; 0; 0];           % [qw; qx; qy; qz]
R_eci_tb = eye(3);              % Initial rotation matrix

% --- Derived chaser initial ECI state ---
r_rel_eci = R_eci_tb * r0_tb;
% Transport theorem: v_eci_rel = v_tb + omega x r  (in body), then rotate to ECI
v_rel_body = v0_tb + cross(omega0_body, r0_tb);
v_rel_eci  = R_eci_tb * v_rel_body;
r_chs0_eci = r_tgt0_eci + r_rel_eci;
v_chs0_eci = v_tgt0_eci + v_rel_eci;

% ODE options for truth propagation
ode_opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

fprintf('  omega_mag = %.2f deg/s,  r_sync = %.1f m\n', ...
    rad2deg(omega_mag), r_sync);
fprintf('  Initial range = %.1f m,  a_max = %.2f m/s^2\n', ...
    norm(r0_tb), a_max);

%% ====================================================================
%  SECTION 2 — HELPER FUNCTIONS (all defined at end of file)
%  ====================================================================
%  See bottom of file for:
%    quat_mult, quat_conj, quat_rotate_vec, rotm_from_quat_fn
%    euler_eom, propagate_euler_quat
%    cwh_stm_fn, eom_eci_fn, propagate_truth_fn
%    build_reference, check_los_fn
%    build_mpc_qp, solve_mpc_qp

%% ====================================================================
%  SECTION 3 — MAIN SIMULATION LOOP
%  ====================================================================

% --- Pre-allocate logging arrays ---
log_r_tb    = zeros(3, N_steps+1);   % Position in body frame
log_v_tb    = zeros(3, N_steps+1);   % Velocity in body frame
log_r_lvlh  = zeros(3, N_steps+1);   % Position in LVLH
log_u_tb    = zeros(3, N_steps);     % Applied control (body frame)
log_omega   = zeros(3, N_steps+1);   % Target angular velocity
log_q       = zeros(4, N_steps+1);   % Target quaternion
log_R_eci_tb = zeros(3, 3, N_steps+1); % Rotation matrices
log_los_margin = zeros(1, N_steps+1);% LOS cone margin
log_dist    = zeros(1, N_steps+1);   % Distance to target
log_status  = cell(1, N_steps);      % Solver status
log_cost    = zeros(1, N_steps);     % MPC stage cost
log_y_ref   = zeros(1, N_steps);     % Reference y_B at each step
log_dock_err_deg = zeros(1, N_steps+1); % Docking axis pointing error [deg]

% --- Initial conditions ---
x_tb = [r0_tb; v0_tb];              % 6x1 state in body frame
omega_body = omega0_body;            % Current angular velocity
r_tgt_eci  = r_tgt0_eci;
v_tgt_eci  = v_tgt0_eci;
r_chs_eci  = r_chs0_eci;
v_chs_eci  = v_chs0_eci;
u_prev     = zeros(nu, 1);          % Previous control input

% Log initial state
log_r_tb(:,1)      = r0_tb;
log_v_tb(:,1)      = v0_tb;
log_omega(:,1)     = omega_body;
log_q(:,1)         = q_tb;
log_R_eci_tb(:,:,1) = R_eci_tb;
log_dist(1)        = norm(r0_tb);

% Compute initial LVLH position
R_eci_lvlh = lvlh_from_rv_fn(r_tgt_eci, v_tgt_eci);
r_rel_eci_0 = r_chs_eci - r_tgt_eci;
log_r_lvlh(:,1) = R_eci_lvlh' * r_rel_eci_0;

% Compute initial LOS margin
[~, margin0] = check_los_fn(r0_tb, cone_k, y_min, n_faces);
log_los_margin(1) = margin0;

% Compute initial docking-axis error: angle between r_tb and +y_B axis
dock_axis = [0; 1; 0];  % Docking axis in body frame
cos_ang0 = dot(r0_tb, dock_axis) / (norm(r0_tb) + 1e-12);
log_dock_err_deg(1) = rad2deg(acos(max(-1, min(1, cos_ang0))));

fprintf('\n--- Starting simulation loop (%d steps) ---\n', N_steps);
sim_tic = tic;
n_los_viol = 0;

for k = 1:N_steps
    t_now = (k-1) * dt;

    % --- 3a. Build reference trajectory ---
    x_ref = build_reference(t_now, dt, Np, y_start, y_end, tau_ref);

    % --- 3b. Linearise dynamics (CWH + body-frame rotation) ---
    [Ad, Bd] = linearise_3d_body(x_tb, omega_body, n_orb, dt);

    % --- 3c. Build and solve MPC QP ---
    [u_opt, qp_status] = solve_mpc_qp(Ad, Bd, x_tb, x_ref, u_prev, ...
        Q, QN, Ru, Rdu, a_max, cone_k, y_min, n_faces, Np, nx, nu);
    log_status{k} = qp_status;

    % Log MPC stage cost: (x - x_ref)' Q (x - x_ref) + u' Ru u
    x_err = x_tb - x_ref(:,1);
    log_cost(k) = x_err' * Q * x_err + u_opt' * Ru * u_opt;
    log_y_ref(k) = x_ref(2,1);

    % Clamp control
    u_tb = max(-a_max, min(a_max, u_opt));

    % --- 3d. Propagate truth dynamics ---
    % Convert control to ECI
    a_ctrl_eci = R_eci_tb * u_tb;

    % Propagate chaser in ECI (two-body + J2 + control)
    x_chs_eci = [r_chs_eci; v_chs_eci];
    x_chs_next = propagate_truth_fn(x_chs_eci, dt, mu, Re, J2, a_ctrl_eci, ode_opts);
    r_chs_eci = x_chs_next(1:3);
    v_chs_eci = x_chs_next(4:6);

    % Propagate target in ECI (two-body + J2, no control)
    x_tgt_eci = [r_tgt_eci; v_tgt_eci];
    x_tgt_next = propagate_truth_fn(x_tgt_eci, dt, mu, Re, J2, [0;0;0], ode_opts);
    r_tgt_eci = x_tgt_next(1:3);
    v_tgt_eci = x_tgt_next(4:6);

    % --- 3e. Propagate target attitude (Euler equations + quaternion) ---
    [q_tb, R_eci_tb, omega_body] = propagate_euler_quat( ...
        q_tb, omega_body, I_body, dt);

    % --- 3f. Transform to body frame ---
    R_tb_eci = R_eci_tb';
    r_rel_eci = r_chs_eci - r_tgt_eci;
    v_rel_eci = v_chs_eci - v_tgt_eci;
    r_tb_new = R_tb_eci * r_rel_eci;
    v_tb_new = R_tb_eci * v_rel_eci - cross(omega_body, r_tb_new);
    x_tb = [r_tb_new; v_tb_new];

    % --- 3g. Log ---
    log_r_tb(:,k+1)      = r_tb_new;
    log_v_tb(:,k+1)      = v_tb_new;
    log_u_tb(:,k)         = u_tb;
    log_omega(:,k+1)      = omega_body;
    log_q(:,k+1)          = q_tb;
    log_R_eci_tb(:,:,k+1) = R_eci_tb;
    log_dist(k+1)         = norm(r_tb_new);

    % LVLH relative position
    R_eci_lvlh_k = lvlh_from_rv_fn(r_tgt_eci, v_tgt_eci);
    log_r_lvlh(:,k+1) = R_eci_lvlh_k' * r_rel_eci;

    % LOS check
    [los_ok, margin] = check_los_fn(r_tb_new, cone_k, y_min, n_faces);
    log_los_margin(k+1) = margin;
    if ~los_ok
        n_los_viol = n_los_viol + 1;
    end

    % Docking-axis pointing error
    cos_ang = dot(r_tb_new, dock_axis) / (norm(r_tb_new) + 1e-12);
    log_dock_err_deg(k+1) = rad2deg(acos(max(-1, min(1, cos_ang))));

    u_prev = u_tb;

    % Progress indicator
    if mod(k, 50) == 0
        fprintf('  Step %d/%d  (t=%.0fs)  dist=%.1fm  margin=%.2f  status=%s\n', ...
            k, N_steps, t_now+dt, norm(r_tb_new), margin, qp_status);
    end
end

sim_time = toc(sim_tic);
fprintf('\n--- Simulation complete (%.1f s wall-clock) ---\n', sim_time);
fprintf('  Final distance: %.2f m\n', log_dist(end));
fprintf('  LOS violations: %d / %d steps (%.1f%%)\n', ...
    n_los_viol, N_steps, 100*n_los_viol/N_steps);
fprintf('  Final omega: [%.4f, %.4f, %.4f] rad/s (|omega|=%.4f)\n', ...
    omega_body(1), omega_body(2), omega_body(3), norm(omega_body));

%% ====================================================================
%  SECTION 4 — OUTPUT 1: STATIC 3D PLOT (BODY FRAME)
%  ====================================================================
fprintf('\nGenerating static 3D plot...\n');

fig1 = figure('Position', [100 100 1000 800], 'Color', 'w');

% Time vector for colouring
t_vec = (0:N_steps) * dt;

% Plot trajectory coloured by time
scatter3(log_r_tb(1,:), log_r_tb(2,:), log_r_tb(3,:), ...
    15, t_vec, 'filled');
colormap(fig1, cool);
cb = colorbar; ylabel(cb, 'Time [s]');
hold on;

% Draw target as cuboid at origin
draw_cuboid([0,0,0], [4, 6, 3], [0.6 0.6 0.6], 0.5);

% Draw LOS cone (transparent)
draw_los_cone_3d(cone_k, y_min, 120, n_faces, [1 0.5 0], 0.15);

% Start and end markers
plot3(log_r_tb(1,1), log_r_tb(2,1), log_r_tb(3,1), ...
    'bs', 'MarkerSize', 14, 'MarkerFaceColor', 'b', 'DisplayName', 'Start');
plot3(log_r_tb(1,end), log_r_tb(2,end), log_r_tb(3,end), ...
    'r^', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'End');

% Docking axis
plot3([0 0], [y_min 120], [0 0], 'g--', 'LineWidth', 1.5, 'DisplayName', 'Docking axis (+y_B)');

xlabel('x_B [m]'); ylabel('y_B [m]'); zlabel('z_B [m]');
title(sprintf('3D Tumble: Chaser Trajectory in Target Body Frame\n|\\omega| = %.1f deg/s, a_{max} = %.2f m/s^2, r_{sync} = %.1f m', ...
    rad2deg(omega_mag), a_max, r_sync));
legend('Location', 'best');
grid on; axis equal;
view(35, 25);
set(gca, 'FontSize', 12);
hold off;

saveas(fig1, fullfile(out_dir, 'fig_3d_tumble_trajectory.png'));
fprintf('  Saved: fig_3d_tumble_trajectory.png\n');

%% ====================================================================
%  SECTION 5 — OUTPUT 2: GIF IN TARGET BODY FRAME
%  ====================================================================
fprintf('Generating body-frame GIF...\n');

gif_file_body = fullfile(out_dir, 'anim_3d_body_frame.gif');
fig2 = figure('Position', [100 100 900 700], 'Color', 'w', ...
    'Visible', 'off', 'Renderer', 'painters');
tmp_png = fullfile(tempdir, 'gif_frame_tmp.png');

% Frame sampling (every 2 seconds)
frame_skip = 2;
frame_idx = 1:frame_skip:N_steps+1;
n_frames = length(frame_idx);

for fi = 1:n_frames
    idx = frame_idx(fi);
    t_now = (idx-1) * dt;

    clf(fig2);

    % Trajectory up to current time (coloured)
    scatter3(log_r_tb(1,1:idx), log_r_tb(2,1:idx), log_r_tb(3,1:idx), ...
        10, t_vec(1:idx), 'filled');
    colormap(fig2, cool);
    hold on;

    % Current position (large marker)
    cp = log_r_tb(:, idx);
    plot3(cp(1), cp(2), cp(3), ...
        'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'LineWidth', 2);

    % Chaser body-frame axes at chaser position (in body frame = identity)
    ax_len_chs = max(5, log_dist(idx) * 0.12);
    quiver3(cp(1),cp(2),cp(3), ax_len_chs,0,0, 0, 'r','LineWidth',2.5,'MaxHeadSize',0.6);
    quiver3(cp(1),cp(2),cp(3), 0,ax_len_chs,0, 0, 'g','LineWidth',2.5,'MaxHeadSize',0.6);
    quiver3(cp(1),cp(2),cp(3), 0,0,ax_len_chs, 0, 'b','LineWidth',2.5,'MaxHeadSize',0.6);
    text(cp(1)+ax_len_chs*1.1, cp(2), cp(3), 'x_B', 'Color','r','FontSize',9,'FontWeight','bold');
    text(cp(1), cp(2)+ax_len_chs*1.1, cp(3), 'y_B', 'Color',[0 0.6 0],'FontSize',9,'FontWeight','bold');
    text(cp(1), cp(2), cp(3)+ax_len_chs*1.1, 'z_B', 'Color','b','FontSize',9,'FontWeight','bold');

    % Target cuboid
    draw_cuboid([0,0,0], [4, 6, 3], [0.6 0.6 0.6], 0.6);

    % LOS cone
    draw_los_cone_3d(cone_k, y_min, max(log_dist(idx)+20, 50), n_faces, [1 0.5 0], 0.1);

    % Docking axis
    plot3([0 0], [y_min max(log_dist(idx)+20, 50)], [0 0], 'g--', 'LineWidth', 1.5);

    xlabel('x_B [m]'); ylabel('y_B [m]'); zlabel('z_B [m]');
    title(sprintf('Body Frame  |  t = %.0f s  |  dist = %.1f m  |  \\omega = [%.2f, %.2f, %.2f] deg/s', ...
        t_now, log_dist(idx), rad2deg(log_omega(1,idx)), ...
        rad2deg(log_omega(2,idx)), rad2deg(log_omega(3,idx))));
    grid on; axis equal;

    % Dynamic axis limits
    rng_max = max(120, max(abs(log_r_tb(:,1:idx)), [], 'all') * 1.2);
    xlim([-rng_max rng_max]); ylim([-10 rng_max]); zlim([-rng_max rng_max]);
    view(35, 25);
    set(gca, 'FontSize', 11);
    hold off;

    % Capture frame
    drawnow;
    exportgraphics(fig2, tmp_png, 'Resolution', 100);
    im = imread(tmp_png);
    [imind, cm] = rgb2ind(im, 256);

    if fi == 1
        imwrite(imind, cm, gif_file_body, 'gif', 'Loopcount', inf, 'DelayTime', 0.08);
    else
        imwrite(imind, cm, gif_file_body, 'gif', 'WriteMode', 'append', 'DelayTime', 0.08);
    end
end
fprintf('  Saved: anim_3d_body_frame.gif (%d frames)\n', n_frames);
close(fig2);

%% ====================================================================
%  SECTION 6 — OUTPUT 3: GIF IN LVLH FRAME (TUMBLING TARGET)
%  ====================================================================
fprintf('Generating LVLH-frame GIF...\n');

gif_file_lvlh = fullfile(out_dir, 'anim_3d_lvlh_frame.gif');
fig3 = figure('Position', [100 100 900 700], 'Color', 'w', ...
    'Visible', 'off', 'Renderer', 'painters');

% Cuboid vertices in body frame (for rotation)
[cub_verts_body, cub_faces] = cuboid_mesh([4, 6, 3]);

% Cone mesh in body frame
[cone_verts_body, cone_faces_mesh] = cone_mesh(cone_k, y_min, 80, n_faces);

for fi = 1:n_frames
    idx = frame_idx(fi);
    t_now = (idx-1) * dt;

    clf(fig3);

    R_now = squeeze(log_R_eci_tb(:,:,idx));

    % --- Transform cuboid to LVLH ---
    R_eci_lvlh_now = lvlh_from_rv_fn(r_tgt0_eci, v_tgt0_eci);  % Approximate LVLH
    R_lvlh_body = R_eci_lvlh_now' * R_now;   % Body-to-LVLH rotation

    cub_verts_lvlh = (R_lvlh_body * cub_verts_body')';
    cone_verts_lvlh = (R_lvlh_body * cone_verts_body')';

    % Draw target (rotated cuboid)
    patch('Vertices', cub_verts_lvlh, 'Faces', cub_faces, ...
        'FaceColor', [0.5 0.5 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'k');
    hold on;

    % Draw rotated LOS cone
    patch('Vertices', cone_verts_lvlh, 'Faces', cone_faces_mesh, ...
        'FaceColor', [1 0.5 0], 'FaceAlpha', 0.12, 'EdgeColor', [1 0.6 0.2], ...
        'EdgeAlpha', 0.3);

    % Draw rotated docking axis
    dock_axis_body = [0 0; y_min 80; 0 0];
    dock_axis_lvlh = R_lvlh_body * dock_axis_body;
    plot3(dock_axis_lvlh(1,:), dock_axis_lvlh(2,:), dock_axis_lvlh(3,:), ...
        'g-', 'LineWidth', 2);

    % Rotated body axes (short, at origin)
    ax_len = 15;
    for ax_i = 1:3
        ax_dir = R_lvlh_body(:,ax_i) * ax_len;
        cols = {'r','g','b'};
        plot3([0 ax_dir(1)], [0 ax_dir(2)], [0 ax_dir(3)], ...
            cols{ax_i}, 'LineWidth', 2);
    end

    % Chaser trajectory in LVLH (up to now)
    scatter3(log_r_lvlh(1,1:idx), log_r_lvlh(2,1:idx), log_r_lvlh(3,1:idx), ...
        8, t_vec(1:idx), 'filled');
    colormap(fig3, cool);

    % Current chaser position
    cp_lvlh = log_r_lvlh(:, idx);
    plot3(cp_lvlh(1), cp_lvlh(2), cp_lvlh(3), ...
        'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'LineWidth', 2);

    % Chaser body-frame axes at chaser position (rotated to LVLH)
    ax_len_chs = max(5, log_dist(idx) * 0.12);
    for ax_i = 1:3
        ax_dir = R_lvlh_body(:, ax_i) * ax_len_chs;
        cols_rgb = [1 0 0; 0 0.7 0; 0 0 1];
        quiver3(cp_lvlh(1), cp_lvlh(2), cp_lvlh(3), ...
            ax_dir(1), ax_dir(2), ax_dir(3), 0, ...
            'Color', cols_rgb(ax_i,:), 'LineWidth', 2.5, 'MaxHeadSize', 0.6);
    end
    % Labels
    lbl = {'x_B','y_B','z_B'}; lbl_col = {[1 0 0],[0 0.6 0],[0 0 1]};
    for ax_i = 1:3
        tip = cp_lvlh + R_lvlh_body(:,ax_i) * ax_len_chs * 1.15;
        text(tip(1),tip(2),tip(3), lbl{ax_i}, 'Color',lbl_col{ax_i}, ...
            'FontSize',8,'FontWeight','bold');
    end

    xlabel('x_{LVLH} [m]'); ylabel('y_{LVLH} [m]'); zlabel('z_{LVLH} [m]');
    title(sprintf('LVLH Frame  |  t = %.0f s  |  dist = %.1f m  |  Target tumbling in 3D', ...
        t_now, log_dist(idx)));
    grid on; axis equal;
    rng_max = max(120, max(abs(log_r_lvlh(:,1:idx)), [], 'all') * 1.3);
    xlim([-rng_max rng_max]); ylim([-rng_max rng_max]); zlim([-rng_max rng_max]);
    view(40, 20);
    set(gca, 'FontSize', 11);
    hold off;

    drawnow;
    exportgraphics(fig3, tmp_png, 'Resolution', 100);
    im = imread(tmp_png);
    [imind, cm] = rgb2ind(im, 256);

    if fi == 1
        imwrite(imind, cm, gif_file_lvlh, 'gif', 'Loopcount', inf, 'DelayTime', 0.08);
    else
        imwrite(imind, cm, gif_file_lvlh, 'gif', 'WriteMode', 'append', 'DelayTime', 0.08);
    end
end
fprintf('  Saved: anim_3d_lvlh_frame.gif (%d frames)\n', n_frames);
close(fig3);

%% ====================================================================
%  SECTION 7 — SUMMARY METRICS PLOT
%  ====================================================================
fprintf('Generating summary plots...\n');

fig4 = figure('Position', [100 100 1200 800], 'Color', 'w');

subplot(2,3,1);
plot(t_vec, log_dist, 'b-', 'LineWidth', 1.5);
hold on;
yline(r_sync, 'r--', 'LineWidth', 1.5, 'Label', 'r_{sync}');
xlabel('Time [s]'); ylabel('Distance [m]');
title('Range to Target'); grid on;

subplot(2,3,2);
plot(t_vec, log_los_margin, 'k-', 'LineWidth', 1.5);
hold on;
yline(0, 'r--', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('Margin');
title('LOS Cone Margin'); grid on;
if any(log_los_margin < 0)
    fill_idx = find(log_los_margin < 0);
    stem(t_vec(fill_idx), log_los_margin(fill_idx), 'r.', 'MarkerSize', 4);
end

subplot(2,3,3);
plot(t_vec, rad2deg(log_omega(1,:)), 'r-', 'LineWidth', 1.2); hold on;
plot(t_vec, rad2deg(log_omega(2,:)), 'g-', 'LineWidth', 1.2);
plot(t_vec, rad2deg(log_omega(3,:)), 'b-', 'LineWidth', 1.2);
plot(t_vec, rad2deg(vecnorm(log_omega)), 'k--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('\omega [deg/s]');
title('Target Angular Velocity'); grid on;
legend('\omega_x','\omega_y','\omega_z','|\omega|','Location','best');

subplot(2,3,4);
plot(t_vec, log_r_tb(1,:), 'r-', 'LineWidth', 1.2); hold on;
plot(t_vec, log_r_tb(2,:), 'g-', 'LineWidth', 1.2);
plot(t_vec, log_r_tb(3,:), 'b-', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Position [m]');
title('Body-Frame Position'); grid on;
legend('x_B','y_B','z_B','Location','best');

subplot(2,3,5);
t_u = (0:N_steps-1) * dt;
plot(t_u, log_u_tb(1,:), 'r-', 'LineWidth', 1.2); hold on;
plot(t_u, log_u_tb(2,:), 'g-', 'LineWidth', 1.2);
plot(t_u, log_u_tb(3,:), 'b-', 'LineWidth', 1.2);
yline(a_max, 'k--'); yline(-a_max, 'k--');
xlabel('Time [s]'); ylabel('Accel [m/s^2]');
title('Control Input (Body Frame)'); grid on;
legend('u_x','u_y','u_z','Location','best');

subplot(2,3,6);
dv_cumul = cumsum(vecnorm(log_u_tb) * dt);
plot(t_u, dv_cumul, 'k-', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('\Delta v [m/s]');
title(sprintf('Cumulative \\Delta v (total = %.1f m/s)', dv_cumul(end)));
grid on;

sgtitle(sprintf('3D Tumble MPC Summary  |  |\\omega_0| = %.1f deg/s  |  LOS violations: %d/%d', ...
    rad2deg(omega_mag), n_los_viol, N_steps), 'FontSize', 14);

saveas(fig4, fullfile(out_dir, 'fig_3d_tumble_summary.png'));
fprintf('  Saved: fig_3d_tumble_summary.png\n');

%% ====================================================================
%  SECTION 8 — ADDITIONAL STATIC PLOTS (PDF + PNG)
%  ====================================================================

% --- 8a. MPC Cost History ---
fprintf('Generating cost history plot...\n');
fig5 = figure('Position', [100 100 800 400], 'Color', 'w');
semilogy(t_u, log_cost, 'b-', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Stage Cost J_k');
title('MPC Stage Cost History');
grid on; set(gca, 'FontSize', 12);
saveas(fig5, fullfile(out_dir, 'fig_cost_history.png'));
exportgraphics(fig5, fullfile(out_dir, 'fig_cost_history.pdf'), 'ContentType', 'vector');
fprintf('  Saved: fig_cost_history.png / .pdf\n');

% --- 8b. Angular Velocity Components ---
fprintf('Generating angular velocity plot...\n');
fig6 = figure('Position', [100 100 800 450], 'Color', 'w');
plot(t_vec, rad2deg(log_omega(1,:)), 'r-', 'LineWidth', 1.5); hold on;
plot(t_vec, rad2deg(log_omega(2,:)), 'g-', 'LineWidth', 1.5);
plot(t_vec, rad2deg(log_omega(3,:)), 'b-', 'LineWidth', 1.5);
plot(t_vec, rad2deg(vecnorm(log_omega)), 'k--', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Angular Velocity [deg/s]');
title('Target Angular Velocity (Euler Torque-Free Motion)');
legend('\omega_x', '\omega_y', '\omega_z', '|\omega|', 'Location', 'best');
grid on; set(gca, 'FontSize', 12);
saveas(fig6, fullfile(out_dir, 'fig_omega_components.png'));
exportgraphics(fig6, fullfile(out_dir, 'fig_omega_components.pdf'), 'ContentType', 'vector');
fprintf('  Saved: fig_omega_components.png / .pdf\n');

% --- 8c. Relative Position in Body Frame (with reference) ---
fprintf('Generating body-frame position plot...\n');
fig7 = figure('Position', [100 100 900 500], 'Color', 'w');

subplot(3,1,1);
plot(t_vec, log_r_tb(1,:), 'r-', 'LineWidth', 1.5);
hold on; yline(0, 'k--', 'LineWidth', 0.8);
ylabel('x_B [m]'); title('Relative Position in Target Body Frame');
grid on; set(gca, 'FontSize', 11, 'XTickLabel', []);

subplot(3,1,2);
plot(t_vec, log_r_tb(2,:), 'g-', 'LineWidth', 1.5); hold on;
% Plot reference
y_ref_full = y_end + (y_start - y_end) * exp(-t_vec / tau_ref);
plot(t_vec, y_ref_full, 'k--', 'LineWidth', 1.2);
ylabel('y_B [m]');
legend('y_B (actual)', 'y_{ref} (exponential)', 'Location', 'best');
grid on; set(gca, 'FontSize', 11, 'XTickLabel', []);

subplot(3,1,3);
plot(t_vec, log_r_tb(3,:), 'b-', 'LineWidth', 1.5);
hold on; yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Time [s]'); ylabel('z_B [m]');
grid on; set(gca, 'FontSize', 11);

saveas(fig7, fullfile(out_dir, 'fig_position_body_frame.png'));
exportgraphics(fig7, fullfile(out_dir, 'fig_position_body_frame.pdf'), 'ContentType', 'vector');
fprintf('  Saved: fig_position_body_frame.png / .pdf\n');

% --- 8d. Docking Axis Tracking Error ---
fprintf('Generating docking-axis error plot...\n');
fig8 = figure('Position', [100 100 800 450], 'Color', 'w');

yyaxis left;
plot(t_vec, log_dock_err_deg, 'b-', 'LineWidth', 1.5);
ylabel('Docking Axis Error [deg]');
ylim([0, max(log_dock_err_deg)*1.1 + 1]);

yyaxis right;
plot(t_vec, log_dist, 'r-', 'LineWidth', 1.2);
ylabel('Range [m]');

hold on;
yline(r_sync, 'r--', 'LineWidth', 1, 'Label', 'r_{sync}');
xlabel('Time [s]');
title(sprintf('Docking Axis Pointing Error  (cone half-angle = %d°)', cone_half_angle));
legend('Pointing error', 'Range', 'Location', 'best');
grid on; set(gca, 'FontSize', 12);

saveas(fig8, fullfile(out_dir, 'fig_docking_axis_error.png'));
exportgraphics(fig8, fullfile(out_dir, 'fig_docking_axis_error.pdf'), 'ContentType', 'vector');
fprintf('  Saved: fig_docking_axis_error.png / .pdf\n');

fprintf('\n=== All outputs saved to: %s ===\n', out_dir);

%% ====================================================================
%  SECTION 9 — LOCAL FUNCTIONS
%  ====================================================================

% ---- Quaternion operations (scalar-first: q = [qw; qx; qy; qz]) ----

function qr = quat_mult_fn(q1, q2)
    % Hamilton product of two scalar-first quaternions
    w1=q1(1); x1=q1(2); y1=q1(3); z1=q1(4);
    w2=q2(1); x2=q2(2); y2=q2(3); z2=q2(4);
    qr = [w1*w2 - x1*x2 - y1*y2 - z1*z2;
          w1*x2 + x1*w2 + y1*z2 - z1*y2;
          w1*y2 - x1*z2 + y1*w2 + z1*x2;
          w1*z2 + x1*y2 - y1*x2 + z1*w2];
end

function qc = quat_conj_fn(q)
    qc = [q(1); -q(2); -q(3); -q(4)];
end

function R = rotm_from_quat_fn(q)
    % Rotation matrix from scalar-first quaternion
    % R rotates vectors from body to reference: v_ref = R * v_body
    q = q / norm(q);
    w=q(1); x=q(2); y=q(3); z=q(4);
    R = [1-2*(y^2+z^2),  2*(x*y - w*z),  2*(x*z + w*y);
         2*(x*y + w*z),  1-2*(x^2+z^2),  2*(y*z - w*x);
         2*(x*z - w*y),  2*(y*z + w*x),  1-2*(x^2+y^2)];
end

% ---- Euler's equations for torque-free rigid body ----

function [q_next, R_next, omega_next] = propagate_euler_quat(q, omega, I_body, dt)
    % Propagate attitude using Euler's equations + quaternion kinematics
    % Uses RK4 for both omega and quaternion.
    %
    % Euler's equations (torque-free):
    %   I1 * dw1/dt = (I2 - I3) * w2 * w3
    %   I2 * dw2/dt = (I3 - I1) * w3 * w1
    %   I3 * dw3/dt = (I1 - I2) * w1 * w2
    %
    % Quaternion kinematics:
    %   dq/dt = 0.5 * q ⊗ [0; omega]

    I1 = I_body(1,1); I2 = I_body(2,2); I3 = I_body(3,3);

    % Combined state: y = [omega(3); q(4)]
    function dy = eom(y)
        w = y(1:3);
        qv = y(4:7);
        % Euler's equations
        dw = [(I2-I3)/I1 * w(2)*w(3);
              (I3-I1)/I2 * w(3)*w(1);
              (I1-I2)/I3 * w(1)*w(2)];
        % Quaternion kinematics: dq/dt = 0.5 * q ⊗ [0; w]
        omega_quat = [0; w];
        dq = 0.5 * quat_mult_fn(qv, omega_quat);
        dy = [dw; dq];
    end

    % RK4 integration
    y0 = [omega; q];
    k1 = eom(y0);
    k2 = eom(y0 + 0.5*dt*k1);
    k3 = eom(y0 + 0.5*dt*k2);
    k4 = eom(y0 + dt*k3);
    y_next = y0 + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

    omega_next = y_next(1:3);
    q_next = y_next(4:7);
    q_next = q_next / norm(q_next);  % Renormalise

    R_next = rotm_from_quat_fn(q_next);
end

% ---- CWH State Transition Matrix ----

function [Ad, Bd] = cwh_stm_fn(dt, n)
    % Clohessy-Wiltshire-Hill discrete dynamics
    % State: [x; y; z; vx; vy; vz] in LVLH
    nt = n * dt;
    snt = sin(nt); cnt = cos(nt);

    Ad = [4-3*cnt,    0, 0,  snt/n,      2*(1-cnt)/n,    0;
          6*(snt-nt), 1, 0, -2*(1-cnt)/n, (4*snt-3*nt)/n, 0;
          0,          0, cnt, 0,          0,              snt/n;
          3*n*snt,    0, 0,  cnt,         2*snt,          0;
         -6*n*(1-cnt),0, 0, -2*snt,       4*cnt-3,        0;
          0,          0, -n*snt, 0,       0,              cnt];

    % Zero-order hold input matrix
    Bd = zeros(6, 3);
    Bd(1,1) = (1-cnt)/n^2;
    Bd(1,2) = 2*(nt - snt)/n^2;
    Bd(2,1) = -2*(nt - snt)/n^2;
    Bd(2,2) = (4*(1-cnt) - 1.5*n^2*dt^2)/n^2;
    Bd(3,3) = (1-cnt)/n^2;
    Bd(4,1) = snt/n;
    Bd(4,2) = 2*(1-cnt)/n;
    Bd(5,1) = -2*(1-cnt)/n;
    Bd(5,2) = (4*snt - 3*nt)/n;
    Bd(6,3) = snt/n;
end

% ---- Linearisation in rotating body frame ----

function [Ad, Bd] = linearise_3d_body(x_tb, omega, n_orb, dt)
    % Approximate linearised discrete dynamics in the rotating body frame.
    % Uses CWH for orbital coupling + rotation transport terms.
    %
    % The body-frame state equation incorporates:
    %   dr_B/dt = v_B + omega × r_B   (transport theorem, already removed)
    %   dv_B/dt = a_CWH - 2*omega × v_B - omega × (omega × r_B)
    %             - domega/dt × r_B + u_B
    %
    % We linearise numerically via finite differences through the
    % nonlinear body-frame propagation.

    eps_x = 1e-6;
    eps_u = 1e-6;

    % Nominal propagation
    x0 = propagate_body_step(x_tb, [0;0;0], omega, n_orb, dt);

    % State Jacobian (6 columns)
    Ad = zeros(6, 6);
    for j = 1:6
        x_pert = x_tb;
        x_pert(j) = x_pert(j) + eps_x;
        x_fwd = propagate_body_step(x_pert, [0;0;0], omega, n_orb, dt);
        Ad(:,j) = (x_fwd - x0) / eps_x;
    end

    % Input Jacobian (3 columns)
    Bd = zeros(6, 3);
    for j = 1:3
        u_pert = zeros(3,1);
        u_pert(j) = eps_u;
        x_fwd = propagate_body_step(x_tb, u_pert, omega, n_orb, dt);
        Bd(:,j) = (x_fwd - x0) / eps_u;
    end
end

function x_next = propagate_body_step(x, u, omega, n_orb, dt)
    % Propagate one step in body frame using RK4.
    % Includes CWH orbital coupling + Coriolis + centrifugal terms.
    %
    % State: x = [r_B; v_B] where v_B is velocity in rotating body frame.
    %
    % Equations of motion in rotating body frame:
    %   dr_B/dt = v_B
    %   dv_B/dt = a_grav(r_B) - 2*omega × v_B - omega × (omega × r_B) + u
    %
    % where a_grav approximates CWH coupling (gravity gradient in LVLH,
    % rotated to body frame — for linearisation purposes we keep the
    % CWH structure directly in body frame).

    function dx = eom_body(x_in)
        r = x_in(1:3);
        v = x_in(4:6);

        % CWH-like gravity gradient in body frame
        % (approximate: treat body frame as quasi-LVLH for gravity gradient)
        n2 = n_orb^2;
        a_grav = [3*n2*r(1); 0; -n2*r(3)];  % CWH tidal + z-axis restoring

        % Coriolis and centrifugal
        a_coriolis    = -2 * cross(omega, v);
        a_centrifugal = -cross(omega, cross(omega, r));

        dv = a_grav + a_coriolis + a_centrifugal + u;
        dx = [v; dv];
    end

    % RK4
    k1 = eom_body(x);
    k2 = eom_body(x + 0.5*dt*k1);
    k3 = eom_body(x + 0.5*dt*k2);
    k4 = eom_body(x + dt*k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% ---- ECI dynamics (truth model) ----

function x_next = propagate_truth_fn(x, dt, mu, Re, J2, a_ctrl_eci, ode_opts)
    % Propagate ECI state with two-body + J2 + control
    function dxdt = eom(~, xv)
        r = xv(1:3); v = xv(4:6);
        r_mag = norm(r);
        % Two-body
        a_2b = -mu / r_mag^3 * r;
        % J2
        z2 = r(3)^2 / r_mag^2;
        fac = -1.5 * J2 * mu * Re^2 / r_mag^5;
        a_j2 = fac * [(1 - 5*z2)*r(1);
                       (1 - 5*z2)*r(2);
                       (3 - 5*z2)*r(3)];
        dxdt = [v; a_2b + a_j2 + a_ctrl_eci];
    end
    [~, X] = ode113(@eom, [0 dt], x, ode_opts);
    x_next = X(end,:)';
end

% ---- LVLH frame rotation matrix ----

function R = lvlh_from_rv_fn(r_eci, v_eci)
    % R_eci_lvlh: columns are LVLH unit vectors in ECI
    x_hat = r_eci / norm(r_eci);         % Radial
    h = cross(r_eci, v_eci);
    z_hat = h / norm(h);                 % Orbit normal
    y_hat = cross(z_hat, x_hat);         % Along-track
    R = [x_hat, y_hat, z_hat];
end

% ---- Reference trajectory builder ----

function x_ref = build_reference(t_now, dt, Np, y_start, y_end, tau_ref)
    % Build Np+1 reference states (exponential approach along +y_B)
    x_ref = zeros(6, Np+1);
    for j = 0:Np
        t_j = t_now + j * dt;
        y_ref_j = y_end + (y_start - y_end) * exp(-t_j / tau_ref);
        x_ref(:, j+1) = [0; y_ref_j; 0; 0; 0; 0];
    end
end

% ---- LOS cone check ----

function [los_ok, margin] = check_los_fn(r_tb, cone_k, y_min, n_faces)
    % Check if position is inside the polyhedral LOS cone
    x = r_tb(1); y = r_tb(2); z = r_tb(3);

    if y < y_min
        los_ok = false;
        margin = y - y_min;
        return;
    end

    % Check each face of the polyhedral cone
    margin = Inf;
    for f = 1:n_faces
        theta = 2*pi*(f-1)/n_faces;
        val = cos(theta)*x + sin(theta)*z - cone_k * y;
        margin = min(margin, -val);  % margin > 0 means inside
    end

    los_ok = (margin >= 0);
end

% ---- MPC QP Solver (dense formulation, no OSQP dependency) ----

function [u_opt, status] = solve_mpc_qp(Ad, Bd, x0, x_ref, u_prev, ...
        Q, QN, Ru, Rdu, a_max, cone_k, y_min, n_faces, Np, nx, nu)
    % Build and solve the MPC QP using MATLAB's quadprog.
    % Decision variables: z = [x_0; x_1; ...; x_Np; u_0; u_1; ...; u_{Np-1}]

    n_x_total = (Np+1) * nx;
    n_u_total = Np * nu;
    n_vars    = n_x_total + n_u_total;

    % --- Cost: 0.5 * z' * H * z + f' * z ---
    H = sparse(n_vars, n_vars);
    f = zeros(n_vars, 1);

    for j = 0:Np
        ix = j*nx + (1:nx);
        if j < Np
            H(ix, ix) = Q;
            f(ix) = -Q * x_ref(:, j+1);
        else
            H(ix, ix) = QN;
            f(ix) = -QN * x_ref(:, j+1);
        end
    end

    for j = 0:Np-1
        iu = n_x_total + j*nu + (1:nu);
        H(iu, iu) = H(iu, iu) + Ru;

        % Delta-u penalty
        if j == 0
            % (u_0 - u_prev)
            H(iu, iu) = H(iu, iu) + Rdu;
            f(iu) = f(iu) - Rdu * u_prev;
        else
            iu_prev = n_x_total + (j-1)*nu + (1:nu);
            H(iu, iu) = H(iu, iu) + Rdu;
            H(iu_prev, iu_prev) = H(iu_prev, iu_prev) + Rdu;
            H(iu, iu_prev) = H(iu, iu_prev) - Rdu;
            H(iu_prev, iu) = H(iu_prev, iu) - Rdu;
        end
    end

    H = (H + H') / 2;  % Ensure symmetry

    % --- Equality constraints: dynamics + initial condition ---
    n_eq = (Np+1) * nx;
    Aeq = sparse(n_eq, n_vars);
    beq = zeros(n_eq, 1);

    % x_0 = x0
    Aeq(1:nx, 1:nx) = eye(nx);
    beq(1:nx) = x0;

    % x_{k+1} = Ad * x_k + Bd * u_k
    for j = 0:Np-1
        row = (j+1)*nx + (1:nx);
        ix_k   = j*nx + (1:nx);
        ix_kp1 = (j+1)*nx + (1:nx);
        iu_k   = n_x_total + j*nu + (1:nu);

        Aeq(row, ix_kp1) = eye(nx);
        Aeq(row, ix_k)   = -Ad;
        Aeq(row, iu_k)   = -Bd;
    end

    % --- Inequality constraints ---
    % 1) Input box bounds: -a_max <= u_k <= a_max
    n_input = 2 * Np * nu;
    A_input = sparse(n_input, n_vars);
    b_input = zeros(n_input, 1);
    row_cnt = 0;
    for j = 0:Np-1
        iu = n_x_total + j*nu + (1:nu);
        % u_k <= a_max
        A_input(row_cnt + (1:nu), iu) = eye(nu);
        b_input(row_cnt + (1:nu)) = a_max;
        row_cnt = row_cnt + nu;
        % -u_k <= a_max
        A_input(row_cnt + (1:nu), iu) = -eye(nu);
        b_input(row_cnt + (1:nu)) = a_max;
        row_cnt = row_cnt + nu;
    end

    % 2) LOS cone constraints: cos(th)*x + sin(th)*z - cone_k*y <= 0, y >= y_min
    n_los = (n_faces + 1) * (Np + 1);
    A_los = sparse(n_los, n_vars);
    b_los = zeros(n_los, 1);
    row_cnt = 0;
    for j = 0:Np
        ix = j*nx + (1:nx);
        for ffi = 1:n_faces
            theta = 2*pi*(ffi-1)/n_faces;
            row_cnt = row_cnt + 1;
            % cos(th)*x + sin(th)*z - cone_k*y <= 0
            A_los(row_cnt, ix(1)) = cos(theta);      % x_B
            A_los(row_cnt, ix(2)) = -cone_k;         % y_B
            A_los(row_cnt, ix(3)) = sin(theta);       % z_B
            b_los(row_cnt) = 0;
        end
        % -y <= -y_min  =>  y >= y_min
        row_cnt = row_cnt + 1;
        A_los(row_cnt, ix(2)) = -1;
        b_los(row_cnt) = -y_min;
    end

    A_ineq = [A_input; A_los];
    b_ineq = [b_input; b_los];

    % --- Solve with quadprog ---
    options = optimoptions('quadprog', 'Display', 'off', ...
        'MaxIterations', 2000, ...
        'OptimalityTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6);

    [z_opt, ~, exitflag] = quadprog(H, f, A_ineq, b_ineq, Aeq, beq, ...
        [], [], [], options);

    if exitflag > 0
        status = 'solved';
        u_opt = z_opt(n_x_total + (1:nu));
    else
        status = sprintf('failed(%d)', exitflag);
        u_opt = zeros(nu, 1);
    end
end

% ---- Visualisation helpers ----

function draw_cuboid(center, dims, color, alpha)
    % Draw a cuboid at 'center' with half-dimensions 'dims'
    cx = center(1); cy = center(2); cz = center(3);
    dx = dims(1)/2; dy = dims(2)/2; dz = dims(3)/2;

    verts = [cx-dx cy-dy cz-dz;
             cx+dx cy-dy cz-dz;
             cx+dx cy+dy cz-dz;
             cx-dx cy+dy cz-dz;
             cx-dx cy-dy cz+dz;
             cx+dx cy-dy cz+dz;
             cx+dx cy+dy cz+dz;
             cx-dx cy+dy cz+dz];

    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

    patch('Vertices', verts, 'Faces', faces, ...
        'FaceColor', color, 'FaceAlpha', alpha, 'EdgeColor', 'k');
end

function [verts, faces] = cuboid_mesh(dims)
    % Returns cuboid vertices and faces for later transformation
    dx = dims(1)/2; dy = dims(2)/2; dz = dims(3)/2;
    verts = [-dx -dy -dz;
              dx -dy -dz;
              dx  dy -dz;
             -dx  dy -dz;
             -dx -dy  dz;
              dx -dy  dz;
              dx  dy  dz;
             -dx  dy  dz];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
end

function draw_los_cone_3d(cone_k, y_min, y_max, n_faces, color, alpha)
    % Draw 3D polyhedral LOS cone along +y_B axis
    n_pts = max(n_faces, 24);  % Use more points for smooth visualisation
    theta = linspace(0, 2*pi, n_pts+1);

    % Cone opening: radius = cone_k * y at each y
    y_vals = linspace(y_min, y_max, 20);

    for yi = 1:length(y_vals)-1
        y1 = y_vals(yi);   r1 = cone_k * y1;
        y2 = y_vals(yi+1); r2 = cone_k * y2;

        for ti = 1:n_pts
            x_ring = [r1*cos(theta(ti)), r2*cos(theta(ti+1)), ...
                      r2*cos(theta(ti+1)), r1*cos(theta(ti))];
            y_ring = [y1, y2, y2, y1];
            z_ring = [r1*sin(theta(ti)), r2*sin(theta(ti+1)), ...
                      r2*sin(theta(ti+1)), r1*sin(theta(ti))];
            fill3(x_ring, y_ring, z_ring, color, ...
                'FaceAlpha', alpha, 'EdgeColor', 'none');
        end
    end
end

function [verts, faces] = cone_mesh(cone_k, y_min, y_max, ~)
    % Generate cone mesh vertices and triangular faces for transformation.
    % All faces are triangles (3 columns) for uniform patch rendering.
    n_ring = 10;  % Number of rings along y
    n_circ = 16;  % Points around circumference

    y_vals = linspace(y_min, y_max, n_ring);
    theta = linspace(0, 2*pi, n_circ+1); theta = theta(1:end-1);

    % Apex vertex (index 1)
    verts = [0, 0, 0];

    % Ring vertices
    for ri = 1:n_ring
        r = cone_k * y_vals(ri);
        for ti = 1:n_circ
            verts = [verts; r*cos(theta(ti)), y_vals(ri), r*sin(theta(ti))]; %#ok<AGROW>
        end
    end

    % Faces: triangles from apex to first ring
    faces = zeros(n_circ, 3);
    for ti = 1:n_circ
        ti_next = mod(ti, n_circ) + 1;
        faces(ti,:) = [1, 1+ti, 1+ti_next];
    end

    % Triangulated quads between rings (2 triangles per quad)
    for ri = 1:n_ring-1
        base1 = 1 + (ri-1)*n_circ;
        base2 = 1 + ri*n_circ;
        for ti = 1:n_circ
            ti_next = mod(ti, n_circ) + 1;
            v1 = base1 + ti;
            v2 = base1 + ti_next;
            v3 = base2 + ti_next;
            v4 = base2 + ti;
            faces = [faces; v1, v2, v4; v2, v3, v4]; %#ok<AGROW>
        end
    end
end
