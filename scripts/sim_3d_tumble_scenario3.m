%% SIM_3D_TUMBLE_SCENARIO3  —  MPC approach to a 3-axis tumbling target
%
%  Scenario 3: Equal angular velocity on all three principal axes (~3 deg/s
%  total), so the tumbling motion does not favour any single axis.
%  This produces a chaotic-looking polhode and visually dramatic 3D tumbling
%  since no axis is "dominant".
%
%  Relative to Scenario 2:
%    - omega0 ≈ [0.030, 0.030, 0.030] rad/s  (|omega| ≈ 3.0 deg/s)
%    - Approach from 200 m  (r0 = [40, 200, -50])
%    - LOS cone rendered as prominent wireframe in GIFs (always visible)
%    - omega_target x,y,z displayed in GIF annotations
%
%  Outputs (saved to results/3d_tumble_scenario3/):
%    Static PNG+PDF:  trajectory, approach profile, summary, cost,
%                     omega, position, control, delta-v, LOS margin,
%                     docking error, polhode
%    Animated GIF:    body frame, LVLH frame
%
%  Author : O. B. Iskender, NTU Singapore
%  Date   : Mar 2026
%  -----------------------------------------------------------------------

clear; clc; close all;
fprintf('=== 3D Tumble MPC — Scenario 3 (Equal 3-Axis Tumble, 3 deg/s) ===\n');

if ~exist('quadprog', 'file')
    error('sim_3d_tumble_scenario3:noToolbox', ...
        'Requires MATLAB Optimization Toolbox (quadprog).');
end

%% ====================================================================
%  SECTION 1 — PARAMETERS
%  ====================================================================

out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', '3d_tumble_scenario3');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% --- Orbit (LEO 500 km circular) ---
mu    = 3.986004418e14;
Re    = 6378137.0;
J2    = 1.08263e-3;
alt   = 500e3;
a_orb = Re + alt;
n_orb = sqrt(mu / a_orb^3);

% --- Target inertia (asymmetric: I1 < I3 < I2) ---
I_body = diag([500, 800, 600]);

% --- Target initial angular velocity [rad/s] ---
%  Equal on all three axes for truly 3-axis tumble.
%  |omega| = sqrt(3)*0.030 = 0.05196 rad/s ≈ 2.98 deg/s ≈ 3 deg/s
omega0_body = [0.030; 0.030; 0.030];
omega_mag   = norm(omega0_body);

% --- Chaser ---
a_max = 0.30;                   % Max thrust [m/s²] (increased for 3-axis tumble)

% --- Timing ---
T_sim   = 600;                  % Longer sim time for 200m approach with aggressive tumble
dt      = 1.0;
N_steps = round(T_sim / dt);

% --- MPC ---
Np = 40;   nx = 6;   nu = 3;
Q     = diag([10, 20, 10, 0.5, 0.5, 0.5]);
QN    = 30 * Q;
Ru    = 5e-3 * eye(nu);         % Lighter input penalty → more thrust authority
Rdu   = 5e3  * eye(nu);         % Softer rate penalty

% --- Reference trajectory ---
%  Linear ramp that reaches docking by T_dock, then holds.
%  Unlike exponential, the chaser has a hard deadline to arrive.
%  y_ref(t) = y_end + (y_start-y_end)*max(0, 1 - t/T_dock)
y_start = 150;   y_end = 5;   T_dock = 250;   % Arrive by 250 s

% --- LOS cone ---
cone_half_angle = 30;
cone_k  = tan(deg2rad(cone_half_angle));
y_min   = 1.0;
n_faces = 8;

% --- Synchronisation range ---
r_sync = 2 * a_max / omega_mag^2;

% --- Initial relative state (body frame) ---
r0_tb = [-20; 150; 20];
v0_tb = [0; 0; 0];

% --- Target initial ECI state ---
r_tgt0_eci = [a_orb; 0; 0];
v_tgt0_eci = [0; sqrt(mu/a_orb); 0];

% --- Initial attitude (identity: body = ECI at t=0) ---
q_tb     = [1; 0; 0; 0];
R_eci_tb = eye(3);

% --- Chaser initial ECI state ---
r_rel_eci  = R_eci_tb * r0_tb;
v_rel_body = v0_tb + cross(omega0_body, r0_tb);
v_rel_eci  = R_eci_tb * v_rel_body;
r_chs0_eci = r_tgt0_eci + r_rel_eci;
v_chs0_eci = v_tgt0_eci + v_rel_eci;

% --- Chaser box dimensions [m] (for visualisation) ---
chs_box = [3, 4, 2.5];   % [x_chs, y_chs, z_chs] — y is along-track
tgt_box = [4, 6, 3];

ode_opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

fprintf('  |omega_0| = %.2f deg/s (equal 3-axis tumble)\n', rad2deg(omega_mag));
fprintf('  omega_0   = [%.3f, %.3f, %.3f] rad/s = [%.1f, %.1f, %.1f] deg/s\n', ...
    omega0_body(1), omega0_body(2), omega0_body(3), ...
    rad2deg(omega0_body(1)), rad2deg(omega0_body(2)), rad2deg(omega0_body(3)));
fprintf('  r_sync    = %.1f m\n', r_sync);
fprintf('  r_0       = [%.0f, %.0f, %.0f] m  (|r_0| = %.1f m)\n', ...
    r0_tb(1), r0_tb(2), r0_tb(3), norm(r0_tb));

%% ====================================================================
%  SECTION 2 — MAIN SIMULATION LOOP
%  ====================================================================

% Pre-allocate
log_r_tb      = zeros(3, N_steps+1);
log_v_tb      = zeros(3, N_steps+1);
log_r_lvlh    = zeros(3, N_steps+1);
log_u_tb      = zeros(3, N_steps);
log_omega     = zeros(3, N_steps+1);
log_q         = zeros(4, N_steps+1);
log_R_eci_tb  = zeros(3,3, N_steps+1);
log_los_margin= zeros(1, N_steps+1);
log_dist      = zeros(1, N_steps+1);
log_cost      = zeros(1, N_steps);
log_y_ref     = zeros(1, N_steps+1);
log_dock_err  = zeros(1, N_steps+1);    % docking-axis pointing error [deg]
log_omega_dir = zeros(1, N_steps+1);    % omega direction change [deg]
log_status    = cell(1, N_steps);

% Initial state
x_tb       = [r0_tb; v0_tb];
omega_body = omega0_body;
r_tgt_eci  = r_tgt0_eci;   v_tgt_eci = v_tgt0_eci;
r_chs_eci  = r_chs0_eci;   v_chs_eci = v_chs0_eci;
u_prev     = zeros(nu,1);
dock_axis  = [0;1;0];

% Log t=0
log_r_tb(:,1)      = r0_tb;
log_v_tb(:,1)      = v0_tb;
log_omega(:,1)     = omega_body;
log_q(:,1)         = q_tb;
log_R_eci_tb(:,:,1)= R_eci_tb;
log_dist(1)        = norm(r0_tb);
log_y_ref(1)       = y_start;

R_eci_lvlh = lvlh_from_rv_fn(r_tgt_eci, v_tgt_eci);
log_r_lvlh(:,1) = R_eci_lvlh' * (r_chs_eci - r_tgt_eci);

[~, m0] = check_los_fn(r0_tb, cone_k, y_min, n_faces);
log_los_margin(1) = m0;

ca0 = dot(r0_tb, dock_axis) / (norm(r0_tb)+1e-12);
log_dock_err(1) = rad2deg(acos(max(-1,min(1,ca0))));
log_omega_dir(1) = 0;

fprintf('\n--- Simulation loop (%d steps) ---\n', N_steps);
tic_sim = tic;
n_los_viol = 0;

for k = 1:N_steps
    t_now = (k-1)*dt;

    % Reference
    x_ref = build_reference(t_now, dt, Np, y_start, y_end, T_dock);

    % Linearise
    [Ad, Bd] = linearise_3d_body(x_tb, omega_body, n_orb, dt);

    % Solve MPC
    [u_opt, qp_status] = solve_mpc_qp(Ad, Bd, x_tb, x_ref, u_prev, ...
        Q, QN, Ru, Rdu, a_max, cone_k, y_min, n_faces, Np, nx, nu);
    log_status{k} = qp_status;

    % Cost
    x_err = x_tb - x_ref(:,1);
    log_cost(k) = x_err' * Q * x_err + u_opt' * Ru * u_opt;

    % Clamp
    u_tb_k = max(-a_max, min(a_max, u_opt));

    % Propagate chaser (truth: 2-body + J2 + control)
    a_eci = R_eci_tb * u_tb_k;
    x_chs_next = propagate_truth_fn([r_chs_eci;v_chs_eci], dt, mu, Re, J2, a_eci, ode_opts);
    r_chs_eci = x_chs_next(1:3);  v_chs_eci = x_chs_next(4:6);

    % Propagate target (truth: 2-body + J2, no control)
    x_tgt_next = propagate_truth_fn([r_tgt_eci;v_tgt_eci], dt, mu, Re, J2, [0;0;0], ode_opts);
    r_tgt_eci = x_tgt_next(1:3);  v_tgt_eci = x_tgt_next(4:6);

    % Attitude (Euler + quaternion)
    [q_tb, R_eci_tb, omega_body] = propagate_euler_quat(q_tb, omega_body, I_body, dt);

    % Body-frame state
    R_tb_eci  = R_eci_tb';
    r_rel_eci = r_chs_eci - r_tgt_eci;
    v_rel_eci = v_chs_eci - v_tgt_eci;
    r_tb_new  = R_tb_eci * r_rel_eci;
    v_tb_new  = R_tb_eci * v_rel_eci - cross(omega_body, r_tb_new);
    x_tb      = [r_tb_new; v_tb_new];

    % Log
    log_r_tb(:,k+1)      = r_tb_new;
    log_v_tb(:,k+1)      = v_tb_new;
    log_u_tb(:,k)         = u_tb_k;
    log_omega(:,k+1)      = omega_body;
    log_q(:,k+1)          = q_tb;
    log_R_eci_tb(:,:,k+1) = R_eci_tb;
    log_dist(k+1)         = norm(r_tb_new);
    log_y_ref(k+1)        = y_end + (y_start-y_end)*max(0, 1-(t_now+dt)/T_dock);

    R_eci_lvlh_k = lvlh_from_rv_fn(r_tgt_eci, v_tgt_eci);
    log_r_lvlh(:,k+1) = R_eci_lvlh_k' * r_rel_eci;

    [lok, mrg] = check_los_fn(r_tb_new, cone_k, y_min, n_faces);
    log_los_margin(k+1) = mrg;
    if ~lok, n_los_viol = n_los_viol + 1; end

    % Docking-axis error
    ca = dot(r_tb_new, dock_axis)/(norm(r_tb_new)+1e-12);
    log_dock_err(k+1) = rad2deg(acos(max(-1,min(1,ca))));

    % Omega direction change from initial
    cd = dot(omega_body, omega0_body)/(norm(omega_body)*norm(omega0_body)+1e-15);
    log_omega_dir(k+1) = rad2deg(acos(max(-1,min(1,cd))));

    u_prev = u_tb_k;

    if mod(k,50)==0
        fprintf('  Step %3d/%d  t=%3.0fs  dist=%6.1fm  margin=%6.2f  omega_dir=%5.1f°  w=[%.2f,%.2f,%.2f]°/s  %s\n',...
            k,N_steps,t_now+dt,norm(r_tb_new),mrg,log_omega_dir(k+1),...
            rad2deg(omega_body(1)),rad2deg(omega_body(2)),rad2deg(omega_body(3)),qp_status);
    end
end

wall_time = toc(tic_sim);
t_vec = (0:N_steps)*dt;
t_u   = (0:N_steps-1)*dt;

fprintf('\n--- Done (%.1f s) ---\n', wall_time);
fprintf('  Final dist : %.2f m\n', log_dist(end));
fprintf('  LOS viol   : %d / %d (%.1f%%)\n', n_los_viol, N_steps, 100*n_los_viol/N_steps);
fprintf('  omega dir  : %.1f° precession from initial\n', log_omega_dir(end));
fprintf('  Delta-v    : %.1f m/s\n', sum(vecnorm(log_u_tb)*dt));

%% ====================================================================
%  SECTION 3 — STATIC 3D TRAJECTORY (BODY FRAME)
%  ====================================================================
fprintf('\nGenerating plots...\n');

fig1 = figure('Position',[50 50 1000 800],'Color','w');
scatter3(log_r_tb(1,:), log_r_tb(2,:), log_r_tb(3,:), 12, t_vec, 'filled');
colormap(fig1, cool); cb=colorbar; ylabel(cb,'Time [s]');
hold on;
draw_cuboid_fn([0,0,0], tgt_box, [0.5 0.5 0.8], 0.6);
draw_los_cone_3d(cone_k, y_min, 220, n_faces, [1 0.5 0], 0.08);
plot3([0 0],[y_min 220],[0 0],'g--','LineWidth',1.5);

% Start / end markers
plot3(r0_tb(1),r0_tb(2),r0_tb(3),'bs','MarkerSize',14,'MarkerFaceColor','b');
plot3(log_r_tb(1,end),log_r_tb(2,end),log_r_tb(3,end),...
    'r^','MarkerSize',14,'MarkerFaceColor','r');

% Draw chaser box at start & end
for pidx = [1, N_steps+1]
    rp = log_r_tb(:,pidx);
    R_chs = chaser_attitude(rp);
    draw_oriented_box(rp, chs_box, R_chs, [0.2 0.8 0.9], 0.5);
    draw_axes_at(rp, R_chs, 6);
end

xlabel('x_B [m]'); ylabel('y_B [m]'); zlabel('z_B [m]');
title(sprintf(['Scenario 3: 3D Trajectory in Target Body Frame\n' ...
    '|\\omega_0| = %.1f deg/s (equal 3-axis), ' ...
    'r_0 = [%.0f,%.0f,%.0f] m'], ...
    rad2deg(omega_mag), r0_tb(1), r0_tb(2), r0_tb(3)));
grid on; axis equal; view(35,25); set(gca,'FontSize',12);
hold off;

saveas(fig1, fullfile(out_dir,'fig_3d_trajectory_body.png'));
exportgraphics(fig1, fullfile(out_dir,'fig_3d_trajectory_body.pdf'),'ContentType','vector');
fprintf('  fig_3d_trajectory_body\n');

%% ====================================================================
%  SECTION 4 — STATIC 3D TRAJECTORY (LVLH)
%  ====================================================================
fig1b = figure('Position',[50 50 1000 800],'Color','w');
scatter3(log_r_lvlh(1,:), log_r_lvlh(2,:), log_r_lvlh(3,:), 12, t_vec, 'filled');
colormap(fig1b, cool); cb=colorbar; ylabel(cb,'Time [s]');
hold on;
plot3(log_r_lvlh(1,1),log_r_lvlh(2,1),log_r_lvlh(3,1),'bs','MarkerSize',14,'MarkerFaceColor','b');
plot3(log_r_lvlh(1,end),log_r_lvlh(2,end),log_r_lvlh(3,end),'r^','MarkerSize',14,'MarkerFaceColor','r');
plot3(0,0,0,'kp','MarkerSize',18,'MarkerFaceColor','k');
xlabel('x_{LVLH} [m]'); ylabel('y_{LVLH} [m]'); zlabel('z_{LVLH} [m]');
title('Chaser Trajectory in LVLH Frame (Scenario 3)');
grid on; axis equal; view(40,20); set(gca,'FontSize',12);
hold off;
saveas(fig1b, fullfile(out_dir,'fig_3d_trajectory_lvlh.png'));
exportgraphics(fig1b, fullfile(out_dir,'fig_3d_trajectory_lvlh.pdf'),'ContentType','vector');
fprintf('  fig_3d_trajectory_lvlh\n');

%% ====================================================================
%  SECTION 5 — GIF: BODY FRAME (with persistent cone + omega text)
%  ====================================================================
fprintf('  Generating body-frame GIF...\n');
gif_body = fullfile(out_dir,'anim_body_frame.gif');
fig2 = figure('Position',[50 50 950 750],'Color','w','Visible','off','Renderer','painters');
tmp_png = fullfile(tempdir,'s3_gif_tmp.png');

frame_skip = 2;
frame_idx  = 1:frame_skip:N_steps+1;
n_frames   = length(frame_idx);

for fi = 1:n_frames
    idx = frame_idx(fi);  t_now = (idx-1)*dt;
    clf(fig2);

    % Trail
    plot3(log_r_tb(1,1:idx), log_r_tb(2,1:idx), log_r_tb(3,1:idx),...
        '-','Color',[0.3 0.6 1],'LineWidth',1);
    hold on;

    % Target cuboid at origin
    draw_cuboid_fn([0,0,0], tgt_box, [0.5 0.5 0.8], 0.6);
    % Target body axes at origin
    draw_axes_at([0;0;0], eye(3), 10);

    % LOS cone — ALWAYS draw with thick, prominent wireframe
    cone_range = max(log_dist(idx)+30, 80);
    draw_los_wireframe_thick(cone_k, y_min, cone_range, n_faces);

    % Docking axis
    plot3([0 0],[y_min cone_range],[0 0],'g--','LineWidth',1.5);

    % Chaser box with pointing attitude
    rp = log_r_tb(:,idx);
    R_chs = chaser_attitude(rp);
    draw_oriented_box(rp, chs_box, R_chs, [0.2 0.8 0.9], 0.7);
    % Axes at tip (front face)
    tip = rp + R_chs(:,2) * (chs_box(2)/2);
    draw_axes_at(tip, R_chs, max(4, log_dist(idx)*0.08));

    xlabel('x_B [m]'); ylabel('y_B [m]'); zlabel('z_B [m]');
    % Title with omega_target components
    title(sprintf(['Body Frame | t=%3.0fs | dist=%.1fm\n' ...
        '\\omega_{tgt} = [%.2f, %.2f, %.2f] deg/s  |\\omega|=%.2f deg/s'],...
        t_now, log_dist(idx), ...
        rad2deg(log_omega(1,idx)), rad2deg(log_omega(2,idx)), rad2deg(log_omega(3,idx)),...
        rad2deg(norm(log_omega(:,idx)))));
    grid on; axis equal;
    rm = max(220, max(abs(log_r_tb(:,1:idx)),[],'all')*1.2);
    xlim([-rm rm]); ylim([-10 rm]); zlim([-rm rm]);
    view(35,25); set(gca,'FontSize',10); hold off;

    drawnow;
    exportgraphics(fig2, tmp_png, 'Resolution', 100);
    im = imread(tmp_png);
    [imind,cm] = rgb2ind(im,256);
    if fi==1
        imwrite(imind,cm,gif_body,'gif','Loopcount',inf,'DelayTime',0.08);
    else
        imwrite(imind,cm,gif_body,'gif','WriteMode','append','DelayTime',0.08);
    end
end
close(fig2);
fprintf('    anim_body_frame.gif (%d frames)\n', n_frames);

%% ====================================================================
%  SECTION 6 — GIF: LVLH FRAME (TUMBLING TARGET + cone + omega text)
%  ====================================================================
fprintf('  Generating LVLH-frame GIF...\n');
gif_lvlh = fullfile(out_dir,'anim_lvlh_frame.gif');
fig3 = figure('Position',[50 50 950 750],'Color','w','Visible','off','Renderer','painters');

[cub_v_body, cub_f] = cuboid_mesh(tgt_box);
R_eci_lvlh_ref = lvlh_from_rv_fn(r_tgt0_eci, v_tgt0_eci);

for fi = 1:n_frames
    idx = frame_idx(fi);  t_now = (idx-1)*dt;
    clf(fig3);

    R_now = squeeze(log_R_eci_tb(:,:,idx));
    R_lvlh_body = R_eci_lvlh_ref' * R_now;

    % Rotated target cuboid
    cv = (R_lvlh_body * cub_v_body')';
    patch('Vertices',cv,'Faces',cub_f,...
        'FaceColor',[0.5 0.5 0.8],'FaceAlpha',0.7,'EdgeColor','k');
    hold on;

    % Target body axes at origin (rotated to LVLH)
    draw_axes_at([0;0;0], R_lvlh_body, 15);

    % Rotated docking axis
    da_len = max(log_dist(idx)+30, 80);
    da = R_lvlh_body * [0 0; y_min da_len; 0 0];
    plot3(da(1,:),da(2,:),da(3,:),'g-','LineWidth',2);

    % Rotated LOS cone in LVLH
    draw_los_wireframe_rotated(cone_k, y_min, da_len, n_faces, R_lvlh_body);

    % Chaser trail
    plot3(log_r_lvlh(1,1:idx),log_r_lvlh(2,1:idx),log_r_lvlh(3,1:idx),...
        '-','Color',[0.3 0.6 1],'LineWidth',1);

    % Chaser box (attitude in LVLH)
    cp = log_r_lvlh(:,idx);
    rp_body = log_r_tb(:,idx);
    R_chs_body = chaser_attitude(rp_body);
    R_chs_lvlh = R_lvlh_body * R_chs_body;
    draw_oriented_box(cp, chs_box, R_chs_lvlh, [0.2 0.8 0.9], 0.7);
    tip_lvlh = cp + R_chs_lvlh(:,2)*(chs_box(2)/2);
    draw_axes_at(tip_lvlh, R_chs_lvlh, max(4, log_dist(idx)*0.08));

    xlabel('x_{LVLH} [m]'); ylabel('y_{LVLH} [m]'); zlabel('z_{LVLH} [m]');
    title(sprintf(['LVLH Frame | t=%3.0fs | dist=%.1fm | axis shift=%.1f°\n'...
        '\\omega_{tgt} = [%.2f, %.2f, %.2f] deg/s  |\\omega|=%.2f deg/s'],...
        t_now, log_dist(idx), log_omega_dir(idx),...
        rad2deg(log_omega(1,idx)), rad2deg(log_omega(2,idx)), rad2deg(log_omega(3,idx)),...
        rad2deg(norm(log_omega(:,idx)))));
    grid on; axis equal;
    rm = max(220, max(abs(log_r_lvlh(:,1:idx)),[],'all')*1.3);
    xlim([-rm rm]); ylim([-rm rm]); zlim([-rm rm]);
    view(40,20); set(gca,'FontSize',10); hold off;

    drawnow;
    exportgraphics(fig3, tmp_png, 'Resolution', 100);
    im = imread(tmp_png);
    [imind,cm] = rgb2ind(im,256);
    if fi==1
        imwrite(imind,cm,gif_lvlh,'gif','Loopcount',inf,'DelayTime',0.08);
    else
        imwrite(imind,cm,gif_lvlh,'gif','WriteMode','append','DelayTime',0.08);
    end
end
close(fig3);
fprintf('    anim_lvlh_frame.gif (%d frames)\n', n_frames);

%% ====================================================================
%  SECTION 7 — SUMMARY 6-PANEL PLOT
%  ====================================================================
fig4 = figure('Position',[50 50 1200 800],'Color','w');
dv = cumsum(vecnorm(log_u_tb)*dt);

subplot(2,3,1);
plot(t_vec, log_dist, 'b-','LineWidth',1.5); hold on;
plot(t_vec, log_y_ref, 'k--','LineWidth',1);
yline(r_sync,'r--','LineWidth',1,'Label','r_{sync}');
xlabel('Time [s]'); ylabel('Distance [m]');
title('Range & Reference'); grid on;
legend('|r|','y_{ref}','Location','best');

subplot(2,3,2);
plot(t_vec, log_los_margin, 'k-','LineWidth',1.5); hold on;
yline(0,'r--','LineWidth',1);
if any(log_los_margin<0)
    vi=find(log_los_margin<0);
    stem(t_vec(vi),log_los_margin(vi),'r.','MarkerSize',4);
end
xlabel('Time [s]'); ylabel('Margin');
title('LOS Cone Margin'); grid on;

subplot(2,3,3);
yyaxis left;
plot(t_vec, rad2deg(log_omega(1,:)),'r-','LineWidth',1.2); hold on;
plot(t_vec, rad2deg(log_omega(2,:)),'g-','LineWidth',1.2);
plot(t_vec, rad2deg(log_omega(3,:)),'b-','LineWidth',1.2);
ylabel('\omega [deg/s]');
yyaxis right;
plot(t_vec, log_omega_dir, 'm-','LineWidth',1.5);
ylabel('Axis shift [deg]');
xlabel('Time [s]'); title('\omega Components & Axis Precession'); grid on;

subplot(2,3,4);
plot(t_vec, log_r_tb(1,:),'r-','LineWidth',1.2); hold on;
plot(t_vec, log_r_tb(2,:),'g-','LineWidth',1.2);
plot(t_vec, log_r_tb(3,:),'b-','LineWidth',1.2);
xlabel('Time [s]'); ylabel('[m]');
title('Body-Frame Position'); grid on;
legend('x_B','y_B','z_B','Location','best');

subplot(2,3,5);
plot(t_u, log_u_tb(1,:),'r-','LineWidth',1.2); hold on;
plot(t_u, log_u_tb(2,:),'g-','LineWidth',1.2);
plot(t_u, log_u_tb(3,:),'b-','LineWidth',1.2);
yline(a_max,'k--'); yline(-a_max,'k--');
xlabel('Time [s]'); ylabel('[m/s^2]');
title('Control Acceleration'); grid on;
legend('u_x','u_y','u_z','Location','best');

subplot(2,3,6);
plot(t_u, dv, 'k-','LineWidth',1.5);
xlabel('Time [s]'); ylabel('\Delta v [m/s]');
title(sprintf('Cumulative \\Delta v = %.1f m/s', dv(end))); grid on;

sgtitle(sprintf(['Scenario 3 Summary | |\\omega_0|=%.1f°/s (equal 3-axis) | ' ...
    'LOS viol: %d/%d | \\Deltav=%.1f m/s'], ...
    rad2deg(omega_mag), n_los_viol, N_steps, dv(end)), 'FontSize',13);

saveas(fig4, fullfile(out_dir,'fig_summary.png'));
exportgraphics(fig4, fullfile(out_dir,'fig_summary.pdf'),'ContentType','vector');
fprintf('  fig_summary\n');

%% ====================================================================
%  SECTION 8 — INDIVIDUAL JOURNAL-QUALITY PLOTS (PNG + PDF)
%  ====================================================================

% --- 8a. MPC Cost History ---
fig5 = figure('Position',[100 100 800 400],'Color','w');
semilogy(t_u, log_cost, 'b-','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Stage Cost J_k');
title('MPC Stage Cost History (Scenario 3 — Equal 3-Axis Tumble)');
grid on; set(gca,'FontSize',12);
saveas(fig5, fullfile(out_dir,'fig_cost_history.png'));
exportgraphics(fig5, fullfile(out_dir,'fig_cost_history.pdf'),'ContentType','vector');
fprintf('  fig_cost_history\n');

% --- 8b. Angular Velocity with Axis-Shift Error ---
fig6 = figure('Position',[100 100 900 450],'Color','w');
yyaxis left;
plot(t_vec, rad2deg(log_omega(1,:)),'r-','LineWidth',1.5); hold on;
plot(t_vec, rad2deg(log_omega(2,:)),'g-','LineWidth',1.5);
plot(t_vec, rad2deg(log_omega(3,:)),'b-','LineWidth',1.5);
plot(t_vec, rad2deg(vecnorm(log_omega)),'k--','LineWidth',2);
ylabel('Angular Velocity [deg/s]');
yyaxis right;
plot(t_vec, log_omega_dir, 'm-','LineWidth',2);
ylabel('Spin-Axis Precession from \omega_0 [deg]');
xlabel('Time [s]');
title('Target Angular Velocity & Spin-Axis Precession (Equal 3-Axis Tumble)');
legend('\omega_x','\omega_y','\omega_z','|\omega|','Axis shift','Location','best');
grid on; set(gca,'FontSize',12);
saveas(fig6, fullfile(out_dir,'fig_omega_components.png'));
exportgraphics(fig6, fullfile(out_dir,'fig_omega_components.pdf'),'ContentType','vector');
fprintf('  fig_omega_components\n');

% --- 8c. Position Components (Body Frame) ---
fig7 = figure('Position',[100 100 900 550],'Color','w');
subplot(3,1,1);
plot(t_vec, log_r_tb(1,:),'r-','LineWidth',1.5);
hold on; yline(0,'k--','LineWidth',0.8);
ylabel('x_B [m]'); title('Relative Position in Target Body Frame (Scenario 3)');
grid on; set(gca,'FontSize',11,'XTickLabel',[]);

subplot(3,1,2);
plot(t_vec, log_r_tb(2,:),'g-','LineWidth',1.5); hold on;
plot(t_vec, log_y_ref,'k--','LineWidth',1.2);
ylabel('y_B [m]');
legend('y_B (actual)','y_{ref} (linear ramp)','Location','best');
grid on; set(gca,'FontSize',11,'XTickLabel',[]);

subplot(3,1,3);
plot(t_vec, log_r_tb(3,:),'b-','LineWidth',1.5);
hold on; yline(0,'k--','LineWidth',0.8);
xlabel('Time [s]'); ylabel('z_B [m]');
grid on; set(gca,'FontSize',11);
saveas(fig7, fullfile(out_dir,'fig_position_body.png'));
exportgraphics(fig7, fullfile(out_dir,'fig_position_body.pdf'),'ContentType','vector');
fprintf('  fig_position_body\n');

% --- 8d. Control Acceleration Components ---
fig8 = figure('Position',[100 100 900 550],'Color','w');
subplot(3,1,1);
plot(t_u, log_u_tb(1,:),'r-','LineWidth',1.5);
hold on; yline(a_max,'k--'); yline(-a_max,'k--');
ylabel('u_{x,B} [m/s^2]');
title('Control Acceleration in Target Body Frame (Scenario 3)');
grid on; set(gca,'FontSize',11,'XTickLabel',[]);

subplot(3,1,2);
plot(t_u, log_u_tb(2,:),'g-','LineWidth',1.5);
hold on; yline(a_max,'k--'); yline(-a_max,'k--');
ylabel('u_{y,B} [m/s^2]');
grid on; set(gca,'FontSize',11,'XTickLabel',[]);

subplot(3,1,3);
plot(t_u, log_u_tb(3,:),'b-','LineWidth',1.5);
hold on; yline(a_max,'k--'); yline(-a_max,'k--');
xlabel('Time [s]'); ylabel('u_{z,B} [m/s^2]');
grid on; set(gca,'FontSize',11);
saveas(fig8, fullfile(out_dir,'fig_control_input.png'));
exportgraphics(fig8, fullfile(out_dir,'fig_control_input.pdf'),'ContentType','vector');
fprintf('  fig_control_input\n');

% --- 8e. Cumulative Delta-V ---
fig9 = figure('Position',[100 100 800 400],'Color','w');
dv_x = cumsum(abs(log_u_tb(1,:))*dt);
dv_y = cumsum(abs(log_u_tb(2,:))*dt);
dv_z = cumsum(abs(log_u_tb(3,:))*dt);
area(t_u, [dv_x; dv_y; dv_z]', 'EdgeColor','none');
colororder([1 0 0; 0 0.7 0; 0 0 1]);
hold on;
plot(t_u, dv, 'k-','LineWidth',2);
xlabel('Time [s]'); ylabel('\Delta v [m/s]');
title(sprintf('Cumulative \\Delta v Budget (total = %.1f m/s)', dv(end)));
legend('\Delta v_x','\Delta v_y','\Delta v_z','Total','Location','best');
grid on; set(gca,'FontSize',12);
saveas(fig9, fullfile(out_dir,'fig_deltav.png'));
exportgraphics(fig9, fullfile(out_dir,'fig_deltav.pdf'),'ContentType','vector');
fprintf('  fig_deltav\n');

% --- 8f. LOS Cone Margin ---
fig10 = figure('Position',[100 100 800 400],'Color','w');
plot(t_vec, log_los_margin, 'k-','LineWidth',1.5); hold on;
yline(0,'r-','LineWidth',1.5);
if any(log_los_margin<0)
    vi=find(log_los_margin<0);
    stem(t_vec(vi), log_los_margin(vi), 'r.','MarkerSize',6);
    fill([t_vec(vi) fliplr(t_vec(vi))], ...
         [log_los_margin(vi) zeros(size(vi))], 'r','FaceAlpha',0.2,'EdgeColor','none');
end
xlabel('Time [s]'); ylabel('Min. Cone Margin');
title(sprintf('LOS Corridor Margin (violations: %d/%d steps)', n_los_viol, N_steps));
grid on; set(gca,'FontSize',12);
saveas(fig10, fullfile(out_dir,'fig_los_margin.png'));
exportgraphics(fig10, fullfile(out_dir,'fig_los_margin.pdf'),'ContentType','vector');
fprintf('  fig_los_margin\n');

% --- 8g. Docking Axis Tracking Error ---
fig11 = figure('Position',[100 100 800 450],'Color','w');
yyaxis left;
plot(t_vec, log_dock_err, 'b-','LineWidth',1.5);
ylabel('Pointing Error [deg]');
ylim([0 max(log_dock_err)*1.1+1]);
yyaxis right;
plot(t_vec, log_dist, 'r-','LineWidth',1.2);
ylabel('Range [m]');
xlabel('Time [s]');
title(sprintf('Docking-Axis Pointing Error (cone half-angle = %d°)', cone_half_angle));
legend('Pointing error','Range','Location','best');
grid on; set(gca,'FontSize',12);
saveas(fig11, fullfile(out_dir,'fig_docking_error.png'));
exportgraphics(fig11, fullfile(out_dir,'fig_docking_error.pdf'),'ContentType','vector');
fprintf('  fig_docking_error\n');

% --- 8h. Approach Profile (journal figure style) ---
fig12 = figure('Position',[100 100 900 650],'Color','w');
subplot(3,1,1);
plot(t_vec, log_r_tb(2,:),'g-','LineWidth',1.5); hold on;
plot(t_vec, log_y_ref,'k--','LineWidth',1.2);
ylabel('y_B [m]'); title('Approach Profile — Scenario 3 (Equal 3-Axis Tumble)');
legend('Docking axis range','y_{ref} (linear ramp)','Location','best');
grid on; set(gca,'FontSize',11,'XTickLabel',[]);

subplot(3,1,2);
lat_dev = sqrt(log_r_tb(1,:).^2 + log_r_tb(3,:).^2);
plot(t_vec, lat_dev, 'b-','LineWidth',1.5);
ylabel('Lateral dev. [m]');
grid on; set(gca,'FontSize',11,'XTickLabel',[]);

subplot(3,1,3);
plot(t_vec, log_los_margin,'k-','LineWidth',1.5); hold on;
yline(0,'r--','LineWidth',1);
xlabel('Time [s]'); ylabel('LOS margin');
grid on; set(gca,'FontSize',11);
saveas(fig12, fullfile(out_dir,'fig_approach_profile.png'));
exportgraphics(fig12, fullfile(out_dir,'fig_approach_profile.pdf'),'ContentType','vector');
fprintf('  fig_approach_profile\n');

% --- 8i. Polhode plot (omega trajectory on Poinsot ellipsoid) ---
fig13 = figure('Position',[100 100 700 600],'Color','w');
plot3(rad2deg(log_omega(1,:)), rad2deg(log_omega(2,:)), rad2deg(log_omega(3,:)),...
    'b-','LineWidth',1.5);
hold on;
plot3(rad2deg(log_omega(1,1)), rad2deg(log_omega(2,1)), rad2deg(log_omega(3,1)),...
    'gs','MarkerSize',12,'MarkerFaceColor','g');
plot3(rad2deg(log_omega(1,end)), rad2deg(log_omega(2,end)), rad2deg(log_omega(3,end)),...
    'r^','MarkerSize',12,'MarkerFaceColor','r');
xlabel('\omega_x [deg/s]'); ylabel('\omega_y [deg/s]'); zlabel('\omega_z [deg/s]');
title('Polhode: Angular Velocity Trajectory (Equal 3-Axis Tumble)');
legend('Polhode','Start','End','Location','best');
grid on; axis equal; view(35,25); set(gca,'FontSize',12);
hold off;
saveas(fig13, fullfile(out_dir,'fig_polhode.png'));
exportgraphics(fig13, fullfile(out_dir,'fig_polhode.pdf'),'ContentType','vector');
fprintf('  fig_polhode\n');

% --- 8j. Omega-target 3-axis time history (dedicated) ---
fig14 = figure('Position',[100 100 900 550],'Color','w');
subplot(3,1,1);
plot(t_vec, rad2deg(log_omega(1,:)),'r-','LineWidth',1.5);
hold on; yline(rad2deg(omega0_body(1)),'k--','LineWidth',0.8);
ylabel('\omega_x [deg/s]'); title('Target Angular Velocity Components (Scenario 3)');
grid on; set(gca,'FontSize',11,'XTickLabel',[]);

subplot(3,1,2);
plot(t_vec, rad2deg(log_omega(2,:)),'g-','LineWidth',1.5);
hold on; yline(rad2deg(omega0_body(2)),'k--','LineWidth',0.8);
ylabel('\omega_y [deg/s]');
grid on; set(gca,'FontSize',11,'XTickLabel',[]);

subplot(3,1,3);
plot(t_vec, rad2deg(log_omega(3,:)),'b-','LineWidth',1.5);
hold on; yline(rad2deg(omega0_body(3)),'k--','LineWidth',0.8);
xlabel('Time [s]'); ylabel('\omega_z [deg/s]');
grid on; set(gca,'FontSize',11);
saveas(fig14, fullfile(out_dir,'fig_omega_target_3axis.png'));
exportgraphics(fig14, fullfile(out_dir,'fig_omega_target_3axis.pdf'),'ContentType','vector');
fprintf('  fig_omega_target_3axis\n');

% --- Summary table ---
fid = fopen(fullfile(out_dir,'summary.txt'),'w');
fprintf(fid,'Scenario 3 — Equal 3-Axis Tumble at 3 deg/s\n');
fprintf(fid,'=============================================\n\n');
fprintf(fid,'Target MOI   : I = diag([%d, %d, %d]) kg.m^2\n', I_body(1,1), I_body(2,2), I_body(3,3));
fprintf(fid,'omega_0      : [%.4f, %.4f, %.4f] rad/s  = [%.2f, %.2f, %.2f] deg/s\n',...
    omega0_body(1),omega0_body(2),omega0_body(3),...
    rad2deg(omega0_body(1)),rad2deg(omega0_body(2)),rad2deg(omega0_body(3)));
fprintf(fid,'|omega_0|    : %.4f rad/s = %.2f deg/s\n', omega_mag, rad2deg(omega_mag));
fprintf(fid,'a_max        : %.2f m/s^2\n', a_max);
fprintf(fid,'r_sync       : %.1f m\n', r_sync);
fprintf(fid,'r_0          : [%.0f, %.0f, %.0f] m  (|r_0| = %.1f m)\n',...
    r0_tb(1), r0_tb(2), r0_tb(3), norm(r0_tb));
fprintf(fid,'T_sim        : %d s  (dt = %.1f s, %d steps)\n', T_sim, dt, N_steps);
fprintf(fid,'MPC horizon  : Np = %d\n', Np);
fprintf(fid,'Cone half-ang: %d deg\n', cone_half_angle);
fprintf(fid,'\nRESULTS:\n');
fprintf(fid,'  Final dist   : %.2f m\n', log_dist(end));
fprintf(fid,'  LOS violations: %d / %d (%.1f%%)\n', n_los_viol, N_steps, 100*n_los_viol/N_steps);
fprintf(fid,'  Omega axis shift: %.1f deg\n', log_omega_dir(end));
fprintf(fid,'  Delta-v total: %.1f m/s\n', sum(vecnorm(log_u_tb)*dt));
fprintf(fid,'  Wall time    : %.1f s\n', wall_time);
fprintf(fid,'  QP solver    : quadprog\n');
fclose(fid);
fprintf('  summary.txt\n');

fprintf('\n=== All outputs saved to: %s ===\n', out_dir);
fprintf('  Total files: 2 GIFs + %d static figure pairs (PNG+PDF) + summary.txt\n', 10);

%% ====================================================================
%  SECTION 9 — LOCAL FUNCTIONS
%  ====================================================================

% ---- Quaternion (scalar-first) ----
function qr = quat_mult_fn(q1, q2)
    w1=q1(1); x1=q1(2); y1=q1(3); z1=q1(4);
    w2=q2(1); x2=q2(2); y2=q2(3); z2=q2(4);
    qr = [w1*w2-x1*x2-y1*y2-z1*z2;
          w1*x2+x1*w2+y1*z2-z1*y2;
          w1*y2-x1*z2+y1*w2+z1*x2;
          w1*z2+x1*y2-y1*x2+z1*w2];
end

function R = rotm_from_quat_fn(q)
    q=q/norm(q); w=q(1);x=q(2);y=q(3);z=q(4);
    R=[1-2*(y^2+z^2),2*(x*y-w*z),2*(x*z+w*y);
       2*(x*y+w*z),1-2*(x^2+z^2),2*(y*z-w*x);
       2*(x*z-w*y),2*(y*z+w*x),1-2*(x^2+y^2)];
end

% ---- Euler + quaternion propagation (RK4) ----
function [q_n, R_n, w_n] = propagate_euler_quat(q, w, Ib, dt)
    I1=Ib(1,1); I2=Ib(2,2); I3=Ib(3,3);
    function dy = eom(y)
        wv=y(1:3); qv=y(4:7);
        dw=[(I2-I3)/I1*wv(2)*wv(3);
            (I3-I1)/I2*wv(3)*wv(1);
            (I1-I2)/I3*wv(1)*wv(2)];
        dq=0.5*quat_mult_fn(qv,[0;wv]);
        dy=[dw;dq];
    end
    y0=[w;q];
    k1=eom(y0); k2=eom(y0+0.5*dt*k1);
    k3=eom(y0+0.5*dt*k2); k4=eom(y0+dt*k3);
    yn=y0+(dt/6)*(k1+2*k2+2*k3+k4);
    w_n=yn(1:3); q_n=yn(4:7); q_n=q_n/norm(q_n);
    R_n=rotm_from_quat_fn(q_n);
end

% ---- Body-frame linearisation ----
function [Ad, Bd] = linearise_3d_body(x, w, n, dt)
    eps_x=1e-6; eps_u=1e-6;
    x0=prop_body(x,[0;0;0],w,n,dt);
    Ad=zeros(6); Bd=zeros(6,3);
    for j=1:6
        xp=x; xp(j)=xp(j)+eps_x;
        Ad(:,j)=(prop_body(xp,[0;0;0],w,n,dt)-x0)/eps_x;
    end
    for j=1:3
        up=zeros(3,1); up(j)=eps_u;
        Bd(:,j)=(prop_body(x,up,w,n,dt)-x0)/eps_u;
    end
end

function xn = prop_body(x, u, w, n, dt)
    function dx = f(xv)
        r=xv(1:3); v=xv(4:6); n2=n^2;
        ag=[3*n2*r(1);0;-n2*r(3)];
        dv=ag-2*cross(w,v)-cross(w,cross(w,r))+u;
        dx=[v;dv];
    end
    k1=f(x); k2=f(x+.5*dt*k1); k3=f(x+.5*dt*k2); k4=f(x+dt*k3);
    xn=x+(dt/6)*(k1+2*k2+2*k3+k4);
end

% ---- ECI truth propagation ----
function xn = propagate_truth_fn(x, dt, mu, Re, J2, a, opts)
    function dxdt = eom(~, xv)
        r=xv(1:3);v=xv(4:6); rm=norm(r);
        a2=-mu/rm^3*r;
        z2=r(3)^2/rm^2; fac=-1.5*J2*mu*Re^2/rm^5;
        aj=fac*[(1-5*z2)*r(1);(1-5*z2)*r(2);(3-5*z2)*r(3)];
        dxdt=[v;a2+aj+a];
    end
    [~,X]=ode113(@eom,[0 dt],x,opts);
    xn=X(end,:)';
end

% ---- LVLH rotation ----
function R = lvlh_from_rv_fn(r, v)
    xh=r/norm(r); h=cross(r,v); zh=h/norm(h); yh=cross(zh,xh);
    R=[xh,yh,zh];
end

% ---- Reference (linear ramp to T_dock, then hold) ----
function xr = build_reference(t, dt, Np, ys, ye, T_dock)
    xr=zeros(6,Np+1);
    for j=0:Np
        tj = t + j*dt;
        xr(2,j+1) = ye + (ys-ye)*max(0, 1 - tj/T_dock);
    end
end

% ---- LOS check ----
function [ok, mrg] = check_los_fn(r, ck, ym, nf)
    x=r(1);y=r(2);z=r(3);
    if y<ym, ok=false; mrg=y-ym; return; end
    mrg=Inf;
    for f=1:nf
        th=2*pi*(f-1)/nf;
        mrg=min(mrg, -(cos(th)*x+sin(th)*z-ck*y));
    end
    ok=(mrg>=0);
end

% ---- MPC QP (quadprog) ----
function [u, st] = solve_mpc_qp(Ad,Bd,x0,xr,up,Q,QN,Ru,Rdu,am,ck,ym,nf,Np,nx,nu)
    nxt=(Np+1)*nx; nut=Np*nu; nv=nxt+nut;
    H=sparse(nv,nv); f=zeros(nv,1);
    for j=0:Np
        ix=j*nx+(1:nx);
        if j<Np, H(ix,ix)=Q; f(ix)=-Q*xr(:,j+1);
        else,    H(ix,ix)=QN; f(ix)=-QN*xr(:,j+1); end
    end
    for j=0:Np-1
        iu=nxt+j*nu+(1:nu);
        H(iu,iu)=H(iu,iu)+Ru;
        if j==0
            H(iu,iu)=H(iu,iu)+Rdu; f(iu)=f(iu)-Rdu*up;
        else
            iup=nxt+(j-1)*nu+(1:nu);
            H(iu,iu)=H(iu,iu)+Rdu; H(iup,iup)=H(iup,iup)+Rdu;
            H(iu,iup)=H(iu,iup)-Rdu; H(iup,iu)=H(iup,iu)-Rdu;
        end
    end
    H=(H+H')/2;
    % Dynamics equality
    neq=(Np+1)*nx; Ae=sparse(neq,nv); be=zeros(neq,1);
    Ae(1:nx,1:nx)=eye(nx); be(1:nx)=x0;
    for j=0:Np-1
        rr=(j+1)*nx+(1:nx);
        Ae(rr,(j+1)*nx+(1:nx))=eye(nx);
        Ae(rr,j*nx+(1:nx))=-Ad;
        Ae(rr,nxt+j*nu+(1:nu))=-Bd;
    end
    % Input bounds
    ni=2*Np*nu; Ai=sparse(ni,nv); bi=am*ones(ni,1); rc=0;
    for j=0:Np-1
        iu=nxt+j*nu+(1:nu);
        Ai(rc+(1:nu),iu)=eye(nu); rc=rc+nu;
        Ai(rc+(1:nu),iu)=-eye(nu); rc=rc+nu;
    end
    % LOS cone
    nl=(nf+1)*(Np+1); Al=sparse(nl,nv); bl=zeros(nl,1); rc=0;
    for j=0:Np
        ix=j*nx+(1:nx);
        for ff=1:nf
            th=2*pi*(ff-1)/nf; rc=rc+1;
            Al(rc,ix(1))=cos(th); Al(rc,ix(2))=-ck; Al(rc,ix(3))=sin(th);
        end
        rc=rc+1; Al(rc,ix(2))=-1; bl(rc)=-ym;
    end
    opts=optimoptions('quadprog','Display','off','MaxIter',2000,...
        'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-6);
    [z,~,ef]=quadprog(H,f,[Ai;Al],[bi;bl],Ae,be,[],[],[],opts);
    if ef>0, st='solved'; u=z(nxt+(1:nu));
    else,    st=sprintf('fail(%d)',ef); u=zeros(nu,1); end
end

% ---- Chaser attitude: +y toward target (origin) ----
function R = chaser_attitude(r_tb)
    y_chs = -r_tb / (norm(r_tb)+1e-12);   % +y toward target
    z_ref = [0; 0; 1];
    x_chs = cross(y_chs, z_ref);
    if norm(x_chs) < 1e-6
        z_ref = [1; 0; 0];
        x_chs = cross(y_chs, z_ref);
    end
    x_chs = x_chs / norm(x_chs);
    z_chs = cross(x_chs, y_chs);
    R = [x_chs, y_chs, z_chs];
end

% ---- Drawing helpers ----
function draw_cuboid_fn(c, d, col, al)
    cx=c(1);cy=c(2);cz=c(3); dx=d(1)/2;dy=d(2)/2;dz=d(3)/2;
    v=[cx-dx cy-dy cz-dz;cx+dx cy-dy cz-dz;cx+dx cy+dy cz-dz;cx-dx cy+dy cz-dz;
       cx-dx cy-dy cz+dz;cx+dx cy-dy cz+dz;cx+dx cy+dy cz+dz;cx-dx cy+dy cz+dz];
    f=[1 2 3 4;5 6 7 8;1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8];
    patch('Vertices',v,'Faces',f,'FaceColor',col,'FaceAlpha',al,'EdgeColor','k');
end

function [v, f] = cuboid_mesh(d)
    dx=d(1)/2;dy=d(2)/2;dz=d(3)/2;
    v=[-dx -dy -dz;dx -dy -dz;dx dy -dz;-dx dy -dz;
       -dx -dy  dz;dx -dy  dz;dx  dy  dz;-dx  dy  dz];
    f=[1 2 3 4;5 6 7 8;1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8];
end

function draw_oriented_box(pos, dims, R, col, al)
    [v0, f] = cuboid_mesh(dims);
    v = (R * v0')' + pos';
    patch('Vertices',v,'Faces',f,'FaceColor',col,'FaceAlpha',al,'EdgeColor','k');
end

function draw_axes_at(pos, R, len)
    cols = [1 0 0; 0 0.7 0; 0 0 1];
    lbls = {'x','y','z'};
    for i = 1:3
        d = R(:,i)*len;
        quiver3(pos(1),pos(2),pos(3), d(1),d(2),d(3), 0,...
            'Color',cols(i,:),'LineWidth',2.5,'MaxHeadSize',0.5);
        tp = pos + d*1.15;
        text(tp(1),tp(2),tp(3), lbls{i},'Color',cols(i,:),...
            'FontSize',8,'FontWeight','bold');
    end
end

function draw_los_cone_3d(ck, ym, ymax, ~, col, al)
    np=24; th=linspace(0,2*pi,np+1);
    yv=linspace(ym,ymax,15);
    for yi=1:length(yv)-1
        y1=yv(yi);r1=ck*y1; y2=yv(yi+1);r2=ck*y2;
        for ti=1:np
            xr=[r1*cos(th(ti)),r2*cos(th(ti+1)),r2*cos(th(ti+1)),r1*cos(th(ti))];
            yr=[y1,y2,y2,y1];
            zr=[r1*sin(th(ti)),r2*sin(th(ti+1)),r2*sin(th(ti+1)),r1*sin(th(ti))];
            fill3(xr,yr,zr,col,'FaceAlpha',al,'EdgeColor','none');
        end
    end
end

function draw_los_wireframe_thick(ck, ym, ymax, nf)
    % Prominent wireframe cone — thicker lines, more rings, always visible
    n_rings = 8;  % More rings for visibility
    n_lon   = max(nf, 24);  % Many longitudinal divisions
    th = linspace(0, 2*pi, n_lon+1);
    yv = linspace(ym, ymax, n_rings);

    % Draw rings (circular cross-sections)
    for yi = 1:n_rings
        r = ck * yv(yi);
        plot3(r*cos(th), yv(yi)*ones(size(th)), r*sin(th), ...
            '-','Color',[1 0.5 0],'LineWidth',1.5);
    end

    % Draw longitudinal lines (apex to base)
    n_longs = 8;  % 8 longitudinal spines
    th_long = linspace(0, 2*pi, n_longs+1);
    for ti = 1:n_longs
        plot3([0, ck*ymax*cos(th_long(ti))], [ym, ymax], [0, ck*ymax*sin(th_long(ti))],...
            '-','Color',[1 0.5 0],'LineWidth',1.2);
    end
end

function draw_los_wireframe_rotated(ck, ym, ymax, nf, R)
    % Wireframe cone rotated by R — for LVLH frame GIF
    n_rings = 6;
    n_lon   = max(nf, 24);
    th = linspace(0, 2*pi, n_lon+1);
    yv = linspace(ym, ymax, n_rings);

    % Rings
    for yi = 1:n_rings
        r = ck * yv(yi);
        pts_body = [r*cos(th); yv(yi)*ones(size(th)); r*sin(th)];
        pts_rot  = R * pts_body;
        plot3(pts_rot(1,:), pts_rot(2,:), pts_rot(3,:), ...
            '-','Color',[1 0.5 0],'LineWidth',1.2);
    end

    % Longitudinal spines
    n_longs = 8;
    th_long = linspace(0, 2*pi, n_longs+1);
    for ti = 1:n_longs
        pts_body = [[0; ym; 0], [ck*ymax*cos(th_long(ti)); ymax; ck*ymax*sin(th_long(ti))]];
        pts_rot  = R * pts_body;
        plot3(pts_rot(1,:), pts_rot(2,:), pts_rot(3,:), ...
            '-','Color',[1 0.5 0],'LineWidth',1.0);
    end
end
