%RUN_RDU_COMPARISON  Compare 4 Rdu scenarios and plot together.
%
%   Rdu values: 1e5, 1e6, 1e7, 1e8
%   y-component of Rdu is 10x lower than x,z to allow approach maneuvers.
%
%   Run:  >> run_rdu_comparison

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'dynamics'));
addpath(fullfile(this_dir, 'frames'));
addpath(fullfile(this_dir, 'mpc'));
addpath(fullfile(this_dir, 'viz'));
addpath(fullfile(this_dir, 'utils'));

%% ===== Define scenarios =====
rdu_vals = [1e5, 1e6, 1e7, 1e8];
n_cases  = length(rdu_vals);
labels   = {'R_{\Delta u}=10^5', 'R_{\Delta u}=10^6', ...
            'R_{\Delta u}=10^7', 'R_{\Delta u}=10^8'};
colors   = {'b', 'r', [0 0.6 0], [0.8 0.4 0]};
styles   = {'-', '--', '-.', ':'};

results  = cell(n_cases, 1);
p_all    = cell(n_cases, 1);
reasons  = cell(n_cases, 1);

%% ===== Run scenarios =====
for ii = 1:n_cases
    p = params();

    % Rdu: y-component 10x lower than x,z
    rdu_xz = rdu_vals(ii);
    rdu_y  = rdu_xz / 10;
    p.Rdu  = diag([rdu_xz, rdu_y, rdu_xz]);

    fprintf('\n====== Scenario %d: Rdu_xz=%.0e, Rdu_y=%.0e ======\n', ...
        ii, rdu_xz, rdu_y);

    [lg, p_out, sim_term, reason] = run_sim_headless(p);

    results{ii}  = lg;
    p_all{ii}    = p_out;
    reasons{ii}  = reason;

    fprintf('  => %s  (t_final=%.0f s, y_final=%.1f m)\n', ...
        reason, lg.t_hist(end), lg.r_tb_hist(2, end));
end

%% ===== Save results =====
results_dir = fullfile(this_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
save(fullfile(results_dir, 'rdu_comparison.mat'), ...
     'results', 'p_all', 'reasons', 'rdu_vals', 'labels');
fprintf('\nComparison data saved.\n');

%% ===== Shared geometry =====
p0 = params();
cone_k = p0.cone_k;
y_max = p0.dr_lvlh0(2) * 1.2;
yy = linspace(0, y_max, 100);

%% ===== Figure 1: 6-panel time history comparison =====
fig1 = figure('Name','Rdu Comparison - Time Histories','Position',[50 50 1100 800]);

subplot(3,2,1); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    plot(lg.t_hist, lg.r_tb_hist(2,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.5);
end
ylabel('y_{TB} [m]');
title('Along-Axis Distance (y_T)');
legend(labels, 'Location', 'northeast', 'FontSize', 8);

subplot(3,2,2); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    rad_dev = sqrt(lg.r_tb_hist(1,:).^2 + lg.r_tb_hist(3,:).^2);
    plot(lg.t_hist, rad_dev, styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.5);
end
ylabel('\surd(x_T^2+z_T^2) [m]');
title('Radial Deviation');

subplot(3,2,3); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    Nc = length(lg.cost_hist);
    plot(lg.t_hist(1:Nc), lg.cost_hist, styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.5);
end
ylabel('J_k'); set(gca, 'YScale', 'log');
title('Stage Cost (log scale)');

subplot(3,2,4); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    u_mag = sqrt(sum(lg.u_hist.^2, 1));
    Nu = length(u_mag);
    plot(lg.t_hist(1:Nu), u_mag, styles{ii}, 'Color', colors{ii}, 'LineWidth', 1);
end
ylabel('||u|| [m/s^2]');
title('Control Magnitude');

subplot(3,2,5); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    v_mag = sqrt(sum(lg.v_tb_hist.^2, 1));
    plot(lg.t_hist, v_mag, styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.5);
end
xlabel('Time [s]'); ylabel('|v_{TB}| [m/s]');
title('Relative Velocity');

subplot(3,2,6); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    Nu = size(lg.u_hist, 2);
    if Nu > 1
        du = diff(lg.u_hist, 1, 2);
        du_mag = sqrt(sum(du.^2, 1));
        plot(lg.t_hist(2:Nu), du_mag, styles{ii}, 'Color', colors{ii}, 'LineWidth', 1);
    end
end
xlabel('Time [s]'); ylabel('|\Delta u| [m/s^2]');
title('Input Rate (\Delta u)');

sgtitle('R_{\Delta u} Comparison Study', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(results_dir, 'fig_rdu_comparison.png'));
fprintf('  Fig 1: time history comparison saved.\n');

%% ===== Figure 2: 3D trajectory overlay -- TB frame =====
fig2 = figure('Name','Rdu Comparison - 3D TB','Position',[100 50 900 750]);
hold on; grid on; axis equal;
draw_los_tetra([0;0;0], eye(3), cone_k, y_max, p0.cone_nfaces, [0.8 0.5 0]);
plot3([0 0], [0 y_max], [0 0], 'k--', 'LineWidth', 1.5);
plot3(0, 0, 0, 'kp', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
h_traj = gobjects(n_cases, 1);
for ii = 1:n_cases
    lg = results{ii};
    h_traj(ii) = plot3(lg.r_tb_hist(1,:), lg.r_tb_hist(2,:), lg.r_tb_hist(3,:), ...
        styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.5);
    plot3(lg.r_tb_hist(1,1), lg.r_tb_hist(2,1), lg.r_tb_hist(3,1), ...
        'o', 'Color', colors{ii}, 'MarkerSize', 6, 'MarkerFaceColor', colors{ii}, ...
        'HandleVisibility', 'off');
    plot3(lg.r_tb_hist(1,end), lg.r_tb_hist(2,end), lg.r_tb_hist(3,end), ...
        's', 'Color', colors{ii}, 'MarkerSize', 8, 'MarkerFaceColor', colors{ii}, ...
        'HandleVisibility', 'off');
end
xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]'); zlabel('z_{TB} [m]');
title('3D Trajectory Comparison -- Target Body Frame');
legend(h_traj, labels, 'Location', 'best', 'FontSize', 8);
view(135, 25);
saveas(fig2, fullfile(results_dir, 'fig_rdu_3d_tb.png'));
fprintf('  Fig 2: 3D TB comparison saved.\n');

%% ===== Figure 3: 2D XY -- TB frame (with cone) =====
fig3 = figure('Name','Rdu Comparison - TB XY','Position',[100 100 900 700]);
hold on; grid on; axis equal;
fill([cone_k*yy, fliplr(-cone_k*yy)], [yy, fliplr(yy)], ...
    [1 0.95 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot( cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2, 'HandleVisibility', 'off');
plot(-cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2, 'HandleVisibility', 'off');
plot([0 0], [0 y_max], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
h_traj = gobjects(n_cases, 1);
for ii = 1:n_cases
    lg = results{ii};
    h_traj(ii) = plot(lg.r_tb_hist(1,:), lg.r_tb_hist(2,:), styles{ii}, ...
        'Color', colors{ii}, 'LineWidth', 1.5);
    plot(lg.r_tb_hist(1,end), lg.r_tb_hist(2,end), 'o', ...
        'Color', colors{ii}, 'MarkerSize', 8, 'MarkerFaceColor', colors{ii}, ...
        'HandleVisibility', 'off');
end
plot(0, p0.dr_lvlh0(2), 'g^', 'MarkerSize', 12, 'MarkerFaceColor', 'g', ...
    'HandleVisibility', 'off');
xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]');
title('2D XY Trajectory Comparison -- Target Body Frame');
legend(h_traj, labels, 'Location', 'best', 'FontSize', 8);
saveas(fig3, fullfile(results_dir, 'fig_rdu_trajectories_xy.png'));
fprintf('  Fig 3: 2D XY TB comparison saved.\n');

%% ===== Figure 4: Control input comparison -- TB frame =====
fig4 = figure('Name','Rdu Comparison - Control TB','Position',[50 50 1100 850]);
for ax = 1:3
    subplot(3,1,ax); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        Nu = size(lg.u_hist, 2);
        plot(lg.t_hist(1:Nu), lg.u_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1);
    end
    ax_lab = {'x','y','z'};
    ylabel(sprintf('u_{%s,TB} [m/s^2]', ax_lab{ax}));
    if ax == 1
        title('Control Input -- Target Body Frame');
        legend(labels, 'Location', 'best', 'FontSize', 8);
    end
    if ax == 3
        xlabel('Time [s]');
        ylim([-0.2 0.2]);  % fix z-axis scale for readability
    end
end
saveas(fig4, fullfile(results_dir, 'fig_rdu_control_tb.png'));
fprintf('  Fig 4: control TB comparison saved.\n');

%% ===== Figure 5: Control input comparison -- Chaser Body frame =====
fig5 = figure('Name','Rdu Comparison - Control CB','Position',[100 50 1100 850]);
for ax = 1:3
    subplot(3,1,ax); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        Nu = size(lg.u_hist, 2);
        plot(lg.t_hist(1:Nu), lg.u_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1);
    end
    ax_lab = {'x','y','z'};
    ylabel(sprintf('u_{%s,CB} [m/s^2]', ax_lab{ax}));
    if ax == 1
        title('Control Input -- Chaser Body Frame');
        legend(labels, 'Location', 'best', 'FontSize', 8);
    end
    if ax == 3, xlabel('Time [s]'); end
end
saveas(fig5, fullfile(results_dir, 'fig_rdu_control_cb.png'));
fprintf('  Fig 5: control CB comparison saved.\n');

%% ===== Figure 6: Control input comparison -- LVLH frame =====
fig6 = figure('Name','Rdu Comparison - Control LVLH','Position',[150 50 1100 850]);
for ax = 1:3
    subplot(3,1,ax); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        Nu = size(lg.u_lvlh_hist, 2);
        plot(lg.t_hist(1:Nu), lg.u_lvlh_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1);
    end
    ax_lab = {'x','y','z'};
    ylabel(sprintf('u_{%s,LVLH} [m/s^2]', ax_lab{ax}));
    if ax == 1
        title('Control Input -- LVLH Frame');
        legend(labels, 'Location', 'best', 'FontSize', 8);
    end
    if ax == 3, xlabel('Time [s]'); end
end
saveas(fig6, fullfile(results_dir, 'fig_rdu_control_lvlh.png'));
fprintf('  Fig 6: control LVLH comparison saved.\n');

%% ===== Figure 7: Delta-u comparison (components) =====
fig7 = figure('Name','Rdu Comparison - DeltaU','Position',[100 50 1100 850]);
for ax = 1:3
    subplot(3,1,ax); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        Nu = size(lg.u_hist, 2);
        if Nu > 1
            du = diff(lg.u_hist, 1, 2);
            plot(lg.t_hist(2:Nu), du(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1);
        end
    end
    ax_lab = {'x','y','z'};
    ylabel(sprintf('\\Delta u_{%s} [m/s^2]', ax_lab{ax}));
    if ax == 1
        title('\Delta u Components -- All Scenarios');
        legend(labels, 'Location', 'best', 'FontSize', 8);
    end
    if ax == 3, xlabel('Time [s]'); end
end
saveas(fig7, fullfile(results_dir, 'fig_rdu_deltau_components.png'));
fprintf('  Fig 7: delta-u comparison saved.\n');

%% ===== Figure 8: Cumulative delta-V comparison =====
fig8 = figure('Name','Rdu Comparison - DeltaV','Position',[150 50 800 500]);
hold on; grid on;
h_traj = gobjects(n_cases, 1);
for ii = 1:n_cases
    lg = results{ii};
    Nu = size(lg.u_hist, 2);
    dv_cum = cumsum(sqrt(sum(lg.u_hist.^2, 1))) * p_all{ii}.dt;
    h_traj(ii) = plot(lg.t_hist(1:Nu), dv_cum, styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.5);
end
xlabel('Time [s]'); ylabel('\Delta V_{cum} [m/s]');
title('Cumulative \Delta V -- All Scenarios');
legend(h_traj, labels, 'Location', 'best', 'FontSize', 8);
saveas(fig8, fullfile(results_dir, 'fig_rdu_deltav.png'));
fprintf('  Fig 8: cumulative delta-V comparison saved.\n');

%% ===== Figure 9: Summary table =====
fig9 = figure('Name','Rdu Summary','Position',[150 150 700 250]);
axes('Visible','off');
txt = sprintf('%-20s | %8s | %8s | %8s | %10s | %s\n', ...
    'Scenario', 't_end[s]', 'y_end[m]', 'rdev[m]', 'DeltaV[m/s]', 'Status');
txt = [txt, repmat('-',1,80), newline];
for ii = 1:n_cases
    lg = results{ii};
    t_end = lg.t_hist(end);
    y_end = lg.r_tb_hist(2, end);
    rd_end = sqrt(lg.r_tb_hist(1,end)^2 + lg.r_tb_hist(3,end)^2);
    dv_total = sum(sqrt(sum(lg.u_hist.^2, 1))) * p_all{ii}.dt;
    txt = [txt, sprintf('Rdu=%.0e (y/10)     | %8.0f | %8.1f | %8.2f | %10.2f | %s\n', ...
        rdu_vals(ii), t_end, y_end, rd_end, dv_total, reasons{ii})];
end
text(0.05, 0.5, txt, 'FontName', 'Courier', 'FontSize', 10, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
saveas(fig9, fullfile(results_dir, 'fig_rdu_summary.png'));
fprintf('  Fig 9: summary table saved.\n');

%% ===== Figure 10: 3D trajectory overlay -- Chaser Body frame =====
fig10 = figure('Name','Rdu Comparison - 3D CB','Position',[200 50 900 750]);
hold on; grid on; axis equal;
draw_los_tetra([0;0;0], eye(3), cone_k, y_max, p0.cone_nfaces, [0.8 0.5 0]);
plot3([0 0], [0 y_max], [0 0], 'k--', 'LineWidth', 1.5);
plot3(0, 0, 0, 'kp', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
h_traj = gobjects(n_cases, 1);
for ii = 1:n_cases
    lg = results{ii};
    h_traj(ii) = plot3(lg.r_tb_hist(1,:), lg.r_tb_hist(2,:), lg.r_tb_hist(3,:), ...
        styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.5);
    plot3(lg.r_tb_hist(1,1), lg.r_tb_hist(2,1), lg.r_tb_hist(3,1), ...
        'o', 'Color', colors{ii}, 'MarkerSize', 6, 'MarkerFaceColor', colors{ii}, ...
        'HandleVisibility', 'off');
    plot3(lg.r_tb_hist(1,end), lg.r_tb_hist(2,end), lg.r_tb_hist(3,end), ...
        's', 'Color', colors{ii}, 'MarkerSize', 8, 'MarkerFaceColor', colors{ii}, ...
        'HandleVisibility', 'off');
end
xlabel('x_{CB} [m]'); ylabel('y_{CB} [m]'); zlabel('z_{CB} [m]');
title('3D Trajectory Comparison -- Chaser Body Frame');
legend(h_traj, labels, 'Location', 'best', 'FontSize', 8);
view(135, 25);
saveas(fig10, fullfile(results_dir, 'fig_rdu_3d_cb.png'));
fprintf('  Fig 10: 3D CB comparison saved.\n');

%% ===== Figure 11: Position time histories -- TB + LVLH =====
fig11 = figure('Name','Rdu Comparison - Positions','Position',[50 50 1100 900]);
labels_tb   = {'x_{TB}','y_{TB}','z_{TB}'};
labels_lvlh = {'x_{LVLH}','y_{LVLH}','z_{LVLH}'};
for ax = 1:3
    % Target body
    subplot(3,2,(ax-1)*2+1); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        plot(lg.t_hist, lg.r_tb_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.2);
    end
    ylabel([labels_tb{ax} ' [m]']);
    title([labels_tb{ax} ' -- Target Body']);
    if ax == 1, legend(labels, 'Location', 'best', 'FontSize', 7); end
    if ax == 3
        xlabel('Time [s]');
        ylim([-50 50]);  % fix z-axis scale for readability
    end

    % LVLH
    subplot(3,2,(ax-1)*2+2); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        plot(lg.t_hist, lg.r_lvlh_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.2);
    end
    ylabel([labels_lvlh{ax} ' [m]']);
    title([labels_lvlh{ax} ' -- LVLH']);
    if ax == 3
        xlabel('Time [s]');
        ylim([-50 50]);  % fix z-axis scale for readability
    end
end
sgtitle('Position Components -- All Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig11, fullfile(results_dir, 'fig_rdu_positions.png'));
fprintf('  Fig 11: position time histories saved.\n');

%% ===== Figure 12: Velocity time histories -- TB + LVLH =====
fig12 = figure('Name','Rdu Comparison - Velocities','Position',[100 50 1100 900]);
vlabels_tb   = {'vx_{TB}','vy_{TB}','vz_{TB}'};
vlabels_lvlh = {'vx_{LVLH}','vy_{LVLH}','vz_{LVLH}'};
for ax = 1:3
    % Target body
    subplot(3,2,(ax-1)*2+1); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        plot(lg.t_hist, lg.v_tb_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.2);
    end
    ylabel([vlabels_tb{ax} ' [m/s]']);
    title([vlabels_tb{ax} ' -- Target Body']);
    if ax == 1, legend(labels, 'Location', 'best', 'FontSize', 7); end
    if ax == 3, xlabel('Time [s]'); end

    % LVLH
    subplot(3,2,(ax-1)*2+2); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        plot(lg.t_hist, lg.v_lvlh_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.2);
    end
    ylabel([vlabels_lvlh{ax} ' [m/s]']);
    title([vlabels_lvlh{ax} ' -- LVLH']);
    if ax == 3, xlabel('Time [s]'); end
end
sgtitle('Velocity Components -- All Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig12, fullfile(results_dir, 'fig_rdu_velocities.png'));
fprintf('  Fig 12: velocity time histories saved.\n');

fprintf('\nAll static comparison plots saved to results/\n');

%% ===== Comparison GIFs =====
fprintf('Generating comparison GIFs...\n');
generate_comparison_gifs(results, p_all, labels, colors, styles, results_dir);

fprintf('Done.\n');
