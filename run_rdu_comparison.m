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

%% ===== Comparison plots =====

% --- Figure 1: y_T vs time ---
fig1 = figure('Name','Rdu Comparison - y_T','Position',[50 50 1100 800]);

subplot(3,2,1); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    plot(lg.t_hist, lg.r_tb_hist(2,:), '-', 'Color', colors{ii}, 'LineWidth', 1.5);
end
ylabel('y_{TB} [m]');
title('Along-Axis Distance (y_T)');
legend(labels, 'Location', 'northeast', 'FontSize', 8);

% --- Radial deviation ---
subplot(3,2,2); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    rad_dev = sqrt(lg.r_tb_hist(1,:).^2 + lg.r_tb_hist(3,:).^2);
    plot(lg.t_hist, rad_dev, '-', 'Color', colors{ii}, 'LineWidth', 1.5);
end
ylabel('\surd(x_T^2+z_T^2) [m]');
title('Radial Deviation');

% --- Cost vs time ---
subplot(3,2,3); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    Nc = length(lg.cost_hist);
    plot(lg.t_hist(1:Nc), lg.cost_hist, '-', 'Color', colors{ii}, 'LineWidth', 1.5);
end
ylabel('J_k'); set(gca, 'YScale', 'log');
title('Stage Cost (log scale)');

% --- Control magnitude ---
subplot(3,2,4); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    u_mag = sqrt(sum(lg.u_hist.^2, 1));
    Nu = length(u_mag);
    plot(lg.t_hist(1:Nu), u_mag, '-', 'Color', colors{ii}, 'LineWidth', 1);
end
ylabel('||u|| [m/s^2]');
title('Control Magnitude');

% --- Velocity magnitude ---
subplot(3,2,5); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    v_mag = sqrt(sum(lg.v_tb_hist.^2, 1));
    plot(lg.t_hist, v_mag, '-', 'Color', colors{ii}, 'LineWidth', 1.5);
end
xlabel('Time [s]'); ylabel('|v_{TB}| [m/s]');
title('Relative Velocity');

% --- Delta-u magnitude ---
subplot(3,2,6); hold on; grid on;
for ii = 1:n_cases
    lg = results{ii};
    Nu = size(lg.u_hist, 2);
    if Nu > 1
        du = diff(lg.u_hist, 1, 2);
        du_mag = sqrt(sum(du.^2, 1));
        plot(lg.t_hist(2:Nu), du_mag, '-', 'Color', colors{ii}, 'LineWidth', 1);
    end
end
xlabel('Time [s]'); ylabel('|\Delta u| [m/s^2]');
title('Input Rate (\Delta u)');

sgtitle('R_{\Delta u} Comparison Study', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(results_dir, 'fig_rdu_comparison.png'));
fprintf('Comparison plot saved.\n');

% --- Figure 2: 2D trajectory overlay in TB xy ---
fig2 = figure('Name','Rdu Comparison - TB XY','Position',[100 100 900 700]);
hold on; grid on; axis equal;

% Cone
p0 = params();
y_max = p0.dr_lvlh0(2) * 1.2;
yy = linspace(0, y_max, 100);
fill([p0.cone_k*yy, fliplr(-p0.cone_k*yy)], [yy, fliplr(yy)], ...
    [1 0.95 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot( p0.cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
plot(-p0.cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
plot([0 0], [0 y_max], 'k--', 'LineWidth', 1.5);

for ii = 1:n_cases
    lg = results{ii};
    plot(lg.r_tb_hist(1,:), lg.r_tb_hist(2,:), '-', ...
        'Color', colors{ii}, 'LineWidth', 1.5);
    plot(lg.r_tb_hist(1,end), lg.r_tb_hist(2,end), 'o', ...
        'Color', colors{ii}, 'MarkerSize', 8, 'MarkerFaceColor', colors{ii});
end
plot(0, p0.dr_lvlh0(2), 'g^', 'MarkerSize', 12, 'MarkerFaceColor', 'g');

xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]');
title('TB Frame XY -- All Scenarios');
legend([{'Cone','','','Axis'}, labels], 'Location', 'best', 'FontSize', 8);
saveas(fig2, fullfile(results_dir, 'fig_rdu_trajectories_xy.png'));

% --- Figure 3: Summary table as text ---
fig3 = figure('Name','Rdu Summary','Position',[150 150 700 250]);
ax = axes('Visible','off');
txt = sprintf('%-20s | %8s | %8s | %8s | %10s | %s\n', ...
    'Scenario', 't_end[s]', 'y_end[m]', 'rdev[m]', 'DeltaV[m/s]', 'Status');
txt = [txt, repmat('-',1,80), newline];
for ii = 1:n_cases
    lg = results{ii};
    t_end = lg.t_hist(end);
    y_end = lg.r_tb_hist(2, end);
    rd_end = sqrt(lg.r_tb_hist(1,end)^2 + lg.r_tb_hist(3,end)^2);
    % Total delta-v = sum(||u|| * dt)
    dv_total = sum(sqrt(sum(lg.u_hist.^2, 1))) * p_all{ii}.dt;
    txt = [txt, sprintf('Rdu=%.0e (y/10)     | %8.0f | %8.1f | %8.2f | %10.2f | %s\n', ...
        rdu_vals(ii), t_end, y_end, rd_end, dv_total, reasons{ii})];
end
text(0.05, 0.5, txt, 'FontName', 'Courier', 'FontSize', 10, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
saveas(fig3, fullfile(results_dir, 'fig_rdu_summary.png'));

fprintf('\nAll comparison plots saved to results/\n');
fprintf('Done.\n');
