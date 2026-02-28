%RUN_MONTE_CARLO  Feasibility region analysis via brute-force MC.
%
%   Samples ~100 initial positions (x_TB, y_TB) across the LOS cone and
%   checks whether the MPC controller keeps the chaser feasible throughout.
%
%   Rdu = diag([1e5, 1e4, 1e5])  (the 1e5 case with y/10)
%
%   Outputs 2 figures:
%     1) Heatmap of initial conditions colored by outcome
%     2) Trajectories with red X at failure points
%
%   Run:  >> run_monte_carlo

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'dynamics'));
addpath(fullfile(this_dir, 'frames'));
addpath(fullfile(this_dir, 'mpc'));
addpath(fullfile(this_dir, 'viz'));
addpath(fullfile(this_dir, 'utils'));

results_dir = fullfile(this_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% ===== Sample initial conditions =====
% Sample (x_TB, y_TB) across the LOS cone.
% z_TB = 0 for all (2D slice in the x-y plane of TB frame).
% Cone constraint: |x| <= cone_k * y,  y >= y_min.

p0 = params();
cone_k = p0.cone_k;   % tan(30 deg) ~ 0.577

% Grid: y from 20 to 300 m, x from -cone_k*y to +cone_k*y
y_vals = linspace(20, 300, 12);     % 12 y-levels
n_x_per_y = 9;                       % 9 x-points per y-level => ~108 samples
% Also add some points near the cone wall
ic_list = [];
for iy = 1:length(y_vals)
    y0 = y_vals(iy);
    x_max = cone_k * y0 * 0.95;  % stay 5% inside to be initially feasible
    x_pts = linspace(-x_max, x_max, n_x_per_y);
    for ix = 1:length(x_pts)
        ic_list = [ic_list; x_pts(ix), y0]; %#ok<AGROW>
    end
end
n_mc = size(ic_list, 1);
fprintf('Monte Carlo: %d initial conditions\n', n_mc);

%% ===== Speed-optimized params =====
p_mc = params();
% Rdu = 1e5 case (y/10)
p_mc.Rdu = diag([1e5, 1e4, 1e5]);
% Speed optimizations
p_mc.Np   = 20;                     % shorter horizon
p_mc.Tsim = 400;                    % shorter sim
p_mc.ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
p_mc.osqp_max_iter = 5000;

%% ===== Run MC (parallel) =====
outcomes  = zeros(n_mc, 1);   % 1=success, 0=failed
fail_x    = cell(n_mc, 1);    % failure point [x; y] in TB
traj_x    = cell(n_mc, 1);    % full trajectory in TB (x, y)
traj_t    = cell(n_mc, 1);
term_reasons = cell(n_mc, 1);

% Start parallel pool if not already running
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
    fprintf('Started parallel pool with %d workers\n', pool.NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n', pool.NumWorkers);
end

% Progress tracking via DataQueue
progress_queue = parallel.pool.DataQueue;
completed_count = 0;
t_start = tic;
afterEach(progress_queue, @(data) report_progress(data, n_mc, t_start));

fprintf('Running %d MC cases on %d workers...\n', n_mc, pool.NumWorkers);

parfor ii = 1:n_mc
    x0 = ic_list(ii, 1);
    y0 = ic_list(ii, 2);

    p_run = p_mc;
    p_run.dr_lvlh0 = [x0; y0; 0];
    p_run.dv_lvlh0 = [0; 0; 0];

    % Adjust reference to start from this y
    p_run.y_hold_start = y0;
    p_run.y_hold_end   = 5;
    p_run.y_hold_tau   = y0;  % scale tau with distance

    % Reset cone draw length
    p_run.cone_draw_L = 0;

    try
        [lg, ~, ~, reason] = run_sim_headless(p_run);

        traj_x{ii} = lg.r_tb_hist([1 2], :);
        traj_t{ii} = lg.t_hist;
        term_reasons{ii} = reason;

        if strcmp(reason, 'completed') || strcmp(reason, 'docked')
            outcomes(ii) = 1;
        else
            outcomes(ii) = 0;
            fail_x{ii} = lg.r_tb_hist([1 2], end);
        end
    catch ME
        outcomes(ii) = 0;
        term_reasons{ii} = sprintf('error: %s', ME.message);
        traj_x{ii} = [x0; y0];
        traj_t{ii} = 0;
        fail_x{ii} = [x0; y0];
    end

    % Send progress update
    send(progress_queue, ii);
end
elapsed_total = toc(t_start);
fprintf('\nMonte Carlo done: %d/%d succeeded (%.0f s total, %d workers)\n', ...
    sum(outcomes), n_mc, elapsed_total, pool.NumWorkers);

%% ===== Figure 1: Heatmap of outcomes =====
fig1 = figure('Position', [100 100 700 550]);
hold on; grid on;

% Draw LOS cone boundary (full cone)
y_cone = linspace(0, 350, 200);
x_cone_pos = cone_k * y_cone;
x_cone_neg = -cone_k * y_cone;
fill([x_cone_neg, fliplr(x_cone_pos)], [y_cone, fliplr(y_cone)], ...
    [0.9 0.9 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.4, ...
    'HandleVisibility', 'off');
plot(x_cone_pos, y_cone, 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(x_cone_neg, y_cone, 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Plot initial conditions colored by outcome
idx_ok   = outcomes == 1;
idx_fail = outcomes == 0;
scatter(ic_list(idx_ok, 1), ic_list(idx_ok, 2), 60, ...
    [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
scatter(ic_list(idx_fail, 1), ic_list(idx_fail, 2), 60, ...
    'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

% Add origin (target)
plot(0, 0, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', [0.3 0.3 0.3], ...
    'HandleVisibility', 'off');
text(2, -8, 'Target', 'FontSize', 10);

legend({'Feasible (all time)', 'Failed'}, ...
    'Location', 'northeast', 'FontSize', 11);
xlabel('x_{TB} [m]', 'FontSize', 12);
ylabel('y_{TB} [m]', 'FontSize', 12);
title(sprintf('MC Feasibility Map  (R_{\\Delta u}=10^5,  \\omega=%.1f deg/s)', ...
    rad2deg(norm(p_mc.omega_body))), 'FontSize', 13);
set(gca, 'FontSize', 11);
xlim([-200 200]);
ylim([-10 320]);

exportgraphics(fig1, fullfile(results_dir, 'fig_mc_feasibility_map.png'), ...
    'Resolution', 200);
fprintf('Saved: fig_mc_feasibility_map.png\n');

%% ===== Figure 2: Trajectories with failure markers =====
fig2 = figure('Position', [100 100 700 550]);
hold on; grid on;

% Draw LOS cone
fill([x_cone_neg, fliplr(x_cone_pos)], [y_cone, fliplr(y_cone)], ...
    [0.9 0.9 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
    'HandleVisibility', 'off');
plot(x_cone_pos, y_cone, 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(x_cone_neg, y_cone, 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Plot trajectories
h_ok = gobjects(0);
h_fail = gobjects(0);
h_failx = gobjects(0);

for ii = 1:n_mc
    xy = traj_x{ii};
    if isempty(xy), continue; end
    if outcomes(ii) == 1
        hh = plot(xy(1,:), xy(2,:), '-', 'Color', [0 0.7 0 0.4], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
        if isempty(h_ok), h_ok = hh; set(hh, 'HandleVisibility', 'on'); end
    else
        hh = plot(xy(1,:), xy(2,:), '-', 'Color', [1 0 0 0.5], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
        if isempty(h_fail), h_fail = hh; set(hh, 'HandleVisibility', 'on'); end
        % Red X at failure point
        hx = plot(xy(1,end), xy(2,end), 'rx', 'MarkerSize', 8, ...
            'LineWidth', 2, 'HandleVisibility', 'off');
        if isempty(h_failx), h_failx = hx; set(hx, 'HandleVisibility', 'on'); end
    end
end

% Add origin (target)
plot(0, 0, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', [0.3 0.3 0.3], ...
    'HandleVisibility', 'off');

% Build legend
leg_h = [];
leg_str = {};
if ~isempty(h_ok)
    leg_h(end+1) = h_ok; leg_str{end+1} = 'Feasible trajectory';
end
if ~isempty(h_fail)
    leg_h(end+1) = h_fail; leg_str{end+1} = 'Failed trajectory';
end
if ~isempty(h_failx)
    leg_h(end+1) = h_failx; leg_str{end+1} = 'Failure point';
end
if ~isempty(leg_h)
    legend(leg_h, leg_str, 'Location', 'northeast', 'FontSize', 11);
end

xlabel('x_{TB} [m]', 'FontSize', 12);
ylabel('y_{TB} [m]', 'FontSize', 12);
title(sprintf('MC Trajectories  (R_{\\Delta u}=10^5,  \\omega=%.1f deg/s)', ...
    rad2deg(norm(p_mc.omega_body))), 'FontSize', 13);
set(gca, 'FontSize', 11);
xlim([-200 200]);
ylim([-10 320]);

exportgraphics(fig2, fullfile(results_dir, 'fig_mc_trajectories.png'), ...
    'Resolution', 200);
fprintf('Saved: fig_mc_trajectories.png\n');

fprintf('\nAll MC figures saved to: %s\n', results_dir);


function report_progress(~, n_total, t_start)
%REPORT_PROGRESS  Callback for parallel progress via DataQueue.
    persistent count
    if isempty(count), count = 0; end
    count = count + 1;
    if mod(count, 5) == 0 || count == n_total
        elapsed = toc(t_start);
        eta = elapsed / count * (n_total - count);
        fprintf('  [%3d/%3d]  %.0fs elapsed, ETA %.0fs\n', ...
            count, n_total, elapsed, eta);
    end
    if count >= n_total, count = 0; end  % reset for next run
end
