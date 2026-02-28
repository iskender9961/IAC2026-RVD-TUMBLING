%RUN_MONTE_CARLO  Feasibility region sweep: tumble rate x chaser a_max.
%
%   Sweeps over:
%     omega_target = [1, 2, 3, 4, 5] deg/s   (5 rows)
%     a_max_chaser = [0.2, 0.1, 0.05, 0.02]  (4 columns)
%
%   For each combo, runs ~195 MC simulations with parfor.
%   Saves data + individual figures + combined 5x4 grid.
%
%   Colors: green = all-time feasible, blue = initially feasible but fails,
%           red = outside cone (infeasible region).
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

%% ===== Sweep parameters =====
omega_vals_deg = [1, 2, 3, 4, 5];      % deg/s
amax_vals      = [0.2, 0.1, 0.05, 0.02]; % m/s^2
n_omega = length(omega_vals_deg);
n_amax  = length(amax_vals);

%% ===== Sample initial conditions =====
p0 = params();
cone_k = p0.cone_k;

n_y = 15;
n_x_per_y = 13;
y_vals = linspace(20, 300, n_y);

ic_list = [];
for iy = 1:length(y_vals)
    y0 = y_vals(iy);
    x_max = cone_k * y0 * 0.95;
    x_pts = linspace(-x_max, x_max, n_x_per_y);
    for ix = 1:length(x_pts)
        ic_list = [ic_list; x_pts(ix), y0]; %#ok<AGROW>
    end
end
n_mc = size(ic_list, 1);
fprintf('Grid: %d omega x %d amax = %d combos, %d MC each\n', ...
    n_omega, n_amax, n_omega*n_amax, n_mc);

%% ===== Start parallel pool =====
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
    fprintf('Started parallel pool with %d workers\n', pool.NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n', pool.NumWorkers);
end

%% ===== Run all combos =====
all_data = cell(n_omega, n_amax);

t_total = tic;
for iw = 1:n_omega
    for ia = 1:n_amax
        omega_deg = omega_vals_deg(iw);
        amax = amax_vals(ia);

        fprintf('\n====== omega=%.0f deg/s, a_max=%.2f m/s^2 ======\n', ...
            omega_deg, amax);

        % Setup params
        p_mc = params();
        p_mc.omega_body = [0; 0; deg2rad(omega_deg)];
        p_mc.u_max = amax;
        p_mc.Rdu = diag([1e5, 1e4, 1e5]);
        p_mc.Np   = 20;
        p_mc.Tsim = 400;
        p_mc.ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
        p_mc.osqp_max_iter = 5000;

        % Run MC
        outcomes = zeros(n_mc, 1);
        traj_x   = cell(n_mc, 1);
        traj_t   = cell(n_mc, 1);
        fail_x   = cell(n_mc, 1);
        term_reasons = cell(n_mc, 1);

        progress_queue = parallel.pool.DataQueue;
        t_start = tic;
        afterEach(progress_queue, @(d) report_progress(d, n_mc, t_start, omega_deg, amax));

        parfor ii = 1:n_mc
            x0 = ic_list(ii, 1);
            y0 = ic_list(ii, 2);

            p_run = p_mc;
            p_run.dr_lvlh0 = [x0; y0; 0];
            p_run.dv_lvlh0 = [0; 0; 0];
            p_run.y_hold_start = y0;
            p_run.y_hold_end   = 5;
            p_run.y_hold_tau   = y0;
            p_run.cone_draw_L  = 0;

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
            send(progress_queue, ii);
        end
        elapsed = toc(t_start);
        fprintf('  => %d/%d feasible (%.0f s)\n', sum(outcomes), n_mc, elapsed);

        mc.ic_list      = ic_list;
        mc.outcomes     = outcomes;
        mc.fail_x       = fail_x;
        mc.traj_x       = traj_x;
        mc.traj_t       = traj_t;
        mc.term_reasons = term_reasons;
        mc.cone_k       = cone_k;
        mc.p_mc         = p_mc;
        mc.n_y          = n_y;
        mc.n_x_per_y    = n_x_per_y;
        mc.y_vals       = y_vals;
        mc.omega_deg    = omega_deg;
        mc.amax         = amax;
        mc.elapsed_s    = elapsed;

        all_data{iw, ia} = mc;

        % Save individual data
        fname = sprintf('mc_data_w%d_a%.2f.mat', omega_deg, amax);
        save(fullfile(results_dir, fname), 'mc');
        fprintf('  Saved: %s\n', fname);

        % Save individual figure
        fig = plot_mc_heatmap(mc, false);
        fname_png = sprintf('fig_mc_w%d_a%.2f.png', omega_deg, amax);
        exportgraphics(fig, fullfile(results_dir, fname_png), 'Resolution', 200);
        close(fig);
        fprintf('  Saved: %s\n', fname_png);
    end
end
elapsed_total = toc(t_total);
fprintf('\n====== All sweeps done in %.0f s ======\n', elapsed_total);

% Save combined data
save(fullfile(results_dir, 'mc_sweep_all.mat'), 'all_data', ...
    'omega_vals_deg', 'amax_vals', 'ic_list', 'cone_k');
fprintf('Saved: mc_sweep_all.mat\n');

%% ===== Combined 5x4 grid figure =====
fig_grid = figure('Position', [50 50 1400 1500]);
for iw = 1:n_omega
    for ia = 1:n_amax
        idx = (iw-1) * n_amax + ia;
        ax = subplot(n_omega, n_amax, idx);
        plot_mc_heatmap_ax(ax, all_data{iw, ia});
    end
end
sgtitle('MC Feasibility Sweep: \omega_{target} vs a_{max,chaser}', ...
    'FontSize', 16, 'FontWeight', 'bold');

exportgraphics(fig_grid, fullfile(results_dir, 'fig_mc_sweep_grid.png'), ...
    'Resolution', 200);
fprintf('Saved: fig_mc_sweep_grid.png\n');

fprintf('\nAll MC outputs saved to: %s\n', results_dir);


%% ========== Helper Functions ==========

function fig = plot_mc_heatmap(mc, show_dots)
%PLOT_MC_HEATMAP  Standalone figure with 3-color feasibility heatmap.
    if nargin < 2, show_dots = true; end

    fig = figure('Position', [100 100 700 550]);
    ax = gca;
    plot_mc_heatmap_ax(ax, mc, show_dots);
end


function plot_mc_heatmap_ax(ax, mc, show_dots)
%PLOT_MC_HEATMAP_AX  Plot 3-color heatmap onto given axes.
%   Green = all-time feasible
%   Blue  = initially feasible but fails
%   Red   = outside cone (infeasible region background)
    if nargin < 3, show_dots = false; end

    axes(ax); hold on;

    ic_list   = mc.ic_list;
    outcomes  = mc.outcomes;
    cone_k    = mc.cone_k;
    omega_deg = mc.omega_deg;
    amax      = mc.amax;
    n_mc      = size(ic_list, 1);

    % Interpolate onto fine grid
    x_fine = linspace(-200, 200, 300);
    y_fine = linspace(0, 320, 300);
    [Xq, Yq] = meshgrid(x_fine, y_fine);

    % Cone mask
    cone_mask = abs(Xq) <= cone_k * Yq & Yq >= 1;

    % Build RGB image
    img = ones(length(y_fine), length(x_fine), 3);  % white background

    % Outside cone = red
    red_r = 0.9; red_g = 0.25; red_b = 0.25;
    for ch = 1:3
        layer = img(:,:,ch);
        c_val = [red_r, red_g, red_b];
        layer(~cone_mask) = c_val(ch);
        img(:,:,ch) = layer;
    end

    % Inside cone: interpolate outcome
    F = scatteredInterpolant(ic_list(:,1), ic_list(:,2), outcomes, ...
        'natural', 'nearest');
    Zq = F(Xq, Yq);
    Zq_bin = Zq >= 0.5;

    % Green = feasible, Blue = fails after start
    green_col = [0.2, 0.75, 0.2];
    blue_col  = [0.3, 0.4, 0.85];

    for ch = 1:3
        layer = img(:,:,ch);
        layer(cone_mask &  Zq_bin) = green_col(ch);
        layer(cone_mask & ~Zq_bin) = blue_col(ch);
        img(:,:,ch) = layer;
    end

    % Display
    image(ax, x_fine, y_fine, img);
    set(ax, 'YDir', 'normal');

    % Cone boundary
    y_cone = linspace(0, 350, 200);
    plot(ax,  cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(ax, -cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Small sample dots
    if show_dots
        idx_ok   = outcomes == 1;
        idx_fail = outcomes == 0;
        plot(ax, ic_list(idx_ok,1),   ic_list(idx_ok,2),   '.', ...
            'Color', [0 0.4 0], 'MarkerSize', 6, 'HandleVisibility', 'off');
        plot(ax, ic_list(idx_fail,1), ic_list(idx_fail,2), '.', ...
            'Color', [0 0 0.5], 'MarkerSize', 6, 'HandleVisibility', 'off');
    end

    % Target marker
    plot(ax, 0, 0, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.3 0.3], ...
        'HandleVisibility', 'off');

    % Legend (dummy patches)
    h1 = fill(ax, NaN, NaN, green_col, 'EdgeColor', 'none');
    h2 = fill(ax, NaN, NaN, blue_col,  'EdgeColor', 'none');
    h3 = fill(ax, NaN, NaN, [red_r red_g red_b], 'EdgeColor', 'none');
    legend(ax, [h1 h2 h3], {'All-time feasible', 'Fails after start', 'Outside LOS cone'}, ...
        'Location', 'northeast', 'FontSize', 9);

    xlabel(ax, 'x_{TB} [m]', 'FontSize', 11);
    ylabel(ax, 'y_{TB} [m]', 'FontSize', 11);
    title(ax, sprintf('\\omega_{tgt}=%d deg/s,  a_{max}=%.2f m/s^2  (%d/%d ok)', ...
        omega_deg, amax, sum(outcomes), n_mc), 'FontSize', 11);
    set(ax, 'FontSize', 10);
    xlim(ax, [-200 200]);
    ylim(ax, [-10 320]);
    grid(ax, 'on');
end


function report_progress(~, n_total, t_start, omega_deg, amax)
%REPORT_PROGRESS  Callback for parallel progress via DataQueue.
    persistent count
    if isempty(count), count = 0; end
    count = count + 1;
    if mod(count, 20) == 0 || count == n_total
        elapsed = toc(t_start);
        eta = elapsed / count * (n_total - count);
        fprintf('  [w=%d a=%.2f] %3d/%3d  %.0fs elapsed, ETA %.0fs\n', ...
            omega_deg, amax, count, n_total, elapsed, eta);
    end
    if count >= n_total, count = 0; end
end
