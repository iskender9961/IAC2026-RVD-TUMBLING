%RUN_NOMINAL_REACHABILITY  Compute nominal deterministic reachability sets.
%
%   Performs the full sweep over omega and a_max combinations, computing:
%     - Forward reachable safe set (erosion-based)
%     - Backward safe-start set (erosion-based + discrete multi-step)
%
%   Saves results to data/ and figures to figures/.
%
%   Run:  >> run_nominal_reachability

if ~exist('CALLED_FROM_MASTER','var'), clear; close all; clc; end

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
root_dir = fileparts(this_dir);
if isempty(root_dir), root_dir = pwd; end
addpath(fullfile(root_dir, 'reachability_common'));
addpath(fullfile(root_dir, 'reachability_nominal'));
addpath(root_dir);  % for params.m

results_dir = fullfile(root_dir, 'data');
fig_dir     = fullfile(root_dir, 'figures');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ===== Parameters =====
rp = reachability_params();

omega_vals = rp.omega_vals_deg;
amax_vals  = rp.amax_vals;
n_omega = length(omega_vals);
n_amax  = length(amax_vals);

fprintf('====== Nominal Deterministic Reachability ======\n');
fprintf('Grid: %d omega x %d amax = %d combos\n', n_omega, n_amax, n_omega*n_amax);
fprintf('Evaluation grid: %d x %d = %d points\n', rp.nx_grid, rp.ny_grid, rp.nx_grid*rp.ny_grid);

%% ===== Forward Reachable Sets =====
fprintf('\n--- Forward Reachable Sets (Erosion-Based) ---\n');
all_forward = cell(n_omega, n_amax);

t_total = tic;
for iw = 1:n_omega
    for ia = 1:n_amax
        omega_deg = omega_vals(iw);
        amax = amax_vals(ia);
        fprintf('  Forward: omega=%d deg/s, a_max=%.2f m/s^2 ... ', omega_deg, amax);

        t_start = tic;
        result = forward_reachable_set(omega_deg, amax, rp);
        elapsed = toc(t_start);

        n_safe = sum(result.safe_mask(:) & result.cone_mask(:));
        n_cone = sum(result.cone_mask(:));
        fprintf('%d/%d safe (%.1f%%), %.1f s\n', n_safe, n_cone, ...
            100*n_safe/max(n_cone,1), elapsed);

        all_forward{iw, ia} = result;
    end
end
forward_time = toc(t_total);
fprintf('Forward total: %.1f s\n', forward_time);

%% ===== Backward Safe-Start Sets =====
fprintf('\n--- Backward Safe-Start Sets ---\n');
all_backward = cell(n_omega, n_amax);

t_total = tic;
for iw = 1:n_omega
    for ia = 1:n_amax
        omega_deg = omega_vals(iw);
        amax = amax_vals(ia);
        fprintf('  Backward: omega=%d deg/s, a_max=%.2f m/s^2 ... ', omega_deg, amax);

        t_start = tic;
        result = backward_safe_start_set(omega_deg, amax, rp);
        elapsed = toc(t_start);

        n_safe_cont = sum(result.safe_mask(:) & result.cone_mask(:));
        n_safe_disc = sum(result.safe_mask_disc(:) & result.cone_mask(:));
        n_cone = sum(result.cone_mask(:));
        fprintf('cont=%d disc=%d /%d (%.1f%%/%.1f%%), %.1f s\n', ...
            n_safe_cont, n_safe_disc, n_cone, ...
            100*n_safe_cont/max(n_cone,1), ...
            100*n_safe_disc/max(n_cone,1), elapsed);

        all_backward{iw, ia} = result;
    end
end
backward_time = toc(t_total);
fprintf('Backward total: %.1f s\n', backward_time);

%% ===== Save Data =====
save(fullfile(results_dir, 'nominal_reachability.mat'), ...
    'all_forward', 'all_backward', 'omega_vals', 'amax_vals', 'rp', ...
    'forward_time', 'backward_time');
fprintf('\nSaved: data/nominal_reachability.mat\n');

%% ===== Generate Figures =====
fprintf('\n--- Generating Figures ---\n');

% --- Individual heatmaps ---
for iw = 1:n_omega
    for ia = 1:n_amax
        fwd = all_forward{iw, ia};
        bwd = all_backward{iw, ia};

        fig = figure('Position', [100 100 1200 500], 'Visible', 'off');

        % Forward
        subplot(1,2,1);
        plot_reachability_heatmap(fwd.x_grid, fwd.y_grid, ...
            fwd.cone_mask, fwd.safe_mask, rp.cone_k, rp.y_min);
        title(sprintf('Forward Nominal\n\\omega_{tgt}=%d deg/s, a_{max}=%.2f m/s^2', ...
            fwd.omega_deg, fwd.a_max), 'FontSize', 11);

        % Backward (discrete)
        subplot(1,2,2);
        plot_reachability_heatmap(bwd.x_grid, bwd.y_grid, ...
            bwd.cone_mask, bwd.safe_mask_disc, rp.cone_k, rp.y_min);
        title(sprintf('Backward Nominal (50-step)\n\\omega_{tgt}=%d deg/s, a_{max}=%.2f m/s^2', ...
            bwd.omega_deg, bwd.a_max), 'FontSize', 11);

        fname = sprintf('nominal_w%d_a%.2f.png', fwd.omega_deg, fwd.a_max);
        exportgraphics(fig, fullfile(fig_dir, fname), 'Resolution', 200);
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end
end

% --- Combined grid: Forward ---
fig_grid = figure('Position', [50 50 1400 1500], 'Visible', 'off');
for iw = 1:n_omega
    for ia = 1:n_amax
        idx = (iw-1) * n_amax + ia;
        ax = subplot(n_omega, n_amax, idx);
        fwd = all_forward{iw, ia};
        plot_reachability_heatmap(fwd.x_grid, fwd.y_grid, ...
            fwd.cone_mask, fwd.safe_mask, rp.cone_k, rp.y_min);
        n_safe = sum(fwd.safe_mask(:) & fwd.cone_mask(:));
        n_cone = sum(fwd.cone_mask(:));
        title(sprintf('\\omega=%d, a_{max}=%.2f (%d/%d)', ...
            fwd.omega_deg, fwd.a_max, n_safe, n_cone), 'FontSize', 9);
        if iw < n_omega, xlabel(''); end
        if ia > 1, ylabel(''); end
    end
end
sgtitle('Nominal Forward Reachable Sets', 'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(fig_grid, fullfile(fig_dir, 'nominal_forward_grid.png'), ...
    'Resolution', 200);
close(fig_grid);
fprintf('  Saved: nominal_forward_grid.png\n');

% --- Combined grid: Backward (discrete) ---
fig_grid = figure('Position', [50 50 1400 1500], 'Visible', 'off');
for iw = 1:n_omega
    for ia = 1:n_amax
        idx = (iw-1) * n_amax + ia;
        ax = subplot(n_omega, n_amax, idx);
        bwd = all_backward{iw, ia};
        plot_reachability_heatmap(bwd.x_grid, bwd.y_grid, ...
            bwd.cone_mask, bwd.safe_mask_disc, rp.cone_k, rp.y_min);
        n_safe = sum(bwd.safe_mask_disc(:) & bwd.cone_mask(:));
        n_cone = sum(bwd.cone_mask(:));
        title(sprintf('\\omega=%d, a_{max}=%.2f (%d/%d)', ...
            bwd.omega_deg, bwd.a_max, n_safe, n_cone), 'FontSize', 9);
        if iw < n_omega, xlabel(''); end
        if ia > 1, ylabel(''); end
    end
end
sgtitle('Nominal Backward Safe-Start Sets (50-step)', 'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(fig_grid, fullfile(fig_dir, 'nominal_backward_grid.png'), ...
    'Resolution', 200);
close(fig_grid);
fprintf('  Saved: nominal_backward_grid.png\n');

fprintf('\n====== Nominal Reachability Complete ======\n');
fprintf('Total time: forward=%.1f s, backward=%.1f s\n', forward_time, backward_time);


%% ========== Helper Functions ==========

function plot_reachability_heatmap(x_grid, y_grid, cone_mask, safe_mask, cone_k, y_min)
%PLOT_REACHABILITY_HEATMAP  3-color heatmap for reachability results.
    hold on;

    ny = length(y_grid);
    nx = length(x_grid);

    % Build RGB image
    img = ones(ny, nx, 3);

    % Red = outside cone
    red_col = [0.9, 0.25, 0.25];
    for ch = 1:3
        layer = img(:,:,ch);
        layer(~cone_mask) = red_col(ch);
        img(:,:,ch) = layer;
    end

    % Blue = inside cone but not safe
    blue_col = [0.3, 0.4, 0.85];
    for ch = 1:3
        layer = img(:,:,ch);
        layer(cone_mask & ~safe_mask) = blue_col(ch);
        img(:,:,ch) = layer;
    end

    % Green = safe
    green_col = [0.2, 0.75, 0.2];
    for ch = 1:3
        layer = img(:,:,ch);
        layer(cone_mask & safe_mask) = green_col(ch);
        img(:,:,ch) = layer;
    end

    image(x_grid, y_grid, img);
    set(gca, 'YDir', 'normal');

    % Cone boundary
    y_cone = linspace(0, max(y_grid)*1.1, 200);
    plot(cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(-cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Target marker
    plot(0, 0, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.3 0.3], ...
        'HandleVisibility', 'off');

    % Legend
    h1 = fill(NaN, NaN, green_col, 'EdgeColor', 'none');
    h2 = fill(NaN, NaN, blue_col, 'EdgeColor', 'none');
    h3 = fill(NaN, NaN, red_col, 'EdgeColor', 'none');
    legend([h1 h2 h3], {'Certified safe', 'Inside cone (unsafe)', 'Outside LOS cone'}, ...
        'Location', 'southeast', 'FontSize', 8);

    xlabel('x_{TB} [m]', 'FontSize', 11);
    ylabel('y_{TB} [m]', 'FontSize', 11);
    xlim([-200 200]);
    ylim([-10 320]);
    grid on;
    set(gca, 'FontSize', 10);
end
