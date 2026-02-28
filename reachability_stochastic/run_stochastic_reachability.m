%RUN_STOCHASTIC_REACHABILITY  Compute stochastic reachability sets.
%
%   Performs the full sweep over omega and a_max combinations, computing
%   chance-constrained safe sets with Gaussian disturbance tightening.
%
%   Run:  >> run_stochastic_reachability

if ~exist('CALLED_FROM_MASTER','var'), clear; close all; clc; end

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
root_dir = fileparts(this_dir);
if isempty(root_dir), root_dir = pwd; end
addpath(fullfile(root_dir, 'reachability_common'));
addpath(fullfile(root_dir, 'reachability_stochastic'));
addpath(root_dir);

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

fprintf('====== Stochastic Reachability ======\n');
fprintf('alpha = %.3f (%.1f%% confidence)\n', rp.alpha, 100*(1-rp.alpha));
fprintf('W_pos_std = %.4f m, W_vel_std = %.6f m/s\n', rp.W_pos_std, rp.W_vel_std);

%% ===== Sweep =====
all_stochastic = cell(n_omega, n_amax);

t_total = tic;
for iw = 1:n_omega
    for ia = 1:n_amax
        omega_deg = omega_vals(iw);
        amax = amax_vals(ia);
        fprintf('  omega=%d deg/s, a_max=%.2f m/s^2 ... ', omega_deg, amax);

        t_start = tic;
        result = stochastic_safe_set(omega_deg, amax, rp);
        elapsed = toc(t_start);

        n_safe = sum(result.safe_mask(:) & result.cone_mask(:));
        n_cone = sum(result.cone_mask(:));
        fprintf('%d/%d safe (%.1f%%), %.1f s\n', n_safe, n_cone, ...
            100*n_safe/max(n_cone,1), elapsed);

        all_stochastic{iw, ia} = result;
    end
end
stoch_time = toc(t_total);

%% ===== Save Data =====
save(fullfile(results_dir, 'stochastic_reachability.mat'), ...
    'all_stochastic', 'omega_vals', 'amax_vals', 'rp', 'stoch_time');
fprintf('\nSaved: data/stochastic_reachability.mat\n');

%% ===== Generate Grid Figure =====
fig_grid = figure('Position', [50 50 1400 1500], 'Visible', 'off');
for iw = 1:n_omega
    for ia = 1:n_amax
        idx = (iw-1) * n_amax + ia;
        subplot(n_omega, n_amax, idx);
        r = all_stochastic{iw, ia};
        plot_reach_heatmap(r.x_grid, r.y_grid, r.cone_mask, r.safe_mask, rp.cone_k);
        n_safe = sum(r.safe_mask(:) & r.cone_mask(:));
        n_cone = sum(r.cone_mask(:));
        title(sprintf('\\omega=%d, a_{max}=%.2f (%d/%d)', ...
            r.omega_deg, r.a_max, n_safe, n_cone), 'FontSize', 9);
        if iw < n_omega, xlabel(''); end
        if ia > 1, ylabel(''); end
    end
end
sgtitle(sprintf('Stochastic Safe Sets (\\alpha=%.2f, %.0f%% confidence)', ...
    rp.alpha, 100*(1-rp.alpha)), 'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(fig_grid, fullfile(fig_dir, 'stochastic_grid.png'), 'Resolution', 200);
close(fig_grid);
fprintf('Saved: stochastic_grid.png\n');

fprintf('\n====== Stochastic Reachability Complete (%.1f s) ======\n', stoch_time);


function plot_reach_heatmap(x_grid, y_grid, cone_mask, safe_mask, cone_k)
    hold on;
    ny = length(y_grid); nx = length(x_grid);
    img = ones(ny, nx, 3);
    red_col = [0.9, 0.25, 0.25]; blue_col = [0.3, 0.4, 0.85]; green_col = [0.2, 0.75, 0.2];
    for ch = 1:3
        layer = img(:,:,ch);
        layer(~cone_mask) = red_col(ch);
        layer(cone_mask & ~safe_mask) = blue_col(ch);
        layer(cone_mask & safe_mask) = green_col(ch);
        img(:,:,ch) = layer;
    end
    image(x_grid, y_grid, img); set(gca, 'YDir', 'normal');
    y_cone = linspace(0, max(y_grid)*1.1, 200);
    plot(cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(-cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]');
    xlim([-200 200]); ylim([-10 320]); grid on; set(gca, 'FontSize', 9);
end
