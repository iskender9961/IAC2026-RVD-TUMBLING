%PLOT_SCALING_COLLAPSE  Safe fraction vs a_max/omega_t^2 collapse plot.
%
%  Loads comparison_results.mat and plots all 20 (omega,amax) combos
%  as safe fraction vs the universal scaling parameter a_max/omega_t^2.
%
%  Run:  >> plot_scaling_collapse

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end

data_dir = fullfile(this_dir, 'data');
results_dir = fullfile(this_dir, 'results');

%% Load data
comp = load(fullfile(data_dir, 'comparison_results.mat'));
c = comp.comparison;

omega_vals = c.omega_vals;   % [1, 2, 3, 4, 5] deg/s
amax_vals  = c.amax_vals;    % [0.02, 0.05, 0.10, 0.20] m/s^2
area_table = c.area_table;   % n_omega x n_amax x 4 [nom, stoch, rob, mc]

n_omega = length(omega_vals);
n_amax  = length(amax_vals);

%% Compute scaling parameter for each combo
% Convert omega from deg/s to rad/s
scaling = zeros(n_omega, n_amax);
safe_nom   = zeros(n_omega, n_amax);
safe_stoch = zeros(n_omega, n_amax);
safe_rob   = zeros(n_omega, n_amax);
safe_mc    = zeros(n_omega, n_amax);

for iw = 1:n_omega
    for ia = 1:n_amax
        omega_rad = omega_vals(iw) * pi / 180;
        scaling(iw, ia) = amax_vals(ia) / omega_rad^2;
        safe_nom(iw, ia)   = 100 * area_table(iw, ia, 1);
        safe_stoch(iw, ia) = 100 * area_table(iw, ia, 2);
        safe_rob(iw, ia)   = 100 * area_table(iw, ia, 3);
        safe_mc(iw, ia)    = 100 * area_table(iw, ia, 4);
    end
end

%% Plot
fig = figure('Position', [100 100 800 500]);

% Flatten arrays
xi = scaling(:);
[xi_sorted, sort_idx] = sort(xi);

markers = {'o', 's', 'd', '^'};
colors  = {[0.2 0.75 0.2], [0.3 0.5 0.85], [0.6 0.2 0.6], [0.7 0.85 0.7]};
labels  = {'Nominal', 'Stochastic (95%)', 'Robust', 'Monte Carlo'};

hold on; grid on;

% Plot each method
datasets = {safe_nom, safe_stoch, safe_rob, safe_mc};
for m = 1:4
    data = datasets{m};
    y = data(:);
    y_sorted = y(sort_idx);
    plot(xi_sorted, y_sorted, '-', 'Color', colors{m}, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    scatter(xi(:), data(:), 60, 'MarkerEdgeColor', colors{m}, ...
        'MarkerFaceColor', colors{m}, 'MarkerFaceAlpha', 0.6, ...
        'Marker', markers{m}, 'DisplayName', labels{m});
end

xlabel('$a_{\max} / \omega_t^2$ [m]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Safe Fraction [%]', 'Interpreter', 'none', 'FontSize', 13);
title('Feasibility Collapse: All Methods vs.\ Universal Scaling Parameter', ...
    'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10);
set(gca, 'FontSize', 11);
xlim([0, max(xi)*1.05]);
ylim([0, 100]);

% Add annotation for r_sync = 2 * scaling
ax2 = axes('Position', get(gca, 'Position'), 'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', 'Color', 'none', 'YTick', []);
xlim(ax2, [0, max(xi)*1.05] * 2);
xlabel(ax2, '$r_{\mathrm{sync}} = 2 a_{\max}/\omega_t^2$ [m]', ...
    'Interpreter', 'latex', 'FontSize', 11);

exportgraphics(fig, fullfile(results_dir, 'fig_scaling_collapse.png'), ...
    'Resolution', 300);
fprintf('Saved: results/fig_scaling_collapse.png\n');
