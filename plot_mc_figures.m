%PLOT_MC_FIGURES  Regenerate all MC figures from saved data.
%
%   Loads mc_sweep_all.mat and produces:
%     1) fig_mc_sweep_grid.png      - 5x4 grid, shared legend at bottom
%     2) fig_mc_omega_X.png (x5)    - fixed omega, overlaid a_max (dark->light)
%     3) fig_mc_amax_X.XX.png (x4)  - fixed a_max, overlaid omega (dark->light)
%     4) fig_mc_wX_aX.XX.png (x20)  - individual heatmaps
%
%   Run:  >> plot_mc_figures

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(this_dir, 'results');

%% ===== Load data =====
S = load(fullfile(results_dir, 'mc_sweep_all.mat'));
all_data       = S.all_data;
omega_vals_deg = S.omega_vals_deg;
amax_vals      = S.amax_vals;
ic_list        = S.ic_list;
cone_k         = S.cone_k;

n_omega = length(omega_vals_deg);
n_amax  = length(amax_vals);

%% ===== 1) Grid figure (5x4) with shared legend =====
fig_grid = figure('Position', [50 50 1400 1600]);

for iw = 1:n_omega
    for ia = 1:n_amax
        idx = (iw-1) * n_amax + ia;
        ax = subplot(n_omega, n_amax, idx);
        plot_mc_heatmap_ax(ax, all_data{iw, ia}, false);  % no legend
    end
end

% Add spacing between subplots
set(fig_grid, 'DefaultAxesLooseInset', [0.05, 0.08, 0.02, 0.02]);
for ii = 1:numel(fig_grid.Children)
    if isa(fig_grid.Children(ii), 'matlab.graphics.axis.Axes')
        fig_grid.Children(ii).Position(4) = fig_grid.Children(ii).Position(4) * 0.88;
    end
end

% Shared legend at the bottom
green_col = [0.2, 0.75, 0.2];
blue_col  = [0.3, 0.4, 0.85];
red_col   = [0.9, 0.25, 0.25];

% Create invisible axes for the legend
ax_leg = axes(fig_grid, 'Position', [0.1 0.01 0.8 0.03], 'Visible', 'off');
hold(ax_leg, 'on');
h1 = fill(ax_leg, NaN, NaN, green_col, 'EdgeColor', 'none');
h2 = fill(ax_leg, NaN, NaN, blue_col,  'EdgeColor', 'none');
h3 = fill(ax_leg, NaN, NaN, red_col,   'EdgeColor', 'none');
legend(ax_leg, [h1 h2 h3], ...
    {'All-time feasible', 'Fails after start', 'Outside LOS cone'}, ...
    'Orientation', 'horizontal', 'Location', 'south', 'FontSize', 12);

sgtitle('MC Feasibility Sweep: \omega_{target} vs a_{max,chaser}', ...
    'FontSize', 16, 'FontWeight', 'bold');

exportgraphics(fig_grid, fullfile(results_dir, 'fig_mc_sweep_grid.png'), ...
    'Resolution', 200);
close(fig_grid);
fprintf('Saved: fig_mc_sweep_grid.png\n');

%% ===== 2) Individual heatmaps (20 figures) =====
for iw = 1:n_omega
    for ia = 1:n_amax
        mc = all_data{iw, ia};
        fig = figure('Position', [100 100 700 550], 'Visible', 'off');
        ax = gca;
        plot_mc_heatmap_ax(ax, mc, true);  % with legend
        fname = sprintf('fig_mc_w%d_a%.2f.png', mc.omega_deg, mc.amax);
        exportgraphics(fig, fullfile(results_dir, fname), 'Resolution', 200);
        close(fig);
        fprintf('Saved: %s\n', fname);
    end
end

%% ===== 3) Overlay figures: fixed omega, varying a_max =====
% For each omega, overlay all a_max levels.
% Plot order: red background, then blue, then a_max from largest to smallest.
% Colors: lightest for a_max=0.2, darkest for a_max=0.02.
amax_greens = [0.75, 0.85, 0.75;   % a_max=0.2  (lightest)
               0.55, 0.70, 0.55;
               0.30, 0.55, 0.30;
               0.10, 0.35, 0.10];  % a_max=0.02 (darkest)

for iw = 1:n_omega
    omega_deg = omega_vals_deg(iw);
    fig = figure('Position', [100 100 750 600], 'Visible', 'off');
    ax = gca; hold(ax, 'on');

    % Red background (outside cone)
    plot_cone_bg(ax, cone_k, red_col);

    % Blue for initially-feasible-but-fails (union of all failures)
    % Draw the full cone as blue first, then overlay green regions
    y_cone = linspace(0, 350, 200);
    fill(ax, [cone_k*y_cone, fliplr(-cone_k*y_cone)], [y_cone, fliplr(y_cone)], ...
        blue_col, 'EdgeColor', 'none', 'FaceAlpha', 1.0, 'HandleVisibility', 'off');

    % Overlay feasible regions from largest (a_max=0.2) to smallest (a_max=0.02)
    % so larger regions are at the bottom
    h_layers = gobjects(n_amax, 1);
    for ia = 1:n_amax
        mc = all_data{iw, ia};
        [boundary_x, boundary_y] = get_feasible_boundary(mc);
        if ~isempty(boundary_x)
            h_layers(ia) = fill(ax, boundary_x, boundary_y, ...
                amax_greens(ia,:), 'EdgeColor', amax_greens(ia,:)*0.7, ...
                'LineWidth', 1.2, 'FaceAlpha', 0.9);
        else
            h_layers(ia) = fill(ax, NaN, NaN, amax_greens(ia,:), 'EdgeColor', 'none');
        end
    end

    % Cone boundary
    plot(ax,  cone_k*y_cone, y_cone, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(ax, -cone_k*y_cone, y_cone, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');

    % Target
    plot(ax, 0, 0, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.3 0.3], ...
        'HandleVisibility', 'off');

    % Legend
    leg_labels = cell(n_amax, 1);
    for ia = 1:n_amax
        n_ok = sum(all_data{iw, ia}.outcomes);
        leg_labels{ia} = sprintf('a_{max}=%.2f m/s^2 (%d/%d)', amax_vals(ia), n_ok, size(ic_list,1));
    end
    legend(ax, h_layers, leg_labels, 'Location', 'southeast', 'FontSize', 10);

    xlabel(ax, 'x_{TB} [m]', 'FontSize', 12);
    ylabel(ax, 'y_{TB} [m]', 'FontSize', 12);
    title(ax, sprintf('Feasible Region  (\\omega_{tgt} = %d deg/s)', omega_deg), ...
        'FontSize', 13);
    set(ax, 'FontSize', 11);
    xlim(ax, [-200 200]); ylim(ax, [-10 320]);
    grid(ax, 'on');

    fname = sprintf('fig_mc_omega_%d.png', omega_deg);
    exportgraphics(fig, fullfile(results_dir, fname), 'Resolution', 200);
    close(fig);
    fprintf('Saved: %s\n', fname);
end

%% ===== 4) Overlay figures: fixed a_max, varying omega =====
% For each a_max, overlay all omega values.
% omega=1 (largest area) at bottom lightest, omega=5 darkest on top.
omega_greens = [0.75, 0.88, 0.70;   % omega=1 (lightest)
                0.55, 0.75, 0.50;
                0.35, 0.60, 0.30;
                0.18, 0.45, 0.15;
                0.05, 0.30, 0.05];  % omega=5 (darkest)

for ia = 1:n_amax
    amax = amax_vals(ia);
    fig = figure('Position', [100 100 750 600], 'Visible', 'off');
    ax = gca; hold(ax, 'on');

    % Red background
    plot_cone_bg(ax, cone_k, red_col);

    % Blue cone interior
    y_cone = linspace(0, 350, 200);
    fill(ax, [cone_k*y_cone, fliplr(-cone_k*y_cone)], [y_cone, fliplr(y_cone)], ...
        blue_col, 'EdgeColor', 'none', 'FaceAlpha', 1.0, 'HandleVisibility', 'off');

    % Overlay from omega=1 (largest, lightest) to omega=5 (smallest, darkest)
    h_layers = gobjects(n_omega, 1);
    for iw = 1:n_omega
        mc = all_data{iw, ia};
        [boundary_x, boundary_y] = get_feasible_boundary(mc);
        if ~isempty(boundary_x)
            h_layers(iw) = fill(ax, boundary_x, boundary_y, ...
                omega_greens(iw,:), 'EdgeColor', omega_greens(iw,:)*0.7, ...
                'LineWidth', 1.2, 'FaceAlpha', 0.9);
        else
            h_layers(iw) = fill(ax, NaN, NaN, omega_greens(iw,:), 'EdgeColor', 'none');
        end
    end

    % Cone boundary
    plot(ax,  cone_k*y_cone, y_cone, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(ax, -cone_k*y_cone, y_cone, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');

    % Target
    plot(ax, 0, 0, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.3 0.3], ...
        'HandleVisibility', 'off');

    % Legend
    leg_labels = cell(n_omega, 1);
    for iw = 1:n_omega
        n_ok = sum(all_data{iw, ia}.outcomes);
        leg_labels{iw} = sprintf('\\omega=%d deg/s (%d/%d)', omega_vals_deg(iw), n_ok, size(ic_list,1));
    end
    legend(ax, h_layers, leg_labels, 'Location', 'southeast', 'FontSize', 10);

    xlabel(ax, 'x_{TB} [m]', 'FontSize', 12);
    ylabel(ax, 'y_{TB} [m]', 'FontSize', 12);
    title(ax, sprintf('Feasible Region  (a_{max} = %.2f m/s^2)', amax), ...
        'FontSize', 13);
    set(ax, 'FontSize', 11);
    xlim(ax, [-200 200]); ylim(ax, [-10 320]);
    grid(ax, 'on');

    fname = sprintf('fig_mc_amax_%.2f.png', amax);
    exportgraphics(fig, fullfile(results_dir, fname), 'Resolution', 200);
    close(fig);
    fprintf('Saved: %s\n', fname);
end

fprintf('\nAll MC figures regenerated from saved data.\n');


%% ========== Helper Functions ==========

function plot_mc_heatmap_ax(ax, mc, with_legend)
%PLOT_MC_HEATMAP_AX  3-color heatmap on given axes.
    if nargin < 3, with_legend = true; end
    axes(ax); hold on;

    ic_list   = mc.ic_list;
    outcomes  = mc.outcomes;
    cone_k    = mc.cone_k;
    omega_deg = mc.omega_deg;
    amax      = mc.amax;
    n_mc      = size(ic_list, 1);

    x_fine = linspace(-200, 200, 300);
    y_fine = linspace(0, 320, 300);
    [Xq, Yq] = meshgrid(x_fine, y_fine);
    cone_mask = abs(Xq) <= cone_k * Yq & Yq >= 1;

    img = ones(length(y_fine), length(x_fine), 3);

    red_col   = [0.9, 0.25, 0.25];
    green_col = [0.2, 0.75, 0.2];
    blue_col  = [0.3, 0.4, 0.85];

    for ch = 1:3
        layer = img(:,:,ch);
        layer(~cone_mask) = red_col(ch);
        img(:,:,ch) = layer;
    end

    F = scatteredInterpolant(ic_list(:,1), ic_list(:,2), outcomes, ...
        'natural', 'nearest');
    Zq = F(Xq, Yq);
    Zq_bin = Zq >= 0.5;

    for ch = 1:3
        layer = img(:,:,ch);
        layer(cone_mask &  Zq_bin) = green_col(ch);
        layer(cone_mask & ~Zq_bin) = blue_col(ch);
        img(:,:,ch) = layer;
    end

    image(ax, x_fine, y_fine, img);
    set(ax, 'YDir', 'normal');

    y_cone = linspace(0, 350, 200);
    plot(ax,  cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(ax, -cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    plot(ax, 0, 0, 'ks', 'MarkerSize', 6, 'MarkerFaceColor', [0.3 0.3 0.3], ...
        'HandleVisibility', 'off');

    if with_legend
        h1 = fill(ax, NaN, NaN, green_col, 'EdgeColor', 'none');
        h2 = fill(ax, NaN, NaN, blue_col,  'EdgeColor', 'none');
        h3 = fill(ax, NaN, NaN, red_col,   'EdgeColor', 'none');
        legend(ax, [h1 h2 h3], {'All-time feasible', 'Fails after start', 'Outside LOS'}, ...
            'Location', 'southeast', 'FontSize', 9);
    end

    xlabel(ax, 'x_{TB} [m]', 'FontSize', 10);
    ylabel(ax, 'y_{TB} [m]', 'FontSize', 10);
    title(ax, sprintf('\\omega_{tgt}=%d deg/s,  a_{max}=%.2f m/s^2  (%d/%d)', ...
        omega_deg, amax, sum(outcomes), n_mc), 'FontSize', 10);
    set(ax, 'FontSize', 9);
    xlim(ax, [-200 200]);
    ylim(ax, [-10 320]);
    grid(ax, 'on');
end


function plot_cone_bg(ax, cone_k, red_col)
%PLOT_CONE_BG  Fill area outside cone with red.
    % Fill entire plot area with red, then we'll overlay cone interior
    fill(ax, [-200 200 200 -200], [-10 -10 320 320], ...
        red_col, 'EdgeColor', 'none', 'FaceAlpha', 1.0, 'HandleVisibility', 'off');
end


function [bx, by] = get_feasible_boundary(mc)
%GET_FEASIBLE_BOUNDARY  Compute boundary polygon of feasible region.
    ic_list  = mc.ic_list;
    outcomes = mc.outcomes;
    cone_k   = mc.cone_k;

    if sum(outcomes) == 0
        bx = []; by = [];
        return;
    end

    % Interpolate onto fine grid
    x_fine = linspace(-200, 200, 400);
    y_fine = linspace(0, 320, 400);
    [Xq, Yq] = meshgrid(x_fine, y_fine);

    F = scatteredInterpolant(ic_list(:,1), ic_list(:,2), outcomes, ...
        'natural', 'nearest');
    Zq = F(Xq, Yq);

    % Also apply cone mask
    cone_mask = abs(Xq) <= cone_k * Yq & Yq >= 1;
    Zq(~cone_mask) = 0;

    % Get contour at 0.5 level
    C = contourc(x_fine, y_fine, Zq, [0.5, 0.5]);

    if isempty(C)
        bx = []; by = [];
        return;
    end

    % Extract the longest contour (may have multiple segments)
    idx = 1;
    best_bx = []; best_by = [];
    best_len = 0;
    while idx < size(C, 2)
        n_pts = C(2, idx);
        cx = C(1, idx+1:idx+n_pts);
        cy = C(2, idx+1:idx+n_pts);
        if n_pts > best_len
            best_bx = cx;
            best_by = cy;
            best_len = n_pts;
        end
        idx = idx + n_pts + 1;
    end

    % Close the polygon
    bx = [best_bx, best_bx(1)];
    by = [best_by, best_by(1)];
end
