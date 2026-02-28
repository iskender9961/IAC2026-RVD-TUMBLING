%RUNSINGLESCENARIO  Entry point: Single scenario comparison overlay.
%
%   Runs one selected (omega, a_max) case and overlays:
%     - MC safe region (if available)
%     - Nominal certified region
%     - Stochastic certified region
%     - Robust certified region
%
%   Edit the parameters below to select the case.
%
%   Usage:  >> RunSingleScenario

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(this_dir));

%% ===== SELECT CASE =====
omega_deg = 2;      % target tumble rate [deg/s]
a_max     = 0.10;   % chaser max acceleration [m/s^2]

fprintf('Single scenario: omega=%d deg/s, a_max=%.2f m/s^2\n', omega_deg, a_max);

%% ===== Compute all three reachability sets =====
rp = reachability_params();

fprintf('\n--- Nominal ---\n');
t1 = tic;
nom = forward_reachable_set(omega_deg, a_max, rp);
t_nom = toc(t1);
fprintf('  %d safe / %d cone (%.1f%%), %.1f s\n', ...
    sum(nom.safe_mask(:) & nom.cone_mask(:)), ...
    sum(nom.cone_mask(:)), ...
    100*sum(nom.safe_mask(:) & nom.cone_mask(:))/max(sum(nom.cone_mask(:)),1), t_nom);

fprintf('\n--- Stochastic ---\n');
t2 = tic;
stoch = stochastic_safe_set(omega_deg, a_max, rp);
t_stoch = toc(t2);
fprintf('  %d safe / %d cone (%.1f%%), %.1f s\n', ...
    sum(stoch.safe_mask(:) & stoch.cone_mask(:)), ...
    sum(stoch.cone_mask(:)), ...
    100*sum(stoch.safe_mask(:) & stoch.cone_mask(:))/max(sum(stoch.cone_mask(:)),1), t_stoch);

fprintf('\n--- Robust ---\n');
t3 = tic;
rob = robust_safe_set(omega_deg, a_max, rp);
t_rob = toc(t3);
fprintf('  %d safe / %d cone (%.1f%%), %.1f s\n', ...
    sum(rob.safe_mask(:) & rob.cone_mask(:)), ...
    sum(rob.cone_mask(:)), ...
    100*sum(rob.safe_mask(:) & rob.cone_mask(:))/max(sum(rob.cone_mask(:)),1), t_rob);

%% ===== Overlay plot =====
fig_dir = fullfile(this_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Colors
mc_col    = [0.7, 0.85, 0.7];
nom_col   = [0.2, 0.75, 0.2];
stoch_col = [0.3, 0.5, 0.85];
rob_col   = [0.6, 0.2, 0.6];
red_col   = [0.9, 0.25, 0.25];
cone_bg   = [0.92, 0.92, 0.92];

x_g = nom.x_grid;
y_g = nom.y_grid;
ny = length(y_g); nx = length(x_g);
cone = nom.cone_mask;

img = ones(ny, nx, 3);
for ch = 1:3
    layer = img(:,:,ch);
    layer(~cone) = red_col(ch);
    layer(cone) = cone_bg(ch);
    img(:,:,ch) = layer;
end

% Try to load MC data
mc_dir = fullfile(this_dir, 'results');
mc_file = fullfile(mc_dir, 'mc_sweep_all.mat');
has_mc = false;
if exist(mc_file, 'file')
    mc_all = load(mc_file);
    % Find matching combo
    iw = find(mc_all.omega_vals_deg == omega_deg);
    ia = find(abs(mc_all.amax_vals - a_max) < 1e-6);
    if ~isempty(iw) && ~isempty(ia) && ~isempty(mc_all.all_data{iw, ia})
        mc = mc_all.all_data{iw, ia};
        F_mc = scatteredInterpolant(mc.ic_list(:,1), mc.ic_list(:,2), ...
            double(mc.outcomes), 'natural', 'nearest');
        [Xq, Yq] = meshgrid(x_g, y_g);
        mc_safe = F_mc(Xq, Yq) >= 0.5 & cone;
        for ch = 1:3
            layer = img(:,:,ch);
            layer(mc_safe) = mc_col(ch);
            img(:,:,ch) = layer;
        end
        has_mc = true;
        fprintf('\nMC data loaded: %d/%d feasible\n', sum(mc.outcomes), length(mc.outcomes));
    end
end

% Layer analytical results
m_nom = nom.safe_mask & cone;
m_stoch = stoch.safe_mask & cone;
m_rob = rob.safe_mask & cone;

for ch = 1:3
    layer = img(:,:,ch);
    layer(m_nom) = nom_col(ch);
    layer(m_stoch) = stoch_col(ch);
    layer(m_rob) = rob_col(ch);
    img(:,:,ch) = layer;
end

fig = figure('Position', [100 100 800 600]);
hold on;
image(x_g, y_g, img);
set(gca, 'YDir', 'normal');

y_cone = linspace(0, max(y_g)*1.1, 200);
plot(rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(-rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(0, 0, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.3 0.3], ...
    'HandleVisibility', 'off');

patches = [];
labels = {};
if has_mc
    patches(end+1) = fill(NaN, NaN, mc_col, 'EdgeColor', 'none');
    labels{end+1} = 'MC safe (empirical)';
end
patches(end+1) = fill(NaN, NaN, nom_col, 'EdgeColor', 'none');
labels{end+1} = 'Nominal (certified)';
patches(end+1) = fill(NaN, NaN, stoch_col, 'EdgeColor', 'none');
labels{end+1} = sprintf('Stochastic (%.0f%%)', 100*(1-rp.alpha));
patches(end+1) = fill(NaN, NaN, rob_col, 'EdgeColor', 'none');
labels{end+1} = 'Robust (worst-case)';
patches(end+1) = fill(NaN, NaN, red_col, 'EdgeColor', 'none');
labels{end+1} = 'Outside LOS cone';
legend(patches, labels, 'Location', 'southeast', 'FontSize', 9);

xlabel('x_{TB} [m]', 'FontSize', 12);
ylabel('y_{TB} [m]', 'FontSize', 12);
title(sprintf('Single Scenario: \\omega_{tgt}=%d deg/s, a_{max}=%.2f m/s^2', ...
    omega_deg, a_max), 'FontSize', 13);
xlim([-200 200]); ylim([-10 320]); grid on;
set(gca, 'FontSize', 11);

fname = sprintf('single_w%d_a%.2f.png', omega_deg, a_max);
exportgraphics(fig, fullfile(fig_dir, fname), 'Resolution', 200);
fprintf('\nSaved: figures/%s\n', fname);

%% ===== Nesting verification =====
fprintf('\n=== Nesting Check ===\n');
v1 = sum(m_rob(:) & ~m_stoch(:));
v2 = sum(m_stoch(:) & ~m_nom(:));
v3 = 0;
if has_mc, v3 = sum(m_nom(:) & ~mc_safe(:)); end
fprintf('  Rob ⊄ Stoch: %d violations\n', v1);
fprintf('  Stoch ⊄ Nom: %d violations\n', v2);
if has_mc, fprintf('  Nom ⊄ MC:    %d violations\n', v3); end
if v1+v2+v3 == 0
    fprintf('  [OK] Hierarchy verified.\n');
end

fprintf('\nTimings: nominal=%.1f s, stochastic=%.1f s, robust=%.1f s\n', ...
    t_nom, t_stoch, t_rob);
