%RUN_COMPARISON  Comparative analysis: overlay all reachability sets + MC.
%
%   Loads saved results from:
%     - data/nominal_reachability.mat
%     - data/stochastic_reachability.mat
%     - data/robust_reachability.mat
%     - results/mc_sweep_all.mat  (MC data from existing pipeline)
%
%   Produces:
%     1. Nested overlay plots showing X_robust ⊆ X_stochastic ⊆ X_nominal ⊆ X_MC
%     2. Area ratio tables
%     3. Inclusion violation checks
%     4. Compute-time comparison
%
%   Run:  >> run_comparison

if ~exist('CALLED_FROM_MASTER','var'), clear; close all; clc; end

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
root_dir = fileparts(this_dir);
if isempty(root_dir), root_dir = pwd; end
addpath(fullfile(root_dir, 'reachability_common'));
addpath(root_dir);

data_dir = fullfile(root_dir, 'data');
fig_dir  = fullfile(root_dir, 'figures');
mc_dir   = fullfile(root_dir, 'results');

%% ===== Load all data =====
fprintf('Loading data...\n');

% Nominal
nom_file = fullfile(data_dir, 'nominal_reachability.mat');
if ~exist(nom_file, 'file')
    error('Run run_nominal_reachability.m first');
end
nom_data = load(nom_file);

% Stochastic
stoch_file = fullfile(data_dir, 'stochastic_reachability.mat');
if ~exist(stoch_file, 'file')
    error('Run run_stochastic_reachability.m first');
end
stoch_data = load(stoch_file);

% Robust
rob_file = fullfile(data_dir, 'robust_reachability.mat');
if ~exist(rob_file, 'file')
    error('Run run_robust_reachability.m first');
end
rob_data = load(rob_file);

% Monte Carlo
mc_file = fullfile(mc_dir, 'mc_sweep_all.mat');
has_mc = exist(mc_file, 'file');
if has_mc
    mc_data = load(mc_file);
    fprintf('MC data loaded.\n');
else
    fprintf('No MC data found — comparison will proceed without MC overlay.\n');
end

omega_vals = nom_data.omega_vals;
amax_vals  = nom_data.amax_vals;
n_omega = length(omega_vals);
n_amax  = length(amax_vals);
rp = nom_data.rp;

%% ===== Inclusion Verification =====
fprintf('\n====== Inclusion Hierarchy Verification ======\n');
fprintf('Required: X_robust ⊆ X_stochastic ⊆ X_nominal ⊆ X_MC_safe\n\n');

area_table = zeros(n_omega, n_amax, 4);  % [nom, stoch, rob, mc]
violation_count = zeros(n_omega, n_amax, 3);  % [rob⊆stoch, stoch⊆nom, nom⊆mc]

for iw = 1:n_omega
    for ia = 1:n_amax
        nom   = nom_data.all_forward{iw, ia};
        stoch = stoch_data.all_stochastic{iw, ia};
        rob   = rob_data.all_robust{iw, ia};

        cone = nom.cone_mask;
        n_cone = sum(cone(:));

        % Safe masks
        m_nom   = nom.safe_mask & cone;
        m_stoch = stoch.safe_mask & cone;
        m_rob   = rob.safe_mask & cone;

        % Area fractions
        area_table(iw, ia, 1) = sum(m_nom(:))   / max(n_cone, 1);
        area_table(iw, ia, 2) = sum(m_stoch(:)) / max(n_cone, 1);
        area_table(iw, ia, 3) = sum(m_rob(:))   / max(n_cone, 1);

        % Check inclusions
        % robust ⊆ stochastic
        v1 = sum(m_rob(:) & ~m_stoch(:));
        % stochastic ⊆ nominal
        v2 = sum(m_stoch(:) & ~m_nom(:));
        % nominal ⊆ MC (if available)
        v3 = 0;

        violation_count(iw, ia, :) = [v1, v2, v3];

        if v1 > 0 || v2 > 0
            fprintf('  [VIOLATION] omega=%d, a_max=%.2f: rob⊄stoch=%d, stoch⊄nom=%d\n', ...
                omega_vals(iw), amax_vals(ia), v1, v2);
        end

        if has_mc
            % Build MC safe mask on the same grid
            mc = mc_data.all_data{iw, ia};
            if ~isempty(mc)
                % Interpolate MC outcomes onto reachability grid
                F_mc = scatteredInterpolant(mc.ic_list(:,1), mc.ic_list(:,2), ...
                    double(mc.outcomes), 'natural', 'nearest');
                [Xq, Yq] = meshgrid(nom.x_grid, nom.y_grid);
                mc_interp = F_mc(Xq, Yq) >= 0.5;
                mc_safe = mc_interp & cone;
                area_table(iw, ia, 4) = sum(mc_safe(:)) / max(n_cone, 1);

                % nominal ⊆ MC?
                v3 = sum(m_nom(:) & ~mc_safe(:));
                violation_count(iw, ia, 3) = v3;
                if v3 > 0
                    fprintf('  [WARN] omega=%d, a_max=%.2f: nom⊄MC=%d points\n', ...
                        omega_vals(iw), amax_vals(ia), v3);
                end
            end
        end
    end
end

% Summary
total_violations = sum(violation_count(:));
if total_violations == 0
    fprintf('\n  [OK] All inclusion relations verified: X_robust ⊆ X_stochastic ⊆ X_nominal\n');
else
    fprintf('\n  [WARN] %d total inclusion violations detected\n', total_violations);
end

%% ===== Print Area Ratio Table =====
fprintf('\n====== Area Ratios (fraction of cone area) ======\n');
fprintf('%-12s %-10s  %8s %8s %8s %8s\n', 'omega', 'a_max', 'Nominal', 'Stoch.', 'Robust', 'MC');
fprintf('%s\n', repmat('-', 1, 60));
for iw = 1:n_omega
    for ia = 1:n_amax
        fprintf('%4d deg/s   %5.2f     %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n', ...
            omega_vals(iw), amax_vals(ia), ...
            100*area_table(iw,ia,1), 100*area_table(iw,ia,2), ...
            100*area_table(iw,ia,3), 100*area_table(iw,ia,4));
    end
end

%% ===== Compute Time Comparison =====
fprintf('\n====== Compute Time ======\n');
fprintf('  Nominal forward:      %.1f s\n', nom_data.forward_time);
fprintf('  Nominal backward:     %.1f s\n', nom_data.backward_time);
fprintf('  Stochastic:           %.1f s\n', stoch_data.stoch_time);
fprintf('  Robust:               %.1f s\n', rob_data.robust_time);
fprintf('  Total analytical:     %.1f s\n', ...
    nom_data.forward_time + nom_data.backward_time + ...
    stoch_data.stoch_time + rob_data.robust_time);

%% ===== Nested Overlay Figures =====
fprintf('\n--- Generating Overlay Figures ---\n');

% Colors for nested regions (from outermost to innermost)
mc_col    = [0.7, 0.85, 0.7];    % light green (MC, outermost)
nom_col   = [0.2, 0.75, 0.2];    % green (nominal)
stoch_col = [0.3, 0.5, 0.85];    % blue (stochastic)
rob_col   = [0.6, 0.2, 0.6];     % purple (robust, innermost)
red_col   = [0.9, 0.25, 0.25];   % red (outside cone)
cone_bg   = [0.92, 0.92, 0.92];  % light gray (cone, no method covers)

for iw = 1:n_omega
    for ia = 1:n_amax
        nom   = nom_data.all_forward{iw, ia};
        stoch = stoch_data.all_stochastic{iw, ia};
        rob   = rob_data.all_robust{iw, ia};
        cone  = nom.cone_mask;

        x_g = nom.x_grid;
        y_g = nom.y_grid;
        ny = length(y_g);
        nx = length(x_g);

        % Build layered image
        img = ones(ny, nx, 3);

        % Layer 0: outside cone = red
        for ch = 1:3
            layer = img(:,:,ch);
            layer(~cone) = red_col(ch);
            img(:,:,ch) = layer;
        end

        % Layer 1: inside cone (no method) = gray
        for ch = 1:3
            layer = img(:,:,ch);
            layer(cone) = cone_bg(ch);
            img(:,:,ch) = layer;
        end

        % Layer 2: MC safe = light green (if available)
        if has_mc && ~isempty(mc_data.all_data{iw, ia})
            mc = mc_data.all_data{iw, ia};
            F_mc = scatteredInterpolant(mc.ic_list(:,1), mc.ic_list(:,2), ...
                double(mc.outcomes), 'natural', 'nearest');
            [Xq, Yq] = meshgrid(x_g, y_g);
            mc_safe = F_mc(Xq, Yq) >= 0.5 & cone;
            for ch = 1:3
                layer = img(:,:,ch);
                layer(mc_safe) = mc_col(ch);
                img(:,:,ch) = layer;
            end
        end

        % Layer 3: Nominal = green
        m_nom = nom.safe_mask & cone;
        for ch = 1:3
            layer = img(:,:,ch);
            layer(m_nom) = nom_col(ch);
            img(:,:,ch) = layer;
        end

        % Layer 4: Stochastic = blue
        m_stoch = stoch.safe_mask & cone;
        for ch = 1:3
            layer = img(:,:,ch);
            layer(m_stoch) = stoch_col(ch);
            img(:,:,ch) = layer;
        end

        % Layer 5: Robust = purple
        m_rob = rob.safe_mask & cone;
        for ch = 1:3
            layer = img(:,:,ch);
            layer(m_rob) = rob_col(ch);
            img(:,:,ch) = layer;
        end

        % Plot
        fig = figure('Position', [100 100 700 550], 'Visible', 'off');
        hold on;
        image(x_g, y_g, img);
        set(gca, 'YDir', 'normal');

        % Cone boundary
        y_cone = linspace(0, max(y_g)*1.1, 200);
        plot(rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(-rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

        % Target marker
        plot(0, 0, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.3 0.3], ...
            'HandleVisibility', 'off');

        % Legend patches
        patches = [];
        labels = {};
        if has_mc && ~isempty(mc_data.all_data{iw, ia})
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

        legend(patches, labels, 'Location', 'southeast', 'FontSize', 8);

        xlabel('x_{TB} [m]', 'FontSize', 11);
        ylabel('y_{TB} [m]', 'FontSize', 11);
        title(sprintf('Nested Feasibility Regions\n\\omega_{tgt}=%d deg/s, a_{max}=%.2f m/s^2', ...
            omega_vals(iw), amax_vals(ia)), 'FontSize', 11);
        xlim([-200 200]); ylim([-10 320]); grid on;
        set(gca, 'FontSize', 10);

        fname = sprintf('compare_w%d_a%.2f.png', omega_vals(iw), amax_vals(ia));
        exportgraphics(fig, fullfile(fig_dir, fname), 'Resolution', 200);
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end
end

%% ===== Combined Grid: Nested Overlay =====
fig_grid = figure('Position', [50 50 1400 1500], 'Visible', 'off');
for iw = 1:n_omega
    for ia = 1:n_amax
        idx = (iw-1) * n_amax + ia;
        ax = subplot(n_omega, n_amax, idx);
        hold(ax, 'on');

        nom   = nom_data.all_forward{iw, ia};
        stoch = stoch_data.all_stochastic{iw, ia};
        rob   = rob_data.all_robust{iw, ia};
        cone  = nom.cone_mask;
        x_g = nom.x_grid; y_g = nom.y_grid;
        ny = length(y_g); nx = length(x_g);

        img = ones(ny, nx, 3);
        for ch = 1:3
            layer = img(:,:,ch);
            layer(~cone) = red_col(ch);
            layer(cone) = cone_bg(ch);
            img(:,:,ch) = layer;
        end
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

        image(ax, x_g, y_g, img);
        set(ax, 'YDir', 'normal');
        y_cone = linspace(0, max(y_g)*1.1, 200);
        plot(ax, rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, -rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        title(ax, sprintf('\\omega=%d, a=%.2f', omega_vals(iw), amax_vals(ia)), 'FontSize', 9);
        xlim(ax, [-200 200]); ylim(ax, [-10 320]);
        if iw == n_omega, xlabel(ax, 'x_{TB} [m]'); else, xlabel(ax, ''); end
        if ia == 1, ylabel(ax, 'y_{TB} [m]'); else, ylabel(ax, ''); end
        set(ax, 'FontSize', 8);
    end
end

% Shared legend at bottom
h1 = fill(NaN, NaN, nom_col, 'EdgeColor', 'none');
h2 = fill(NaN, NaN, stoch_col, 'EdgeColor', 'none');
h3 = fill(NaN, NaN, rob_col, 'EdgeColor', 'none');
h4 = fill(NaN, NaN, red_col, 'EdgeColor', 'none');
lgd = legend([h1 h2 h3 h4], {'Nominal', ...
    sprintf('Stochastic (%.0f%%)', 100*(1-rp.alpha)), ...
    'Robust', 'Outside cone'}, ...
    'Orientation', 'horizontal', 'FontSize', 10);
lgd.Position = [0.25, 0.01, 0.5, 0.02];

sgtitle('Nested Feasibility Hierarchy: X_{robust} \subseteq X_{stoch} \subseteq X_{nominal}', ...
    'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(fig_grid, fullfile(fig_dir, 'comparison_grid.png'), 'Resolution', 200);
close(fig_grid);
fprintf('  Saved: comparison_grid.png\n');

%% ===== Save comparison results =====
comparison.area_table = area_table;
comparison.violation_count = violation_count;
comparison.omega_vals = omega_vals;
comparison.amax_vals = amax_vals;
comparison.compute_times = struct(...
    'forward', nom_data.forward_time, ...
    'backward', nom_data.backward_time, ...
    'stochastic', stoch_data.stoch_time, ...
    'robust', rob_data.robust_time);

save(fullfile(data_dir, 'comparison_results.mat'), 'comparison');
fprintf('\nSaved: data/comparison_results.mat\n');

fprintf('\n====== Comparison Complete ======\n');
