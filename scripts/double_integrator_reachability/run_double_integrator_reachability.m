%RUN_DOUBLE_INTEGRATOR_REACHABILITY
%  Grid-based forward reachable set for a 2D constrained double integrator.
%
%  Approach: discretise the state space on a fixed grid. At each step,
%  a grid point x' is reachable if there exists a predecessor x in the
%  current reachable set and a control u such that x' = Ax + Bu (+ w).
%  For the robust case, require this for ALL disturbance vertices.
%
%  Generates: figures/fig_double_integrator_reachability.png
%
%  Run:  >> run_double_integrator_reachability

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
root_dir = fullfile(this_dir, '..', '..');
fig_dir  = fullfile(root_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% System definition
A = [1 1; 0 1];
B = [0; 1];
u_min = -1;
u_max =  1;
N = 5;
numU = 10;

% State constraints
xmin = -5; xmax = 5;
ymin = -5; ymax = 5;

% Disturbance bounds
dist_bound = [0.2; 0.2];
dist_vertices = [
     dist_bound(1),  dist_bound(2);
     dist_bound(1), -dist_bound(2);
    -dist_bound(1),  dist_bound(2);
    -dist_bound(1), -dist_bound(2);
];
n_dist = size(dist_vertices, 1);

%% Evaluation grid
n_grid = 200;
x1_grid = linspace(xmin, xmax, n_grid);
x2_grid = linspace(ymin, ymax, n_grid);
dx1 = x1_grid(2) - x1_grid(1);
dx2 = x2_grid(2) - x2_grid(1);
[X1G, X2G] = meshgrid(x1_grid, x2_grid);

% Control samples
u_samples = linspace(u_min, u_max, numU);

%% Helper: snap to grid indices
snap = @(val, grid_min, d) round((val - grid_min) / d) + 1;

%% Initial reachable set: small box around origin
nom_reach = false(n_grid, n_grid);
rob_reach = false(n_grid, n_grid);

% Initial set: |x1| <= 0.5, |x2| <= 0.5
for i1 = 1:n_grid
    for i2 = 1:n_grid
        if abs(x1_grid(i1)) <= 0.5 && abs(x2_grid(i2)) <= 0.5
            nom_reach(i2, i1) = true;
            rob_reach(i2, i1) = true;
        end
    end
end

fprintf('Initial set: %d grid points\n', sum(nom_reach(:)));

%% Forward reachable set — Nominal
fprintf('Computing nominal forward reachable set...\n');
nom_current = nom_reach;

for k = 1:N
    nom_next = false(n_grid, n_grid);
    % For each currently reachable state
    [r_idx, c_idx] = find(nom_current);
    for pp = 1:length(r_idx)
        x = [x1_grid(c_idx(pp)); x2_grid(r_idx(pp))];
        for ju = 1:numU
            x_next = A * x + B * u_samples(ju);
            % Check constraints
            if x_next(1) < xmin || x_next(1) > xmax || ...
               x_next(2) < ymin || x_next(2) > ymax
                continue;
            end
            % Snap to grid
            j1 = snap(x_next(1), xmin, dx1);
            j2 = snap(x_next(2), ymin, dx2);
            if j1 >= 1 && j1 <= n_grid && j2 >= 1 && j2 <= n_grid
                nom_next(j2, j1) = true;
            end
        end
    end
    nom_current = nom_next;
    fprintf('  Step %d: %d reachable points (nominal)\n', k, sum(nom_current(:)));
end

nom_final = nom_current;

%% Forward reachable set — Robust
fprintf('Computing robust forward reachable set...\n');
rob_current = rob_reach;

for k = 1:N
    rob_next = false(n_grid, n_grid);
    [r_idx, c_idx] = find(rob_current);
    for pp = 1:length(r_idx)
        x = [x1_grid(c_idx(pp)); x2_grid(r_idx(pp))];
        for ju = 1:numU
            % Check if ALL disturbance vertices keep state in bounds
            all_ok = true;
            for jw = 1:n_dist
                x_next_w = A * x + B * u_samples(ju) + dist_vertices(jw,:)';
                if x_next_w(1) < xmin || x_next_w(1) > xmax || ...
                   x_next_w(2) < ymin || x_next_w(2) > ymax
                    all_ok = false;
                    break;
                end
            end
            if all_ok
                % Use nominal next state for grid mapping
                x_next = A * x + B * u_samples(ju);
                j1 = snap(x_next(1), xmin, dx1);
                j2 = snap(x_next(2), ymin, dx2);
                if j1 >= 1 && j1 <= n_grid && j2 >= 1 && j2 <= n_grid
                    rob_next(j2, j1) = true;
                end
            end
        end
    end
    rob_current = rob_next;
    fprintf('  Step %d: %d reachable points (robust)\n', k, sum(rob_current(:)));
end

rob_final = rob_current;

%% Generate figure with 3 panels
fprintf('Generating figure...\n');

fig = figure('Position', [50 50 1500 450]);

% Colors
nom_col = [0.2 0.75 0.2];
rob_col = [0.6 0.2 0.6];
box_col = [0.5 0.5 0.5];

% ---- Panel 1: Nominal ----
subplot(1,3,1);
hold on; grid on;
imagesc(x1_grid, x2_grid, double(nom_final));
colormap(gca, [1 1 1; nom_col]);
set(gca, 'YDir', 'normal');
rectangle('Position', [xmin, ymin, xmax-xmin, ymax-ymin], ...
    'EdgeColor', box_col, 'LineWidth', 1.5, 'LineStyle', '--');
rectangle('Position', [-0.5, -0.5, 1, 1], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
xlabel('$x_1$ (position)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$ (velocity)', 'Interpreter', 'latex', 'FontSize', 12);
title('(a) Nominal Reachable Set', 'FontSize', 13);
xlim([xmin-0.3, xmax+0.3]); ylim([ymin-0.3, ymax+0.3]);
set(gca, 'FontSize', 11);

% ---- Panel 2: Robust ----
subplot(1,3,2);
hold on; grid on;
imagesc(x1_grid, x2_grid, double(rob_final));
colormap(gca, [1 1 1; rob_col]);
set(gca, 'YDir', 'normal');
rectangle('Position', [xmin, ymin, xmax-xmin, ymax-ymin], ...
    'EdgeColor', box_col, 'LineWidth', 1.5, 'LineStyle', '--');
rectangle('Position', [-0.5, -0.5, 1, 1], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
xlabel('$x_1$ (position)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$ (velocity)', 'Interpreter', 'latex', 'FontSize', 12);
title('(b) Robust Reachable Set', 'FontSize', 13);
xlim([xmin-0.3, xmax+0.3]); ylim([ymin-0.3, ymax+0.3]);
set(gca, 'FontSize', 11);

% ---- Panel 3: Overlay ----
subplot(1,3,3);
hold on; grid on;
% Build RGB overlay image
img = ones(n_grid, n_grid, 3);
for ch = 1:3
    layer = img(:,:,ch);
    layer(nom_final) = nom_col(ch);
    img(:,:,ch) = layer;
end
for ch = 1:3
    layer = img(:,:,ch);
    layer(rob_final) = rob_col(ch);
    img(:,:,ch) = layer;
end
image(x1_grid, x2_grid, img);
set(gca, 'YDir', 'normal');
rectangle('Position', [xmin, ymin, xmax-xmin, ymax-ymin], ...
    'EdgeColor', box_col, 'LineWidth', 1.5, 'LineStyle', '--');
rectangle('Position', [-0.5, -0.5, 1, 1], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
% Legend patches
h1 = fill(NaN, NaN, nom_col, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
h2 = fill(NaN, NaN, rob_col, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
legend([h1, h2], {'Nominal', 'Robust'}, 'Location', 'southeast', 'FontSize', 10);
xlabel('$x_1$ (position)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$ (velocity)', 'Interpreter', 'latex', 'FontSize', 12);
title('(c) Comparison', 'FontSize', 13);
xlim([xmin-0.3, xmax+0.3]); ylim([ymin-0.3, ymax+0.3]);
set(gca, 'FontSize', 11);

% Area stats
nom_area = sum(nom_final(:)) * dx1 * dx2;
rob_area = sum(rob_final(:)) * dx1 * dx2;
state_area = (xmax-xmin)*(ymax-ymin);
fprintf('\nNominal reachable area: %.2f (%.1f%% of state space)\n', nom_area, 100*nom_area/state_area);
fprintf('Robust reachable area:  %.2f (%.1f%% of state space)\n', rob_area, 100*rob_area/state_area);
fprintf('Shrinkage (robust/nominal): %.1f%%\n', 100*rob_area/nom_area);

exportgraphics(fig, fullfile(fig_dir, 'fig_double_integrator_reachability.png'), ...
    'Resolution', 300);
fprintf('\nSaved: figures/fig_double_integrator_reachability.png\n');
