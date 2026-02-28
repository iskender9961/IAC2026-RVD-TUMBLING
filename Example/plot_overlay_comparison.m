%PLOT_OVERLAY_COMPARISON  Overlay all three backward reachable sets.
%
%   Produces a single figure showing:
%     - Nominal backward (green, largest)
%     - Stochastic backward (blue, middle)
%     - Robust backward (purple, smallest)
%     - Target set (gold)
%     - State constraint box (gray)
%
%   Verifies the nesting: X_robust ⊆ X_stochastic ⊆ X_nominal
%
%   Run:  >> plot_overlay_comparison

clear; clc; close all;

cu = common_utils();
A = cu.A;  B = cu.B;
u_max = cu.u_max;
N = cu.N;

H_state = [eye(2); -eye(2)];
h_state = [cu.x_max; -cu.x_min];

% Grid
x1 = linspace(cu.x_min(1), cu.x_max(1), cu.nx_grid);
x2 = linspace(cu.x_min(2), cu.x_max(2), cu.nx_grid);
[X1, X2] = meshgrid(x1, x2);

%% === Nominal backward ===
H = cu.H_target; h = cu.h_target;
for k = N:-1:1
    HA = H * A; HB = H * B;
    h_tight = h - abs(HB) * u_max;
    H = [HA; H_state]; h = [h_tight; h_state];
end
mask_nom = evaluate_polytope(H, h, X1, X2);

%% === Stochastic backward ===
W = cu.w_cov; alpha = cu.alpha;
n_con = size(H_state,1) + size(cu.H_target,1);
z_q = norminv(1 - alpha / n_con);

H = cu.H_target; h = cu.h_target;
for k = N:-1:1
    Sig_k = zeros(2);
    Aj = eye(2);
    for j = 0:N-k-1
        Sig_k = Sig_k + Aj * W * Aj';
        Aj = A * Aj;
    end

    HA = H * A; HB = H * B;
    h_tight = h - abs(HB) * u_max;
    for i = 1:size(HA,1)
        ai = HA(i,:)';
        h_tight(i) = h_tight(i) - z_q * sqrt(max(ai'*Sig_k*ai, 0));
    end
    h_state_t = h_state;
    for i = 1:size(H_state,1)
        ai = H_state(i,:)';
        h_state_t(i) = h_state_t(i) - z_q * sqrt(max(ai'*Sig_k*ai, 0));
    end
    H = [HA; H_state]; h = [h_tight; h_state_t];
end
mask_stoch = evaluate_polytope(H, h, X1, X2);

%% === Robust backward ===
w_max = cu.w_max;
H = cu.H_target; h = cu.h_target;
for k = N:-1:1
    HA = H * A; HB = H * B;
    h_tight = h - abs(HB) * u_max;
    for i = 1:size(HA,1)
        h_tight(i) = h_tight(i) - abs(HA(i,:)) * w_max;
    end
    h_state_t = h_state;
    for i = 1:size(H_state,1)
        h_state_t(i) = h_state_t(i) - abs(H_state(i,:)) * w_max;
    end
    H = [HA; H_state]; h = [h_tight; h_state_t];
end
mask_rob = evaluate_polytope(H, h, X1, X2);

%% === Target set ===
mask_target = evaluate_polytope(cu.H_target, cu.h_target, X1, X2);

%% === Verify nesting ===
fprintf('=== Nesting verification ===\n');
n_nom   = sum(mask_nom(:));
n_stoch = sum(mask_stoch(:));
n_rob   = sum(mask_rob(:));

rob_in_stoch = sum(mask_rob(:) & ~mask_stoch(:));
stoch_in_nom = sum(mask_stoch(:) & ~mask_nom(:));

fprintf('  |X_nominal|   = %d points\n', n_nom);
fprintf('  |X_stochastic| = %d points\n', n_stoch);
fprintf('  |X_robust|     = %d points\n', n_rob);
fprintf('  Robust ⊄ Stoch violations: %d\n', rob_in_stoch);
fprintf('  Stoch ⊄ Nominal violations: %d\n', stoch_in_nom);

if rob_in_stoch == 0 && stoch_in_nom == 0
    fprintf('  [OK] Nesting verified: X_robust ⊆ X_stochastic ⊆ X_nominal\n');
else
    fprintf('  [WARN] Nesting violations detected!\n');
end

%% === Plot ===
fig = figure('Position', [100 100 800 600]);
hold on;

% State constraint box
fill([cu.x_min(1) cu.x_max(1) cu.x_max(1) cu.x_min(1)], ...
     [cu.x_min(2) cu.x_min(2) cu.x_max(2) cu.x_max(2)], ...
     cu.col_constraint, 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.3);

% Nominal (outermost, green)
contourf(X1, X2, double(mask_nom), [0.5 0.5], ...
    'FaceColor', cu.col_nominal, 'FaceAlpha', 0.5);

% Stochastic (middle, blue)
contourf(X1, X2, double(mask_stoch), [0.5 0.5], ...
    'FaceColor', cu.col_stochastic, 'FaceAlpha', 0.6);

% Robust (innermost, purple)
contourf(X1, X2, double(mask_rob), [0.5 0.5], ...
    'FaceColor', cu.col_robust, 'FaceAlpha', 0.7);

% Target set (gold)
contourf(X1, X2, double(mask_target), [0.5 0.5], ...
    'FaceColor', cu.col_target, 'FaceAlpha', 0.9);

% Origin
plot(0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2);

% Legend
p1 = fill(NaN, NaN, cu.col_nominal, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
p2 = fill(NaN, NaN, cu.col_stochastic, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
p3 = fill(NaN, NaN, cu.col_robust, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
p4 = fill(NaN, NaN, cu.col_target, 'FaceAlpha', 0.9, 'EdgeColor', 'none');
p5 = fill(NaN, NaN, cu.col_constraint, 'FaceAlpha', 0.3, 'EdgeColor', 'k');

legend([p1 p2 p3 p4 p5], ...
    {sprintf('Nominal (%d pts)', n_nom), ...
     sprintf('Stochastic, 95%% (%d pts)', n_stoch), ...
     sprintf('Robust (%d pts)', n_rob), ...
     'Target set', 'State constraints'}, ...
    'Location', 'southeast', 'FontSize', 10);

xlabel('Position x_1', 'FontSize', 12);
ylabel('Velocity x_2', 'FontSize', 12);
title(sprintf('Nested Backward Reachable Sets (N=%d steps)\n\\mathcal{X}_{rob} \\subseteq \\mathcal{X}_{stoch} \\subseteq \\mathcal{X}_{nom}', N), ...
    'FontSize', 13);
grid on; axis equal;
xlim([cu.x_min(1) cu.x_max(1)]);
ylim([cu.x_min(2) cu.x_max(2)]);
set(gca, 'FontSize', 11);

saveas(fig, fullfile('Example', 'figures', 'overlay_comparison.png'));
exportgraphics(fig, fullfile('Example', 'figures', 'overlay_comparison.pdf'), 'ContentType', 'vector');
fprintf('Saved: Example/figures/overlay_comparison.png/.pdf\n');


function mask = evaluate_polytope(H, h, X1, X2)
    mask = true(size(X1));
    for i = 1:size(H, 1)
        mask = mask & (H(i,1)*X1 + H(i,2)*X2 <= h(i) + 1e-9);
    end
end
