%ROBUST_FORWARD  Forward reachable set with bounded disturbances.
%
%   Computes the forward reachable tube guaranteeing state constraint
%   satisfaction for ALL disturbance realizations.
%
%   Model: x_{k+1} = A*x_k + B*u_k + w_k,  w_k ∈ W
%
%   The robust forward set is TIGHTER than nominal because at each step
%   the state constraints must be tightened by the worst-case accumulated
%   disturbance, and the Minkowski sum with B*U is REDUCED by the
%   Minkowski sum with W (the disturbance expands uncertainty but we
%   need to guarantee constraints).
%
%   In H-representation:
%       After image through A and Minkowski sum with B*U:
%           h_expanded(i) = h(i) + |H_img(i,:)*B| * u_max
%       Tighten for worst-case disturbance:
%           h_tight(i) = h_expanded(i) - |H_img(i,:)| * w_max (accumulated)
%       Intersect with tightened state constraints.
%
%   Run:  >> robust_forward

clear; clc; close all;

cu = common_utils();
A = cu.A;  B = cu.B;
u_max = cu.u_max;
N = cu.N;
w_max = cu.w_max;

H_state = [eye(2); -eye(2)];
h_state = [cu.x_max; -cu.x_min];

Ainv = inv(A);

H = cu.H_target;
h = cu.h_target;

% Grid
x1 = linspace(cu.x_min(1), cu.x_max(1), cu.nx_grid);
x2 = linspace(cu.x_min(2), cu.x_max(2), cu.nx_grid);
[X1, X2] = meshgrid(x1, x2);

safe_masks = cell(N+1, 1);
safe_masks{1} = evaluate_polytope(H, h, X1, X2);

fprintf('Robust forward reachability (N=%d):\n', N);

% Accumulated worst-case disturbance per step
% At step k, total accumulated = sum of per-step worst-case
Aj = eye(2);
accum_support = zeros(size(H_state, 1), 1);

for k = 1:N
    % Image through A
    H_img = H * Ainv;
    HB = H_img * B;
    h_expanded = h + abs(HB) * u_max;

    % Per-step disturbance tightening on image constraints
    h_tight = h_expanded;
    for i = 1:size(H_img, 1)
        support_w = abs(H_img(i,:)) * w_max;
        h_tight(i) = h_tight(i) - support_w;
    end

    % Accumulate disturbance effect on state constraints
    for i = 1:size(H_state, 1)
        accum_support(i) = accum_support(i) + abs(H_state(i,:) * Aj) * w_max;
    end
    Aj = A * Aj;

    h_state_tight = h_state - accum_support;

    % Intersect
    H = [H_img; H_state];
    h = [h_tight; h_state_tight];

    safe_masks{k+1} = evaluate_polytope(H, h, X1, X2);

    n_safe = sum(safe_masks{k+1}(:));
    fprintf('  Step %d: %d reachable points\n', k, n_safe);
end

% Plot
fig = figure('Position', [100 100 700 550]);
hold on;

fill([cu.x_min(1) cu.x_max(1) cu.x_max(1) cu.x_min(1)], ...
     [cu.x_min(2) cu.x_min(2) cu.x_max(2) cu.x_max(2)], ...
     cu.col_constraint, 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.3);

contourf(X1, X2, double(safe_masks{N+1}), [0.5 0.5], ...
    'FaceColor', cu.col_robust, 'FaceAlpha', 0.6);

contourf(X1, X2, double(safe_masks{1}), [0.5 0.5], ...
    'FaceColor', cu.col_target, 'FaceAlpha', 0.8);

h1 = fill(NaN, NaN, cu.col_robust, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
h2 = fill(NaN, NaN, cu.col_target, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
legend([h1 h2], {sprintf('Robust forward (%d steps)', N), 'Initial set'}, ...
    'Location', 'southeast', 'FontSize', 10);

xlabel('Position x_1', 'FontSize', 12);
ylabel('Velocity x_2', 'FontSize', 12);
title('Robust Forward Reachable Set (Double Integrator)', 'FontSize', 13);
grid on; axis equal;
xlim([cu.x_min(1) cu.x_max(1)]);
ylim([cu.x_min(2) cu.x_max(2)]);
set(gca, 'FontSize', 11);

saveas(fig, fullfile('Example', 'figures', 'robust_forward.png'));
fprintf('Saved: Example/figures/robust_forward.png\n');


function mask = evaluate_polytope(H, h, X1, X2)
    mask = true(size(X1));
    for i = 1:size(H, 1)
        mask = mask & (H(i,1)*X1 + H(i,2)*X2 <= h(i) + 1e-9);
    end
end
