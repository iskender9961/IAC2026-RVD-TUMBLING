%NOMINAL_BACKWARD  Backward reachable set (nominal, deterministic).
%
%   Computes the set of states from which the target can be reached
%   in exactly N steps while satisfying state constraints at every step,
%   assuming NO disturbance.
%
%   Mathematical definition (backward recursion):
%       S_N = Target set
%       S_k = Pre(S_{k+1}) ∩ X_state     for k = N-1, ..., 0
%
%   where Pre(S) = { x : ∃ u ∈ U s.t. Ax + Bu ∈ S }
%
%   For S = {z: H*z ≤ h}, the pre-image in H-representation is:
%       Pre(S) = { x : H*A*x ≤ h - max_{u∈U} H*B*u }
%              = { x : H*A*x ≤ h_tight }
%   where h_tight(i) = h(i) - |H(i,:)*B| * u_max  (for symmetric box U).
%
%   The intersection with state constraints adds rows to [H; h].
%
%   Run:  >> nominal_backward

clear; clc; close all;

cu = common_utils();
A = cu.A;  B = cu.B;
u_max = cu.u_max;
N = cu.N;

% State constraint polytope: H_state * x <= h_state
H_state = [eye(2); -eye(2)];
h_state = [cu.x_max; -cu.x_min];

% Target set
H = cu.H_target;
h = cu.h_target;

% Grid for evaluation
x1 = linspace(cu.x_min(1), cu.x_max(1), cu.nx_grid);
x2 = linspace(cu.x_min(2), cu.x_max(2), cu.nx_grid);
[X1, X2] = meshgrid(x1, x2);

% Store safe sets at each backward step
safe_masks = cell(N+1, 1);
safe_masks{N+1} = evaluate_polytope(H, h, X1, X2);

fprintf('Nominal backward reachability (N=%d):\n', N);

for k = N:-1:1
    % Pre-image: H*A*x <= h - support(H*B, U)
    HA = H * A;
    HB = H * B;  % n_con x 1

    % Tighten: for box U = [-u_max, u_max], support = |HB| * u_max
    h_tight = h - abs(HB) * u_max;

    % Intersect with state constraints
    H_new = [HA; H_state];
    h_new = [h_tight; h_state];

    % Remove redundant constraints (keep all for simplicity)
    H = H_new;
    h = h_new;

    safe_masks{k} = evaluate_polytope(H, h, X1, X2);

    n_safe = sum(safe_masks{k}(:));
    fprintf('  Step %d: %d safe points\n', k, n_safe);
end

% Plot
fig = figure('Position', [100 100 700 550]);
hold on;

% State constraint box
fill([cu.x_min(1) cu.x_max(1) cu.x_max(1) cu.x_min(1)], ...
     [cu.x_min(2) cu.x_min(2) cu.x_max(2) cu.x_max(2)], ...
     cu.col_constraint, 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.3);

% Backward reachable set (step 1 = full N-step backward)
contourf(X1, X2, double(safe_masks{1}), [0.5 0.5], ...
    'FaceColor', cu.col_nominal, 'FaceAlpha', 0.6);

% Target set
contourf(X1, X2, double(safe_masks{N+1}), [0.5 0.5], ...
    'FaceColor', cu.col_target, 'FaceAlpha', 0.8);

% Legend
h1 = fill(NaN, NaN, cu.col_nominal, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
h2 = fill(NaN, NaN, cu.col_target, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
h3 = fill(NaN, NaN, cu.col_constraint, 'FaceAlpha', 0.3, 'EdgeColor', 'k');
legend([h1 h2 h3], {sprintf('Backward reachable (%d steps)', N), ...
    'Target set', 'State constraints'}, ...
    'Location', 'southeast', 'FontSize', 10);

xlabel('Position x_1', 'FontSize', 12);
ylabel('Velocity x_2', 'FontSize', 12);
title('Nominal Backward Reachable Set (Double Integrator)', 'FontSize', 13);
grid on; axis equal;
xlim([cu.x_min(1) cu.x_max(1)]);
ylim([cu.x_min(2) cu.x_max(2)]);
set(gca, 'FontSize', 11);

saveas(fig, fullfile('Example', 'figures', 'nominal_backward.png'));
fprintf('Saved: Example/figures/nominal_backward.png\n');


function mask = evaluate_polytope(H, h, X1, X2)
%EVALUATE_POLYTOPE  Check which grid points satisfy H*x <= h.
    mask = true(size(X1));
    for i = 1:size(H, 1)
        mask = mask & (H(i,1)*X1 + H(i,2)*X2 <= h(i) + 1e-9);
    end
end
