%NOMINAL_FORWARD  Forward reachable set (nominal, deterministic).
%
%   Computes the set of states reachable from the target set in exactly
%   N steps while satisfying state constraints at every step,
%   assuming NO disturbance.
%
%   Mathematical definition (forward recursion):
%       R_0 = Target set
%       R_{k+1} = (A * R_k ⊕ B * U) ∩ X_state
%
%   where A*R ⊕ B*U is the Minkowski sum:
%       { A*x + B*u : x ∈ R, u ∈ U }
%
%   In H-representation {x: H*x ≤ h}, the image A*R has constraints:
%       H*A^{-1}*z ≤ h  =>  (H*A^{-1})*z ≤ h
%   and the Minkowski sum expands:
%       h_new(i) = h(i) + max_{u∈U} (H*A^{-1})(i,:) * B * u
%                = h(i) + |(H*A^{-1})(i,:) * B| * u_max
%
%   Intersect with state constraints by appending rows.
%
%   Run:  >> nominal_forward

clear; clc; close all;

cu = common_utils();
A = cu.A;  B = cu.B;
u_max = cu.u_max;
N = cu.N;

H_state = [eye(2); -eye(2)];
h_state = [cu.x_max; -cu.x_min];

% Start from target set
H = cu.H_target;
h = cu.h_target;

% Grid
x1 = linspace(cu.x_min(1), cu.x_max(1), cu.nx_grid);
x2 = linspace(cu.x_min(2), cu.x_max(2), cu.nx_grid);
[X1, X2] = meshgrid(x1, x2);

safe_masks = cell(N+1, 1);
safe_masks{1} = evaluate_polytope(H, h, X1, X2);

fprintf('Nominal forward reachability (N=%d):\n', N);

Ainv = inv(A);

for k = 1:N
    % Image through A: H_new = H * A^{-1}
    H_img = H * Ainv;

    % Minkowski sum with B*U: expand bounds
    HB = H_img * B;  % n_con x 1
    h_expanded = h + abs(HB) * u_max;

    % Intersect with state constraints
    H_new = [H_img; H_state];
    h_new = [h_expanded; h_state];

    H = H_new;
    h = h_new;

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
    'FaceColor', cu.col_nominal, 'FaceAlpha', 0.6);

contourf(X1, X2, double(safe_masks{1}), [0.5 0.5], ...
    'FaceColor', cu.col_target, 'FaceAlpha', 0.8);

h1 = fill(NaN, NaN, cu.col_nominal, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
h2 = fill(NaN, NaN, cu.col_target, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
h3 = fill(NaN, NaN, cu.col_constraint, 'FaceAlpha', 0.3, 'EdgeColor', 'k');
legend([h1 h2 h3], {sprintf('Forward reachable (%d steps)', N), ...
    'Initial set', 'State constraints'}, ...
    'Location', 'southeast', 'FontSize', 10);

xlabel('Position x_1', 'FontSize', 12);
ylabel('Velocity x_2', 'FontSize', 12);
title('Nominal Forward Reachable Set (Double Integrator)', 'FontSize', 13);
grid on; axis equal;
xlim([cu.x_min(1) cu.x_max(1)]);
ylim([cu.x_min(2) cu.x_max(2)]);
set(gca, 'FontSize', 11);

saveas(fig, fullfile('Example', 'figures', 'nominal_forward.png'));
fprintf('Saved: Example/figures/nominal_forward.png\n');


function mask = evaluate_polytope(H, h, X1, X2)
    mask = true(size(X1));
    for i = 1:size(H, 1)
        mask = mask & (H(i,1)*X1 + H(i,2)*X2 <= h(i) + 1e-9);
    end
end
