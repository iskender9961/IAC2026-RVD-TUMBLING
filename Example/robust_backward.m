%ROBUST_BACKWARD  Backward reachable set with bounded disturbances.
%
%   Computes the backward reachable set (safe-start region) guaranteeing
%   feasibility for ALL disturbance realizations in a bounded set.
%
%   Model: x_{k+1} = A*x_k + B*u_k + w_k,  w_k ∈ W = {w: |w_i| ≤ w_max}
%
%   The robust backward recursion tightens constraints by the WORST-CASE
%   disturbance (support function) at each step:
%
%       S_N = Target set (tightened)
%       S_k = RobustPre(S_{k+1}) ∩ X_state_tight
%
%   where RobustPre accounts for BOTH the "∃ u" (expands via B*U)
%   and the "∀ w" (tightens via W):
%
%       RobustPre(S) = { x : ∀ w ∈ W, ∃ u ∈ U s.t. Ax + Bu + w ∈ S }
%
%   In H-representation { z: H*z ≤ h }:
%       RobustPre(S) = { x : H*A*x ≤ h - max_{u∈U} H*B*u - max_{w∈W} H*w }
%
%   Wait — this is WRONG! The correct robust pre-image is:
%       RobustPre(S) = { x : ∃ u ∈ U s.t. ∀ w ∈ W: Ax + Bu + w ∈ S }
%                    = { x : ∃ u ∈ U: H*A*x + H*B*u ≤ h - max_{w∈W} H*w }
%   Eliminating u:  H*A*x ≤ h - max_{w∈W} H*w - max_{u∈U}(-H*B*u)
%                 = h - support_W(H) + support_U(H*B)
%
%   Actually: since ∃u means we pick the best u:
%       H*(Ax + Bu + w) ≤ h  for all w ∈ W
%       H*A*x ≤ h - H*B*u - max_{w∈W} H*w
%   We need ∃u: H*A*x + H*B*u ≤ h_tight  where h_tight = h - support_W(H)
%   Eliminating u (support over U):
%       H*A*x ≤ h_tight - support_{-U}(H*B)  ... no, this eliminates u by
%       requiring the constraint for the WORST u, but we need the BEST u.
%
%   Correct elimination: ∃u means we subtract the support:
%       h_final(i) = h(i) - support_W(H(i,:)) - support_{-U}(H(i,:)*B)
%                  = h(i) - |H(i,:)| * w_max  - |H(i,:)*B| * u_max
%
%   Wait, that double-subtracts. Let me be precise:
%
%   ∃ u ∈ U: H*(Ax + Bu + w) ≤ h  ∀ w ∈ W
%   ⟺ ∃ u ∈ U: H*A*x + H*B*u ≤ h - max_{w∈W} H*w
%   ⟺ H*A*x ≤ h - max_{w∈W} H*w - min_{u∈U}(-H*B*u)
%           wait, we want ∃u satisfying row-wise, not all rows same u.
%
%   For the H-rep elimination (conservative): the standard approach is
%   h_tight = h - support_W - support_{-U}(HB)
%   = h - |H|*w_max - |HB|*u_max
%
%   But this is OVERLY conservative because it eliminates u row-by-row.
%   The correct approach uses: ∃u: H*B*u ≤ h_shifted - H*A*x
%   which is an LP feasibility check.
%
%   For this pedagogical example, we use the conservative row-by-row
%   approach, which gives a valid INNER approximation.
%
%   Run:  >> robust_backward

clear; clc; close all;

cu = common_utils();
A = cu.A;  B = cu.B;
u_max = cu.u_max;
N = cu.N;
w_max = cu.w_max;

H_state = [eye(2); -eye(2)];
h_state = [cu.x_max; -cu.x_min];

% Target set
H = cu.H_target;
h = cu.h_target;

% Grid
x1 = linspace(cu.x_min(1), cu.x_max(1), cu.nx_grid);
x2 = linspace(cu.x_min(2), cu.x_max(2), cu.nx_grid);
[X1, X2] = meshgrid(x1, x2);

safe_masks = cell(N+1, 1);
safe_masks{N+1} = evaluate_polytope(H, h, X1, X2);

fprintf('Robust backward reachability (N=%d):\n', N);

for k = N:-1:1
    % Pre-image through dynamics
    HA = H * A;
    HB = H * B;

    % Tighten for "∃ u ∈ U": subtract support of U through B
    h_u = h - abs(HB) * u_max;

    % Tighten for "∀ w ∈ W": subtract support function of W
    % support_W(a) = |a_1| * w_max(1) + |a_2| * w_max(2)
    h_uw = h_u;
    for i = 1:size(HA, 1)
        support_w = abs(HA(i,:)) * w_max;
        h_uw(i) = h_uw(i) - support_w;
    end

    % State constraints tightened by worst-case disturbance
    h_state_tight = h_state;
    for i = 1:size(H_state, 1)
        support_w = abs(H_state(i,:)) * w_max;
        h_state_tight(i) = h_state_tight(i) - support_w;
    end

    % Intersect
    H = [HA; H_state];
    h = [h_uw; h_state_tight];

    safe_masks{k} = evaluate_polytope(H, h, X1, X2);

    n_safe = sum(safe_masks{k}(:));
    fprintf('  Step %d: %d safe points\n', k, n_safe);
end

% Plot
fig = figure('Position', [100 100 700 550]);
hold on;

fill([cu.x_min(1) cu.x_max(1) cu.x_max(1) cu.x_min(1)], ...
     [cu.x_min(2) cu.x_min(2) cu.x_max(2) cu.x_max(2)], ...
     cu.col_constraint, 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.3);

contourf(X1, X2, double(safe_masks{1}), [0.5 0.5], ...
    'FaceColor', cu.col_robust, 'FaceAlpha', 0.6);

contourf(X1, X2, double(safe_masks{N+1}), [0.5 0.5], ...
    'FaceColor', cu.col_target, 'FaceAlpha', 0.8);

h1 = fill(NaN, NaN, cu.col_robust, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
h2 = fill(NaN, NaN, cu.col_target, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
legend([h1 h2], {sprintf('Robust backward (%d steps)', N), 'Target set'}, ...
    'Location', 'southeast', 'FontSize', 10);

xlabel('Position x_1', 'FontSize', 12);
ylabel('Velocity x_2', 'FontSize', 12);
title('Robust Backward Reachable Set (Double Integrator)', 'FontSize', 13);
grid on; axis equal;
xlim([cu.x_min(1) cu.x_max(1)]);
ylim([cu.x_min(2) cu.x_max(2)]);
set(gca, 'FontSize', 11);

saveas(fig, fullfile('Example', 'figures', 'robust_backward.png'));
fprintf('Saved: Example/figures/robust_backward.png\n');


function mask = evaluate_polytope(H, h, X1, X2)
    mask = true(size(X1));
    for i = 1:size(H, 1)
        mask = mask & (H(i,1)*X1 + H(i,2)*X2 <= h(i) + 1e-9);
    end
end
