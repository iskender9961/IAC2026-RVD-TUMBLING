%STOCHASTIC_FORWARD  Forward reachable set with chance constraints.
%
%   Computes the forward reachable tube under Gaussian process noise,
%   with probabilistic constraint satisfaction at every step.
%
%   Model: x_{k+1} = A*x_k + B*u_k + w_k,  w_k ~ N(0, W)
%
%   The set at step k represents all states that can be reached
%   from the initial set while satisfying tightened state constraints:
%       h_state_tight(i) = h_state(i) - z_{1-alpha/n_c} * sigma_i(k)
%
%   The Minkowski sum now includes the noise expansion. For the
%   stochastic forward set, we compute the REACHABLE TUBE where
%   at each step the constraints are tightened to account for
%   accumulated noise variance.
%
%   Run:  >> stochastic_forward

clear; clc; close all;

cu = common_utils();
A = cu.A;  B = cu.B;
u_max = cu.u_max;
N = cu.N;
W = cu.w_cov;
alpha = cu.alpha;

H_state = [eye(2); -eye(2)];
h_state = [cu.x_max; -cu.x_min];
n_state_con = size(H_state, 1);
n_target_con = size(cu.H_target, 1);
n_total_con = n_state_con + n_target_con;
z_q = norminv(1 - alpha / n_total_con);

Ainv = inv(A);

% Start from initial set
H = cu.H_target;
h = cu.h_target;

% Grid
x1 = linspace(cu.x_min(1), cu.x_max(1), cu.nx_grid);
x2 = linspace(cu.x_min(2), cu.x_max(2), cu.nx_grid);
[X1, X2] = meshgrid(x1, x2);

safe_masks = cell(N+1, 1);
safe_masks{1} = evaluate_polytope(H, h, X1, X2);

fprintf('Stochastic forward reachability (N=%d, alpha=%.3f):\n', N, alpha);

% Accumulated covariance
Sigma_accum = zeros(2);

for k = 1:N
    % Accumulate noise
    Sigma_accum = A * Sigma_accum * A' + W;

    % Image through A
    H_img = H * Ainv;
    HB = H_img * B;
    h_expanded = h + abs(HB) * u_max;

    % Tighten state constraints by accumulated noise
    h_state_tight = h_state;
    for i = 1:n_state_con
        ai = H_state(i,:)';
        sigma_i = sqrt(max(ai' * Sigma_accum * ai, 0));
        h_state_tight(i) = h_state_tight(i) - z_q * sigma_i;
    end

    % Also tighten the image set constraints
    h_img_tight = h_expanded;
    for i = 1:size(H_img, 1)
        ai = H_img(i,:)';
        sigma_i = sqrt(max(ai' * Sigma_accum * ai, 0));
        h_img_tight(i) = h_img_tight(i) - z_q * sigma_i;
    end

    % Intersect
    H = [H_img; H_state];
    h = [h_img_tight; h_state_tight];

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
    'FaceColor', cu.col_stochastic, 'FaceAlpha', 0.6);

contourf(X1, X2, double(safe_masks{1}), [0.5 0.5], ...
    'FaceColor', cu.col_target, 'FaceAlpha', 0.8);

h1 = fill(NaN, NaN, cu.col_stochastic, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
h2 = fill(NaN, NaN, cu.col_target, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
legend([h1 h2], {sprintf('Stochastic forward (%d steps)', N), 'Initial set'}, ...
    'Location', 'southeast', 'FontSize', 10);

xlabel('Position x_1', 'FontSize', 12);
ylabel('Velocity x_2', 'FontSize', 12);
title(sprintf('Stochastic Forward Reachable Set (\\alpha=%.2f)', alpha), 'FontSize', 13);
grid on; axis equal;
xlim([cu.x_min(1) cu.x_max(1)]);
ylim([cu.x_min(2) cu.x_max(2)]);
set(gca, 'FontSize', 11);

saveas(fig, fullfile('Example', 'figures', 'stochastic_forward.png'));
fprintf('Saved: Example/figures/stochastic_forward.png\n');


function mask = evaluate_polytope(H, h, X1, X2)
    mask = true(size(X1));
    for i = 1:size(H, 1)
        mask = mask & (H(i,1)*X1 + H(i,2)*X2 <= h(i) + 1e-9);
    end
end
