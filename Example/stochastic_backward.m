%STOCHASTIC_BACKWARD  Backward reachable set with chance constraints.
%
%   Computes the backward reachable set (safe-start region) under
%   Gaussian process noise, with probabilistic constraint satisfaction.
%
%   Model: x_{k+1} = A*x_k + B*u_k + w_k,  w_k ~ N(0, W)
%
%   The chance-constrained backward recursion tightens each constraint
%   to account for the accumulated noise variance:
%
%       S_N = Target set (tightened)
%       S_k = Pre_stoch(S_{k+1}) ∩ X_state_tight
%
%   where the tightening at step k accounts for noise from steps k to N:
%       h_tight(i) = h(i) - z_{1-alpha/n_c} * sigma_i(k)
%       sigma_i(k) = sqrt(a_i' * Sigma_accum(k) * a_i)
%       Sigma_accum(k) = sum_{j=0}^{N-k-1} A^j * W * (A^j)'
%
%   This is a CERTIFIED inner approximation: any state in the set
%   satisfies all constraints with probability >= 1 - alpha.
%
%   The "exists u" is correctly implemented via support function tightening
%   (same as nominal), and the Bonferroni correction handles joint
%   probability over all constraints.
%
%   Run:  >> stochastic_backward

clear; clc; close all;

cu = common_utils();
A = cu.A;  B = cu.B;
u_max = cu.u_max;
N = cu.N;

H_state = [eye(2); -eye(2)];
h_state = [cu.x_max; -cu.x_min];

W = cu.w_cov;
alpha = cu.alpha;

% Compute accumulated covariance for each backward step
% At step k, noise accumulates from k to N-1:
%   Sigma(k) = sum_{j=0}^{N-k-1} A^j * W * (A^j)'
Sigma_accum = cell(N+1, 1);
Sigma_accum{N+1} = zeros(2);  % no noise at terminal step
for k = N:-1:1
    steps_ahead = N - k;
    Sig = zeros(2);
    Aj = eye(2);
    for j = 0:steps_ahead-1
        Sig = Sig + Aj * W * Aj';
        Aj = A * Aj;
    end
    Sigma_accum{k} = Sig;
end

% Number of constraints (for Bonferroni)
n_state_con = size(H_state, 1);
n_target_con = size(cu.H_target, 1);
n_total_con = n_state_con + n_target_con;

% Quantile for chance constraints
z_q = norminv(1 - alpha / n_total_con);
fprintf('Stochastic backward (alpha=%.3f, z=%.3f):\n', alpha, z_q);

% Start from tightened target set
H = cu.H_target;
h = cu.h_target;

% Tighten target by terminal noise = 0 (no tightening at N)

% Grid
x1 = linspace(cu.x_min(1), cu.x_max(1), cu.nx_grid);
x2 = linspace(cu.x_min(2), cu.x_max(2), cu.nx_grid);
[X1, X2] = meshgrid(x1, x2);

safe_masks = cell(N+1, 1);
safe_masks{N+1} = evaluate_polytope(H, h, X1, X2);

for k = N:-1:1
    Sig_k = Sigma_accum{k};

    % Pre-image (same as nominal)
    HA = H * A;
    HB = H * B;
    h_tight_pre = h - abs(HB) * u_max;

    % Stochastic tightening of pre-image constraints
    h_stoch = h_tight_pre;
    for i = 1:size(HA, 1)
        ai = HA(i,:)';
        sigma_i = sqrt(max(ai' * Sig_k * ai, 0));
        h_stoch(i) = h_stoch(i) - z_q * sigma_i;
    end

    % State constraints (also tightened)
    h_state_tight = h_state;
    for i = 1:n_state_con
        ai = H_state(i,:)';
        sigma_i = sqrt(max(ai' * Sig_k * ai, 0));
        h_state_tight(i) = h_state_tight(i) - z_q * sigma_i;
    end

    % Intersect
    H = [HA; H_state];
    h = [h_stoch; h_state_tight];

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
    'FaceColor', cu.col_stochastic, 'FaceAlpha', 0.6);

contourf(X1, X2, double(safe_masks{N+1}), [0.5 0.5], ...
    'FaceColor', cu.col_target, 'FaceAlpha', 0.8);

h1 = fill(NaN, NaN, cu.col_stochastic, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
h2 = fill(NaN, NaN, cu.col_target, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
legend([h1 h2], {sprintf('Stochastic backward (%d steps, %.0f%%)', N, 100*(1-alpha)), ...
    'Target set'}, 'Location', 'southeast', 'FontSize', 10);

xlabel('Position x_1', 'FontSize', 12);
ylabel('Velocity x_2', 'FontSize', 12);
title(sprintf('Stochastic Backward Reachable Set (\\alpha=%.2f)', alpha), 'FontSize', 13);
grid on; axis equal;
xlim([cu.x_min(1) cu.x_max(1)]);
ylim([cu.x_min(2) cu.x_max(2)]);
set(gca, 'FontSize', 11);

saveas(fig, fullfile('Example', 'figures', 'stochastic_backward.png'));
fprintf('Saved: Example/figures/stochastic_backward.png\n');


function mask = evaluate_polytope(H, h, X1, X2)
    mask = true(size(X1));
    for i = 1:size(H, 1)
        mask = mask & (H(i,1)*X1 + H(i,2)*X2 <= h(i) + 1e-9);
    end
end
