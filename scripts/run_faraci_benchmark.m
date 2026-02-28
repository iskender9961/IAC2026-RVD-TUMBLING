%RUN_FARACI_BENCHMARK  Benchmark against Faraci & Lampariello (2025).
%
%   Reproduces the tumbling-target interception scenario from:
%     Faraci & Lampariello, "Reachability-Guaranteed Optimal Control for
%     the Interception of Dynamic Targets under Uncertainty," arXiv:2503.16935v2
%
%   Their approach: Reachability-Guaranteed OCP (RG-OCP)
%     - Forward RA on controlled system to guarantee X(T_f) ⊃ Y_target(T_f)
%     - NLP-based, solved in 2658 ms per instance
%     - Guarantees the chaser CAN REACH all target orientations
%
%   Our approach: Hierarchical safe-start certification
%     - Erosion-based + multi-step backward RA
%     - Analytical, solved in <1 s per (omega, a_max) combo
%     - Guarantees the chaser STAYS INSIDE the LOS corridor
%
%   These are complementary: safe-start (ours) + guaranteed reachability
%   (theirs) = complete certification pipeline.
%
%   Benchmark: Match Faraci's parameters and compare:
%     1. Safe-start region coverage vs tumble rate
%     2. Computation time (analytical vs NLP)
%     3. Effect of orientation uncertainty on safe region
%
%   Run:  >> run_faraci_benchmark

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
root_dir = fileparts(this_dir);
if isempty(root_dir), root_dir = pwd; end

addpath(fullfile(root_dir, 'dynamics'));
addpath(fullfile(root_dir, 'reachability_common'));
addpath(fullfile(root_dir, 'reachability_nominal'));
addpath(root_dir);

results_dir = fullfile(root_dir, 'results', 'benchmark');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% ===== Faraci's Scenario Parameters =====
% From Section VII of arXiv:2503.16935v2
fprintf('====== Benchmark: Faraci & Lampariello (2025) ======\n\n');

% Target tumble: omega = [0, 0.0698, 0] rad/s => ~4 deg/s about y-axis
omega_faraci_rad = 0.0698;   % rad/s
omega_faraci_deg = rad2deg(omega_faraci_rad);  % ~4 deg/s

% Inertia: J = diag(29.2, 30, 38.4) kg*m^2
J_target = diag([29.2, 30, 38.4]);

% Chaser: m = 32 kg, grasping point c_target = [5.5, 0, 0] m
m_chaser = 32.0;

% Time horizon: T_f = 30 s, dt = 1 s, N = 30
T_f = 30;
dt_faraci = 1.0;
N_steps_faraci = T_f / dt_faraci;

% Orientation uncertainty: rho_SO3 = 0.17 rad = 10 deg
rho_uncertainty_deg = 10;

% Their RG-OCP solve time: 2658 ms
rgocp_time_ms = 2658;

fprintf('Faraci scenario:\n');
fprintf('  Tumble rate:    %.4f rad/s (%.1f deg/s)\n', omega_faraci_rad, omega_faraci_deg);
fprintf('  Inertia:        diag(%.1f, %.1f, %.1f) kg*m^2\n', J_target(1,1), J_target(2,2), J_target(3,3));
fprintf('  Chaser mass:    %.1f kg\n', m_chaser);
fprintf('  Time horizon:   %.0f s (dt=%.0f s, N=%d)\n', T_f, dt_faraci, N_steps_faraci);
fprintf('  Uncertainty:    %.0f deg\n', rho_uncertainty_deg);
fprintf('  RG-OCP time:    %d ms\n\n', rgocp_time_ms);

%% ===== Our Framework: Safe-Start Analysis =====
% Use our reachability infrastructure with Faraci-matched parameters

rp = reachability_params();
% Override dt to match Faraci
rp.dt = dt_faraci;

% Sweep: Faraci uses ~4 deg/s; we compare across a range
omega_vals = [1, 2, 3, 4, 5, 6];       % deg/s (Faraci's is ~4)
% Thrust authority sweep: derive from chaser mass
% Typical small satellite thruster: 0.5-5 N => a_max = F/m
% Faraci doesn't specify a_max directly; they optimize control
% We sweep representative values
amax_vals = [0.20, 0.10, 0.05, 0.02];  % m/s^2

n_omega = length(omega_vals);
n_amax = length(amax_vals);

% Build LOS constraints
[A_body, b_body] = los_constraints_body(rp.cone_k, rp.y_min, rp.cone_nfaces);

fprintf('--- Safe-Start Region Analysis (our method) ---\n');
fprintf('Grid: %d x %d, dt=%.0f s, N_back=%d\n', rp.nx_grid, rp.ny_grid, rp.dt, N_steps_faraci);

%% ===== Part 1: Erosion-based analysis (fast) =====
safe_fractions_erosion = zeros(n_omega, n_amax);
erosion_times = zeros(n_omega, n_amax);

x_grid = linspace(rp.x_range(1), rp.x_range(2), rp.nx_grid);
y_grid = linspace(rp.y_range(1), rp.y_range(2), rp.ny_grid);
ny = length(y_grid); nx = length(x_grid);

% Cone mask (computed once)
cone_mask = false(ny, nx);
for iy = 1:ny
    for ix = 1:nx
        p_body = [x_grid(ix); y_grid(iy); 0; 0; 0; 0];
        cone_mask(iy, ix) = all(A_body * p_body <= b_body + 1e-9);
    end
end
n_cone = sum(cone_mask(:));

% Store results for the Faraci-matching case (omega=4, various a_max)
faraci_results = cell(1, n_amax);

for iw = 1:n_omega
    for ia = 1:n_amax
        omega_deg = omega_vals(iw);
        a_max = amax_vals(ia);
        omega_rad = deg2rad(omega_deg);

        t0 = tic;

        safe_mask = false(ny, nx);
        w2 = omega_rad^2;
        r_sync_max = 2 * a_max / max(w2, 1e-30);

        for iy = 1:ny
            for ix = 1:nx
                if ~cone_mask(iy, ix), continue; end
                xb = x_grid(ix); yb = y_grid(iy);
                rng = sqrt(xb^2 + yb^2);
                if rng >= r_sync_max, continue; end

                p_body = [xb; yb; 0; 0; 0; 0];
                slacks = b_body - A_body * p_body;
                v_rot = [omega_rad * yb; -omega_rad * xb; 0];

                safe_point = true;
                for ic = 1:size(A_body, 1)
                    a_pos = A_body(ic, 1:3)';
                    slack_rate = -a_pos' * v_rot;
                    if slack_rate < 0
                        erosion = 0.5 * slack_rate^2 / a_max;
                        if slacks(ic) < erosion
                            safe_point = false; break;
                        end
                    end
                end
                safe_mask(iy, ix) = safe_point;
            end
        end

        elapsed = toc(t0);
        n_safe = sum(safe_mask(:) & cone_mask(:));
        frac = 100 * n_safe / max(n_cone, 1);

        safe_fractions_erosion(iw, ia) = frac;
        erosion_times(iw, ia) = elapsed;

        % Store Faraci-matching case
        if omega_deg == 4
            faraci_results{ia}.safe_mask = safe_mask;
            faraci_results{ia}.frac = frac;
            faraci_results{ia}.time = elapsed;
            faraci_results{ia}.a_max = a_max;
        end

        fprintf('  w=%d deg/s, a_max=%.2f: %.1f%% safe, %.3f s\n', ...
            omega_deg, a_max, frac, elapsed);
    end
end

%% ===== Part 2: Multi-step backward at Faraci's omega =====
fprintf('\n--- Multi-step backward (Faraci omega=4 deg/s, N=%d) ---\n', N_steps_faraci);

% Coarser grid for LP-based analysis
nx_disc = 60; ny_disc = 60;
x_disc = linspace(rp.x_range(1), rp.x_range(2), nx_disc);
y_disc = linspace(rp.y_range(1), rp.y_range(2), ny_disc);

omega_faraci_rad_s = deg2rad(4);  % exact 4 deg/s
[Ad, Bd] = cwh_stm(dt_faraci, rp.n);

lp_opts = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');

multistep_fracs = zeros(1, n_amax);
multistep_times = zeros(1, n_amax);

for ia = 1:n_amax
    a_max = amax_vals(ia);
    t0 = tic;

    safe_mask_ms = false(ny_disc, nx_disc);
    cone_mask_disc = false(ny_disc, nx_disc);

    for iy = 1:ny_disc
        for ix = 1:nx_disc
            xb = x_disc(ix); yb = y_disc(iy);
            p_body = [xb; yb; 0; 0; 0; 0];
            cone_mask_disc(iy, ix) = all(A_body * p_body <= b_body + 1e-9);
            if ~cone_mask_disc(iy, ix), continue; end

            x_lvlh = [xb; yb; 0; -omega_faraci_rad_s*yb; omega_faraci_rad_s*xb; 0];
            feasible = true;
            x_k = x_lvlh;

            for k = 1:N_steps_faraci
                t_k = k * dt_faraci;
                [A_los_k, b_los_k] = los_constraints_at_time(...
                    rp.cone_k, rp.y_min, rp.cone_nfaces, omega_faraci_rad_s, t_k);

                x_next_zero = Ad * x_k;
                if all(b_los_k - A_los_k * x_next_zero >= -1e-9)
                    u_k = [0; 0; 0];
                else
                    rhs = b_los_k - A_los_k * (Ad * x_k);
                    A_u = A_los_k * Bd;
                    nu_ctrl = 3;
                    f_lp = [0; 0; 0; 1; 1; 1];
                    A_lp = [A_u, zeros(size(A_u,1), nu_ctrl);
                            eye(nu_ctrl), -eye(nu_ctrl);
                           -eye(nu_ctrl), -eye(nu_ctrl)];
                    b_lp = [rhs; zeros(nu_ctrl,1); zeros(nu_ctrl,1)];
                    lb = [-a_max*ones(nu_ctrl,1); zeros(nu_ctrl,1)];
                    ub = [a_max*ones(nu_ctrl,1); inf(nu_ctrl,1)];
                    [z_opt, ~, exitflag] = linprog(f_lp, A_lp, b_lp, ...
                        [], [], lb, ub, lp_opts);
                    if exitflag == 1
                        u_k = z_opt(1:nu_ctrl);
                    else
                        feasible = false; break;
                    end
                end
                x_k = Ad * x_k + Bd * u_k;
            end
            safe_mask_ms(iy, ix) = feasible;
        end
    end

    elapsed = toc(t0);
    n_safe_ms = sum(safe_mask_ms(:) & cone_mask_disc(:));
    n_cone_disc = sum(cone_mask_disc(:));
    frac_ms = 100 * n_safe_ms / max(n_cone_disc, 1);

    multistep_fracs(ia) = frac_ms;
    multistep_times(ia) = elapsed;
    faraci_results{ia}.safe_mask_ms = safe_mask_ms;
    faraci_results{ia}.frac_ms = frac_ms;
    faraci_results{ia}.time_ms = elapsed;

    fprintf('  a_max=%.2f: erosion=%.1f%%, multi-step=%.1f%%, time=%.1f s\n', ...
        a_max, faraci_results{ia}.frac, frac_ms, elapsed);
end

%% ===== Part 3: Uncertainty effect =====
% Faraci accounts for 10 deg orientation uncertainty via MEGB on SO(3).
% In our framework, orientation uncertainty effectively widens the apparent
% tumble rate envelope. We model this by computing the safe region at
% omega + delta_omega, where delta_omega accounts for the uncertainty.
%
% The angular uncertainty rho = 10 deg over dt=1s corresponds to an
% additional angular rate uncertainty of delta_omega = rho/dt = 10 deg/s.
% This is a worst-case bound; the actual effect depends on the inertia.
%
% More precisely: with principal-axis rotation and J known,
% omega_uncertainty = rho_SO3 / T_f ≈ 10/30 = 0.33 deg/s
% (attitude uncertainty spread over the time horizon)

fprintf('\n--- Uncertainty effect (Faraci rho=10 deg) ---\n');
omega_nominal = 4;  % deg/s
delta_omega_vals = [0, 0.33, 1.0, 3.0];  % deg/s uncertainty
a_max_ref = 0.10;

uncertainty_fracs = zeros(1, length(delta_omega_vals));
for id = 1:length(delta_omega_vals)
    omega_eff = omega_nominal + delta_omega_vals(id);
    omega_eff_rad = deg2rad(omega_eff);
    w2 = omega_eff_rad^2;
    r_sync_max = 2 * a_max_ref / max(w2, 1e-30);

    safe_mask = false(ny, nx);
    for iy = 1:ny
        for ix = 1:nx
            if ~cone_mask(iy, ix), continue; end
            xb = x_grid(ix); yb = y_grid(iy);
            rng = sqrt(xb^2 + yb^2);
            if rng >= r_sync_max, continue; end

            p_body = [xb; yb; 0; 0; 0; 0];
            slacks = b_body - A_body * p_body;
            v_rot = [omega_eff_rad * yb; -omega_eff_rad * xb; 0];

            safe_point = true;
            for ic = 1:size(A_body, 1)
                a_pos = A_body(ic, 1:3)';
                slack_rate = -a_pos' * v_rot;
                if slack_rate < 0
                    erosion = 0.5 * slack_rate^2 / a_max_ref;
                    if slacks(ic) < erosion
                        safe_point = false; break;
                    end
                end
            end
            safe_mask(iy, ix) = safe_point;
        end
    end
    n_safe = sum(safe_mask(:) & cone_mask(:));
    uncertainty_fracs(id) = 100 * n_safe / max(n_cone, 1);
    fprintf('  delta_omega=%.2f deg/s (eff=%.2f): %.1f%% safe\n', ...
        delta_omega_vals(id), omega_eff, uncertainty_fracs(id));
end

%% ===== Summary =====
fprintf('\n====== Benchmark Summary ======\n\n');

% Computation time comparison
our_total_time_ms = sum(erosion_times(:)) * 1000;
our_per_combo_ms = our_total_time_ms / (n_omega * n_amax);
fprintf('Computation time comparison:\n');
fprintf('  Faraci RG-OCP:       %d ms (single trajectory)\n', rgocp_time_ms);
fprintf('  Our erosion (per combo): %.0f ms (full 300x300 grid)\n', our_per_combo_ms);
fprintf('  Our total (%d combos): %.0f ms\n', n_omega*n_amax, our_total_time_ms);
fprintf('  Speedup factor:      %.0fx per combo\n\n', rgocp_time_ms / our_per_combo_ms);

% Coverage at Faraci's omega = 4 deg/s
fprintf('Safe-start coverage at omega=4 deg/s:\n');
fprintf('  %-10s  %-12s  %-12s\n', 'a_max', 'Erosion', 'Multi-step');
for ia = 1:n_amax
    fprintf('  %-10.2f  %-12.1f%%  %-12.1f%%\n', ...
        amax_vals(ia), faraci_results{ia}.frac, faraci_results{ia}.frac_ms);
end

% Scaling parameter
fprintf('\nScaling parameter a_max/omega^2 at omega=4 deg/s:\n');
omega_4_rad = deg2rad(4);
for ia = 1:n_amax
    scaling = amax_vals(ia) / omega_4_rad^2;
    fprintf('  a_max=%.2f: a_max/w^2 = %.1f m, r_sync = %.1f m\n', ...
        amax_vals(ia), scaling, 2*scaling);
end

%% ===== Figures =====
fprintf('\n--- Generating Figures ---\n');

green_col = [0.2, 0.75, 0.2];
blue_col = [0.3, 0.4, 0.85];
red_col = [0.9, 0.25, 0.25];
orange_col = [0.95, 0.55, 0.1];

% Figure 1: Safe-start maps at Faraci's omega=4 for different a_max
fig1 = figure('Position', [50 50 1600 400], 'Visible', 'off');
for ia = 1:n_amax
    ax = subplot(1, n_amax, ia);
    hold(ax, 'on');

    r = faraci_results{ia};
    img = ones(ny, nx, 3);
    for ch = 1:3
        layer = img(:,:,ch);
        layer(~cone_mask) = red_col(ch);
        layer(cone_mask & ~r.safe_mask) = blue_col(ch);
        layer(cone_mask & r.safe_mask) = green_col(ch);
        img(:,:,ch) = layer;
    end
    image(ax, x_grid, y_grid, img);
    set(ax, 'YDir', 'normal');

    y_cone = linspace(0, max(y_grid)*1.1, 200);
    plot(ax, rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5);
    plot(ax, -rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5);

    % Mark r_sync
    omega_4_rad = deg2rad(4);
    r_sync = 2 * amax_vals(ia) / omega_4_rad^2;
    if r_sync < 320
        theta = linspace(0, 2*pi, 200);
        plot(ax, r_sync*cos(theta), r_sync*sin(theta), 'w--', 'LineWidth', 1.5);
    end

    title(ax, sprintf('a_{max}=%.2f m/s^2\n%.1f%% safe, r_{sync}=%.0f m', ...
        amax_vals(ia), r.frac, r_sync), 'FontSize', 10);
    xlim(ax, [-200 200]); ylim(ax, [-10 320]);
    xlabel(ax, 'x_{TB} [m]');
    if ia == 1, ylabel(ax, 'y_{TB} [m]'); end
    set(ax, 'FontSize', 9);
end
sgtitle(sprintf('Safe-Start Regions at \\omega_t = 4 deg/s (Faraci benchmark scenario)'), ...
    'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(fig1, fullfile(results_dir, 'faraci_benchmark_maps.png'), 'Resolution', 200);
close(fig1);
fprintf('  Saved: faraci_benchmark_maps.png\n');

% Figure 2: Safe fraction vs scaling parameter for ALL combos
fig2 = figure('Position', [100 100 800 500], 'Visible', 'off');
hold on;

markers = {'o', 's', 'd', '^', 'v', 'p'};
colors = lines(n_omega);

for iw = 1:n_omega
    scaling_vals = amax_vals / deg2rad(omega_vals(iw))^2;
    plot(scaling_vals, safe_fractions_erosion(iw, :), ...
        ['-' markers{iw}], 'Color', colors(iw,:), 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(iw,:), ...
        'DisplayName', sprintf('\\omega=%d deg/s', omega_vals(iw)));
end

% Mark Faraci's operating point
omega_4_scaling = amax_vals / deg2rad(4)^2;
for ia = 1:n_amax
    if ia == 1
        plot(omega_4_scaling(ia), safe_fractions_erosion(4,ia), 'kp', ...
            'MarkerSize', 18, 'MarkerFaceColor', orange_col, 'LineWidth', 2, ...
            'DisplayName', 'Faraci scenario');
    else
        plot(omega_4_scaling(ia), safe_fractions_erosion(4,ia), 'kp', ...
            'MarkerSize', 18, 'MarkerFaceColor', orange_col, 'LineWidth', 2, ...
            'HandleVisibility', 'off');
    end
end

xlabel('Scaling parameter a_{max}/\omega_t^2 [m]', 'FontSize', 12, 'Interpreter', 'tex');
ylabel('Safe fraction [%]', 'FontSize', 12);
title('Safe-Start Coverage: Erosion-Based Certification', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
grid on;
set(gca, 'FontSize', 11);
exportgraphics(fig2, fullfile(results_dir, 'faraci_scaling_collapse.png'), 'Resolution', 200);
close(fig2);
fprintf('  Saved: faraci_scaling_collapse.png\n');

% Figure 3: Computation time comparison (bar chart)
fig3 = figure('Position', [100 100 700 450], 'Visible', 'off');

time_data = [our_per_combo_ms/1000; ...
             mean(multistep_times); ...
             rgocp_time_ms/1000];
time_labels = {'Our erosion\n(300x300 grid)', ...
               'Our multi-step\n(60x60, 30 steps)', ...
               'Faraci RG-OCP\n(single trajectory)'};

b = bar(time_data, 0.6);
b.FaceColor = 'flat';
b.CData = [green_col; blue_col; orange_col];
set(gca, 'XTickLabel', {'Erosion', 'Multi-step', 'RG-OCP'}, 'FontSize', 11);
ylabel('Computation time [s]', 'FontSize', 12);
title('Computation Time Comparison', 'FontSize', 13);
grid on;

% Add value labels on bars
for i = 1:3
    if time_data(i) < 1
        txt = sprintf('%.0f ms', time_data(i)*1000);
    else
        txt = sprintf('%.1f s', time_data(i));
    end
    text(i, time_data(i) + max(time_data)*0.03, txt, ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
end
ylim([0, max(time_data)*1.15]);

exportgraphics(fig3, fullfile(results_dir, 'faraci_timing_comparison.png'), 'Resolution', 200);
close(fig3);
fprintf('  Saved: faraci_timing_comparison.png\n');

% Figure 4: Uncertainty effect
fig4 = figure('Position', [100 100 700 450], 'Visible', 'off');
plot(delta_omega_vals, uncertainty_fracs, 'b-o', 'LineWidth', 2, ...
    'MarkerSize', 10, 'MarkerFaceColor', blue_col);
hold on;
% Mark Faraci's uncertainty
faraci_delta = rho_uncertainty_deg / T_f;  % ~0.33 deg/s
xline(faraci_delta, 'r--', 'LineWidth', 1.5);
text(faraci_delta + 0.1, max(uncertainty_fracs) - 0.2, ...
    sprintf('Faraci \\rho=10^\\circ\n(\\Delta\\omega=%.2f deg/s)', faraci_delta), ...
    'FontSize', 10, 'Color', 'r');

xlabel('\Delta\omega uncertainty [deg/s]', 'FontSize', 12, 'Interpreter', 'tex');
ylabel('Safe fraction [%]', 'FontSize', 12);
title(sprintf('Effect of Orientation Uncertainty on Safe Region\n(\\omega_t=4 deg/s, a_{max}=0.10 m/s^2)'), ...
    'FontSize', 12);
grid on;
set(gca, 'FontSize', 11);
exportgraphics(fig4, fullfile(results_dir, 'faraci_uncertainty_effect.png'), 'Resolution', 200);
close(fig4);
fprintf('  Saved: faraci_uncertainty_effect.png\n');

%% ===== Save =====
benchmark.safe_fractions_erosion = safe_fractions_erosion;
benchmark.erosion_times = erosion_times;
benchmark.multistep_fracs = multistep_fracs;
benchmark.multistep_times = multistep_times;
benchmark.uncertainty_fracs = uncertainty_fracs;
benchmark.omega_vals = omega_vals;
benchmark.amax_vals = amax_vals;
benchmark.delta_omega_vals = delta_omega_vals;
benchmark.faraci_results = faraci_results;
benchmark.rgocp_time_ms = rgocp_time_ms;

save(fullfile(results_dir, 'faraci_benchmark.mat'), 'benchmark');
fprintf('\nSaved: results/benchmark/faraci_benchmark.mat\n');

fprintf('\n====== Benchmark Complete ======\n');
