%RUN_MC_HIRES  High-resolution Monte Carlo (Part IV: resolution upgrade).
%
%   Increases grid density while keeping the SAME domain bounds as the
%   original run_monte_carlo.m. This script is a NEW file that does NOT
%   modify the original MC code.
%
%   Changes from original:
%     - n_y:        15 -> 25   (y-axis samples)
%     - n_x_per_y:  13 -> 21   (x-axis samples per y-level)
%     - Total:     195 -> 525  initial conditions per combo
%     - Random seed set for reproducibility
%     - Runtime logged
%
%   Saves to results/mc_hires_*.mat and results/mc_hires_sweep_all.mat
%
%   Run:  >> run_mc_hires

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'dynamics'));
addpath(fullfile(this_dir, 'frames'));
addpath(fullfile(this_dir, 'mpc'));
addpath(fullfile(this_dir, 'viz'));
addpath(fullfile(this_dir, 'utils'));

results_dir = fullfile(this_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% ===== Reproducibility =====
rng(2026, 'twister');

%% ===== Sweep parameters (SAME as original) =====
omega_vals_deg = [1, 2, 3, 4, 5];
amax_vals      = [0.2, 0.1, 0.05, 0.02];
n_omega = length(omega_vals_deg);
n_amax  = length(amax_vals);

%% ===== Sample initial conditions (HIGHER density, SAME bounds) =====
p0 = params();
cone_k = p0.cone_k;

% Higher resolution grid
n_y = 25;          % was 15
n_x_per_y = 21;    % was 13
y_vals = linspace(20, 300, n_y);  % SAME bounds

ic_list = [];
for iy = 1:length(y_vals)
    y0 = y_vals(iy);
    x_max = cone_k * y0 * 0.95;  % SAME margin
    x_pts = linspace(-x_max, x_max, n_x_per_y);
    for ix = 1:length(x_pts)
        ic_list = [ic_list; x_pts(ix), y0]; %#ok<AGROW>
    end
end
n_mc = size(ic_list, 1);
fprintf('Hi-Res MC Grid: %d omega x %d amax = %d combos, %d MC each\n', ...
    n_omega, n_amax, n_omega*n_amax, n_mc);

%% ===== Start parallel pool =====
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
    fprintf('Started parallel pool with %d workers\n', pool.NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n', pool.NumWorkers);
end

%% ===== Run all combos =====
all_data_hires = cell(n_omega, n_amax);

t_total = tic;
for iw = 1:n_omega
    for ia = 1:n_amax
        omega_deg = omega_vals_deg(iw);
        amax = amax_vals(ia);

        fprintf('\n====== omega=%.0f deg/s, a_max=%.2f m/s^2 ======\n', ...
            omega_deg, amax);

        p_mc = params();
        p_mc.omega_body = [0; 0; deg2rad(omega_deg)];
        p_mc.u_max = amax;
        p_mc.Rdu = diag([1e5, 1e4, 1e5]);
        p_mc.Np   = 20;
        p_mc.Tsim = 400;
        p_mc.ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
        p_mc.osqp_max_iter = 5000;

        outcomes = zeros(n_mc, 1);
        traj_x   = cell(n_mc, 1);
        traj_t   = cell(n_mc, 1);
        fail_x   = cell(n_mc, 1);
        term_reasons = cell(n_mc, 1);

        progress_queue = parallel.pool.DataQueue;
        t_start = tic;
        afterEach(progress_queue, @(d) report_progress(d, n_mc, t_start, omega_deg, amax));

        parfor ii = 1:n_mc
            x0 = ic_list(ii, 1);
            y0 = ic_list(ii, 2);

            p_run = p_mc;
            p_run.dr_lvlh0 = [x0; y0; 0];
            p_run.dv_lvlh0 = [0; 0; 0];
            p_run.y_hold_start = y0;
            p_run.y_hold_end   = 5;
            p_run.y_hold_tau   = y0;
            p_run.cone_draw_L  = 0;

            try
                [lg, ~, ~, reason] = run_sim_headless(p_run);
                traj_x{ii} = lg.r_tb_hist([1 2], :);
                traj_t{ii} = lg.t_hist;
                term_reasons{ii} = reason;
                if strcmp(reason, 'completed') || strcmp(reason, 'docked')
                    outcomes(ii) = 1;
                else
                    outcomes(ii) = 0;
                    fail_x{ii} = lg.r_tb_hist([1 2], end);
                end
            catch ME
                outcomes(ii) = 0;
                term_reasons{ii} = sprintf('error: %s', ME.message);
                traj_x{ii} = [x0; y0];
                traj_t{ii} = 0;
                fail_x{ii} = [x0; y0];
            end
            send(progress_queue, ii);
        end
        elapsed = toc(t_start);
        fprintf('  => %d/%d feasible (%.0f s)\n', sum(outcomes), n_mc, elapsed);

        mc.ic_list      = ic_list;
        mc.outcomes     = outcomes;
        mc.fail_x       = fail_x;
        mc.traj_x       = traj_x;
        mc.traj_t       = traj_t;
        mc.term_reasons = term_reasons;
        mc.cone_k       = cone_k;
        mc.p_mc         = p_mc;
        mc.n_y          = n_y;
        mc.n_x_per_y    = n_x_per_y;
        mc.y_vals       = y_vals;
        mc.omega_deg    = omega_deg;
        mc.amax         = amax;
        mc.elapsed_s    = elapsed;
        mc.seed         = 2026;

        all_data_hires{iw, ia} = mc;

        fname = sprintf('mc_hires_w%d_a%.2f.mat', omega_deg, amax);
        save(fullfile(results_dir, fname), 'mc');
        fprintf('  Saved: %s\n', fname);
    end
end
elapsed_total = toc(t_total);
fprintf('\n====== Hi-Res MC done in %.0f s ======\n', elapsed_total);

save(fullfile(results_dir, 'mc_hires_sweep_all.mat'), 'all_data_hires', ...
    'omega_vals_deg', 'amax_vals', 'ic_list', 'cone_k', 'elapsed_total');
fprintf('Saved: mc_hires_sweep_all.mat\n');


%% ========== Helper ==========
function report_progress(~, n_total, t_start, omega_deg, amax)
    persistent count
    if isempty(count), count = 0; end
    count = count + 1;
    if mod(count, 50) == 0 || count == n_total
        elapsed = toc(t_start);
        eta = elapsed / count * (n_total - count);
        fprintf('  [w=%d a=%.2f] %3d/%3d  %.0fs elapsed, ETA %.0fs\n', ...
            omega_deg, amax, count, n_total, elapsed, eta);
    end
    if count >= n_total, count = 0; end
end
