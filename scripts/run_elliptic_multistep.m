%RUN_ELLIPTIC_MULTISTEP  Multi-step backward reachability: CWH vs YA.
%
%   Computes multi-step backward safe-start sets using discrete LP-based
%   analysis. This is where elliptic effects actually show up, because
%   the YA STM varies with orbital phase over many steps.
%
%   For elliptic orbits, the STM (Ad, Bd) changes at each step as the
%   true anomaly advances. This creates a time-varying linear system that
%   can produce different safe regions compared to CWH.
%
%   Run:  >> run_elliptic_multistep

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
root_dir = fileparts(this_dir);
if isempty(root_dir), root_dir = pwd; end

addpath(fullfile(root_dir, 'dynamics'));
addpath(fullfile(root_dir, 'reachability_common'));
addpath(fullfile(root_dir, 'reachability_nominal'));
addpath(fullfile(root_dir, 'Elliptic_orbit'));
addpath(root_dir);

results_dir = fullfile(root_dir, 'results', 'elliptic');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% ===== Parameters =====
rp = reachability_params();
a_km = (rp.Re + rp.alt) / 1e3;
mu   = rp.mu;
dt   = rp.dt;

% Cases
ecc_vals = [0, 0.1, 0.3];
ecc_labels = {'CWH (e=0)', 'YA (e=0.1)', 'YA (e=0.3)'};
nu0_start = 0;  % start at periapsis

% Representative scenario
omega_deg = 2;
a_max = 0.10;
omega_rad = deg2rad(omega_deg);

% Backward steps (more steps = more elliptic effect)
N_back = 50;

% Use a coarser grid for the LP-based analysis (speed)
nx_disc = 60;
ny_disc = 60;

fprintf('====== Multi-Step Backward Reachability: CWH vs YA ======\n');
fprintf('omega = %d deg/s, a_max = %.2f m/s^2, N_back = %d\n', ...
    omega_deg, a_max, N_back);
fprintf('Grid: %d x %d\n\n', nx_disc, ny_disc);

%% ===== LOS constraints =====
[A_body, b_body] = los_constraints_body(rp.cone_k, rp.y_min, rp.cone_nfaces);

x_disc = linspace(rp.x_range(1), rp.x_range(2), nx_disc);
y_disc = linspace(rp.y_range(1), rp.y_range(2), ny_disc);

lp_opts = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');

%% ===== Run for each eccentricity =====
all_safe = cell(1, length(ecc_vals));
all_frac = zeros(1, length(ecc_vals));

for ie = 1:length(ecc_vals)
    ecc = ecc_vals(ie);
    fprintf('--- %s ---\n', ecc_labels{ie});

    % Pre-compute STM sequence for all N_back steps
    Ad_seq = zeros(6, 6, N_back);
    Bd_seq = zeros(6, 3, N_back);

    if ecc == 0
        [Ad0, Bd0] = cwh_stm(dt, rp.n);
        for k = 1:N_back
            Ad_seq(:,:,k) = Ad0;
            Bd_seq(:,:,k) = Bd0;
        end
    else
        cfg = elliptic_orbit_config(ecc, a_km, mu);
        nu_current = nu0_start;
        for k = 1:N_back
            [Ad_k, Bd_k, nu_next] = ya_stm(dt, cfg, nu_current);
            Ad_seq(:,:,k) = Ad_k;
            Bd_seq(:,:,k) = Bd_k;
            nu_current = nu_next;
        end
        fprintf('  True anomaly after %d steps: %.1f deg\n', N_back, rad2deg(nu_current));
    end

    % Evaluate each grid point
    safe_mask = false(ny_disc, nx_disc);
    t0 = tic;

    for iy = 1:ny_disc
        for ix = 1:nx_disc
            xb = x_disc(ix);
            yb = y_disc(iy);

            % Quick cone check
            p_body = [xb; yb; 0; 0; 0; 0];
            if ~all(A_body * p_body <= b_body + 1e-9), continue; end

            % Initial LVLH state (body = LVLH at t=0)
            x_lvlh = [xb; yb; 0; -omega_rad*yb; omega_rad*xb; 0];

            feasible = true;
            x_k = x_lvlh;

            for k = 1:N_back
                t_k = k * dt;
                [A_los_k, b_los_k] = los_constraints_at_time(...
                    rp.cone_k, rp.y_min, rp.cone_nfaces, omega_rad, t_k);

                Ad_k = Ad_seq(:,:,k);
                Bd_k = Bd_seq(:,:,k);

                % Check if zero control is feasible
                x_next_zero = Ad_k * x_k;
                if all(b_los_k - A_los_k * x_next_zero >= -1e-9)
                    u_k = [0; 0; 0];
                else
                    % LP: find minimum-norm feasible control
                    rhs = b_los_k - A_los_k * (Ad_k * x_k);
                    A_u = A_los_k * Bd_k;
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
                x_k = Ad_k * x_k + Bd_k * u_k;
            end
            safe_mask(iy, ix) = feasible;
        end
    end
    elapsed = toc(t0);

    % Cone mask
    cone_mask = false(ny_disc, nx_disc);
    for iy = 1:ny_disc
        for ix = 1:nx_disc
            p_body = [x_disc(ix); y_disc(iy); 0; 0; 0; 0];
            cone_mask(iy, ix) = all(A_body * p_body <= b_body + 1e-9);
        end
    end

    n_safe = sum(safe_mask(:) & cone_mask(:));
    n_cone = sum(cone_mask(:));
    frac = 100 * n_safe / max(n_cone, 1);

    fprintf('  Safe: %d/%d = %.1f%%, Time: %.1f s\n', n_safe, n_cone, frac, elapsed);

    all_safe{ie} = safe_mask;
    all_frac(ie) = frac;
end

%% ===== Summary =====
fprintf('\n====== Multi-Step Safe Fraction (N=%d steps) ======\n', N_back);
for ie = 1:length(ecc_vals)
    fprintf('  %-20s: %.1f%%\n', ecc_labels{ie}, all_frac(ie));
end

%% ===== Figure: Side-by-side comparison =====
green_col = [0.2, 0.75, 0.2];
blue_col = [0.3, 0.4, 0.85];
red_col = [0.9, 0.25, 0.25];

fig = figure('Position', [50 50 1500 450], 'Visible', 'off');
for ie = 1:length(ecc_vals)
    ax = subplot(1, 3, ie);
    hold(ax, 'on');

    ny = ny_disc; nx = nx_disc;
    img = ones(ny, nx, 3);
    for ch = 1:3
        layer = img(:,:,ch);
        layer(~cone_mask) = red_col(ch);
        layer(cone_mask & ~all_safe{ie}) = blue_col(ch);
        layer(cone_mask & all_safe{ie}) = green_col(ch);
        img(:,:,ch) = layer;
    end
    image(ax, x_disc, y_disc, img);
    set(ax, 'YDir', 'normal');

    y_cone = linspace(0, max(y_disc)*1.1, 200);
    plot(ax, rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5);
    plot(ax, -rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1.5);

    title(ax, sprintf('%s\n%.1f%% safe', ecc_labels{ie}, all_frac(ie)), 'FontSize', 11);
    xlim(ax, [-200 200]); ylim(ax, [-10 320]);
    xlabel(ax, 'x_{TB} [m]');
    if ie == 1, ylabel(ax, 'y_{TB} [m]'); end
    set(ax, 'FontSize', 10);
end

sgtitle(sprintf('Multi-Step Backward Safe Set (N=%d): \\omega=%d deg/s, a_{max}=%.2f m/s^2', ...
    N_back, omega_deg, a_max), 'FontSize', 13, 'FontWeight', 'bold');

exportgraphics(fig, fullfile(results_dir, 'elliptic_multistep_compare.png'), ...
    'Resolution', 200);
close(fig);
fprintf('\nSaved: elliptic_multistep_compare.png\n');

%% ===== Save =====
save(fullfile(results_dir, 'elliptic_multistep.mat'), ...
    'all_safe', 'all_frac', 'ecc_vals', 'cone_mask', ...
    'x_disc', 'y_disc', 'omega_deg', 'a_max', 'N_back', 'rp');
fprintf('Saved: elliptic_multistep.mat\n');

fprintf('\n====== Multi-Step Analysis Complete ======\n');
