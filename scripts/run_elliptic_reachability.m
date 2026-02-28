%RUN_ELLIPTIC_REACHABILITY  Compare reachability sets: circular vs elliptic orbits.
%
%   Computes forward reachable safe sets using:
%     - CWH (e=0): time-invariant STM (baseline)
%     - YA (e=0.1): mildly elliptic, periapsis vs apoapsis
%     - YA (e=0.3): moderately elliptic, periapsis vs apoapsis
%
%   For elliptic orbits, the safe set depends on the orbital phase (true
%   anomaly nu0) because the YA STM is time-varying. We evaluate at
%   periapsis (nu=0) and apoapsis (nu=pi) to show the range.
%
%   Produces:
%     - Comparison figures: CWH vs YA at different eccentricities
%     - Summary table of safe fractions
%     - Saved to results/elliptic/
%
%   Run:  >> run_elliptic_reachability

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

% Orbit parameters
a_km = (rp.Re + rp.alt) / 1e3;  % semi-major axis [km]
mu   = rp.mu;

% Eccentricity cases
ecc_vals = [0, 0.1, 0.3];
ecc_labels = {'e=0 (CWH)', 'e=0.1 (YA)', 'e=0.3 (YA)'};

% Orbital phase cases for elliptic
nu_vals = [0, pi];
nu_labels = {'periapsis', 'apoapsis'};

% Representative sweep: 2 tumble rates x 2 thrust levels
omega_vals = [2, 4];       % deg/s
amax_vals  = [0.10, 0.05]; % m/s^2

n_omega = length(omega_vals);
n_amax  = length(amax_vals);
n_ecc   = length(ecc_vals);
n_nu    = length(nu_vals);

fprintf('====== Elliptic Orbit Reachability Comparison ======\n');
fprintf('Eccentricities: %s\n', mat2str(ecc_vals));
fprintf('Tumble rates: %s deg/s\n', mat2str(omega_vals));
fprintf('Thrust levels: %s m/s^2\n', mat2str(amax_vals));
fprintf('Grid: %d x %d\n\n', rp.nx_grid, rp.ny_grid);

%% ===== Compute safe sets =====
% Store results: {ecc_idx, nu_idx, omega_idx, amax_idx}
all_results = cell(n_ecc, n_nu, n_omega, n_amax);
safe_fractions = zeros(n_ecc, n_nu, n_omega, n_amax);

for ie = 1:n_ecc
    ecc = ecc_vals(ie);

    for inu = 1:n_nu
        nu0 = nu_vals(inu);

        % Skip apoapsis for circular (identical to periapsis)
        if ecc == 0 && inu > 1
            for iw = 1:n_omega
                for ia = 1:n_amax
                    all_results{ie, inu, iw, ia} = all_results{ie, 1, iw, ia};
                    safe_fractions(ie, inu, iw, ia) = safe_fractions(ie, 1, iw, ia);
                end
            end
            continue;
        end

        fprintf('--- %s, %s ---\n', ecc_labels{ie}, nu_labels{inu});

        for iw = 1:n_omega
            for ia = 1:n_amax
                omega_deg = omega_vals(iw);
                a_max = amax_vals(ia);
                omega_rad = deg2rad(omega_deg);

                fprintf('  omega=%d, a_max=%.2f ... ', omega_deg, a_max);
                t0 = tic;

                % Build evaluation grid
                x_grid = linspace(rp.x_range(1), rp.x_range(2), rp.nx_grid);
                y_grid = linspace(rp.y_range(1), rp.y_range(2), rp.ny_grid);
                ny = length(y_grid);
                nx = length(x_grid);

                % LOS constraints
                [A_body, b_body] = los_constraints_body(rp.cone_k, rp.y_min, rp.cone_nfaces);

                % Cone mask
                cone_mask = false(ny, nx);
                for iy = 1:ny
                    for ix = 1:nx
                        p_body = [x_grid(ix); y_grid(iy); 0; 0; 0; 0];
                        cone_mask(iy, ix) = all(A_body * p_body <= b_body + 1e-9);
                    end
                end

                % Compute safe mask using erosion model
                % For elliptic orbits, the dynamics coupling modifies the
                % effective drift rate. However, the erosion model is a
                % continuous-time approximation that depends primarily on:
                %   1. Target rotation rate (omega)
                %   2. Chaser thrust authority (a_max)
                %   3. Range to target
                % The orbital dynamics enter through the CWH/YA secular drift.
                %
                % For the erosion analysis, we augment the rotation-induced
                % drift with the orbital coupling contribution from the STM.
                %
                % The key quantity is the one-step drift of a zero-velocity
                % body-frame point through the dynamics.

                safe_mask = false(ny, nx);
                w2 = omega_rad^2;
                r_sync_max = 2 * a_max / max(w2, 1e-30);

                if ecc == 0
                    % CWH: time-invariant, use standard erosion
                    [Ad, Bd] = cwh_stm(rp.dt, rp.n);
                else
                    % YA: get STM at this orbital phase
                    cfg = elliptic_orbit_config(ecc, a_km, mu);
                    [Ad, Bd, ~] = ya_stm(rp.dt, cfg, nu0);
                end

                for iy = 1:ny
                    for ix = 1:nx
                        if ~cone_mask(iy, ix), continue; end

                        xb = x_grid(ix);
                        yb = y_grid(iy);
                        rng = sqrt(xb^2 + yb^2);

                        if rng >= r_sync_max, continue; end

                        p_body = [xb; yb; 0; 0; 0; 0];
                        slacks = b_body - A_body * p_body;

                        % Apparent body-frame velocity from target rotation
                        v_rot = [omega_rad * yb; -omega_rad * xb; 0];

                        % For elliptic orbits, also account for orbital
                        % coupling drift. The CWH/YA dynamics cause a
                        % secular drift even without rotation. We add the
                        % one-step position drift from the STM:
                        %   delta_r = (Ad(1:3,1:3) - I) * r
                        % This is typically small for short dt but grows
                        % with eccentricity.
                        if ecc > 0
                            dyn_drift = (Ad(1:3,1:3) - eye(3)) * [xb; yb; 0] / rp.dt;
                            v_eff = v_rot + dyn_drift;
                        else
                            v_eff = v_rot;
                        end

                        safe_point = true;
                        for ic = 1:size(A_body, 1)
                            a_pos = A_body(ic, 1:3)';
                            slack_rate = -a_pos' * v_eff;
                            if slack_rate < 0
                                erosion = 0.5 * slack_rate^2 / a_max;
                                if slacks(ic) < erosion
                                    safe_point = false;
                                    break;
                                end
                            end
                        end

                        safe_mask(iy, ix) = safe_point;
                    end
                end

                n_safe = sum(safe_mask(:) & cone_mask(:));
                n_cone = sum(cone_mask(:));
                frac = 100 * n_safe / max(n_cone, 1);
                elapsed = toc(t0);

                fprintf('%.1f%% safe, %.1f s\n', frac, elapsed);

                result.safe_mask = safe_mask;
                result.cone_mask = cone_mask;
                result.x_grid = x_grid;
                result.y_grid = y_grid;
                result.omega_deg = omega_deg;
                result.a_max = a_max;
                result.ecc = ecc;
                result.nu0 = nu0;

                all_results{ie, inu, iw, ia} = result;
                safe_fractions(ie, inu, iw, ia) = frac;
            end
        end
    end
end

%% ===== Summary Table =====
fprintf('\n====== Safe Fraction Summary (%%  of LOS cone) ======\n');
fprintf('%-20s', 'Scenario');
for iw = 1:n_omega
    for ia = 1:n_amax
        fprintf('  w=%d,a=%.2f', omega_vals(iw), amax_vals(ia));
    end
end
fprintf('\n%s\n', repmat('-', 1, 20 + n_omega*n_amax*12));

for ie = 1:n_ecc
    for inu = 1:n_nu
        if ecc_vals(ie) == 0 && inu > 1, continue; end
        if ecc_vals(ie) == 0
            label = ecc_labels{ie};
        else
            label = sprintf('%s (%s)', ecc_labels{ie}, nu_labels{inu});
        end
        fprintf('%-20s', label);
        for iw = 1:n_omega
            for ia = 1:n_amax
                fprintf('  %10.1f', safe_fractions(ie, inu, iw, ia));
            end
        end
        fprintf('\n');
    end
end

%% ===== Figures =====
fprintf('\n--- Generating Figures ---\n');

% Color palette
green_col = [0.2, 0.75, 0.2];
blue_col = [0.3, 0.4, 0.85];
red_col = [0.9, 0.25, 0.25];
cone_bg = [0.92, 0.92, 0.92];

% Figure 1: Side-by-side comparison for each (omega, amax)
for iw = 1:n_omega
    for ia = 1:n_amax
        fig = figure('Position', [50 50 1800 500], 'Visible', 'off');

        % Subplot indices: CWH | YA(0.1,peri) | YA(0.1,apo) | YA(0.3,peri) | YA(0.3,apo)
        panel_configs = {
            {1, 1}, {2, 1}, {2, 2}, {3, 1}, {3, 2}
        };
        panel_titles = {
            'CWH (e=0)', 'YA e=0.1 (peri)', 'YA e=0.1 (apo)', ...
            'YA e=0.3 (peri)', 'YA e=0.3 (apo)'
        };

        for ip = 1:5
            ie = panel_configs{ip}{1};
            inu = panel_configs{ip}{2};
            r = all_results{ie, inu, iw, ia};

            ax = subplot(1, 5, ip);
            hold(ax, 'on');

            ny = length(r.y_grid); nx = length(r.x_grid);
            img = ones(ny, nx, 3);
            for ch = 1:3
                layer = img(:,:,ch);
                layer(~r.cone_mask) = red_col(ch);
                layer(r.cone_mask & ~r.safe_mask) = blue_col(ch);
                layer(r.cone_mask & r.safe_mask) = green_col(ch);
                img(:,:,ch) = layer;
            end
            image(ax, r.x_grid, r.y_grid, img);
            set(ax, 'YDir', 'normal');

            y_cone = linspace(0, max(r.y_grid)*1.1, 200);
            plot(ax, rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1);
            plot(ax, -rp.cone_k*y_cone, y_cone, 'k-', 'LineWidth', 1);

            frac = safe_fractions(ie, inu, iw, ia);
            title(ax, sprintf('%s\n%.1f%% safe', panel_titles{ip}, frac), 'FontSize', 9);
            xlim(ax, [-200 200]); ylim(ax, [-10 320]);
            if ip == 1, ylabel(ax, 'y_{TB} [m]'); end
            xlabel(ax, 'x_{TB} [m]');
            set(ax, 'FontSize', 8);
        end

        sgtitle(sprintf('Elliptic Orbit Effect: \\omega=%d deg/s, a_{max}=%.2f m/s^2', ...
            omega_vals(iw), amax_vals(ia)), 'FontSize', 12, 'FontWeight', 'bold');

        fname = sprintf('elliptic_compare_w%d_a%.2f.png', omega_vals(iw), amax_vals(ia));
        exportgraphics(fig, fullfile(results_dir, fname), 'Resolution', 200);
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end
end

% Figure 2: Bar chart of safe fractions
fig_bar = figure('Position', [100 100 900 500], 'Visible', 'off');
scenarios = {};
data_matrix = [];
for ie = 1:n_ecc
    for inu = 1:n_nu
        if ecc_vals(ie) == 0 && inu > 1, continue; end
        if ecc_vals(ie) == 0
            scenarios{end+1} = 'CWH';
        else
            scenarios{end+1} = sprintf('e=%.1f %s', ecc_vals(ie), nu_labels{inu}(1:4));
        end
        row = [];
        for iw = 1:n_omega
            for ia = 1:n_amax
                row(end+1) = safe_fractions(ie, inu, iw, ia);
            end
        end
        data_matrix(end+1, :) = row;
    end
end

bar(data_matrix');
set(gca, 'XTickLabel', {});
xtick_labels = {};
for iw = 1:n_omega
    for ia = 1:n_amax
        xtick_labels{end+1} = sprintf('w=%d,a=%.2f', omega_vals(iw), amax_vals(ia));
    end
end
set(gca, 'XTickLabel', xtick_labels, 'FontSize', 9);
xtickangle(30);
ylabel('Safe fraction [%]', 'FontSize', 11);
legend(scenarios, 'Location', 'northeast', 'FontSize', 9);
title('Effect of Orbital Eccentricity on Safe-Start Region', 'FontSize', 12);
grid on;

exportgraphics(fig_bar, fullfile(results_dir, 'elliptic_bar_chart.png'), 'Resolution', 200);
close(fig_bar);
fprintf('  Saved: elliptic_bar_chart.png\n');

%% ===== Save results =====
save(fullfile(results_dir, 'elliptic_reachability.mat'), ...
    'all_results', 'safe_fractions', 'ecc_vals', 'nu_vals', ...
    'omega_vals', 'amax_vals', 'rp');
fprintf('\nSaved: results/elliptic/elliptic_reachability.mat\n');

fprintf('\n====== Elliptic Reachability Complete ======\n');
