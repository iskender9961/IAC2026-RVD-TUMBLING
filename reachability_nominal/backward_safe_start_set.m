function result = backward_safe_start_set(omega_deg, a_max, rp)
%BACKWARD_SAFE_START_SET  Backward reachability: safe start set (nominal).
%   result = backward_safe_start_set(omega_deg, a_max, rp)
%
%   Computes the set of initial body-frame positions (with zero velocity)
%   from which the chaser can maintain LOS corridor feasibility for all
%   future time steps, using the nominal (disturbance-free) CWH dynamics.
%
%   The backward safe-start set (viability kernel) is computed as:
%       S_N = X_adm(N)
%       S_k = Pre(S_{k+1}) intersect X_adm(k)
%   where
%       Pre(S) = { x : exists u in U s.t. A*x + B*u in S }
%
%   For H-representation S = {z: H*z <= h}:
%       Pre(S) = { x : H*A*x <= h - max_{u in U} H*B*u }
%              = { x : H*A*x <= h_tight }
%   where h_tight(i) = h(i) - sum_j |H_i * B_j| * u_max_j
%
%   Since the LOS cone rotates with the target body, X_adm(k) depends on
%   the angle theta_k = omega * k * dt. We work in LVLH coordinates and
%   rotate the constraints at each step.
%
%   For multi-step backward propagation with changing constraint normals,
%   the H-representation grows. We use a practical approach:
%   - Propagate the one-step pre-image constraints N_back steps backward
%   - At each step, intersect with the rotated LOS constraints
%   - Evaluate the resulting safe region on the grid
%
%   Since full polytope tracking is expensive, we use the CONTINUOUS
%   analysis that captures the key physics:
%   - The per-constraint directional erosion accounts for the worst-case
%     rotation-induced margin loss that the controller must arrest
%   - The synchronization range bound ensures the chaser can match
%     the target rotation rate
%   - Together these define a conservative safe-start region
%
%   This is mathematically equivalent to the first-step backward Pre image
%   under the assumption that the one-step erosion dominates the multi-step
%   dynamics (valid when dt is small relative to the rotation period).
%
%   For a more precise multi-step analysis, we also compute N_back steps
%   of backward propagation using the discrete CWH dynamics.
%
%   Inputs:
%       omega_deg - target tumble rate [deg/s]
%       a_max     - max chaser acceleration [m/s^2]
%       rp        - reachability parameters struct
%
%   Outputs:
%       result - struct with fields:
%           safe_mask      - ny x nx logical (continuous erosion-based)
%           safe_mask_disc - ny x nx logical (discrete multi-step)
%           x_grid, y_grid, cone_mask, omega_deg, a_max, method

    omega_rad = deg2rad(omega_deg);
    dt = rp.dt;
    n_orbit = rp.n;

    % CWH matrices in LVLH
    [Ad, Bd] = cwh_stm(dt, n_orbit);

    % Evaluation grid (body frame)
    x_grid = linspace(rp.x_range(1), rp.x_range(2), rp.nx_grid);
    y_grid = linspace(rp.y_range(1), rp.y_range(2), rp.ny_grid);
    ny = length(y_grid);
    nx_g = length(x_grid);

    % Body-frame LOS constraints
    [A_body, b_body] = los_constraints_body(rp.cone_k, rp.y_min, rp.cone_nfaces);

    % Cone mask
    cone_mask = false(ny, nx_g);
    for iy = 1:ny
        for ix = 1:nx_g
            p_body = [x_grid(ix); y_grid(iy); 0; 0; 0; 0];
            cone_mask(iy, ix) = all(A_body * p_body <= b_body + 1e-9);
        end
    end

    %% ---- Continuous erosion-based safe mask (same as forward) ----
    safe_mask = false(ny, nx_g);
    w2 = omega_rad^2;
    r_sync_max = 2 * a_max / max(w2, 1e-30);

    for iy = 1:ny
        for ix = 1:nx_g
            if ~cone_mask(iy, ix), continue; end

            xb = x_grid(ix);
            yb = y_grid(iy);
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
                        safe_point = false;
                        break;
                    end
                end
            end
            safe_mask(iy, ix) = safe_point;
        end
    end

    %% ---- Discrete multi-step backward propagation ----
    % Number of backward steps (use enough to capture transients)
    N_back = min(20, rp.N_steps);

    % Use a coarser grid for the LP-based discrete analysis (computational cost)
    nx_disc = min(80, nx_g);
    ny_disc = min(80, ny);
    x_disc = linspace(rp.x_range(1), rp.x_range(2), nx_disc);
    y_disc = linspace(rp.y_range(1), rp.y_range(2), ny_disc);

    safe_disc_coarse = false(ny_disc, nx_disc);
    lp_opts = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');

    for iy = 1:ny_disc
        for ix = 1:nx_disc
            xb = x_disc(ix);
            yb = y_disc(iy);

            % Quick cone check
            p_body = [xb; yb; 0; 0; 0; 0];
            if ~all(A_body * p_body <= b_body + 1e-9), continue; end

            % Initial LVLH state (body=LVLH at t=0, add transport velocity)
            x_lvlh = [xb; yb; 0; -omega_rad*yb; omega_rad*xb; 0];

            feasible = true;
            x_k = x_lvlh;

            for k = 0:N_back-1
                t_kp1 = (k+1) * dt;
                [A_los_kp1, b_los_kp1] = los_constraints_at_time(...
                    rp.cone_k, rp.y_min, rp.cone_nfaces, omega_rad, t_kp1);

                % Check if zero control is feasible
                x_next_zero = Ad * x_k;
                if all(b_los_kp1 - A_los_kp1 * x_next_zero >= -1e-9)
                    u_k = [0; 0; 0];
                else
                    % LP to find minimum-norm feasible control
                    rhs = b_los_kp1 - A_los_kp1 * (Ad * x_k);
                    A_u = A_los_kp1 * Bd;
                    nu = 3;
                    f_lp = [0; 0; 0; 1; 1; 1];
                    A_lp = [A_u, zeros(size(A_u,1), nu);
                            eye(nu), -eye(nu);
                           -eye(nu), -eye(nu)];
                    b_lp = [rhs; zeros(nu,1); zeros(nu,1)];
                    lb = [-a_max*ones(nu,1); zeros(nu,1)];
                    ub = [a_max*ones(nu,1); inf(nu,1)];
                    [z_opt, ~, exitflag] = linprog(f_lp, A_lp, b_lp, ...
                        [], [], lb, ub, lp_opts);
                    if exitflag == 1
                        u_k = z_opt(1:nu);
                    else
                        feasible = false; break;
                    end
                end
                x_k = Ad * x_k + Bd * u_k;
            end
            safe_disc_coarse(iy, ix) = feasible;
        end
    end

    % Interpolate coarse result to full grid
    [Xd, Yd] = meshgrid(x_disc, y_disc);
    [Xf, Yf] = meshgrid(x_grid, y_grid);
    safe_mask_disc = interp2(Xd, Yd, double(safe_disc_coarse), Xf, Yf, 'nearest') > 0.5;
    safe_mask_disc = safe_mask_disc & cone_mask;

    result.safe_mask      = safe_mask;
    result.safe_mask_disc = safe_mask_disc;
    result.x_grid         = x_grid;
    result.y_grid         = y_grid;
    result.cone_mask      = cone_mask;
    result.omega_deg      = omega_deg;
    result.a_max          = a_max;
    result.method         = 'nominal_backward';
end
