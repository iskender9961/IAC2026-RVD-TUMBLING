function result = robust_safe_set(omega_deg, a_max, rp)
%ROBUST_SAFE_SET  Robust reachability with bounded disturbances.
%   result = robust_safe_set(omega_deg, a_max, rp)
%
%   Extends the nominal erosion-based analysis with worst-case disturbance
%   tightening. The robust safe set is a subset of the stochastic set,
%   providing deterministic guarantees for bounded disturbances.
%
%   Model: x_{k+1} = A*x_k + B*u_k + w_k
%   where ||w_k||_inf <= w_max (bounded disturbance).
%
%   For robust constraint satisfaction for ALL w in W:
%       A_c * (A*x + B*u + w) <= b_c  for all w in W
%   We tighten each constraint by the worst-case disturbance:
%       b_tight(i) = b(i) - max_{w in W} A_c(i,:) * w
%                  = b(i) - sum_j |A_c(i,j)| * w_max_j
%
%   This is the support function of W in the direction A_c(i,:).
%
%   For the erosion model, the robust tightening adds:
%       erosion_robust(i) = erosion_nom(i) + support_W(A_c(i,:)) * N_eff
%   where N_eff accounts for accumulated worst-case disturbance over
%   multiple steps.
%
%   Inputs:
%       omega_deg - target tumble rate [deg/s]
%       a_max     - max chaser acceleration [m/s^2]
%       rp        - reachability parameters struct
%
%   Outputs:
%       result - struct with safe_mask, x_grid, y_grid, cone_mask, etc.

    omega_rad = deg2rad(omega_deg);
    dt = rp.dt;
    n_orbit = rp.n;

    % CWH matrices
    [Ad, Bd] = cwh_stm(dt, n_orbit);

    % Body-frame LOS constraints
    [A_body, b_body] = los_constraints_body(rp.cone_k, rp.y_min, rp.cone_nfaces);
    n_con = size(A_body, 1);

    % Bounded disturbance: w_max per state component
    w_max = [rp.w_pos_max; rp.w_pos_max; rp.w_pos_max; ...
             rp.w_vel_max; rp.w_vel_max; rp.w_vel_max];

    % Accumulated worst-case disturbance over N_eff steps
    % At step k, total accumulated disturbance through A propagation:
    %   d_total = sum_{j=0}^{k-1} A^j * w_j
    % Worst case: max |A_c(i,:) * d_total| = sum_{j=0}^{k-1} support(A_c(i,:) * A^j, W)
    N_eff = min(50, rp.N_steps);

    % Compute per-constraint accumulated worst-case tightening
    robust_tighten = zeros(n_con, 1);
    Aj = eye(6);
    for j = 0:N_eff-1
        for i = 1:n_con
            ai_Aj = A_body(i,:) * Aj;  % 1x6
            support_val = sum(abs(ai_Aj(:)) .* w_max(:));
            robust_tighten(i) = robust_tighten(i) + support_val;
        end
        Aj = Ad * Aj;
    end

    % Evaluation grid
    x_grid = linspace(rp.x_range(1), rp.x_range(2), rp.nx_grid);
    y_grid = linspace(rp.y_range(1), rp.y_range(2), rp.ny_grid);
    ny = length(y_grid);
    nx_g = length(x_grid);

    % Cone mask
    cone_mask = false(ny, nx_g);
    for iy = 1:ny
        for ix = 1:nx_g
            p_body = [x_grid(ix); y_grid(iy); 0; 0; 0; 0];
            cone_mask(iy, ix) = all(A_body * p_body <= b_body + 1e-9);
        end
    end

    % Compute robust safe mask
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

            % Nominal erosion
            v_rot = [omega_rad * yb; -omega_rad * xb; 0];
            safe_point = true;
            for ic = 1:n_con
                a_pos = A_body(ic, 1:3)';
                slack_rate = -a_pos' * v_rot;
                erosion_nom = 0;
                if slack_rate < 0
                    erosion_nom = 0.5 * slack_rate^2 / a_max;
                end

                % Robust tightening: add worst-case accumulated disturbance
                erosion_robust = erosion_nom + robust_tighten(ic);

                if slacks(ic) < erosion_robust
                    safe_point = false;
                    break;
                end
            end

            safe_mask(iy, ix) = safe_point;
        end
    end

    result.safe_mask  = safe_mask;
    result.x_grid    = x_grid;
    result.y_grid    = y_grid;
    result.cone_mask = cone_mask;
    result.omega_deg = omega_deg;
    result.a_max     = a_max;
    result.robust_tighten = robust_tighten;
    result.method    = 'robust';
end
