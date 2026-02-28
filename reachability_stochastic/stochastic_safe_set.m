function result = stochastic_safe_set(omega_deg, a_max, rp)
%STOCHASTIC_SAFE_SET  Stochastic reachability with chance constraints.
%   result = stochastic_safe_set(omega_deg, a_max, rp)
%
%   Extends the nominal erosion-based analysis with Gaussian disturbance
%   tightening. The stochastic safe set is a subset of the nominal set,
%   providing probabilistic guarantees.
%
%   Model: x_{k+1} = A*x_k + B*u_k + w_k
%   where w_k ~ N(0, W) is i.i.d. process noise.
%
%   For chance constraints P(A_c * x_k <= b_c) >= 1 - alpha,
%   we tighten each constraint by the Gaussian quantile:
%       b_tight(i) = b(i) - Phi^{-1}(1-alpha/n_c) * sigma_i
%   where sigma_i = sqrt(a_i' * Sigma_k * a_i) is the standard deviation
%   of the constraint slack, and Sigma_k is the state covariance at step k.
%
%   The one-step disturbance tightening for the erosion model gives:
%       erosion_stoch(i) = erosion_nom(i) + z_{1-alpha/n_c} * sigma_w_i
%   where sigma_w_i is the disturbance's contribution to constraint i.
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

    % Process noise covariance (per step, in body frame)
    W = diag([rp.W_pos_std^2, rp.W_pos_std^2, rp.W_pos_std^2, ...
              rp.W_vel_std^2, rp.W_vel_std^2, rp.W_vel_std^2]);

    % Accumulated covariance after N steps of noise propagation
    % Sigma_{k+1} = A * Sigma_k * A' + W
    % For the erosion model, we use steady-state Sigma (DARE solution)
    % or the N-step accumulated covariance.
    % For a conservative one-step analysis, use W directly.
    %
    % Multi-step: Sigma_N = sum_{j=0}^{N-1} A^j W (A^j)'
    % For the safe region, we use N_eff steps of accumulated noise
    N_eff = min(50, rp.N_steps);
    Sigma_accum = zeros(6);
    Aj = eye(6);
    for j = 0:N_eff-1
        Sigma_accum = Sigma_accum + Aj * W * Aj';
        Aj = Ad * Aj;
    end

    % Quantile for Bonferroni-corrected chance constraints
    % P(all constraints satisfied) >= 1 - alpha
    % => P(constraint i violated) <= alpha / n_con
    alpha = rp.alpha;
    z_quantile = norminv(1 - alpha / n_con);

    % Per-constraint standard deviation from accumulated noise
    sigma_per_con = zeros(n_con, 1);
    for i = 1:n_con
        ai = A_body(i,:)';
        sigma_per_con(i) = sqrt(ai' * Sigma_accum * ai);
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

    % Compute stochastic safe mask
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

                % Stochastic tightening: add z * sigma
                erosion_stoch = erosion_nom + z_quantile * sigma_per_con(ic);

                if slacks(ic) < erosion_stoch
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
    result.alpha     = alpha;
    result.z_quantile = z_quantile;
    result.sigma_per_con = sigma_per_con;
    result.method    = 'stochastic';
end
