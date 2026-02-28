function rp = reachability_params()
%REACHABILITY_PARAMS  Parameters for reachability analysis.
%   rp = reachability_params()
%
%   Returns a struct with all parameters needed by the reachability modules.
%   Pulls orbital/cone parameters from the main params.m and adds
%   reachability-specific settings.

    % Import base simulation parameters (read-only, no modification)
    p = params();

    %% ---- Orbital mechanics ----
    rp.mu  = p.mu;
    rp.Re  = p.Re;
    rp.alt = p.alt;
    rp.a   = p.a;
    rp.n   = p.n;  % mean motion [rad/s]

    %% ---- LOS cone (body-fixed) ----
    rp.cone_half_angle_deg = p.cone_half_angle_deg;
    rp.cone_k = p.cone_k;
    rp.cone_nfaces = p.cone_nfaces;
    rp.y_min = p.y_min;

    %% ---- Control authority ----
    rp.u_max_default = p.u_max;  % from params.m

    %% ---- Time step ----
    rp.dt = p.dt;

    %% ---- Sweep parameters (matching MC) ----
    rp.omega_vals_deg = [1, 2, 3, 4, 5];      % deg/s
    rp.amax_vals      = [0.2, 0.1, 0.05, 0.02]; % m/s^2

    %% ---- Reachability-specific ----
    rp.N_steps   = 400;    % number of time steps to propagate
    rp.n_proj_dirs = 180;  % directions for polytope projection

    %% ---- Grid for point-wise evaluation ----
    rp.x_range = [-200, 200];  % body-frame x range [m]
    rp.y_range = [0, 320];     % body-frame y range [m]
    rp.nx_grid = 300;
    rp.ny_grid = 300;

    %% ---- Stochastic parameters ----
    rp.alpha = 0.05;  % violation probability (1-alpha = 95% confidence)
    % Process noise covariance (per-step disturbance in LVLH)
    rp.W_pos_std = 0.01;   % position disturbance std [m]
    rp.W_vel_std = 1e-4;   % velocity disturbance std [m/s]

    %% ---- Robust parameters ----
    rp.w_pos_max = 0.05;   % worst-case position disturbance bound [m]
    rp.w_vel_max = 5e-4;   % worst-case velocity disturbance bound [m/s]
end
