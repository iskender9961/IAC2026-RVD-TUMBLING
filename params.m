function p = params()
%PARAMS  All simulation parameters in one place.
%   p = params() returns a struct with fields grouped by category.

    %% ---- Earth constants (SI) ----
    p.mu  = 3.986004418e14;      % gravitational parameter  [m^3/s^2]
    p.Re  = 6378137.0;           % mean equatorial radius    [m]
    p.J2  = 1.08263e-3;          % J2 zonal harmonic         [-]

    %% ---- Target orbit (circular LEO) ----
    p.alt  = 500e3;              % altitude [m]
    p.a    = p.Re + p.alt;       % semi-major axis [m]
    p.n    = sqrt(p.mu / p.a^3); % mean motion [rad/s]

    %% ---- Target initial ECI state ----
    p.r_tgt0_eci = [p.a; 0; 0];
    p.v_tgt0_eci = [0; sqrt(p.mu / p.a); 0];

    %% ---- Chaser initial relative state (LVLH) ----
    p.dr_lvlh0 = [0; 200; 0];       % [m]
    p.dv_lvlh0 = [0;   0; 0];       % [m/s]
    %  NOTE: transport term v_TB = dv_LVLH - omega x r_LVLH gives
    %  v_TB_0 ~ [3;0;0] m/s lateral drift.  Controller must sync.

    %% ---- Target tumble model ----
    p.omega_body = [0; 0; 0.015];  % 0.86 deg/s about body-z [rad/s]

    %% ---- Simulation timing ----
    p.Tsim = 400;               % total simulation time [s]
    p.dt   = 1.0;               % control step [s]

    %% ---- MPC parameters ----
    p.Np   = 40;                % prediction horizon (steps)
    p.nx   = 6;                 % states  [r_TB; v_TB]
    p.nu   = 3;                 % inputs  (accel in TB frame)
    p.u_max = 0.15;             % max accel per axis [m/s^2]

    % Cost weights  --  TUNING RATIONALE:
    %   Heavy Q on x_TB, z_TB drives chaser onto docking axis.
    %   Light Q on y_TB lets reference set approach rate.
    %   Moderate velocity penalty damps oscillation.
    %   Rdu moderate: allows maneuver without chatter.
    %   Ru tiny: numerical conditioning only.
    p.Q  = diag([15, 1, 15, 1, 1, 1]);     % stage state
    p.QN = 30 * p.Q;                        % terminal
    p.Rdu = 3e2 * eye(3);                   % delta-u (smoothness)
    p.Ru  = 1e-2 * eye(3);                  % regularization

    %% ---- Reference / approach strategy ----
    p.y_hold_start = 200;      % initial hold distance on +yT [m]
    p.y_hold_end   = 5;        % final approach distance [m]
    p.y_hold_tau   = 200;      % exponential time constant [s]

    %% ---- LOS tetrahedral corridor (body-fixed, +yT axis) ----
    p.cone_half_angle_deg = 30;
    p.cone_k = tan(deg2rad(p.cone_half_angle_deg));
    p.cone_nfaces = 8;            % polyhedral faces
    p.y_min = 1.0;                % corridor floor [m]

    %% ---- Integrator tolerances ----
    p.ode_opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

    %% ---- Visualization ----
    p.save_mp4    = false;
    p.mp4_file    = 'rvd_approach.mp4';
    p.cube_size   = 3;           % half-edge of drawn cubes [m]
    p.triad_len   = 20;          % triad axis length [m]
    p.cone_draw_L = 0;           % 0 = set at runtime from initial y_TB

    %% ---- OSQP settings ----
    p.osqp_max_iter  = 20000;
    p.osqp_eps_abs   = 1e-5;
    p.osqp_eps_rel   = 1e-5;
    p.osqp_warm_start = true;
    p.osqp_verbose   = false;
end
