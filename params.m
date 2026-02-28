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
    %  Place target on x-axis, velocity along y-axis (circular equatorial)
    p.r_tgt0_eci = [p.a; 0; 0];                         % [m]
    p.v_tgt0_eci = [0; sqrt(p.mu / p.a); 0];            % [m/s]

    %% ---- Chaser initial relative state (LVLH) ----
    %  x=radial, y=in-track (along-track), z=cross-track  (meters, m/s)
    p.dr_lvlh0 = [0; 200; 0];       % relative position in LVLH [m]
    p.dv_lvlh0 = [0;   0; 0];       % relative velocity in LVLH [m/s]

    %% ---- Target tumble model ----
    %  Constant body-rate vector in target-body frame [rad/s]
    %  TB is aligned with LVLH at t=0, then tumbles at this rate.
    p.omega_body = [0; 0; 0.02];   % ~1.15 deg/s about body-z

    %% ---- Simulation timing ----
    p.Tsim = 300;               % total simulation time [s]
    p.dt   = 1.0;               % control / MPC step    [s]

    %% ---- MPC parameters ----
    p.Np   = 30;                % prediction horizon (steps)
    p.nx   = 6;                 % states  [r_TB; v_TB]
    p.nu   = 3;                 % inputs  (accel in TB frame)
    p.u_max = 0.1;              % max accel per axis [m/s^2]

    % Cost weights
    p.Q  = diag([1, 1, 1, 0.1, 0.1, 0.1]);    % stage state weight
    p.QN = 10 * p.Q;                            % terminal state weight
    p.Rdu = 1e4 * eye(3);                       % delta-u weight (smoothness)
    p.Ru  = 1e0 * eye(3);                       % small u penalty for regularization

    %% ---- Reference / approach strategy ----
    %  Hold-point on +yT axis. The reference y component shrinks over time
    %  from y_hold_start toward y_hold_end.
    p.y_hold_start = 200;      % initial hold distance on +yT [m]
    p.y_hold_end   = 5;        % final approach distance [m]
    p.y_hold_tau   = 150;      % exponential time constant [s]

    %% ---- LOS cone constraint (body-fixed, +yT axis) ----
    p.cone_half_angle_deg = 30;                          % half-angle [deg]
    p.cone_k = tan(deg2rad(p.cone_half_angle_deg));      % tan(alpha)
    p.cone_nfaces = 8;            % number of half-space faces for polyhedral approx
    p.y_min = 1.0;                % minimum yT (approach corridor floor) [m]

    %% ---- Integrator tolerances ----
    p.ode_opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

    %% ---- Visualization ----
    p.save_mp4    = false;
    p.mp4_file    = 'rvd_approach.mp4';
    p.cube_size   = 3;           % half-edge of drawn cubes [m]
    p.triad_len   = 20;          % length of drawn frame triads [m]
    p.cone_draw_L = 80;          % length of drawn cone [m]

    %% ---- OSQP settings ----
    p.osqp_max_iter  = 4000;
    p.osqp_eps_abs   = 1e-5;
    p.osqp_eps_rel   = 1e-5;
    p.osqp_warm_start = true;
    p.osqp_verbose   = false;
end
