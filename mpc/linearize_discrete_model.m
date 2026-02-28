function [Ad, Bd] = linearize_discrete_model(x_tb, R_eci_tb, ...
                        r_target_eci, v_target_eci, omega_body, ...
                        u_nom_tb, dt, p)
%LINEARIZE_DISCRETE_MODEL  Numerical linearization via finite differences.
%   [Ad, Bd] = linearize_discrete_model(x_tb, R_eci_tb, r_t, v_t, omega_body, u_nom, dt, p)
%
%   Produces discrete-time A_d, B_d such that
%       x_{k+1} ~ Ad * x_k + Bd * u_k    (in target-body frame)
%   by perturbing the nonlinear propagation around the current state and nominal input.
%
%   The nonlinear propagation proceeds as:
%     1) Convert (x_tb, u_tb) to ECI chaser state + control.
%     2) Propagate chaser in ECI for dt.
%     3) Propagate target in ECI for dt (with zero control).
%     4) Propagate target attitude for dt.
%     5) Transform the new states back to TB frame.
%
%   This gives us the nonlinear map f(x_tb, u_tb) -> x_tb_next.

    nx = 6;
    nu = 3;
    eps_x = 1e-6;   % perturbation for state
    eps_u = 1e-6;   % perturbation for input

    % Nominal next state
    x_nom_next = propagate_one_step_TB(x_tb, u_nom_tb, R_eci_tb, ...
                    r_target_eci, v_target_eci, omega_body, dt, p);

    % ---- Compute Ad via forward finite differences ----
    Ad = zeros(nx, nx);
    for j = 1:nx
        x_pert = x_tb;
        x_pert(j) = x_pert(j) + eps_x;
        x_next_pert = propagate_one_step_TB(x_pert, u_nom_tb, R_eci_tb, ...
                        r_target_eci, v_target_eci, omega_body, dt, p);
        Ad(:, j) = (x_next_pert - x_nom_next) / eps_x;
    end

    % ---- Compute Bd via forward finite differences ----
    Bd = zeros(nx, nu);
    for j = 1:nu
        u_pert = u_nom_tb;
        u_pert(j) = u_pert(j) + eps_u;
        x_next_pert = propagate_one_step_TB(x_tb, u_pert, R_eci_tb, ...
                        r_target_eci, v_target_eci, omega_body, dt, p);
        Bd(:, j) = (x_next_pert - x_nom_next) / eps_u;
    end
end


function x_tb_next = propagate_one_step_TB(x_tb, u_tb, R_eci_tb, ...
                        r_target_eci, v_target_eci, omega_body, dt, p)
%PROPAGATE_ONE_STEP_TB  Full nonlinear propagation mapped into TB frame.
    % Recover chaser ECI state from TB relative state + target ECI state
    R_tb_eci = R_eci_tb';
    r_rel_eci = R_eci_tb * x_tb(1:3);
    % v_tb = R_tb_eci * dv_eci - omega x r_tb  =>  dv_eci = R_eci_tb * (v_tb + omega x r_tb)
    v_rel_eci = R_eci_tb * (x_tb(4:6) + cross(omega_body, x_tb(1:3)));
    r_chaser_eci = r_target_eci + r_rel_eci;
    v_chaser_eci = v_target_eci + v_rel_eci;

    % Control in ECI
    a_ctrl_eci = R_eci_tb * u_tb;

    % Propagate both spacecraft
    x_tgt_next = propagate_truth_step([r_target_eci; v_target_eci], ...
                    dt, p.mu, p.Re, p.J2, [0;0;0], p.ode_opts);
    x_chs_next = propagate_truth_step([r_chaser_eci; v_chaser_eci], ...
                    dt, p.mu, p.Re, p.J2, a_ctrl_eci, p.ode_opts);

    % Propagate target attitude
    q_tb = quat_from_rotm(R_eci_tb);
    [~, R_eci_tb_next] = target_attitude_model(q_tb, omega_body, dt);

    % Transform back to TB
    [r_tb_next, v_tb_next] = transform_rel_to_TB( ...
        x_chs_next(1:3), x_chs_next(4:6), ...
        x_tgt_next(1:3), x_tgt_next(4:6), ...
        R_eci_tb_next, omega_body);

    x_tb_next = [r_tb_next; v_tb_next];
end


function q = quat_from_rotm(R)
%QUAT_FROM_ROTM  Extract scalar-first quaternion from rotation matrix.
    tr = trace(R);
    if tr > 0
        s = 0.5 / sqrt(tr + 1);
        qw = 0.25 / s;
        qx = (R(3,2) - R(2,3)) * s;
        qy = (R(1,3) - R(3,1)) * s;
        qz = (R(2,1) - R(1,2)) * s;
    elseif R(1,1) > R(2,2) && R(1,1) > R(3,3)
        s = 2 * sqrt(1 + R(1,1) - R(2,2) - R(3,3));
        qw = (R(3,2) - R(2,3)) / s;
        qx = 0.25 * s;
        qy = (R(1,2) + R(2,1)) / s;
        qz = (R(1,3) + R(3,1)) / s;
    elseif R(2,2) > R(3,3)
        s = 2 * sqrt(1 + R(2,2) - R(1,1) - R(3,3));
        qw = (R(1,3) - R(3,1)) / s;
        qx = (R(1,2) + R(2,1)) / s;
        qy = 0.25 * s;
        qz = (R(2,3) + R(3,2)) / s;
    else
        s = 2 * sqrt(1 + R(3,3) - R(1,1) - R(2,2));
        qw = (R(2,1) - R(1,2)) / s;
        qx = (R(1,3) + R(3,1)) / s;
        qy = (R(2,3) + R(3,2)) / s;
        qz = 0.25 * s;
    end
    q = [qw; qx; qy; qz];
    q = q / norm(q);
end
