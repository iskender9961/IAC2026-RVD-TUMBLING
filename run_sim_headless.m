function [lg, p, sim_terminated, term_reason] = run_sim_headless(p)
%RUN_SIM_HEADLESS  Run simulation with given params, no plots/GIFs.
%   [lg, p, sim_terminated, term_reason] = run_sim_headless(p)

    this_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(this_dir, 'dynamics'));
    addpath(fullfile(this_dir, 'frames'));
    addpath(fullfile(this_dir, 'mpc'));
    addpath(fullfile(this_dir, 'viz'));
    addpath(fullfile(this_dir, 'utils'));

    Nsteps = floor(p.Tsim / p.dt);
    term_reason = 'completed';

    % Target initial ECI state
    r_tgt_eci = p.r_tgt0_eci;
    v_tgt_eci = p.v_tgt0_eci;

    % LVLH frame at t=0
    R_eci_lvlh_0 = lvlh_from_rv(r_tgt_eci, v_tgt_eci);

    % Target body frame starts aligned with LVLH
    R_eci_tb = R_eci_lvlh_0;
    q_tb = quat_from_rotm_local(R_eci_tb);

    % Chaser initial ECI state
    r_chs_eci = r_tgt_eci + R_eci_lvlh_0 * p.dr_lvlh0;
    v_chs_eci = v_tgt_eci + R_eci_lvlh_0 * p.dv_lvlh0;

    % Initial relative state in target-body frame
    [r_tb, v_tb] = transform_rel_to_TB(r_chs_eci, v_chs_eci, ...
                                        r_tgt_eci, v_tgt_eci, ...
                                        R_eci_tb, p.omega_body);
    x_tb = [r_tb; v_tb];
    u_prev = zeros(3, 1);

    % Logging
    lg.r_tb_hist   = zeros(3, Nsteps+1);
    lg.v_tb_hist   = zeros(3, Nsteps+1);
    lg.r_lvlh_hist = zeros(3, Nsteps+1);
    lg.v_lvlh_hist = zeros(3, Nsteps+1);
    lg.u_hist      = zeros(3, Nsteps);
    lg.u_lvlh_hist = zeros(3, Nsteps);
    lg.t_hist      = zeros(1, Nsteps+1);
    lg.cost_hist   = zeros(1, Nsteps);
    lg.status_hist = cell(1, Nsteps);
    lg.R_eci_tb_hist   = zeros(3, 3, Nsteps+1);
    lg.R_eci_lvlh_hist = zeros(3, 3, Nsteps+1);
    lg.dock_axis_lvlh_hist = zeros(3, Nsteps+1);

    R_eci_lvlh = R_eci_lvlh_0;
    dr_eci = r_chs_eci - r_tgt_eci;
    r_lvlh = R_eci_lvlh' * dr_eci;
    v_lvlh = R_eci_lvlh' * (v_chs_eci - v_tgt_eci);

    lg.r_tb_hist(:,1)   = r_tb;
    lg.v_tb_hist(:,1)   = v_tb;
    lg.r_lvlh_hist(:,1) = r_lvlh;
    lg.v_lvlh_hist(:,1) = v_lvlh;
    lg.R_eci_tb_hist(:,:,1)   = R_eci_tb;
    lg.R_eci_lvlh_hist(:,:,1) = R_eci_lvlh;
    lg.dock_axis_lvlh_hist(:,1) = R_eci_lvlh' * R_eci_tb(:,2);
    lg.t_hist(1) = 0;

    % Set cone draw length
    if p.cone_draw_L <= 0
        p.cone_draw_L = r_tb(2) * 1.2;
    end

    % Build initial OSQP
    [Ad, Bd] = linearize_discrete_model(x_tb, R_eci_tb, ...
                    r_tgt_eci, v_tgt_eci, p.omega_body, ...
                    zeros(3,1), p.dt, p);
    x_ref = compute_ref(0, p);
    [prob, qp_data] = build_qp_osqp(Ad, Bd, x_tb, x_ref, u_prev, p);

    sim_terminated = false;
    final_step = Nsteps;

    for k = 1:Nsteps
        t_now = (k-1) * p.dt;

        % Check LOS
        [los_ok, ~] = check_los_fn(r_tb, p.cone_k, p.y_min);
        if ~los_ok
            final_step = k - 1;
            sim_terminated = true;
            term_reason = 'LOS_violation';
            break;
        end

        x_ref = compute_ref(t_now, p);

        % Linearize
        [Ad, Bd] = linearize_discrete_model(x_tb, R_eci_tb, ...
                        r_tgt_eci, v_tgt_eci, p.omega_body, ...
                        u_prev, p.dt, p);

        % Solve
        [u_opt, status, ~] = update_qp_osqp(prob, qp_data, Ad, Bd, ...
                                x_tb, x_ref, u_prev, p);

        lg.status_hist{k} = status;

        if ~contains(status, 'solved')
            final_step = k - 1;
            sim_terminated = true;
            term_reason = sprintf('infeasible:%s', status);
            break;
        end

        % Cost
        du = u_opt - u_prev;
        J_k = (x_tb - x_ref)' * p.Q * (x_tb - x_ref) + ...
              u_opt' * p.Ru * u_opt + du' * p.Rdu * du;
        lg.cost_hist(k) = J_k;

        u_applied = clamp(u_opt, -p.u_max, p.u_max);
        lg.u_hist(:, k) = u_applied;

        % Input in LVLH frame
        a_ctrl_eci = control_TB_to_ECI(u_applied, R_eci_tb);
        lg.u_lvlh_hist(:, k) = R_eci_lvlh' * a_ctrl_eci;

        % Propagate
        x_tgt_next = propagate_truth_step([r_tgt_eci; v_tgt_eci], ...
                        p.dt, p.mu, p.Re, p.J2, [0;0;0], p.ode_opts);
        r_tgt_eci = x_tgt_next(1:3);
        v_tgt_eci = x_tgt_next(4:6);

        x_chs_next = propagate_truth_step([r_chs_eci; v_chs_eci], ...
                        p.dt, p.mu, p.Re, p.J2, a_ctrl_eci, p.ode_opts);
        r_chs_eci = x_chs_next(1:3);
        v_chs_eci = x_chs_next(4:6);

        [q_tb, R_eci_tb] = target_attitude_model(q_tb, p.omega_body, p.dt);
        R_eci_lvlh = lvlh_from_rv(r_tgt_eci, v_tgt_eci);

        [r_tb, v_tb] = transform_rel_to_TB(r_chs_eci, v_chs_eci, ...
                                            r_tgt_eci, v_tgt_eci, ...
                                            R_eci_tb, p.omega_body);
        x_tb = [r_tb; v_tb];

        dr_eci = r_chs_eci - r_tgt_eci;
        r_lvlh = R_eci_lvlh' * dr_eci;
        v_lvlh = R_eci_lvlh' * (v_chs_eci - v_tgt_eci);

        % Check docking
        rad_dev = sqrt(r_tb(1)^2 + r_tb(3)^2);
        if r_tb(2) <= p.dock_y_thresh && ...
           rad_dev <= p.dock_rad_thresh && ...
           norm(v_tb) <= p.dock_v_thresh
            final_step = k;
            sim_terminated = true;
            term_reason = 'docked';
        end

        % Log
        lg.r_tb_hist(:, k+1)   = r_tb;
        lg.v_tb_hist(:, k+1)   = v_tb;
        lg.r_lvlh_hist(:, k+1) = r_lvlh;
        lg.v_lvlh_hist(:, k+1) = v_lvlh;
        lg.R_eci_tb_hist(:,:,k+1)   = R_eci_tb;
        lg.R_eci_lvlh_hist(:,:,k+1) = R_eci_lvlh;
        lg.dock_axis_lvlh_hist(:,k+1) = R_eci_lvlh' * R_eci_tb(:,2);
        lg.t_hist(k+1) = k * p.dt;

        u_prev = u_applied;

        if sim_terminated, break; end

        % Print progress
        if mod(k, 50) == 0
            fprintf('  t=%5.0fs y=%.1f rdev=%.2f |u|=%.3f J=%.1f %s\n', ...
                t_now, r_tb(2), rad_dev, norm(u_applied), J_k, status);
        end
    end

    % Trim
    N = final_step;
    lg.r_tb_hist   = lg.r_tb_hist(:, 1:N+1);
    lg.v_tb_hist   = lg.v_tb_hist(:, 1:N+1);
    lg.r_lvlh_hist = lg.r_lvlh_hist(:, 1:N+1);
    lg.v_lvlh_hist = lg.v_lvlh_hist(:, 1:N+1);
    lg.u_hist      = lg.u_hist(:, 1:N);
    lg.u_lvlh_hist = lg.u_lvlh_hist(:, 1:N);
    lg.t_hist      = lg.t_hist(1:N+1);
    lg.cost_hist   = lg.cost_hist(1:N);
    lg.status_hist = lg.status_hist(1:N);
    lg.R_eci_tb_hist   = lg.R_eci_tb_hist(:,:,1:N+1);
    lg.R_eci_lvlh_hist = lg.R_eci_lvlh_hist(:,:,1:N+1);
    lg.dock_axis_lvlh_hist = lg.dock_axis_lvlh_hist(:,1:N+1);
end


function x_ref = compute_ref(t, p)
    y_hold = p.y_hold_end + (p.y_hold_start - p.y_hold_end) * exp(-t / p.y_hold_tau);
    x_ref = [0; y_hold; 0; 0; 0; 0];
end

function [ok, margin] = check_los_fn(r_tb, cone_k, y_min)
    xT = r_tb(1); yT = r_tb(2); zT = r_tb(3);
    rad = sqrt(xT^2 + zT^2);
    if yT < y_min
        ok = false; margin = yT - y_min;
    elseif rad > cone_k * yT
        ok = false; margin = cone_k * yT - rad;
    else
        ok = true; margin = cone_k * yT - rad;
    end
end

function q = quat_from_rotm_local(R)
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
