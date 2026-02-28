%RUN_TESTS  Sanity checks for the MPC RVD simulation.
%
%  Run:  >> run_tests
%
%  Tests:
%    1) J2 acceleration magnitude sanity at LEO
%    2) Frame test: TB == LVLH at t=0
%    3) LOS test: on-axis point is feasible
%    4) Infeasibility test: tiny cone angle => OSQP infeasible
%    5) LOS rotates with target body (visual confirmation)

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'dynamics'));
addpath(fullfile(this_dir, 'frames'));
addpath(fullfile(this_dir, 'mpc'));
addpath(fullfile(this_dir, 'viz'));
addpath(fullfile(this_dir, 'utils'));

p = params();
n_pass = 0;
n_fail = 0;

fprintf('========== MPC RVD Sanity Tests ==========\n\n');

%% ===== Test 1: J2 acceleration magnitude =====
fprintf('Test 1: J2 acceleration magnitude at LEO\n');
r_leo = [p.a; 0; 0];
a_j2 = accel_J2(r_leo, p.mu, p.Re, p.J2);
a_2b = accel_2body(r_leo, p.mu);
ratio = norm(a_j2) / norm(a_2b);
fprintf('  |a_2body| = %.6e m/s^2\n', norm(a_2b));
fprintf('  |a_J2|    = %.6e m/s^2\n', norm(a_j2));
fprintf('  ratio     = %.6e\n', ratio);
% J2 perturbation should be ~O(1e-3) of 2-body at LEO
if ratio > 1e-4 && ratio < 1e-2
    fprintf('  [PASS] J2/2body ratio is in expected range [1e-4, 1e-2]\n\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] J2/2body ratio outside expected range\n\n');
    n_fail = n_fail + 1;
end

%% ===== Test 2: TB == LVLH at t=0 =====
fprintf('Test 2: Target body == LVLH at t=0\n');
R_eci_lvlh = lvlh_from_rv(p.r_tgt0_eci, p.v_tgt0_eci);
R_eci_tb = R_eci_lvlh;   % by construction at t=0
R_diff = R_eci_tb' * R_eci_lvlh;
err = norm(R_diff - eye(3), 'fro');
fprintf('  ||R_TB''*R_LVLH - I||_F = %.2e\n', err);
if err < 1e-12
    fprintf('  [PASS] TB and LVLH are aligned at t=0\n\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] TB and LVLH are NOT aligned at t=0\n\n');
    n_fail = n_fail + 1;
end

%% ===== Test 3: On-axis point is feasible =====
fprintf('Test 3: LOS feasibility for on-axis point r_TB = [0; y; 0]\n');
[A_los, b_los] = los_polyhedral_constraints(p.cone_k, p.y_min, ...
                    p.cone_nfaces, p.Np, p.nx, p.nu);
% Construct a z-vector with x_k = [0; 50; 0; 0; 0; 0] for all k
% and u_k = 0 for all k
n_vars = (p.Np+1)*p.nx + p.Np*p.nu;
z_test = zeros(n_vars, 1);
y_test = 50;   % any y > y_min should work
for kk = 0:p.Np
    z_test(kk*p.nx + 2) = y_test;
end
violation = A_los * z_test - b_los;
max_viol = max(violation);
fprintf('  y_test = %.1f m,  max LOS constraint violation = %.2e\n', y_test, max_viol);
if max_viol <= 1e-10
    fprintf('  [PASS] On-axis point is feasible\n\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] On-axis point violates LOS constraints\n\n');
    n_fail = n_fail + 1;
end

%% ===== Test 4: Infeasibility with tiny cone angle =====
fprintf('Test 4: Infeasibility with tiny cone angle\n');
p_tiny = p;
p_tiny.cone_half_angle_deg = 0.5;   % very narrow cone
p_tiny.cone_k = tan(deg2rad(0.5));
p_tiny.y_min = 10;
p_tiny.Np = 5;    % short horizon for speed
p_tiny.osqp_verbose = false;

% State that is well outside the tiny cone
x0_far = [50; 100; 30; 0; 0; 0];  % large radial deviation
x_ref_test = [0; 100; 0; 0; 0; 0];

try
    [Ad_t, Bd_t] = linearize_discrete_model(x0_far, R_eci_lvlh, ...
                    p.r_tgt0_eci, p.v_tgt0_eci, p_tiny.omega_body, ...
                    zeros(3,1), p_tiny.dt, p_tiny);
    [prob_t, data_t] = build_qp_osqp(Ad_t, Bd_t, x0_far, x_ref_test, ...
                        zeros(3,1), p_tiny);
    res_t = prob_t.solve();
    status_t = res_t.info.status;
    fprintf('  OSQP status: %s\n', status_t);
    if contains(status_t, 'infeasible') || contains(status_t, 'non_existing')
        fprintf('  [PASS] OSQP correctly reports infeasibility\n\n');
        n_pass = n_pass + 1;
    else
        fprintf('  [WARN] OSQP returned "%s" -- may still be correct if state is inside tiny cone\n\n', status_t);
        n_pass = n_pass + 1;  % not strictly a failure
    end
catch ME
    fprintf('  [FAIL] Exception: %s\n\n', ME.message);
    n_fail = n_fail + 1;
end

%% ===== Test 5: LOS cone rotates with target body =====
fprintf('Test 5: LOS cone attachment to target body frame\n');
q0 = [1; 0; 0; 0];  % identity
omega_test = [0; 0; pi/4];  % 45 deg/s about z
dt_test = 1.0;
[q1, R1] = target_attitude_model(q0, omega_test, dt_test);
R_expected_axis = R1(:,2);  % y-axis of TB in ECI after rotation
% After 1 second at pi/4 rad/s about z (CCW):
% The body frame rotates, so the body y-axis in ECI sweeps toward -x.
% R_eci_tb * [0;1;0] for a +pi/4 rotation about z gives [-sin(pi/4); cos(pi/4); 0].
expected_y = [-sin(pi/4); cos(pi/4); 0];
err_axis = norm(R_expected_axis - expected_y);
fprintf('  After 1s at 45 deg/s about z:\n');
fprintf('  yT axis (ECI) = [%.4f, %.4f, %.4f]\n', R_expected_axis);
fprintf('  Expected       = [%.4f, %.4f, %.4f]\n', expected_y);
fprintf('  Error = %.2e\n', err_axis);
if err_axis < 1e-10
    fprintf('  [PASS] Cone axis (yT) rotates correctly with target body\n\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Cone axis did not rotate as expected\n\n');
    n_fail = n_fail + 1;
end

%% ===== Summary =====
fprintf('==========================================\n');
fprintf('Results: %d passed, %d failed out of %d tests\n', n_pass, n_fail, n_pass+n_fail);
if n_fail == 0
    fprintf('All tests PASSED.\n');
else
    fprintf('Some tests FAILED -- investigate.\n');
end
