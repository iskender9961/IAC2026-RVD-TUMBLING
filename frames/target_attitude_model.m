function [q_tb, R_eci_tb] = target_attitude_model(q_tb_prev, omega_body, dt)
%TARGET_ATTITUDE_MODEL  Update target body attitude for one step.
%   [q_tb, R_eci_tb] = target_attitude_model(q_tb_prev, omega_body, dt)
%   q_tb_prev  : 4x1 previous quaternion (ECI-from-TB, scalar-first)
%   omega_body : 3x1 constant body angular velocity [rad/s]
%   dt         : time step [s]
%   q_tb       : 4x1 updated quaternion
%   R_eci_tb   : 3x3 rotation matrix  v_eci = R_eci_tb * v_tb
    q_tb = quat_propagate(q_tb_prev, omega_body, dt);
    R_eci_tb = rotm_from_quat(q_tb);
end
