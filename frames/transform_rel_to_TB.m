function [r_tb, v_tb] = transform_rel_to_TB(r_chaser_eci, v_chaser_eci, ...
                                              r_target_eci, v_target_eci, ...
                                              R_eci_tb, omega_body)
%TRANSFORM_REL_TO_TB  Compute relative position/velocity in target body frame.
%   [r_tb, v_tb] = transform_rel_to_TB(r_c, v_c, r_t, v_t, R_eci_tb, omega_body)
%
%   r_tb = R_tb_eci * (r_c - r_t)
%   v_tb = R_tb_eci * (v_c - v_t) - omega_body x r_tb   (transport theorem)
%
%   omega_body : 3x1 angular velocity of TB in TB frame [rad/s]
    R_tb_eci = R_eci_tb';    % transpose = inverse for rotation matrix
    dr_eci = r_chaser_eci - r_target_eci;
    dv_eci = v_chaser_eci - v_target_eci;
    r_tb = R_tb_eci * dr_eci;
    v_tb = R_tb_eci * dv_eci - cross(omega_body, r_tb);
end
