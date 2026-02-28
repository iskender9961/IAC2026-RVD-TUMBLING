function [A_body, B_body] = body_frame_transform(Ad_lvlh, Bd_lvlh, omega_z, dt)
%BODY_FRAME_TRANSFORM  Transform LVLH discrete dynamics to body frame.
%   [A_body, B_body] = body_frame_transform(Ad_lvlh, Bd_lvlh, omega_z, dt)
%
%   For a target rotating at omega_z about z-axis:
%       r_body = Rz(-theta) * r_lvlh
%       v_body = Rz(-theta) * v_lvlh - omega x r_body
%
%   The state in body frame x_B = T(theta) * x_L where T is the
%   6x6 transformation. The body-frame discrete dynamics become:
%       x_B(k+1) = T(theta_{k+1}) * Ad * T(theta_k)^{-1} * x_B(k)
%                + T(theta_{k+1}) * Bd * u
%
%   Since theta_{k+1} = theta_k + omega_z*dt, we compute:
%       A_body = T(omega_z*dt) * [T(0)^{-1} composed with Ad]
%
%   For reachability analysis with rotating constraints, we work in LVLH
%   and rotate the constraints instead. This function returns the
%   body-to-LVLH and LVLH-to-body transformations for a given angle.
%
%   Inputs:
%       Ad_lvlh  - 6x6 LVLH state transition matrix
%       Bd_lvlh  - 6x3 LVLH input matrix
%       omega_z  - body rotation rate about z [rad/s]
%       dt       - time step [s]
%
%   Outputs:
%       A_body   - 6x6 body-frame state transition matrix
%       B_body   - 6x3 body-frame input matrix

    dtheta = omega_z * dt;

    % T(theta) maps LVLH to body frame:
    %   r_B = Rz(-theta) * r_L
    %   v_B = Rz(-theta) * v_L - omega x r_B
    % T_inv(theta) maps body to LVLH:
    %   r_L = Rz(theta) * r_B
    %   v_L = Rz(theta) * (v_B + omega x r_B)

    T_0     = build_T(0, omega_z);
    T_dt    = build_T(dtheta, omega_z);
    T_0_inv = build_T_inv(0, omega_z);

    A_body = T_dt * Ad_lvlh * T_0_inv;
    B_body = T_dt * Bd_lvlh;
end


function T = build_T(theta, omega_z)
%BUILD_T  LVLH-to-body transformation matrix (6x6).
    Rz = rotz_mat(-theta);
    Omega_cross = [0, -omega_z, 0; omega_z, 0, 0; 0, 0, 0];
    T = [Rz, zeros(3); -Omega_cross*Rz, Rz];
end


function Ti = build_T_inv(theta, omega_z)
%BUILD_T_INV  Body-to-LVLH transformation matrix (6x6).
    Rz = rotz_mat(theta);
    Omega_cross = [0, -omega_z, 0; omega_z, 0, 0; 0, 0, 0];
    Ti = [Rz, zeros(3); Rz*Omega_cross, Rz];
end


function R = rotz_mat(theta)
%ROTZ_MAT  Rotation about z-axis.
    c = cos(theta);
    s = sin(theta);
    R = [c, -s, 0; s, c, 0; 0, 0, 1];
end
