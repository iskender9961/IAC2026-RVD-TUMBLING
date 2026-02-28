function [A_c_lvlh, b_c] = los_constraints_at_time(cone_k, y_min, nfaces, omega_z, t)
%LOS_CONSTRAINTS_AT_TIME  LOS cone constraints rotated to LVLH at time t.
%   [A_c_lvlh, b_c] = los_constraints_at_time(cone_k, y_min, nfaces, omega_z, t)
%
%   The LOS cone is defined in body frame. At time t, the body frame is
%   rotated by theta = omega_z * t from LVLH. The constraints in LVLH are:
%       A_body * T(theta) * x_lvlh <= b_c
%   where T(theta) maps LVLH position to body frame (rotation Rz(-theta)).
%
%   Since constraints only involve position (first 3 states), we construct:
%       A_c_lvlh = A_c_body * T_pos
%   where T_pos is the 6x6 transformation that rotates positions.
%
%   Inputs:
%       cone_k   - tan(half_angle)
%       y_min    - minimum y_B distance
%       nfaces   - number of polyhedral faces
%       omega_z  - body rotation rate [rad/s]
%       t        - current time [s]
%
%   Outputs:
%       A_c_lvlh - (nfaces+1) x 6 constraint matrix in LVLH
%       b_c      - (nfaces+1) x 1 constraint RHS

    % Get body-frame constraints
    [A_c_body, b_c] = los_constraints_body(cone_k, y_min, nfaces);

    % Rotation angle
    theta = omega_z * t;

    % Rotation matrix: r_body = Rz(-theta) * r_lvlh
    c = cos(theta);
    s = sin(theta);
    Rz_neg = [c, s, 0; -s, c, 0; 0, 0, 1];  % Rz(-theta)

    % Build full 6x6 position-only rotation (velocity cols are zero in A_c_body)
    T_full = blkdiag(Rz_neg, Rz_neg);

    % Transform: A_c_lvlh * x_lvlh <= b_c
    A_c_lvlh = A_c_body * T_full;
end
