function [A_c, b_c] = los_constraints_body(cone_k, y_min, nfaces)
%LOS_CONSTRAINTS_BODY  LOS polyhedral cone constraints in body frame.
%   [A_c, b_c] = los_constraints_body(cone_k, y_min, nfaces)
%
%   Returns A_c, b_c such that x is inside the LOS cone iff A_c*x <= b_c,
%   where x = [x_B; y_B; z_B; vx; vy; vz] is the body-frame state.
%
%   The LOS cone (body-fixed, axis = +yB):
%       sqrt(x_B^2 + z_B^2) <= cone_k * y_B   and   y_B >= y_min
%
%   Polyhedral approximation with nfaces half-spaces:
%       cos(theta_i)*x_B + sin(theta_i)*z_B - cone_k*y_B <= 0
%       -y_B <= -y_min
%
%   Inputs:
%       cone_k  - tan(half_angle), slope of LOS cone
%       y_min   - minimum y_B distance (floor)
%       nfaces  - number of polyhedral faces
%
%   Outputs:
%       A_c - (nfaces+1) x 6 constraint matrix
%       b_c - (nfaces+1) x 1 constraint RHS

    n_con = nfaces + 1;
    A_c = zeros(n_con, 6);
    b_c = zeros(n_con, 1);

    thetas = 2*pi*(0:nfaces-1)'/nfaces;

    for f = 1:nfaces
        ct = cos(thetas(f));
        st = sin(thetas(f));
        A_c(f, 1) = ct;       % x_B coefficient
        A_c(f, 2) = -cone_k;  % y_B coefficient
        A_c(f, 3) = st;       % z_B coefficient
        b_c(f) = 0;
    end

    % Floor constraint: -y_B <= -y_min  =>  y_B >= y_min
    A_c(n_con, 2) = -1;
    b_c(n_con) = -y_min;
end
