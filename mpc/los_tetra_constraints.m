function [A_los, b_los] = los_tetra_constraints(cone_k, y_min, nfaces, Np, nx, nu)
%LOS_TETRA_CONSTRAINTS  Build polyhedral LOS corridor constraints for QP.
%   [A_los, b_los] = los_tetra_constraints(cone_k, y_min, nfaces, Np, nx, nu)
%
%   The LOS cone (body-fixed, axis = +yT) is:
%       sqrt(xT^2 + zT^2) <= cone_k * yT   and   yT >= y_min
%
%   Polyhedral approximation with nfaces half-spaces:
%       For i = 1..nfaces:
%           cos(theta_i)*xT + sin(theta_i)*zT <= cone_k * yT
%       where theta_i = 2*pi*(i-1)/nfaces
%       Plus:  -yT <= -y_min   (ensure yT >= y_min)
%
%   The constraints are applied at each of the (Np+1) state positions
%   in the QP decision variable vector.
%
%   QP decision layout:  z = [x_0; x_1; ... x_Np; u_0; u_1; ... u_{Np-1}]
%   Total states: (Np+1)*nx,  total inputs: Np*nu
%
%   Returns A_los, b_los such that  A_los * z <= b_los.

    n_constraints_per_step = nfaces + 1;   % nfaces cone + 1 y_min
    total_constraints = n_constraints_per_step * (Np + 1);
    total_vars = (Np+1)*nx + Np*nu;

    % Pre-compute cone face normals
    thetas = 2*pi*(0:nfaces-1)'/nfaces;

    rows = zeros(total_constraints * 3, 1);  % sparse triplets (over-allocate)
    cols = zeros(total_constraints * 3, 1);
    vals = zeros(total_constraints * 3, 1);
    b_los = zeros(total_constraints, 1);

    idx = 0;   % triplet counter
    row = 0;   % constraint row counter

    for k = 0:Np
        x_offset = k * nx;  % offset to x_k in z vector

        % --- Cone face constraints:  cos(th)*xT + sin(th)*zT - cone_k*yT <= 0 ---
        for f = 1:nfaces
            row = row + 1;
            ct = cos(thetas(f));
            st = sin(thetas(f));

            % xT component (index 1 in x_k)
            idx = idx + 1;
            rows(idx) = row; cols(idx) = x_offset + 1; vals(idx) = ct;

            % yT component (index 2 in x_k) with -cone_k
            idx = idx + 1;
            rows(idx) = row; cols(idx) = x_offset + 2; vals(idx) = -cone_k;

            % zT component (index 3 in x_k)
            idx = idx + 1;
            rows(idx) = row; cols(idx) = x_offset + 3; vals(idx) = st;

            b_los(row) = 0;
        end

        % --- y_min constraint:  -yT <= -y_min  =>  yT >= y_min ---
        row = row + 1;
        idx = idx + 1;
        rows(idx) = row; cols(idx) = x_offset + 2; vals(idx) = -1;
        b_los(row) = -y_min;
    end

    % Trim and build sparse matrix
    rows = rows(1:idx);
    cols = cols(1:idx);
    vals = vals(1:idx);
    A_los = sparse(rows, cols, vals, total_constraints, total_vars);
end
