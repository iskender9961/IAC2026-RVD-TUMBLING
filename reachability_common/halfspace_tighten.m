function b_tight = halfspace_tighten(A_c, b_c, B, u_max)
%HALFSPACE_TIGHTEN  Tighten half-space constraints for backward reachability.
%   b_tight = halfspace_tighten(A_c, b_c, B, u_max)
%
%   For backward reachability (pre-image computation), we need:
%       Pre(X_next) = { x : exists u in U s.t. Ax + Bu in X_next } cap X_adm
%
%   If X_next = { z : H*z <= h }, then:
%       H * (Ax + Bu) <= h
%       H*A*x <= h - max_{u in U} H*B*u
%
%   For box input constraints |u_i| <= u_max:
%       max_{u in U} n^T * B * u = sum_j |n^T * B_j| * u_max
%   where n = H(i,:) is a constraint normal and B_j is the j-th column of B.
%
%   This function computes the tightened RHS:
%       b_tight(i) = b_c(i) - max_{u in U} A_c(i,:) * B * u
%
%   Inputs:
%       A_c   - n_con x nx constraint matrix (H)
%       b_c   - n_con x 1 constraint RHS (h)
%       B     - nx x nu input matrix
%       u_max - scalar or nu x 1 input bound
%
%   Outputs:
%       b_tight - n_con x 1 tightened RHS

    n_con = size(A_c, 1);
    nu = size(B, 2);

    if isscalar(u_max)
        u_max = u_max * ones(nu, 1);
    end

    b_tight = b_c;
    for i = 1:n_con
        % n^T * B where n = A_c(i,:)
        nTB = A_c(i,:) * B;  % 1 x nu
        % max over box: sum of |nTB_j| * u_max_j
        support = sum(abs(nTB(:)) .* u_max(:));
        b_tight(i) = b_c(i) - support;
    end
end
