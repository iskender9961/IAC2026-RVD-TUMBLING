function b_expanded = minkowski_sum_halfspace(A_c, b_c, B, u_max)
%MINKOWSKI_SUM_HALFSPACE  Expand half-space constraints for forward reachability.
%   b_expanded = minkowski_sum_halfspace(A_c, b_c, B, u_max)
%
%   For forward reachability, the reachable set is:
%       X_{k+1} = A * X_k  oplus  B * U
%
%   If X_k = { x : H*x <= h } and U = { u : |u_i| <= u_max }, then
%   the Minkowski sum A*X oplus B*U in the H-representation adds the
%   support function of B*U to each constraint:
%
%       h_new(i) = h(i) + max_{u in U} H_i * B * u
%                = h(i) + sum_j |H_i * B_j| * u_max_j
%
%   Inputs:
%       A_c   - n_con x nx constraint matrix (H) — only used for dimensions
%       b_c   - n_con x 1 constraint RHS (h)
%       B     - nx x nu input matrix
%       u_max - scalar or nu x 1 input bound
%
%   Outputs:
%       b_expanded - n_con x 1 expanded RHS

    n_con = size(A_c, 1);
    nu = size(B, 2);

    if isscalar(u_max)
        u_max = u_max * ones(nu, 1);
    end

    b_expanded = b_c;
    for i = 1:n_con
        nTB = A_c(i,:) * B;  % 1 x nu
        support = sum(abs(nTB(:)) .* u_max(:));
        b_expanded(i) = b_c(i) + support;
    end
end
