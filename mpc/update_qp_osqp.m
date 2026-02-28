function [u_opt, status, x_pred] = update_qp_osqp(prob, data, Ad, Bd, x0, x_ref, u_prev, p)
%UPDATE_QP_OSQP  Rebuild and solve the OSQP QP for the current MPC step.
%   [u_opt, status, x_pred] = update_qp_osqp(prob, data, Ad, Bd, x0, x_ref, u_prev, p)
%
%   Because Ad/Bd change sparsity pattern each step (finite-difference
%   linearization introduces/removes numerical zeros), we rebuild and re-setup
%   the OSQP problem each step. This is fast enough for the problem size
%   (~300 variables) and avoids sparsity pattern mismatch issues.
%
%   The 'prob' argument is ignored (kept for interface compatibility).
%   A fresh OSQP object is created each call.

    nx = data.nx;
    nu = data.nu;
    Np = data.Np;
    n_x_vars = data.n_x_vars;
    n_vars   = data.n_vars;

    % Rebuild full problem from scratch with current Ad, Bd, x0, x_ref, u_prev
    [prob_new, ~] = build_qp_osqp(Ad, Bd, x0, x_ref, u_prev, p);

    % Solve
    res = prob_new.solve();
    status = res.info.status;

    if ~contains(status, 'solved')
        u_opt = zeros(nu, 1);
        x_pred = [];
        return;
    end

    z_opt = res.x;

    % Extract first input
    u_opt = z_opt(n_x_vars + (1:nu));

    % Extract predicted states
    x_pred = zeros(nx, Np+1);
    for k = 0:Np
        x_pred(:, k+1) = z_opt(k*nx + (1:nx));
    end
end
