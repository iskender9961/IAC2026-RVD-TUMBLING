function [prob, data] = build_qp_osqp(Ad, Bd, x0, x_ref, u_prev, p)
%BUILD_QP_OSQP  Construct and initialize the OSQP problem for the MPC QP.
%
%   Decision variable layout:
%     z = [x_0; x_1; ...; x_Np; u_0; u_1; ...; u_{Np-1}]
%         |<-- (Np+1)*nx -->|   |<--    Np*nu         -->|
%
%   Cost:  min 0.5 z' P z + q' z
%   Constraints:  l <= A z <= u
%     1) Dynamics equalities (including initial condition)
%     2) Input box bounds
%     3) LOS polyhedral cone

    nx = p.nx;
    nu = p.nu;
    Np = p.Np;
    n_x_vars = (Np+1) * nx;
    n_u_vars = Np * nu;
    n_vars   = n_x_vars + n_u_vars;

    %% ===== Cost matrix P (constant, only depends on weights) =====
    P = sparse(n_vars, n_vars);

    % Stage state cost: Q for k=0..Np-1
    for k = 0:Np-1
        idx = k*nx + (1:nx);
        P(idx, idx) = p.Q;
    end
    % Terminal state cost: QN for k=Np
    idx = Np*nx + (1:nx);
    P(idx, idx) = p.QN;

    % Input cost Ru
    for k = 0:Np-1
        u_idx = n_x_vars + k*nu + (1:nu);
        P(u_idx, u_idx) = P(u_idx, u_idx) + p.Ru;
    end

    % Delta-u cost Rdu:  sum_{k=0}^{Np-1} (u_k - u_{k-1})' Rdu (u_k - u_{k-1})
    % u_{-1} = u_prev (parameter => only affects q for k=0)
    % Quadratic expansion:
    %   (u_k - u_{k-1})' Rdu (u_k - u_{k-1})
    %   = u_k' Rdu u_k - 2 u_k' Rdu u_{k-1} + u_{k-1}' Rdu u_{k-1}
    for k = 0:Np-1
        u_idx = n_x_vars + k*nu + (1:nu);
        P(u_idx, u_idx) = P(u_idx, u_idx) + p.Rdu;
        if k > 0
            u_prev_idx = n_x_vars + (k-1)*nu + (1:nu);
            P(u_prev_idx, u_prev_idx) = P(u_prev_idx, u_prev_idx) + p.Rdu;
            P(u_idx, u_prev_idx) = P(u_idx, u_prev_idx) - p.Rdu;
            P(u_prev_idx, u_idx) = P(u_prev_idx, u_idx) - p.Rdu;
        end
    end

    % OSQP wants 0.5 z' P_osqp z + q' z
    % Our cost is z' P z + q' z  =>  P_osqp = 2*P
    P = (P + P') / 2;    % ensure exact symmetry
    P_osqp = 2 * P;

    %% ===== Cost vector q =====
    q = zeros(n_vars, 1);
    for k = 0:Np-1
        idx = k*nx + (1:nx);
        q(idx) = -2 * p.Q * x_ref;
    end
    idx = Np*nx + (1:nx);
    q(idx) = -2 * p.QN * x_ref;
    % Delta-u linear term at k=0:  -2 * u_0' Rdu u_prev  =>  q += -2*Rdu*u_prev
    u0_idx = n_x_vars + (1:nu);
    q(u0_idx) = q(u0_idx) - 2 * p.Rdu * u_prev;

    %% ===== Constraint matrix A =====
    n_dyn = (Np+1) * nx;
    n_input = Np * nu;

    % --- Dynamics + initial condition ---
    Adyn = sparse(n_dyn, n_vars);
    ldyn = zeros(n_dyn, 1);
    udyn = zeros(n_dyn, 1);

    Adyn(1:nx, 1:nx) = speye(nx);
    ldyn(1:nx) = x0;
    udyn(1:nx) = x0;

    for k = 0:Np-1
        row_idx   = (k+1)*nx + (1:nx);
        x_k_idx   = k*nx + (1:nx);
        x_kp1_idx = (k+1)*nx + (1:nx);
        u_k_idx   = n_x_vars + k*nu + (1:nu);
        Adyn(row_idx, x_kp1_idx) = speye(nx);
        Adyn(row_idx, x_k_idx)   = -Ad;
        Adyn(row_idx, u_k_idx)   = -Bd;
    end

    % --- Input bounds ---
    Ainput = sparse(n_input, n_vars);
    for k = 0:Np-1
        row_idx = k*nu + (1:nu);
        u_k_idx = n_x_vars + k*nu + (1:nu);
        Ainput(row_idx, u_k_idx) = speye(nu);
    end
    linput = -p.u_max * ones(n_input, 1);
    uinput =  p.u_max * ones(n_input, 1);

    % --- LOS cone (polyhedral) ---
    [A_los, b_los] = los_polyhedral_constraints(p.cone_k, p.y_min, ...
                        p.cone_nfaces, Np, nx, nu);
    n_los = size(A_los, 1);
    l_los = -inf(n_los, 1);
    u_los = b_los;

    % --- Stack ---
    A_all = [Adyn; Ainput; A_los];
    l_all = [ldyn; linput; l_los];
    u_all = [udyn; uinput; u_los];

    %% ===== Create OSQP problem =====
    prob = osqp;
    prob.setup(P_osqp, q, A_all, l_all, u_all, ...
        'max_iter', p.osqp_max_iter, ...
        'eps_abs', p.osqp_eps_abs, ...
        'eps_rel', p.osqp_eps_rel, ...
        'warm_start', p.osqp_warm_start, ...
        'verbose', p.osqp_verbose, ...
        'polish', true, ...
        'adaptive_rho', true);

    %% ===== Data struct for updates =====
    data.P_osqp    = P_osqp;
    data.n_vars    = n_vars;
    data.n_x_vars  = n_x_vars;
    data.n_u_vars  = n_u_vars;
    data.n_dyn     = n_dyn;
    data.n_input   = n_input;
    data.n_los     = n_los;
    data.nx        = nx;
    data.nu        = nu;
    data.Np        = Np;
end
