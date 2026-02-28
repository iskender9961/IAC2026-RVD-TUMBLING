function result = forward_reachable_set(omega_deg, a_max, rp)
%FORWARD_REACHABLE_SET  Compute forward reachable set (nominal, deterministic).
%   result = forward_reachable_set(omega_deg, a_max, rp)
%
%   Propagates the admissible set forward in time:
%       X_{k+1} = (A * X_k) intersect X_adm(k+1)
%   where the control input expands the set via Minkowski sum with B*U.
%
%   The admissible set X_adm(k) is the LOS cone rotated to LVLH at time k*dt.
%   We work with the H-representation {x: H*x <= h} throughout.
%
%   For the initial set, we use the LOS cone at t=0 as X_0 = X_adm(0).
%
%   At each step:
%     1. Transform X_k through dynamics: H_k * A^{-1} * x <= h_k
%     2. Expand by Minkowski sum with B*U: h += support(H, B*U)
%     3. Intersect with X_adm(k+1): add constraints from LOS at t=(k+1)*dt
%
%   Since tracking full polytope evolution in H-rep is complex with
%   increasing constraint count, we use a practical approach:
%     - At each step, evaluate the safe mask on a grid of body-frame points
%     - A point x_B is forward-reachable if at step k there exists a
%       sequence of controls u_0,...,u_{k-1} that keeps the trajectory
%       inside the LOS cone from step 0 to k.
%
%   We use the one-step backward reachable set characterization:
%     A point x_B at step k is reachable from X_{k-1} if:
%       exists x_{k-1} in X_{k-1} and u in U s.t. x_k = A*x_{k-1} + B*u
%     AND x_k in X_adm(k)
%
%   Equivalently, x_k is reachable if x_k in X_adm(k) AND
%   there exists u in U s.t. A^{-1}*(x_k - B*u) in X_{k-1}.
%
%   For the grid evaluation, we check each point for membership at each step.
%
%   Inputs:
%       omega_deg - target tumble rate [deg/s]
%       a_max     - max chaser acceleration [m/s^2]
%       rp        - reachability parameters struct
%
%   Outputs:
%       result - struct with fields:
%           safe_mask    - ny x nx logical, true if point is forward-reachable
%           x_grid       - 1 x nx body-frame x values
%           y_grid       - 1 x ny body-frame y values
%           cone_mask    - ny x nx logical, true if inside LOS cone
%           omega_deg    - target tumble rate
%           a_max        - max acceleration
%           method       - 'nominal_forward'

    omega_rad = deg2rad(omega_deg);
    dt = rp.dt;
    n_orbit = rp.n;

    % Get CWH matrices in LVLH
    [Ad, Bd] = cwh_stm(dt, n_orbit);

    % Build body-frame evaluation grid
    x_grid = linspace(rp.x_range(1), rp.x_range(2), rp.nx_grid);
    y_grid = linspace(rp.y_range(1), rp.y_range(2), rp.ny_grid);
    ny = length(y_grid);
    nx = length(x_grid);

    % Body-frame LOS constraints (time-invariant in body frame)
    [A_body, b_body] = los_constraints_body(rp.cone_k, rp.y_min, rp.cone_nfaces);

    % Cone mask: which grid points are inside the LOS cone
    cone_mask = false(ny, nx);
    for iy = 1:ny
        for ix = 1:nx
            p_body = [x_grid(ix); y_grid(iy); 0; 0; 0; 0];
            cone_mask(iy, ix) = all(A_body * p_body <= b_body + 1e-9);
        end
    end

    % For the nominal forward analysis, we compute the directional erosion
    % at each grid point. This tells us whether a point starting at rest
    % in the body frame can maintain LOS feasibility as the cone rotates.
    %
    % The forward reachable set from a starting point x_0 in body frame
    % with v_0=0 is safe if the controller can counteract the apparent
    % rotation-induced drift before any constraint is violated.
    %
    % Using the one-step dynamic analysis:
    % At t=0, state in body frame is [x_B, y_B, 0, 0, 0, 0]
    % In LVLH, the state is obtained by inverse rotation.
    % After dt, LVLH state evolves by CWH + control.
    % In new body frame (rotated by omega*dt), the new body state must
    % still satisfy LOS constraints.
    %
    % For multi-step analysis, we propagate the safe region backward
    % from the final step. See backward_safe_start_set.m for that approach.
    %
    % Here we use the directional per-constraint erosion model which gives
    % a conservative (inner) approximation of the forward-reachable safe set.

    safe_mask = false(ny, nx);

    % Synchronization range limit: r_max = 2 * a_max / omega^2
    w2 = omega_rad^2;
    r_sync_max = 2 * a_max / max(w2, 1e-30);

    for iy = 1:ny
        for ix = 1:nx
            if ~cone_mask(iy, ix)
                continue;
            end

            xb = x_grid(ix);
            yb = y_grid(iy);
            rng = sqrt(xb^2 + yb^2);

            % Check synchronization range
            if rng >= r_sync_max
                continue;
            end

            % Body-frame state (zero initial velocity)
            p_body = [xb; yb; 0; 0; 0; 0];

            % Per-constraint slack
            slacks = b_body - A_body * p_body;

            % Directional per-constraint erosion
            % Apparent body-frame velocity for inertially-stationary chaser:
            %   v_rot = [omega*y_B, -omega*x_B, 0]
            v_rot = [omega_rad * yb; -omega_rad * xb; 0];

            safe_point = true;
            for ic = 1:size(A_body, 1)
                a_pos = A_body(ic, 1:3)';
                slack_rate = -a_pos' * v_rot;
                if slack_rate < 0
                    erosion = 0.5 * slack_rate^2 / a_max;
                    if slacks(ic) < erosion
                        safe_point = false;
                        break;
                    end
                end
            end

            safe_mask(iy, ix) = safe_point;
        end
    end

    result.safe_mask  = safe_mask;
    result.x_grid    = x_grid;
    result.y_grid    = y_grid;
    result.cone_mask = cone_mask;
    result.omega_deg = omega_deg;
    result.a_max     = a_max;
    result.method    = 'nominal_forward';
end
