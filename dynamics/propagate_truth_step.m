function x_next = propagate_truth_step(x, dt, mu, Re, J2, a_ctrl_eci, ode_opts)
%PROPAGATE_TRUTH_STEP  Propagate one spacecraft state by dt using nonlinear ECI dynamics.
%   x_next = propagate_truth_step(x, dt, mu, Re, J2, a_ctrl_eci, ode_opts)
%   x           : 6x1 ECI state [r; v] [m; m/s]
%   dt          : time step [s]
%   a_ctrl_eci  : 3x1 constant control acceleration over this step [m/s^2]
%   ode_opts    : odeset options struct
%   x_next      : 6x1 ECI state after dt
    odefun = @(t, s) eom_eci(t, s, mu, Re, J2, a_ctrl_eci);
    [~, X] = ode113(odefun, [0 dt], x, ode_opts);
    x_next = X(end, :)';
end
