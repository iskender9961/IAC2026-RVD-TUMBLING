function dxdt = eom_eci(~, x, mu, Re, J2, a_ctrl_eci)
%EOM_ECI  Equations of motion for a spacecraft in ECI (2-body + J2 + control).
%   dxdt = eom_eci(t, x, mu, Re, J2, a_ctrl_eci)
%   x            : 6x1 state [r; v] in ECI [m; m/s]
%   a_ctrl_eci   : 3x1 control acceleration in ECI [m/s^2]
%   dxdt         : 6x1 state derivative
    r = x(1:3);
    v = x(4:6);
    a_grav = accel_2body(r, mu);
    a_j2   = accel_J2(r, mu, Re, J2);
    dxdt = [v; a_grav + a_j2 + a_ctrl_eci];
end
