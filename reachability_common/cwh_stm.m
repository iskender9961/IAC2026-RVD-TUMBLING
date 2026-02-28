function [Ad, Bd] = cwh_stm(dt, n_orbit)
%CWH_STM  Exact discrete-time CWH/HCW state-transition and input matrices.
%   [Ad, Bd] = cwh_stm(dt, n_orbit)
%
%   State: x = [x; y; z; vx; vy; vz]  in LVLH
%       x  = radial (outward)
%       y  = along-track
%       z  = cross-track
%
%   Continuous CWH:
%       x''  = 3n^2*x + 2n*y'  + a_x
%       y''  =        - 2n*x'  + a_y
%       z''  = -n^2*z           + a_z
%
%   Returns exact matrix exponential (closed-form) and ZOH input matrix.
%
%   Inputs:
%       dt      - time step [s]
%       n_orbit - mean motion [rad/s]
%
%   Outputs:
%       Ad - 6x6 state transition matrix Phi(dt)
%       Bd - 6x3 input matrix (zero-order hold)

    nt = n_orbit * dt;
    c  = cos(nt);
    s  = sin(nt);

    % State transition Phi(dt) - exact closed-form
    Ad = [
        4 - 3*c,          0, 0,  s/n_orbit,           2*(1 - c)/n_orbit,    0;
        6*(s - nt),        1, 0, -2*(1 - c)/n_orbit,   (4*s - 3*nt)/n_orbit, 0;
        0,                 0, c,  0,                    0,                    s/n_orbit;
        3*n_orbit*s,       0, 0,  c,                    2*s,                  0;
       -6*n_orbit*(1 - c), 0, 0, -2*s,                  4*c - 3,              0;
        0,                 0, -n_orbit*s, 0,             0,                    c
    ];

    % ZOH input matrix Bd = integral_0^dt Phi(tau) dtau * B_c
    Bd = [
        (1 - c)/n_orbit^2,              2*(nt - s)/n_orbit^2,              0;
       -2*(nt - s)/n_orbit^2,           (4*(1 - c) - 1.5*nt^2)/n_orbit^2, 0;
        0,                               0,                                (1 - c)/n_orbit^2;
        s/n_orbit,                       2*(1 - c)/n_orbit,                0;
       -2*(1 - c)/n_orbit,              (4*s - 3*nt)/n_orbit,             0;
        0,                               0,                                s/n_orbit
    ];
end
