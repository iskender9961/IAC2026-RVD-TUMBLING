function [Ad, Bd, nu1] = ya_stm(dt, cfg, nu0)
%YA_STM  Discrete-time YA state-transition and input matrices for elliptic orbits.
%
%   [Ad, Bd, nu1] = ya_stm(dt, cfg, nu0)
%
%   Wraps the Yamanaka-Ankersen STM from Elliptic_orbit/ into the same
%   (Ad, Bd) interface used by cwh_stm.m.
%
%   Inputs:
%       dt   - time step [s]
%       cfg  - orbital config struct from elliptic_orbit_config.m
%              (must contain: ecc, ecc2, mu, a, p, ya_k2)
%       nu0  - current true anomaly [rad]
%
%   Outputs:
%       Ad   - 6x6 state transition matrix
%       Bd   - 6x3 input matrix
%       nu1  - true anomaly at end of step [rad]
%
%   State ordering: x = [x; y; z; vx; vy; vz] in LVLH-like frame
%
%   The input matrix follows the standard YA convention:
%       Bd = Ad * [0_{3x3}; I_{3x3}]
%   as used in Hartley's FORCES Pro implementation and consistent with
%   the Yamanaka-Ankersen formulation where inputs enter as velocity
%   impulses propagated through the STM.
%
%   Reference:
%       K. Yamanaka and F. Ankersen, "New state transition matrix for
%       relative motion on an arbitrary elliptical orbit," J. Guidance,
%       Control, and Dynamics, vol. 25, no. 1, pp. 60-66, 2002.

% --- Ensure Elliptic_orbit/ is on the path ---
persistent ya_path_added
if isempty(ya_path_added)
    ya_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Elliptic_orbit');
    addpath(ya_dir);
    ya_path_added = true;
end

% --- Get the 6x6 STM from the YA implementation ---
[Ad, nu1] = YA_A_matrix(cfg, nu0, dt);

% --- Input matrix: standard YA formulation ---
% Bd = Ad * [0; I3] maps velocity impulses through the STM.
% This is equivalent to a zero-order hold on the input with the
% impulse applied at the beginning of the interval.
Bd = Ad * [zeros(3); eye(3)];
end
