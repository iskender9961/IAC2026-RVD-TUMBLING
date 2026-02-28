function cu = common_utils()
%COMMON_UTILS  Shared definitions for the toy reachability example.
%
%   cu = common_utils() returns a struct with:
%       A, B       - double-integrator discrete dynamics
%       u_min/max  - input bounds
%       x_min/max  - state constraint box
%       N          - reachability horizon
%       nx_grid    - grid resolution per axis
%       target     - target set (polytope H-rep)
%       colors     - plot colors for methods
%
%   System:  x_{k+1} = A x_k + B u_k  (+w_k for stochastic/robust)
%   State:   x = [position; velocity]  (double integrator)

    % Double integrator (dt = 1)
    cu.A = [1 1; 0 1];
    cu.B = [0.5; 1];  % exact ZOH: [0.5*dt^2; dt] with dt=1

    % Input bounds
    cu.u_min = -1;
    cu.u_max =  1;

    % State constraints (box)
    cu.x_min = [-5; -3];
    cu.x_max = [ 5;  3];

    % Horizon
    cu.N = 5;

    % Grid resolution
    cu.nx_grid = 200;

    % Target set: small box around origin
    % {x: H_t * x <= h_t}
    cu.H_target = [eye(2); -eye(2)];
    cu.h_target = [0.5; 0.3; 0.5; 0.3];

    % Stochastic parameters
    cu.w_cov = diag([0.01, 0.005]);  % per-step Gaussian noise covariance
    cu.alpha = 0.05;                  % chance constraint violation probability

    % Robust parameters
    cu.w_max = [0.1; 0.05];  % per-step bounded disturbance (box)

    % Plot colors
    cu.col_nominal   = [0.2, 0.75, 0.2];   % green
    cu.col_stochastic = [0.3, 0.5, 0.85];  % blue
    cu.col_robust    = [0.6, 0.2, 0.6];    % purple
    cu.col_target    = [1.0, 0.85, 0.0];   % gold
    cu.col_constraint = [0.9, 0.9, 0.9];   % light gray
end
