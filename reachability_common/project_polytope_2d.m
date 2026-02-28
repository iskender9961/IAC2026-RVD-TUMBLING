function [boundary_x, boundary_y] = project_polytope_2d(A, b, proj_dims, n_dirs)
%PROJECT_POLYTOPE_2D  Project a polytope {x: Ax<=b} onto 2D plane.
%   [bx, by] = project_polytope_2d(A, b, proj_dims, n_dirs)
%
%   Uses support function evaluation: for each direction d in the 2D plane,
%   solve max d'*x_proj subject to Ax<=b, where x_proj are the projected
%   dimensions.
%
%   This is done via LP for each direction, yielding the boundary of the
%   projected polytope.
%
%   Inputs:
%       A         - n_con x nx constraint matrix
%       b         - n_con x 1 constraint RHS
%       proj_dims - [d1, d2] indices of dimensions to project onto (1-based)
%       n_dirs    - number of directions to sample (default: 120)
%
%   Outputs:
%       boundary_x - n_dirs x 1 x-coordinates of projected boundary
%       boundary_y - n_dirs x 1 y-coordinates of projected boundary

    if nargin < 4, n_dirs = 120; end

    nx = size(A, 2);
    d1 = proj_dims(1);
    d2 = proj_dims(2);

    angles = linspace(0, 2*pi, n_dirs + 1);
    angles = angles(1:end-1);

    boundary_x = zeros(n_dirs, 1);
    boundary_y = zeros(n_dirs, 1);

    opts = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');

    valid = true(n_dirs, 1);

    for i = 1:n_dirs
        % Direction in projected plane
        dir_x = cos(angles(i));
        dir_y = sin(angles(i));

        % Objective: maximize dir' * [x_{d1}; x_{d2}]
        % = minimize -dir' * [x_{d1}; x_{d2}]
        f = zeros(nx, 1);
        f(d1) = -dir_x;
        f(d2) = -dir_y;

        [x_opt, fval, exitflag] = linprog(f, A, b, [], [], [], [], opts);

        if exitflag == 1
            boundary_x(i) = x_opt(d1);
            boundary_y(i) = x_opt(d2);
        else
            valid(i) = false;
        end
    end

    boundary_x = boundary_x(valid);
    boundary_y = boundary_y(valid);

    % Order points by angle for clean polygon
    cx = mean(boundary_x);
    cy = mean(boundary_y);
    ang = atan2(boundary_y - cy, boundary_x - cx);
    [~, idx] = sort(ang);
    boundary_x = boundary_x(idx);
    boundary_y = boundary_y(idx);
end
