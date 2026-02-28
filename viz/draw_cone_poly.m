function h = draw_cone_poly(origin, R, cone_k, cone_L, nfaces, color)
%DRAW_CONE_POLY  Draw a polyhedral cone (wireframe) attached to a body frame.
%   h = draw_cone_poly(origin, R, cone_k, cone_L, nfaces, color)
%   origin  : 3x1 apex position (target center)
%   R       : 3x3 rotation matrix (body frame to plot frame; cone axis = R*[0;1;0])
%   cone_k  : tan(half_angle)
%   cone_L  : length of cone to draw
%   nfaces  : number of polygon faces
%   color   : edge color

    thetas = linspace(0, 2*pi, nfaces+1);
    thetas = thetas(1:end-1);

    % Tip at origin, base circle at y=cone_L in body frame
    % radius at base = cone_k * cone_L
    radius = cone_k * cone_L;

    % Base points in body frame
    base_body = zeros(3, nfaces);
    for i = 1:nfaces
        base_body(:, i) = [radius*cos(thetas(i)); cone_L; radius*sin(thetas(i))];
    end

    % Transform to plot frame
    base_plot = R * base_body + origin(:);
    o = origin(:);

    h = gobjects(nfaces*2, 1);
    idx = 0;
    % Draw cone edges (apex to each base vertex)
    for i = 1:nfaces
        idx = idx + 1;
        h(idx) = plot3([o(1) base_plot(1,i)], [o(2) base_plot(2,i)], ...
                       [o(3) base_plot(3,i)], '--', 'Color', color, 'LineWidth', 0.8);
    end
    % Draw base polygon
    for i = 1:nfaces
        j = mod(i, nfaces) + 1;
        idx = idx + 1;
        h(idx) = plot3([base_plot(1,i) base_plot(1,j)], ...
                       [base_plot(2,i) base_plot(2,j)], ...
                       [base_plot(3,i) base_plot(3,j)], '-', 'Color', color, 'LineWidth', 1);
    end
end
