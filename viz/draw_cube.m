function h = draw_cube(center, half_edge, R, color, alpha_val)
%DRAW_CUBE  Draw a 3D cube centered at 'center' with orientation R.
%   h = draw_cube(center, half_edge, R, color, alpha_val)
%   center    : 3x1 position
%   half_edge : scalar half-edge length
%   R         : 3x3 rotation matrix (body-to-plot frame)
%   color     : color spec (e.g. 'r', [1 0 0])
%   alpha_val : face alpha (transparency)

    s = half_edge;
    % Unit cube vertices (centered at origin)
    v = s * [ -1 -1 -1;  1 -1 -1;  1  1 -1; -1  1 -1;
              -1 -1  1;  1 -1  1;  1  1  1; -1  1  1 ];
    % Rotate and translate
    v = (R * v')' + center(:)';

    % 6 faces (indices into v)
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
    h = patch('Vertices', v, 'Faces', faces, ...
              'FaceColor', color, 'FaceAlpha', alpha_val, ...
              'EdgeColor', 'k', 'LineWidth', 0.5);
end
