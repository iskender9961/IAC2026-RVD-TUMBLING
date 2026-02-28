function h = draw_triad(origin, R, len, labels, line_styles)
%DRAW_TRIAD  Draw a coordinate frame triad (3 axes).
%   h = draw_triad(origin, R, len, labels, line_styles)
%   origin      : 3x1 position of the frame origin
%   R           : 3x3 rotation matrix (frame axes as columns, in plot frame)
%   len         : scalar axis length
%   labels      : 3x1 cell array of axis label strings, e.g. {'xT','yT','zT'}
%   line_styles : 3x1 cell array of line style strings, e.g. {'-','--','-.'}
%
%   Returns h as a struct with fields: lines(3), texts(3).

    colors = {'r', [0 0.7 0], 'b'};   % x=red, y=green, z=blue
    o = origin(:);
    h.lines = gobjects(3,1);
    h.texts = gobjects(3,1);
    for i = 1:3
        tip = o + len * R(:, i);
        h.lines(i) = plot3([o(1) tip(1)], [o(2) tip(2)], [o(3) tip(3)], ...
            line_styles{i}, 'Color', colors{i}, 'LineWidth', 2);
        h.texts(i) = text(tip(1), tip(2), tip(3), ['  ' labels{i}], ...
            'Color', colors{i}, 'FontSize', 10, 'FontWeight', 'bold');
    end
end
