function generate_gifs(lg, p, results_dir)
%GENERATE_GIFS  Create animated GIFs from simulation results.
%
%   GIF 1: gif_tb_topview     -- TB xy with fixed cone, cubes, triads, HUD
%   GIF 2: gif_lvlh_rotating_cone -- LVLH 3D with rotating cone, cubes, triads, HUD
%   GIF 3: gif_tb_3d          -- TB 3D perspective with cubes, triads, cone, HUD
%   GIF 4: gif_2d_xy_tb       -- TB xy 2D with cone
%   GIF 5: gif_2d_xy_lvlh     -- LVLH xy 2D with rotating cone projection

    if nargin < 3, results_dir = 'results'; end

    Nt = length(lg.t_hist);
    cone_k = p.cone_k;
    frame_skip = max(1, floor(Nt / 150));

    % Fixed cone length = initial y distance (never shrinks)
    cone_L_fixed = lg.r_tb_hist(2, 1) * 1.2;
    y_max = cone_L_fixed;
    x_max = max(max(abs(lg.r_tb_hist(1,:))), cone_k * y_max) * 1.3;
    r_max = max(vecnorm(lg.r_lvlh_hist, 2, 1)) * 1.3;

    % Pre-compute derived quantities for HUD
    Nu = size(lg.u_hist, 2);

    %% ================================================================
    %  GIF 1: Top view TB (x_TB vs y_TB) -- cubes, triads, fixed cone, HUD
    %  ================================================================
    fprintf('  GIF 1: TB top view...\n');
    gif_file = fullfile(results_dir, 'gif_tb_topview.gif');
    fig = figure('Name','GIF TB Top','Position',[50 50 900 700],'Visible','off');

    for ii = 1:frame_skip:Nt
        clf;
        hold on; grid on; axis equal;
        t_now = lg.t_hist(ii);
        r_now = lg.r_tb_hist(:, ii);

        % Cone (fixed size)
        yy = linspace(0, y_max, 100);
        fill([ cone_k*yy, fliplr(-cone_k*yy)], [yy, fliplr(yy)], ...
            [1 0.95 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot( cone_k*yy, yy, 'Color',[0.8 0.5 0], 'LineWidth', 1.5);
        plot(-cone_k*yy, yy, 'Color',[0.8 0.5 0], 'LineWidth', 1.5);
        % Docking axis
        plot([0 0], [0 y_max], 'k--', 'LineWidth', 2);

        % Target cube (2D square) + triad
        draw_square_2d([0,0], p.cube_size, [0.3 0.3 0.8]);
        draw_triad_2d([0,0], eye(2), p.triad_len*0.6, {'x_T','y_T'}, '-');
        % Chaser cube + triad
        draw_square_2d([r_now(1),r_now(2)], p.cube_size, [0.8 0.3 0.3]);
        draw_triad_2d([r_now(1),r_now(2)], eye(2), p.triad_len*0.4, {'x_C','y_C'}, '--');

        % Trajectory
        plot(lg.r_tb_hist(1,1:ii), lg.r_tb_hist(2,1:ii), 'b-', 'LineWidth', 1.5);
        plot(lg.r_tb_hist(1,1), lg.r_tb_hist(2,1), 'g^','MarkerSize',10,'MarkerFaceColor','g');

        xlim([-x_max x_max]); ylim([-10 y_max*1.05]);
        xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]');
        title('Target Body Frame -- Top View (x_{TB} vs y_{TB})');

        % HUD
        add_hud(lg, p, ii, Nu);

        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, ii == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 2: LVLH 3D with rotating cone, cubes, triads, HUD
    %  ================================================================
    fprintf('  GIF 2: LVLH 3D rotating cone...\n');
    gif_file = fullfile(results_dir, 'gif_lvlh_rotating_cone.gif');
    fig = figure('Name','GIF LVLH','Position',[50 50 1000 750],'Visible','off');

    for ii = 1:frame_skip:Nt
        clf;
        hold on; grid on; axis equal;
        r_now_l = lg.r_lvlh_hist(:, ii);

        % Trajectory
        plot3(lg.r_lvlh_hist(1,1:ii), lg.r_lvlh_hist(2,1:ii), ...
              lg.r_lvlh_hist(3,1:ii), 'b-', 'LineWidth', 1.5);
        plot3(lg.r_lvlh_hist(1,1), lg.r_lvlh_hist(2,1), ...
              lg.r_lvlh_hist(3,1), 'g^','MarkerSize',10,'MarkerFaceColor','g');

        % Cubes
        draw_cube([0;0;0], p.cube_size, eye(3), [0.3 0.3 0.8], 0.4);
        draw_cube(r_now_l, p.cube_size, eye(3), [0.8 0.3 0.3], 0.6);

        % Triads at cubes
        draw_triad([0;0;0], eye(3), p.triad_len*0.6, {'x_L','y_L','z_L'}, {'-','-','-'});
        draw_triad(r_now_l, eye(3), p.triad_len*0.4, {'x_C','y_C','z_C'}, {'--','--','--'});

        % Docking axis in LVLH
        d = lg.dock_axis_lvlh_hist(:, ii);
        plot3([0 d(1)*cone_L_fixed], [0 d(2)*cone_L_fixed], [0 d(3)*cone_L_fixed], ...
            'k--', 'LineWidth', 2.5);

        % Rotating cone (fixed length)
        R_eci_tb   = lg.R_eci_tb_hist(:,:,ii);
        R_eci_lvlh = lg.R_eci_lvlh_hist(:,:,ii);
        R_lvlh_tb  = R_eci_lvlh' * R_eci_tb;
        draw_los_tetra([0;0;0], R_lvlh_tb, cone_k, cone_L_fixed, p.cone_nfaces, [0.8 0.5 0]);

        % Target body triad at origin (rotated into LVLH)
        draw_triad([0;0;0], R_lvlh_tb, p.triad_len*0.5, {'x_T','y_T','z_T'}, {'-.','-.','-.'});

        xlim([-r_max r_max]); ylim([-r_max r_max]); zlim([-r_max r_max]);
        xlabel('x_{LVLH} [m]'); ylabel('y_{LVLH} [m]'); zlabel('z_{LVLH} [m]');
        title('LVLH Frame -- Rotating Corridor');
        view(135, 25);

        add_hud(lg, p, ii, Nu);
        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, ii == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 3: 3D TB frame -- cubes, triads, fixed cone, HUD
    %  ================================================================
    fprintf('  GIF 3: TB 3D perspective...\n');
    gif_file = fullfile(results_dir, 'gif_tb_3d.gif');
    fig = figure('Name','GIF TB 3D','Position',[50 50 1000 750],'Visible','off');

    for ii = 1:frame_skip:Nt
        clf;
        hold on; grid on; axis equal;
        r_now = lg.r_tb_hist(:, ii);

        % Fixed-size cone + docking axis
        draw_los_tetra([0;0;0], eye(3), cone_k, cone_L_fixed, p.cone_nfaces, [0.8 0.5 0]);
        plot3([0 0], [0 cone_L_fixed], [0 0], 'k--', 'LineWidth', 2);

        % Trajectory
        plot3(lg.r_tb_hist(1,1:ii), lg.r_tb_hist(2,1:ii), ...
              lg.r_tb_hist(3,1:ii), 'b-', 'LineWidth', 1.5);

        % Cubes
        draw_cube([0;0;0], p.cube_size, eye(3), [0.3 0.3 0.8], 0.4);
        draw_cube(r_now, p.cube_size, eye(3), [0.8 0.3 0.3], 0.6);

        % Triads
        draw_triad([0;0;0], eye(3), p.triad_len, {'x_T','y_T','z_T'}, {'-','-','-'});
        draw_triad(r_now, eye(3), p.triad_len*0.6, {'x_C','y_C','z_C'}, {'--','--','--'});

        % LVLH triad at origin (TB coords)
        R_eci_tb   = lg.R_eci_tb_hist(:,:,ii);
        R_eci_lvlh = lg.R_eci_lvlh_hist(:,:,ii);
        R_tb_lvlh  = R_eci_tb' * R_eci_lvlh;
        draw_triad([0;0;0], R_tb_lvlh, p.triad_len*0.5, {'x_L','y_L','z_L'}, {'-.','-.','-.'});

        xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]'); zlabel('z_{TB} [m]');
        title('Target Body Frame -- 3D');
        view(135, 25);
        lim = max(r_max, cone_L_fixed);
        xlim([-lim lim]); ylim([-10 lim*1.1]); zlim([-lim lim]);

        add_hud(lg, p, ii, Nu);
        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, ii == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 4: 2D XY -- TB frame with fixed cone
    %  ================================================================
    fprintf('  GIF 4: 2D xy TB...\n');
    gif_file = fullfile(results_dir, 'gif_2d_xy_tb.gif');
    fig = figure('Visible','off','Position',[50 50 900 700]);

    for ii = 1:frame_skip:Nt
        clf;
        hold on; grid on; axis equal;
        r_now = lg.r_tb_hist(:, ii);

        % Fixed cone
        yy = linspace(0, y_max, 100);
        fill([cone_k*yy, fliplr(-cone_k*yy)], [yy, fliplr(yy)], ...
            [1 0.95 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot( cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
        plot(-cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
        plot([0 0], [0 y_max], 'k--', 'LineWidth', 2);

        % Cubes + triads
        draw_square_2d([0,0], p.cube_size, [0.3 0.3 0.8]);
        draw_square_2d([r_now(1),r_now(2)], p.cube_size, [0.8 0.3 0.3]);
        draw_triad_2d([0,0], eye(2), p.triad_len*0.5, {'x_T','y_T'}, '-');
        draw_triad_2d([r_now(1),r_now(2)], eye(2), p.triad_len*0.3, {'x_C','y_C'}, '--');

        % Trajectory
        plot(lg.r_tb_hist(1,1:ii), lg.r_tb_hist(2,1:ii), 'b-', 'LineWidth', 1.5);
        plot(lg.r_tb_hist(1,1), lg.r_tb_hist(2,1), 'g^','MarkerSize',10,'MarkerFaceColor','g');

        xlim([-x_max x_max]); ylim([-10 y_max*1.05]);
        xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]');
        title(sprintf('TB Frame -- XY  t=%.1f s', lg.t_hist(ii)));
        add_hud(lg, p, ii, Nu);
        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, ii == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 5: 2D XY -- LVLH frame with rotating cone projection
    %  ================================================================
    fprintf('  GIF 5: 2D xy LVLH...\n');
    gif_file = fullfile(results_dir, 'gif_2d_xy_lvlh.gif');
    fig = figure('Visible','off','Position',[50 50 900 700]);

    lim_l = max([max(abs(lg.r_lvlh_hist(1,:))), max(abs(lg.r_lvlh_hist(2,:))), ...
                 cone_L_fixed]) * 1.3;

    for ii = 1:frame_skip:Nt
        clf;
        hold on; grid on; axis equal;
        r_now_l = lg.r_lvlh_hist(:, ii);

        % Rotating docking axis projected to xy
        d = lg.dock_axis_lvlh_hist(:, ii);
        plot([0 d(1)*cone_L_fixed], [0 d(2)*cone_L_fixed], 'k--', 'LineWidth', 2.5);

        % Rotating cone edges projected to xy (2D cross-section)
        % Project the 3D cone edges onto the LVLH xy plane
        R_eci_tb   = lg.R_eci_tb_hist(:,:,ii);
        R_eci_lvlh = lg.R_eci_lvlh_hist(:,:,ii);
        R_lt       = R_eci_lvlh' * R_eci_tb;
        thetas = linspace(0, 2*pi, p.cone_nfaces+1);
        for ff = 1:p.cone_nfaces
            edge_tb = [cone_k*cos(thetas(ff))*cone_L_fixed; cone_L_fixed; cone_k*sin(thetas(ff))*cone_L_fixed];
            edge_l  = R_lt * edge_tb;
            plot([0 edge_l(1)], [0 edge_l(2)], '-', 'Color', [0.8 0.5 0 0.4], 'LineWidth', 0.8);
        end

        % Cubes + triads
        draw_square_2d([0,0], p.cube_size, [0.3 0.3 0.8]);
        draw_square_2d([r_now_l(1),r_now_l(2)], p.cube_size, [0.8 0.3 0.3]);
        draw_triad_2d([0,0], eye(2), p.triad_len*0.5, {'x_L','y_L'}, '-');
        draw_triad_2d([r_now_l(1),r_now_l(2)], eye(2), p.triad_len*0.3, {'x_C','y_C'}, '--');

        % Trajectory
        plot(lg.r_lvlh_hist(1,1:ii), lg.r_lvlh_hist(2,1:ii), 'b-', 'LineWidth', 1.5);
        plot(lg.r_lvlh_hist(1,1), lg.r_lvlh_hist(2,1), 'g^','MarkerSize',10,'MarkerFaceColor','g');

        xlim([-lim_l lim_l]); ylim([-lim_l lim_l]);
        xlabel('x_{LVLH} [m]'); ylabel('y_{LVLH} [m]');
        title(sprintf('LVLH Frame -- XY  t=%.1f s', lg.t_hist(ii)));
        add_hud(lg, p, ii, Nu);
        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, ii == 1, 0.05);
    end
    close(fig);

    fprintf('  All GIFs saved to %s/\n', results_dir);
end


%% ===== HUD overlay =====
function add_hud(lg, p, ii, Nu)
    r = lg.r_tb_hist(:, ii);
    rad_dev = sqrt(r(1)^2 + r(3)^2);
    y_dist = r(2);
    if y_dist > 0
        los_m = p.cone_k * y_dist - rad_dev;
    else
        los_m = -inf;
    end
    if ii <= Nu
        u_mag = norm(lg.u_hist(:, ii));
    else
        u_mag = 0;
    end
    hud = sprintf(['t=%.0fs  \\omega=[%.3f,%.3f,%.3f]r/s\n' ...
        'u_{max}=%.2f  ||u||=%.3f m/s^2\n' ...
        'y_T=%.1fm  r_{dev}=%.1fm  LOS=%.1fm\n' ...
        'x_0=[%.0f,%.0f,%.0f]m'], ...
        lg.t_hist(ii), p.omega_body, ...
        p.u_max, u_mag, y_dist, rad_dev, los_m, ...
        lg.r_tb_hist(1,1), lg.r_tb_hist(2,1), lg.r_tb_hist(3,1));
    annotation('textbox', [0.01 0.72 0.30 0.26], 'String', hud, ...
        'FitBoxToText','on', 'BackgroundColor',[1 1 1 0.85], ...
        'EdgeColor','k', 'FontSize',8, 'FontName','Courier', ...
        'Interpreter','tex');
end


%% ===== 2D drawing helpers =====
function draw_square_2d(center, half_edge, color)
    s = half_edge;
    cx = center(1); cy = center(2);
    x = cx + s*[-1 1 1 -1 -1];
    y = cy + s*[-1 -1 1 1 -1];
    fill(x, y, color, 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 0.5);
end

function draw_triad_2d(origin, R, len, labels, style)
    colors = {'r', [0 0.7 0]};
    for i = 1:2
        tip = origin(:)' + len * R(:,i)';
        plot([origin(1) tip(1)], [origin(2) tip(2)], style, ...
            'Color', colors{i}, 'LineWidth', 1.5);
        text(tip(1), tip(2), ['  ' labels{i}], 'Color', colors{i}, ...
            'FontSize', 8, 'FontWeight', 'bold');
    end
end


%% ===== GIF frame writer =====
function write_gif_frame(filename, frame, is_first, delay)
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);
    if is_first
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', delay);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
    end
end
