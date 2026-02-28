function generate_gifs(lg, p, results_dir)
%GENERATE_GIFS  Create animated GIFs from simulation results.
%   generate_gifs(lg, p, results_dir)
%
%   GIF 1: Top view in TB frame (x_TB vs y_TB) with LOS cone, trajectory, docking axis
%   GIF 2: LVLH frame with rotating cone, trajectory, docking axis
%   GIF 3: 3D TB frame rotation view
%   GIF 4-6: 2D animated projections (xy, xz, yz) in TB frame
%   GIF 7-9: 2D animated projections (xy, xz, yz) in LVLH frame

    if nargin < 3, results_dir = 'results'; end

    Nt = length(lg.t_hist);
    cone_k = p.cone_k;
    frame_skip = max(1, floor(Nt / 150));  % ~150 frames max for GIF size

    %% ================================================================
    %  GIF 1: Top view TB frame (x_TB vs y_TB)
    %  ================================================================
    fprintf('  GIF 1: TB top view (x_TB vs y_TB)...\n');
    gif_file = fullfile(results_dir, 'gif_tb_topview.gif');
    fig = figure('Name','GIF TB Top View','Position',[50 50 800 600],'Visible','off');

    y_max = max(lg.r_tb_hist(2,:)) * 1.3;
    x_max = max(abs(lg.r_tb_hist(1,:))) * 1.5;
    x_max = max(x_max, cone_k * y_max * 1.1);

    for ii = 1:frame_skip:Nt
        clf;
        hold on; grid on; axis equal;
        % Cone edges
        yy = linspace(0, y_max, 100);
        fill([ cone_k*yy, fliplr(-cone_k*yy)], [yy, fliplr(yy)], ...
            [1 0.95 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot( cone_k*yy, yy, 'Color',[0.8 0.5 0], 'LineWidth', 1.5);
        plot(-cone_k*yy, yy, 'Color',[0.8 0.5 0], 'LineWidth', 1.5);
        % Docking axis
        plot([0 0], [0 y_max], 'k--', 'LineWidth', 2);
        % Trajectory history
        plot(lg.r_tb_hist(1,1:ii), lg.r_tb_hist(2,1:ii), 'b-', 'LineWidth', 1.5);
        % Current position
        plot(lg.r_tb_hist(1,ii), lg.r_tb_hist(2,ii), 'ro', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'r');
        % Start
        plot(lg.r_tb_hist(1,1), lg.r_tb_hist(2,1), 'g^', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'g');
        % Target
        plot(0, 0, 'kp', 'MarkerSize', 15, 'MarkerFaceColor', 'k');

        xlim([-x_max x_max]); ylim([-10 y_max]);
        xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]');
        title(sprintf('Target Body Frame -- Top View  t = %.1f s', lg.t_hist(ii)));

        drawnow;
        frame = getframe(fig);
        write_gif_frame(gif_file, frame, ii == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 2: LVLH frame with rotating cone + trajectory
    %  ================================================================
    fprintf('  GIF 2: LVLH 3D with rotating cone...\n');
    gif_file = fullfile(results_dir, 'gif_lvlh_rotating_cone.gif');
    fig = figure('Name','GIF LVLH Rotating Cone','Position',[50 50 900 700],'Visible','off');

    r_max = max(vecnorm(lg.r_lvlh_hist, 2, 1)) * 1.3;

    for ii = 1:frame_skip:Nt
        clf;
        hold on; grid on; axis equal;

        % Trajectory history in LVLH
        plot3(lg.r_lvlh_hist(1,1:ii), lg.r_lvlh_hist(2,1:ii), ...
              lg.r_lvlh_hist(3,1:ii), 'b-', 'LineWidth', 1.5);
        plot3(lg.r_lvlh_hist(1,ii), lg.r_lvlh_hist(2,ii), ...
              lg.r_lvlh_hist(3,ii), 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
        plot3(lg.r_lvlh_hist(1,1), lg.r_lvlh_hist(2,1), ...
              lg.r_lvlh_hist(3,1), 'g^', 'MarkerSize', 10, 'MarkerFaceColor','g');
        plot3(0,0,0, 'kp', 'MarkerSize', 15, 'MarkerFaceColor','k');

        % Instantaneous docking axis in LVLH
        d = lg.dock_axis_lvlh_hist(:, ii);
        dock_len = max(50, lg.r_tb_hist(2,ii) * 1.5);
        plot3([0 d(1)*dock_len], [0 d(2)*dock_len], [0 d(3)*dock_len], ...
            'k--', 'LineWidth', 2.5);

        % Rotating cone in LVLH
        R_eci_tb   = lg.R_eci_tb_hist(:,:,ii);
        R_eci_lvlh = lg.R_eci_lvlh_hist(:,:,ii);
        R_lvlh_tb  = R_eci_lvlh' * R_eci_tb;
        draw_cone_poly([0;0;0], R_lvlh_tb, cone_k, dock_len, p.cone_nfaces, [0.8 0.5 0]);

        xlim([-r_max r_max]); ylim([-r_max r_max]); zlim([-r_max r_max]);
        xlabel('x_{LVLH} [m]'); ylabel('y_{LVLH} [m]'); zlabel('z_{LVLH} [m]');
        title(sprintf('LVLH Frame -- Rotating Cone  t = %.1f s', lg.t_hist(ii)));
        view(135, 25);

        drawnow;
        frame = getframe(fig);
        write_gif_frame(gif_file, frame, ii == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 3: 3D TB frame with cone, trajectory building up
    %  ================================================================
    fprintf('  GIF 3: TB 3D perspective...\n');
    gif_file = fullfile(results_dir, 'gif_tb_3d.gif');
    fig = figure('Name','GIF TB 3D','Position',[50 50 900 700],'Visible','off');

    for ii = 1:frame_skip:Nt
        clf;
        hold on; grid on; axis equal;

        % Cone
        y_max_dyn = max(50, lg.r_tb_hist(2,ii)*1.5);
        draw_cone_poly([0;0;0], eye(3), cone_k, y_max_dyn, p.cone_nfaces, [0.8 0.5 0]);
        % Docking axis
        plot3([0 0], [0 y_max_dyn], [0 0], 'k--', 'LineWidth', 2);
        % Trajectory
        plot3(lg.r_tb_hist(1,1:ii), lg.r_tb_hist(2,1:ii), ...
              lg.r_tb_hist(3,1:ii), 'b-', 'LineWidth', 1.5);
        plot3(lg.r_tb_hist(1,ii), lg.r_tb_hist(2,ii), ...
              lg.r_tb_hist(3,ii), 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
        plot3(0,0,0,'kp','MarkerSize',15,'MarkerFaceColor','k');
        % Target and chaser cubes
        draw_cube([0;0;0], p.cube_size, eye(3), [0.3 0.3 0.8], 0.4);
        r_now = lg.r_tb_hist(:,ii);
        draw_cube(r_now, p.cube_size, eye(3), [0.8 0.3 0.3], 0.6);

        xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]'); zlabel('z_{TB} [m]');
        title(sprintf('Target Body Frame  t = %.1f s', lg.t_hist(ii)));
        view(135, 25);
        r_lim = max(r_max, 50);
        xlim([-r_lim r_lim]); ylim([-10 r_lim*1.5]); zlim([-r_lim r_lim]);

        drawnow;
        frame = getframe(fig);
        write_gif_frame(gif_file, frame, ii == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 4-6: 2D animated projections in TB frame
    %  ================================================================
    generate_2d_gifs(lg.r_tb_hist, lg.t_hist, 'TB', cone_k, p, results_dir, frame_skip);

    %% ================================================================
    %  GIF 7-9: 2D animated projections in LVLH frame
    %  ================================================================
    generate_2d_gifs_lvlh(lg.r_lvlh_hist, lg.t_hist, lg.dock_axis_lvlh_hist, ...
        p, results_dir, frame_skip);

    fprintf('  All GIFs saved to %s/\n', results_dir);
end


%% ===== 2D GIF helper (TB / CB frames with static cone) =====
function generate_2d_gifs(r_hist, t_hist, tag, cone_k, p, results_dir, frame_skip)
    planes   = {'xy','xz','yz'};
    idx_map  = {[1 2],[1 3],[2 3]};
    ax_labs  = {{'x','y'},{'x','z'},{'y','z'}};
    Nt = length(t_hist);

    y_max = max(r_hist(2,:)) * 1.3;

    for pp = 1:3
        ix = idx_map{pp};
        al = ax_labs{pp};
        gif_file = fullfile(results_dir, sprintf('gif_2d_%s_%s.gif', planes{pp}, lower(tag)));
        fig = figure('Visible','off','Position',[50 50 750 600]);
        fprintf('  GIF 2D %s %s...\n', planes{pp}, tag);

        x_lim = max(abs(r_hist(ix(1),:))) * 1.5;
        y_lim = max(abs(r_hist(ix(2),:))) * 1.3;
        if pp == 1, x_lim = max(x_lim, cone_k*y_max*1.1); end

        for ii = 1:frame_skip:Nt
            clf;
            hold on; grid on; axis equal;

            % Cone / docking axis projections
            if pp == 1  % xy
                yy = linspace(0, y_max, 100);
                fill([cone_k*yy, fliplr(-cone_k*yy)], [yy, fliplr(yy)], ...
                    [1 0.95 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                plot( cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
                plot(-cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
                plot([0 0], [0 y_max], 'k--', 'LineWidth', 2);
            elseif pp == 3  % yz
                yy = linspace(0, y_max, 100);
                fill([yy, fliplr(yy)], [cone_k*yy, fliplr(-cone_k*yy)], ...
                    [1 0.95 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                plot(yy,  cone_k*yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
                plot(yy, -cone_k*yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
                plot([0 y_max], [0 0], 'k--', 'LineWidth', 2);
            elseif pp == 2  % xz
                for yslice = linspace(20, y_max*0.8, 4)
                    th = linspace(0,2*pi,60);
                    rr = cone_k*yslice;
                    plot(rr*cos(th), rr*sin(th), ':', 'Color', [0.8 0.5 0], 'LineWidth', 0.6);
                end
            end

            % Trajectory up to current step
            plot(r_hist(ix(1),1:ii), r_hist(ix(2),1:ii), 'b-', 'LineWidth', 1.5);
            plot(r_hist(ix(1),ii), r_hist(ix(2),ii), 'ro', 'MarkerSize', 10, ...
                'MarkerFaceColor','r');
            plot(r_hist(ix(1),1), r_hist(ix(2),1), 'g^', 'MarkerSize', 10, ...
                'MarkerFaceColor','g');
            plot(0, 0, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'k');

            if pp == 1
                xlim([-x_lim x_lim]); ylim([-10 y_lim]);
            else
                xlim([-x_lim x_lim]); ylim([-y_lim y_lim]);
            end
            xlabel(sprintf('%s_{%s} [m]', al{1}, tag));
            ylabel(sprintf('%s_{%s} [m]', al{2}, tag));
            title(sprintf('%s Frame -- %s  t=%.1f s', tag, upper(planes{pp}), t_hist(ii)));

            drawnow;
            frame = getframe(fig);
            write_gif_frame(gif_file, frame, ii == 1, 0.05);
        end
        close(fig);
    end
end


%% ===== 2D GIF helper (LVLH frame with rotating docking axis) =====
function generate_2d_gifs_lvlh(r_hist, t_hist, dock_hist, p, results_dir, frame_skip)
    planes   = {'xy','xz','yz'};
    idx_map  = {[1 2],[1 3],[2 3]};
    ax_labs  = {{'x','y'},{'x','z'},{'y','z'}};
    Nt = length(t_hist);

    for pp = 1:3
        ix = idx_map{pp};
        al = ax_labs{pp};
        gif_file = fullfile(results_dir, sprintf('gif_2d_%s_lvlh.gif', planes{pp}));
        fig = figure('Visible','off','Position',[50 50 750 600]);
        fprintf('  GIF 2D %s LVLH...\n', planes{pp});

        x_lim = max(abs(r_hist(ix(1),:))) * 1.5;
        y_lim = max(abs(r_hist(ix(2),:))) * 1.3;
        dock_len = max(abs(r_hist(2,:))) * 1.2;

        for ii = 1:frame_skip:Nt
            clf;
            hold on; grid on; axis equal;

            % Instantaneous docking axis (rotating)
            d = dock_hist(:, ii);
            plot([0 d(ix(1))*dock_len], [0 d(ix(2))*dock_len], ...
                'k--', 'LineWidth', 2);

            % Trajectory
            plot(r_hist(ix(1),1:ii), r_hist(ix(2),1:ii), 'b-', 'LineWidth', 1.5);
            plot(r_hist(ix(1),ii), r_hist(ix(2),ii), 'ro', 'MarkerSize', 10, ...
                'MarkerFaceColor','r');
            plot(r_hist(ix(1),1), r_hist(ix(2),1), 'g^', 'MarkerSize', 10, ...
                'MarkerFaceColor','g');
            plot(0, 0, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'k');

            xlim([-max(x_lim,50) max(x_lim,50)]);
            ylim([-max(y_lim,50) max(y_lim,50)]);
            xlabel(sprintf('%s_{LVLH} [m]', al{1}));
            ylabel(sprintf('%s_{LVLH} [m]', al{2}));
            title(sprintf('LVLH -- %s  t=%.1f s', upper(planes{pp}), t_hist(ii)));

            drawnow;
            frame = getframe(fig);
            write_gif_frame(gif_file, frame, ii == 1, 0.05);
        end
        close(fig);
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
