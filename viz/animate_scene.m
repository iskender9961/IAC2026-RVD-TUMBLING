function animate_scene(log, p)
%ANIMATE_SCENE  Animate the RVD approach with cubes, triads, and LOS cone.
%   animate_scene(log, p)
%   log : struct with fields:
%         .r_tb_hist     (3 x Nsteps+1)  relative pos in TB
%         .v_tb_hist     (3 x Nsteps+1)  relative vel in TB
%         .u_hist        (3 x Nsteps)    applied inputs in TB
%         .R_eci_tb_hist (3x3 x Nsteps+1) target body rotation matrices
%         .R_eci_lvlh_hist (3x3 x Nsteps+1) LVLH rotation matrices
%         .t_hist        (1 x Nsteps+1)  time vector
%         .status_hist   cell array of solver status strings
%   p   : params struct

    Nsteps = length(log.t_hist) - 1;

    fig = figure('Name', 'RVD Approach Animation', 'Position', [100 100 1200 800]);

    if p.save_mp4
        vid = VideoWriter(p.mp4_file, 'MPEG-4');
        vid.FrameRate = 15;
        open(vid);
    end

    for step = 1:Nsteps+1
        clf;
        ax = axes('Parent', fig);
        hold(ax, 'on');
        grid(ax, 'on');
        axis(ax, 'equal');

        r_tb = log.r_tb_hist(:, step);
        t_now = log.t_hist(step);

        % We draw everything in "relative" coordinates where target is at origin.
        % Chaser is at r_tb (in target body frame).
        % We plot in the target body frame.

        R_eci_tb   = log.R_eci_tb_hist(:,:,step);
        R_eci_lvlh = log.R_eci_lvlh_hist(:,:,step);
        R_tb_eci   = R_eci_tb';

        % Rotation from LVLH to TB (for drawing LVLH triad in TB coords)
        R_tb_lvlh = R_tb_eci * R_eci_lvlh;

        % Chaser body frame: assume chaser is attitude-aligned with TB for simplicity
        R_chaser_in_tb = eye(3);

        % --- Draw target cube at origin ---
        draw_cube([0;0;0], p.cube_size, eye(3), [0.3 0.3 0.8], 0.5);

        % --- Draw chaser cube at r_tb ---
        draw_cube(r_tb, p.cube_size, R_chaser_in_tb, [0.8 0.3 0.3], 0.7);

        % --- Draw target body frame triad (solid lines) ---
        draw_triad([0;0;0], eye(3), p.triad_len, ...
            {'x_T','y_T','z_T'}, {'-','-','-'});

        % --- Draw chaser body frame triad (dashed lines) ---
        draw_triad(r_tb, R_chaser_in_tb, p.triad_len * 0.7, ...
            {'x_C','y_C','z_C'}, {'--','--','--'});

        % --- Draw LVLH frame triad at origin (dash-dot lines) ---
        draw_triad([0;0;0], R_tb_lvlh, p.triad_len * 0.8, ...
            {'x_L','y_L','z_L'}, {'-.','-.','-.'});

        % --- Draw body-fixed LOS cone ---
        draw_los_tetra([0;0;0], eye(3), p.cone_k, p.cone_draw_L, ...
            p.cone_nfaces, [0.8 0.6 0]);

        % --- Draw trajectory so far (in TB frame) ---
        if step > 1
            plot3(log.r_tb_hist(1,1:step), log.r_tb_hist(2,1:step), ...
                  log.r_tb_hist(3,1:step), 'k-', 'LineWidth', 1.5);
        end

        % --- Plot chaser marker ---
        plot3(r_tb(1), r_tb(2), r_tb(3), 'r.', 'MarkerSize', 20);

        % --- Labels and view ---
        xlabel(ax, 'x_{TB} [m]');
        ylabel(ax, 'y_{TB} [m]');
        zlabel(ax, 'z_{TB} [m]');
        title(ax, sprintf('RVD Approach (Target Body Frame)  t = %.1f s', t_now));

        % Compute overlay info
        rad_dev = sqrt(r_tb(1)^2 + r_tb(3)^2);
        y_dist = r_tb(2);
        if y_dist > 0
            los_margin = p.cone_k * y_dist - rad_dev;
        else
            los_margin = -inf;
        end

        if step <= Nsteps
            u_now = log.u_hist(:, step);
            u_mag = norm(u_now);
            solver_stat = log.status_hist{step};
        else
            u_mag = 0;
            solver_stat = 'N/A';
        end

        % --- Text overlay ---
        info_str = sprintf(['Time: %.1f s\n' ...
            'Tumble rate: [%.3f, %.3f, %.3f] rad/s\n' ...
            'y_T dist: %.2f m\n' ...
            'Radial dev: %.2f m\n' ...
            'LOS margin: %.2f m\n' ...
            '||u||: %.4f m/s^2\n' ...
            'Solver: %s'], ...
            t_now, p.omega_body(1), p.omega_body(2), p.omega_body(3), ...
            y_dist, rad_dev, los_margin, u_mag, solver_stat);

        annotation('textbox', [0.01 0.70 0.25 0.28], 'String', info_str, ...
            'FitBoxToText', 'on', 'BackgroundColor', [1 1 1 0.8], ...
            'EdgeColor', 'k', 'FontSize', 9, 'FontName', 'Courier');

        % Set view
        view(ax, 135, 25);
        max_range = max(300, max(abs(r_tb)) * 1.5);
        xlim(ax, [-max_range max_range]);
        ylim(ax, [-50 max_range*1.5]);
        zlim(ax, [-max_range max_range]);

        drawnow;

        if p.save_mp4
            frame = getframe(fig);
            writeVideo(vid, frame);
        end

        pause(0.02);
    end

    if p.save_mp4
        close(vid);
        fprintf('Animation saved to %s\n', p.mp4_file);
    end
end
