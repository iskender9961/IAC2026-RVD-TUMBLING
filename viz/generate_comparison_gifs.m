function generate_comparison_gifs(results, p_all, labels, colors, results_dir)
%GENERATE_COMPARISON_GIFS  Animated GIFs overlaying all Rdu scenarios.
%
%   GIF 1: gif_rdu_tb_topview       -- TB xy with cone, all trajectories
%   GIF 2: gif_rdu_lvlh_rotating    -- LVLH 3D with rotating cone
%   GIF 3: gif_rdu_tb_3d            -- TB 3D perspective
%   GIF 4: gif_rdu_2d_xy_tb         -- TB xy 2D
%   GIF 5: gif_rdu_2d_xy_lvlh       -- LVLH xy 2D
%   GIF 6: gif_rdu_combined         -- Big 2x2 combined overview

    if nargin < 5, results_dir = 'results'; end

    n_cases = length(results);
    p0 = p_all{1};
    cone_k = p0.cone_k;

    % Determine common time grid (use shortest scenario)
    Nt_min = inf;
    for ii = 1:n_cases
        Nt_min = min(Nt_min, length(results{ii}.t_hist));
    end
    frame_skip = max(1, floor(Nt_min / 150));

    % Fixed cone length
    cone_L_fixed = results{1}.r_tb_hist(2, 1) * 1.2;
    y_max = cone_L_fixed;

    % Compute axis limits across all scenarios
    x_max_tb = 0; r_max_lvlh = 0;
    for ii = 1:n_cases
        lg = results{ii};
        x_max_tb = max(x_max_tb, max(abs(lg.r_tb_hist(1,:))));
        x_max_tb = max(x_max_tb, max(abs(lg.r_tb_hist(3,:))));
        r_max_lvlh = max(r_max_lvlh, max(vecnorm(lg.r_lvlh_hist, 2, 1)));
    end
    x_max_tb = max(x_max_tb, cone_k * y_max) * 1.3;
    r_max_lvlh = r_max_lvlh * 1.3;

    %% ================================================================
    %  GIF 1: Top view TB (x_TB vs y_TB)
    %  ================================================================
    fprintf('  Comp GIF 1: TB top view...\n');
    gif_file = fullfile(results_dir, 'gif_rdu_tb_topview.gif');
    fig = figure('Visible','off','Position',[50 50 900 700]);

    for jj = 1:frame_skip:Nt_min
        clf; hold on; grid on; axis equal;

        % Cone (fixed)
        yy = linspace(0, y_max, 100);
        fill([cone_k*yy, fliplr(-cone_k*yy)], [yy, fliplr(yy)], ...
            [1 0.95 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot( cone_k*yy, yy, 'Color',[0.8 0.5 0], 'LineWidth', 1.5);
        plot(-cone_k*yy, yy, 'Color',[0.8 0.5 0], 'LineWidth', 1.5);
        plot([0 0], [0 y_max], 'k--', 'LineWidth', 2);

        h = gobjects(n_cases, 1);
        for ii = 1:n_cases
            lg = results{ii};
            idx = min(jj, length(lg.t_hist));
            h(ii) = plot(lg.r_tb_hist(1,1:idx), lg.r_tb_hist(2,1:idx), '-', ...
                'Color', colors{ii}, 'LineWidth', 1.5);
            plot(lg.r_tb_hist(1,idx), lg.r_tb_hist(2,idx), 'o', ...
                'Color', colors{ii}, 'MarkerSize', 8, 'MarkerFaceColor', colors{ii}, ...
                'HandleVisibility', 'off');
        end

        xlim([-x_max_tb x_max_tb]); ylim([-10 y_max*1.05]);
        xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]');
        title(sprintf('TB Top View -- All Scenarios  t=%.0f s', results{1}.t_hist(jj)));
        legend(h, labels, 'Location', 'northeast', 'FontSize', 7);
        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, jj == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 2: LVLH 3D with rotating cone
    %  ================================================================
    fprintf('  Comp GIF 2: LVLH 3D rotating cone...\n');
    gif_file = fullfile(results_dir, 'gif_rdu_lvlh_rotating.gif');
    fig = figure('Visible','off','Position',[50 50 1000 750]);

    for jj = 1:frame_skip:Nt_min
        clf; hold on; grid on; axis equal;

        % Use scenario 1 for cone/axis (all share same target)
        lg1 = results{1};
        R_eci_tb1   = lg1.R_eci_tb_hist(:,:,jj);
        R_eci_lvlh1 = lg1.R_eci_lvlh_hist(:,:,jj);
        R_lvlh_tb   = R_eci_lvlh1' * R_eci_tb1;

        d = lg1.dock_axis_lvlh_hist(:, jj);
        plot3([0 d(1)*cone_L_fixed], [0 d(2)*cone_L_fixed], [0 d(3)*cone_L_fixed], ...
            'k--', 'LineWidth', 2.5);
        draw_los_tetra([0;0;0], R_lvlh_tb, cone_k, cone_L_fixed, p0.cone_nfaces, [0.8 0.5 0]);
        plot3(0, 0, 0, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'k');

        h = gobjects(n_cases, 1);
        for ii = 1:n_cases
            lg = results{ii};
            idx = min(jj, length(lg.t_hist));
            h(ii) = plot3(lg.r_lvlh_hist(1,1:idx), lg.r_lvlh_hist(2,1:idx), ...
                lg.r_lvlh_hist(3,1:idx), '-', 'Color', colors{ii}, 'LineWidth', 1.5);
            plot3(lg.r_lvlh_hist(1,idx), lg.r_lvlh_hist(2,idx), ...
                lg.r_lvlh_hist(3,idx), 'o', 'Color', colors{ii}, ...
                'MarkerSize', 8, 'MarkerFaceColor', colors{ii}, 'HandleVisibility', 'off');
        end

        xlim([-r_max_lvlh r_max_lvlh]);
        ylim([-r_max_lvlh r_max_lvlh]);
        zlim([-r_max_lvlh r_max_lvlh]);
        xlabel('x_{LVLH} [m]'); ylabel('y_{LVLH} [m]'); zlabel('z_{LVLH} [m]');
        title(sprintf('LVLH -- Rotating Corridor  t=%.0f s', lg1.t_hist(jj)));
        legend(h, labels, 'Location', 'best', 'FontSize', 7);
        view(135, 25);
        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, jj == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 3: TB 3D perspective
    %  ================================================================
    fprintf('  Comp GIF 3: TB 3D...\n');
    gif_file = fullfile(results_dir, 'gif_rdu_tb_3d.gif');
    fig = figure('Visible','off','Position',[50 50 1000 750]);
    lim3d = max(r_max_lvlh, cone_L_fixed);

    for jj = 1:frame_skip:Nt_min
        clf; hold on; grid on; axis equal;

        draw_los_tetra([0;0;0], eye(3), cone_k, cone_L_fixed, p0.cone_nfaces, [0.8 0.5 0]);
        plot3([0 0], [0 cone_L_fixed], [0 0], 'k--', 'LineWidth', 2);
        plot3(0, 0, 0, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'k');

        h = gobjects(n_cases, 1);
        for ii = 1:n_cases
            lg = results{ii};
            idx = min(jj, length(lg.t_hist));
            h(ii) = plot3(lg.r_tb_hist(1,1:idx), lg.r_tb_hist(2,1:idx), ...
                lg.r_tb_hist(3,1:idx), '-', 'Color', colors{ii}, 'LineWidth', 1.5);
            plot3(lg.r_tb_hist(1,idx), lg.r_tb_hist(2,idx), ...
                lg.r_tb_hist(3,idx), 'o', 'Color', colors{ii}, ...
                'MarkerSize', 8, 'MarkerFaceColor', colors{ii}, 'HandleVisibility', 'off');
        end

        xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]'); zlabel('z_{TB} [m]');
        title(sprintf('TB 3D -- All Scenarios  t=%.0f s', results{1}.t_hist(jj)));
        legend(h, labels, 'Location', 'best', 'FontSize', 7);
        view(135, 25);
        xlim([-lim3d lim3d]); ylim([-10 lim3d*1.1]); zlim([-lim3d lim3d]);
        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, jj == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 4: 2D XY -- TB frame
    %  ================================================================
    fprintf('  Comp GIF 4: 2D xy TB...\n');
    gif_file = fullfile(results_dir, 'gif_rdu_2d_xy_tb.gif');
    fig = figure('Visible','off','Position',[50 50 900 700]);

    for jj = 1:frame_skip:Nt_min
        clf; hold on; grid on; axis equal;

        yy = linspace(0, y_max, 100);
        fill([cone_k*yy, fliplr(-cone_k*yy)], [yy, fliplr(yy)], ...
            [1 0.95 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot( cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
        plot(-cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
        plot([0 0], [0 y_max], 'k--', 'LineWidth', 2);

        h = gobjects(n_cases, 1);
        for ii = 1:n_cases
            lg = results{ii};
            idx = min(jj, length(lg.t_hist));
            h(ii) = plot(lg.r_tb_hist(1,1:idx), lg.r_tb_hist(2,1:idx), '-', ...
                'Color', colors{ii}, 'LineWidth', 1.5);
            plot(lg.r_tb_hist(1,idx), lg.r_tb_hist(2,idx), 'o', ...
                'Color', colors{ii}, 'MarkerSize', 8, 'MarkerFaceColor', colors{ii}, ...
                'HandleVisibility', 'off');
        end

        xlim([-x_max_tb x_max_tb]); ylim([-10 y_max*1.05]);
        xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]');
        title(sprintf('TB XY -- All Scenarios  t=%.0f s', results{1}.t_hist(jj)));
        legend(h, labels, 'Location', 'northeast', 'FontSize', 7);
        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, jj == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 5: 2D XY -- LVLH frame with rotating cone
    %  ================================================================
    fprintf('  Comp GIF 5: 2D xy LVLH...\n');
    gif_file = fullfile(results_dir, 'gif_rdu_2d_xy_lvlh.gif');
    fig = figure('Visible','off','Position',[50 50 900 700]);

    lim_l = 0;
    for ii = 1:n_cases
        lg = results{ii};
        lim_l = max(lim_l, max(abs(lg.r_lvlh_hist(1,:))));
        lim_l = max(lim_l, max(abs(lg.r_lvlh_hist(2,:))));
    end
    lim_l = max(lim_l, cone_L_fixed) * 1.3;

    for jj = 1:frame_skip:Nt_min
        clf; hold on; grid on; axis equal;

        % Rotating docking axis (from scenario 1, shared target)
        lg1 = results{1};
        d = lg1.dock_axis_lvlh_hist(:, jj);
        plot([0 d(1)*cone_L_fixed], [0 d(2)*cone_L_fixed], 'k--', 'LineWidth', 2.5);

        % Rotating cone edges
        R_eci_tb1   = lg1.R_eci_tb_hist(:,:,jj);
        R_eci_lvlh1 = lg1.R_eci_lvlh_hist(:,:,jj);
        R_lt = R_eci_lvlh1' * R_eci_tb1;
        thetas = linspace(0, 2*pi, p0.cone_nfaces+1);
        for ff = 1:p0.cone_nfaces
            edge_tb = [cone_k*cos(thetas(ff))*cone_L_fixed; cone_L_fixed; ...
                       cone_k*sin(thetas(ff))*cone_L_fixed];
            edge_l = R_lt * edge_tb;
            plot([0 edge_l(1)], [0 edge_l(2)], '-', 'Color', [0.8 0.5 0 0.4], 'LineWidth', 0.8);
        end

        h = gobjects(n_cases, 1);
        for ii = 1:n_cases
            lg = results{ii};
            idx = min(jj, length(lg.t_hist));
            h(ii) = plot(lg.r_lvlh_hist(1,1:idx), lg.r_lvlh_hist(2,1:idx), '-', ...
                'Color', colors{ii}, 'LineWidth', 1.5);
            plot(lg.r_lvlh_hist(1,idx), lg.r_lvlh_hist(2,idx), 'o', ...
                'Color', colors{ii}, 'MarkerSize', 8, 'MarkerFaceColor', colors{ii}, ...
                'HandleVisibility', 'off');
        end

        xlim([-lim_l lim_l]); ylim([-lim_l lim_l]);
        xlabel('x_{LVLH} [m]'); ylabel('y_{LVLH} [m]');
        title(sprintf('LVLH XY -- All Scenarios  t=%.0f s', lg1.t_hist(jj)));
        legend(h, labels, 'Location', 'best', 'FontSize', 7);
        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, jj == 1, 0.05);
    end
    close(fig);

    %% ================================================================
    %  GIF 6: Big combined 2x2 overview
    %  ================================================================
    fprintf('  Comp GIF 6: Combined overview...\n');
    gif_file = fullfile(results_dir, 'gif_rdu_combined.gif');
    fig = figure('Visible','off','Position',[50 50 1400 1000]);

    for jj = 1:frame_skip:Nt_min
        clf;
        t_now = results{1}.t_hist(jj);

        % --- Panel 1: TB XY (top-left) ---
        subplot(2,2,1); hold on; grid on; axis equal;
        yy = linspace(0, y_max, 100);
        fill([cone_k*yy, fliplr(-cone_k*yy)], [yy, fliplr(yy)], ...
            [1 0.95 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot( cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
        plot(-cone_k*yy, yy, 'Color',[0.8 0.5 0],'LineWidth',1.2);
        plot([0 0], [0 y_max], 'k--', 'LineWidth', 1.5);
        h = gobjects(n_cases, 1);
        for ii = 1:n_cases
            lg = results{ii};
            idx = min(jj, length(lg.t_hist));
            h(ii) = plot(lg.r_tb_hist(1,1:idx), lg.r_tb_hist(2,1:idx), '-', ...
                'Color', colors{ii}, 'LineWidth', 1.5);
            plot(lg.r_tb_hist(1,idx), lg.r_tb_hist(2,idx), 'o', ...
                'Color', colors{ii}, 'MarkerSize', 6, 'MarkerFaceColor', colors{ii}, ...
                'HandleVisibility', 'off');
        end
        xlim([-x_max_tb x_max_tb]); ylim([-10 y_max*1.05]);
        xlabel('x_{TB}'); ylabel('y_{TB}');
        title('TB XY');
        legend(h, labels, 'Location', 'northeast', 'FontSize', 6);

        % --- Panel 2: TB 3D (top-right) ---
        subplot(2,2,2); hold on; grid on; axis equal;
        draw_los_tetra([0;0;0], eye(3), cone_k, cone_L_fixed, p0.cone_nfaces, [0.8 0.5 0]);
        plot3([0 0], [0 cone_L_fixed], [0 0], 'k--', 'LineWidth', 1.5);
        for ii = 1:n_cases
            lg = results{ii};
            idx = min(jj, length(lg.t_hist));
            plot3(lg.r_tb_hist(1,1:idx), lg.r_tb_hist(2,1:idx), ...
                lg.r_tb_hist(3,1:idx), '-', 'Color', colors{ii}, 'LineWidth', 1.5);
            plot3(lg.r_tb_hist(1,idx), lg.r_tb_hist(2,idx), ...
                lg.r_tb_hist(3,idx), 'o', 'Color', colors{ii}, ...
                'MarkerSize', 6, 'MarkerFaceColor', colors{ii});
        end
        xlabel('x_{TB}'); ylabel('y_{TB}'); zlabel('z_{TB}');
        title('TB 3D');
        view(135, 25);
        xlim([-lim3d lim3d]); ylim([-10 lim3d*1.1]); zlim([-lim3d lim3d]);

        % --- Panel 3: LVLH XY (bottom-left) ---
        subplot(2,2,3); hold on; grid on; axis equal;
        lg1 = results{1};
        d = lg1.dock_axis_lvlh_hist(:, jj);
        plot([0 d(1)*cone_L_fixed], [0 d(2)*cone_L_fixed], 'k--', 'LineWidth', 2);
        R_eci_tb1   = lg1.R_eci_tb_hist(:,:,jj);
        R_eci_lvlh1 = lg1.R_eci_lvlh_hist(:,:,jj);
        R_lt = R_eci_lvlh1' * R_eci_tb1;
        thetas = linspace(0, 2*pi, p0.cone_nfaces+1);
        for ff = 1:p0.cone_nfaces
            edge_tb = [cone_k*cos(thetas(ff))*cone_L_fixed; cone_L_fixed; ...
                       cone_k*sin(thetas(ff))*cone_L_fixed];
            edge_l = R_lt * edge_tb;
            plot([0 edge_l(1)], [0 edge_l(2)], '-', 'Color', [0.8 0.5 0 0.4], 'LineWidth', 0.8);
        end
        for ii = 1:n_cases
            lg = results{ii};
            idx = min(jj, length(lg.t_hist));
            plot(lg.r_lvlh_hist(1,1:idx), lg.r_lvlh_hist(2,1:idx), '-', ...
                'Color', colors{ii}, 'LineWidth', 1.5);
            plot(lg.r_lvlh_hist(1,idx), lg.r_lvlh_hist(2,idx), 'o', ...
                'Color', colors{ii}, 'MarkerSize', 6, 'MarkerFaceColor', colors{ii});
        end
        xlim([-lim_l lim_l]); ylim([-lim_l lim_l]);
        xlabel('x_{LVLH}'); ylabel('y_{LVLH}');
        title('LVLH XY');

        % --- Panel 4: Time histories (bottom-right) ---
        subplot(2,2,4); hold on; grid on;
        for ii = 1:n_cases
            lg = results{ii};
            idx = min(jj, length(lg.t_hist));
            plot(lg.t_hist(1:idx), lg.r_tb_hist(2,1:idx), '-', ...
                'Color', colors{ii}, 'LineWidth', 1.5);
        end
        xlabel('Time [s]'); ylabel('y_{TB} [m]');
        title('Along-Axis Distance');
        % vertical time marker
        xline(t_now, 'k--', 'LineWidth', 1);
        xlim([0 results{1}.t_hist(end)]);

        sgtitle(sprintf('R_{\\Delta u} Comparison  --  t = %.0f s', t_now), ...
            'FontSize', 14, 'FontWeight', 'bold');

        drawnow;
        fr = getframe(fig);
        write_gif_frame(gif_file, fr, jj == 1, 0.05);
    end
    close(fig);

    fprintf('  All comparison GIFs saved to %s/\n', results_dir);
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
