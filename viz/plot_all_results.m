function plot_all_results(lg, p, results_dir)
%PLOT_ALL_RESULTS  Generate all static plots and save as PNG.
%   plot_all_results(lg, p, results_dir)

    if nargin < 3, results_dir = 'results'; end

    t = lg.t_hist;
    Nt = length(t);
    Nu = size(lg.u_hist, 2);

    % ---- Cone boundary helpers ----
    cone_k = p.cone_k;

    %% ================================================================
    %  1. 3D TRAJECTORY -- TARGET BODY FRAME
    %  ================================================================
    fig = figure('Name','3D Target Body Frame','Position',[50 50 900 750],'Visible','on');
    hold on; grid on; axis equal;
    % Trajectory
    plot3(lg.r_tb_hist(1,:), lg.r_tb_hist(2,:), lg.r_tb_hist(3,:), ...
        'b-', 'LineWidth', 1.8);
    plot3(lg.r_tb_hist(1,1), lg.r_tb_hist(2,1), lg.r_tb_hist(3,1), ...
        'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot3(lg.r_tb_hist(1,end), lg.r_tb_hist(2,end), lg.r_tb_hist(3,end), ...
        'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    plot3(0,0,0,'kp','MarkerSize',15,'MarkerFaceColor','k');
    % Instantaneous docking axis (+yT) as black dashed line
    y_max_tb = max(lg.r_tb_hist(2,:)) * 1.2;
    plot3([0 0], [0 y_max_tb], [0 0], 'k--', 'LineWidth', 2);
    % LOS cone wireframe
    draw_cone_poly([0;0;0], eye(3), cone_k, y_max_tb, p.cone_nfaces, [0.8 0.5 0]);
    xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]'); zlabel('z_{TB} [m]');
    title('3D Trajectory -- Target Body Frame');
    legend('Trajectory','Start','End','Target','Docking axis (+y_{TB})', ...
           'Location','best');
    view(135, 25);
    saveas(fig, fullfile(results_dir, 'fig_3d_target_body.png'));

    %% ================================================================
    %  2. 3D TRAJECTORY -- LVLH FRAME
    %  ================================================================
    fig = figure('Name','3D LVLH Frame','Position',[100 50 900 750],'Visible','on');
    hold on; grid on; axis equal;
    plot3(lg.r_lvlh_hist(1,:), lg.r_lvlh_hist(2,:), lg.r_lvlh_hist(3,:), ...
        'b-', 'LineWidth', 1.8);
    plot3(lg.r_lvlh_hist(1,1), lg.r_lvlh_hist(2,1), lg.r_lvlh_hist(3,1), ...
        'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot3(lg.r_lvlh_hist(1,end), lg.r_lvlh_hist(2,end), lg.r_lvlh_hist(3,end), ...
        'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    plot3(0,0,0,'kp','MarkerSize',15,'MarkerFaceColor','k');
    % Instantaneous docking axis at each sample (show every 10th as short arrows)
    step = max(1, floor(Nt/20));
    for ii = 1:step:Nt
        d = lg.dock_axis_lvlh_hist(:,ii);
        base = [0;0;0];
        tip  = d * y_max_tb * 0.3;
        plot3([base(1) tip(1)], [base(2) tip(2)], [base(3) tip(3)], ...
            'k--', 'LineWidth', 0.8);
    end
    xlabel('x_{LVLH} [m]'); ylabel('y_{LVLH} [m]'); zlabel('z_{LVLH} [m]');
    title('3D Trajectory -- LVLH Frame');
    legend('Trajectory','Start','End','Target','Location','best');
    view(135, 25);
    saveas(fig, fullfile(results_dir, 'fig_3d_lvlh.png'));

    %% ================================================================
    %  3. 3D TRAJECTORY -- CHASER BODY FRAME
    %  ================================================================
    fig = figure('Name','3D Chaser Body Frame','Position',[150 50 900 750],'Visible','on');
    hold on; grid on; axis equal;
    plot3(lg.r_cb_hist(1,:), lg.r_cb_hist(2,:), lg.r_cb_hist(3,:), ...
        'b-', 'LineWidth', 1.8);
    plot3(lg.r_cb_hist(1,1), lg.r_cb_hist(2,1), lg.r_cb_hist(3,1), ...
        'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot3(lg.r_cb_hist(1,end), lg.r_cb_hist(2,end), lg.r_cb_hist(3,end), ...
        'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    plot3(0,0,0,'kp','MarkerSize',15,'MarkerFaceColor','k');
    plot3([0 0], [0 y_max_tb], [0 0], 'k--', 'LineWidth', 2);
    draw_cone_poly([0;0;0], eye(3), cone_k, y_max_tb, p.cone_nfaces, [0.8 0.5 0]);
    xlabel('x_{CB} [m]'); ylabel('y_{CB} [m]'); zlabel('z_{CB} [m]');
    title('3D Trajectory -- Chaser Body Frame');
    legend('Trajectory','Start','End','Target','Docking axis','Location','best');
    view(135, 25);
    saveas(fig, fullfile(results_dir, 'fig_3d_chaser_body.png'));

    %% ================================================================
    %  4-6. 2D PROJECTIONS -- TARGET BODY FRAME
    %  ================================================================
    plot_2d_projections(lg.r_tb_hist, 'Target Body', 'TB', cone_k, ...
        y_max_tb, p, results_dir, lg);

    %% ================================================================
    %  7-9. 2D PROJECTIONS -- LVLH FRAME
    %  ================================================================
    plot_2d_projections(lg.r_lvlh_hist, 'LVLH', 'LVLH', cone_k, ...
        y_max_tb, p, results_dir, lg);

    %% ================================================================
    %  10-12. 2D PROJECTIONS -- CHASER BODY FRAME
    %  ================================================================
    plot_2d_projections(lg.r_cb_hist, 'Chaser Body', 'CB', cone_k, ...
        y_max_tb, p, results_dir, lg);

    %% ================================================================
    %  13. TIME VS COST
    %  ================================================================
    fig = figure('Name','Time vs Cost','Position',[50 400 800 400],'Visible','on');
    plot(t(1:Nu), lg.cost_hist(1:Nu), 'b-', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Stage cost J_k');
    title('MPC Stage Cost vs Time'); grid on;
    saveas(fig, fullfile(results_dir, 'fig_cost_vs_time.png'));

    %% ================================================================
    %  14. TIME HISTORIES -- POSITION IN ALL FRAMES (subplots)
    %  ================================================================
    labels_tb   = {'x_{TB}','y_{TB}','z_{TB}'};
    labels_lvlh = {'x_{LVLH}','y_{LVLH}','z_{LVLH}'};
    labels_cb   = {'x_{CB}','y_{CB}','z_{CB}'};

    fig = figure('Name','Position Time Histories','Position',[50 50 1100 900],'Visible','on');
    for ax = 1:3
        % Target body
        subplot(3,3,(ax-1)*3+1);
        plot(t, lg.r_tb_hist(ax,:), 'b-', 'LineWidth', 1.2); hold on;
        if ax == 2  % y_TB: add docking axis ref and reference
            plot(t, lg.ref_hist(2,:), 'r--', 'LineWidth', 1);
            yline(0, 'k--', 'LineWidth', 1.5);
            legend('y_{TB}','ref','dock axis','Location','best');
        end
        xlabel('Time [s]'); ylabel([labels_tb{ax} ' [m]']);
        title([labels_tb{ax} ' -- Target Body']); grid on;

        % LVLH
        subplot(3,3,(ax-1)*3+2);
        plot(t, lg.r_lvlh_hist(ax,:), 'b-', 'LineWidth', 1.2);
        xlabel('Time [s]'); ylabel([labels_lvlh{ax} ' [m]']);
        title([labels_lvlh{ax} ' -- LVLH']); grid on;

        % Chaser body
        subplot(3,3,(ax-1)*3+3);
        plot(t, lg.r_cb_hist(ax,:), 'b-', 'LineWidth', 1.2);
        xlabel('Time [s]'); ylabel([labels_cb{ax} ' [m]']);
        title([labels_cb{ax} ' -- Chaser Body']); grid on;
    end
    sgtitle('Position Components vs Time in All Frames');
    saveas(fig, fullfile(results_dir, 'fig_time_positions_all.png'));

    %% ================================================================
    %  15. TIME HISTORIES -- VELOCITY IN ALL FRAMES (subplots)
    %  ================================================================
    fig = figure('Name','Velocity Time Histories','Position',[100 50 1100 900],'Visible','on');
    vlabels_tb   = {'vx_{TB}','vy_{TB}','vz_{TB}'};
    vlabels_lvlh = {'vx_{LVLH}','vy_{LVLH}','vz_{LVLH}'};
    vlabels_cb   = {'vx_{CB}','vy_{CB}','vz_{CB}'};
    for ax = 1:3
        subplot(3,3,(ax-1)*3+1);
        plot(t, lg.v_tb_hist(ax,:), 'b-', 'LineWidth', 1.2);
        xlabel('Time [s]'); ylabel([vlabels_tb{ax} ' [m/s]']);
        title([vlabels_tb{ax} ' -- Target Body']); grid on;

        subplot(3,3,(ax-1)*3+2);
        plot(t, lg.v_lvlh_hist(ax,:), 'b-', 'LineWidth', 1.2);
        xlabel('Time [s]'); ylabel([vlabels_lvlh{ax} ' [m/s]']);
        title([vlabels_lvlh{ax} ' -- LVLH']); grid on;

        subplot(3,3,(ax-1)*3+3);
        plot(t, lg.v_cb_hist(ax,:), 'b-', 'LineWidth', 1.2);
        xlabel('Time [s]'); ylabel([vlabels_cb{ax} ' [m/s]']);
        title([vlabels_cb{ax} ' -- Chaser Body']); grid on;
    end
    sgtitle('Velocity Components vs Time in All Frames');
    saveas(fig, fullfile(results_dir, 'fig_time_velocities_all.png'));

    %% ================================================================
    %  16. CONTROL INPUTS (TB and LVLH)
    %  ================================================================
    fig = figure('Name','Control Inputs','Position',[50 50 1000 700],'Visible','on');
    t_u = t(1:Nu);
    subplot(2,1,1);
    hold on;
    plot(t_u, lg.u_hist(1,:), 'r-', 'LineWidth', 1.2);
    plot(t_u, lg.u_hist(2,:), 'g-', 'LineWidth', 1.2);
    plot(t_u, lg.u_hist(3,:), 'b-', 'LineWidth', 1.2);
    plot(t_u, sqrt(sum(lg.u_hist.^2,1)), 'k-', 'LineWidth', 1.5);
    ylabel('u [m/s^2]'); xlabel('Time [s]');
    title('Control Input -- Target Body Frame');
    legend('u_{xTB}','u_{yTB}','u_{zTB}','||u||','Location','best');
    grid on;

    subplot(2,1,2);
    hold on;
    plot(t_u, lg.u_lvlh_hist(1,:), 'r-', 'LineWidth', 1.2);
    plot(t_u, lg.u_lvlh_hist(2,:), 'g-', 'LineWidth', 1.2);
    plot(t_u, lg.u_lvlh_hist(3,:), 'b-', 'LineWidth', 1.2);
    plot(t_u, sqrt(sum(lg.u_lvlh_hist.^2,1)), 'k-', 'LineWidth', 1.5);
    ylabel('u [m/s^2]'); xlabel('Time [s]');
    title('Control Input -- LVLH Frame');
    legend('u_{xL}','u_{yL}','u_{zL}','||u||','Location','best');
    grid on;
    saveas(fig, fullfile(results_dir, 'fig_control_inputs.png'));

    %% ================================================================
    %  17. APPROACH PROFILE (y_T, radial deviation, LOS margin)
    %  ================================================================
    fig = figure('Name','Approach Profile','Position',[50 50 900 700],'Visible','on');
    rad_dev = sqrt(lg.r_tb_hist(1,:).^2 + lg.r_tb_hist(3,:).^2);
    y_T = lg.r_tb_hist(2,:);
    los_margin = cone_k * y_T - rad_dev;

    subplot(3,1,1);
    plot(t, y_T, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, lg.ref_hist(2,:), 'r--', 'LineWidth', 1);
    ylabel('y_{TB} [m]'); title('Distance Along Docking Axis'); grid on;
    legend('y_{TB}','reference','Location','best');

    subplot(3,1,2);
    plot(t, rad_dev, 'r-', 'LineWidth', 1.5);
    ylabel('\surd(x_{TB}^2+z_{TB}^2) [m]');
    title('Radial Deviation from Docking Axis'); grid on;

    subplot(3,1,3);
    plot(t, los_margin, 'm-', 'LineWidth', 1.5); hold on;
    yline(0, 'k--', 'LineWidth', 1);
    ylabel('LOS margin [m]'); xlabel('Time [s]');
    title('LOS Cone Margin (>0 = inside cone)'); grid on;
    saveas(fig, fullfile(results_dir, 'fig_approach_profile.png'));

    fprintf('  All static plots saved to %s/\n', results_dir);
end


%% ===== Helper: 2D projections for a given frame =====
function plot_2d_projections(r_hist, frame_name, tag, cone_k, y_max, p, results_dir, lg)
    planes = {'xy', 'xz', 'yz'};
    idx_map = {[1 2], [1 3], [2 3]};
    ax_labels = {{'x','y'}, {'x','z'}, {'y','z'}};

    for pp = 1:3
        ix = idx_map{pp};
        al = ax_labels{pp};

        fig = figure('Name', sprintf('2D %s -- %s', planes{pp}, frame_name), ...
                     'Position', [50+pp*50 50+pp*50 750 600], 'Visible', 'on');
        hold on; grid on; axis equal;

        % Trajectory
        plot(r_hist(ix(1),:), r_hist(ix(2),:), 'b-', 'LineWidth', 1.5);
        plot(r_hist(ix(1),1), r_hist(ix(2),1), 'go', 'MarkerSize', 10, 'MarkerFaceColor','g');
        plot(r_hist(ix(1),end), r_hist(ix(2),end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor','r');
        plot(0, 0, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'k');

        % Docking axis and cone projections (only for TB/CB frames)
        if strcmp(tag,'TB') || strcmp(tag,'CB')
            if pp == 1  % xy plane: docking axis is along +y, cone edges in x
                plot([0 0], [0 y_max], 'k--', 'LineWidth', 2);
                yy = linspace(0, y_max, 100);
                plot( cone_k*yy, yy, 'Color', [0.8 0.5 0], 'LineWidth', 1.2);
                plot(-cone_k*yy, yy, 'Color', [0.8 0.5 0], 'LineWidth', 1.2);
            elseif pp == 3  % yz plane: docking axis along +y, cone edges in z
                plot([0 y_max], [0 0], 'k--', 'LineWidth', 2);
                yy = linspace(0, y_max, 100);
                plot(yy,  cone_k*yy, 'Color', [0.8 0.5 0], 'LineWidth', 1.2);
                plot(yy, -cone_k*yy, 'Color', [0.8 0.5 0], 'LineWidth', 1.2);
            elseif pp == 2  % xz plane: show cone circle at some y-slices
                for yslice = [50, 100, 150, 200]
                    th = linspace(0, 2*pi, 60);
                    rr = cone_k * yslice;
                    plot(rr*cos(th), rr*sin(th), ':', 'Color', [0.8 0.5 0], 'LineWidth', 0.8);
                end
            end
        end

        % For LVLH: show docking axis at start and end
        if strcmp(tag, 'LVLH')
            if isfield(lg, 'dock_axis_lvlh_hist')
                d0 = lg.dock_axis_lvlh_hist(:, 1);
                dN = lg.dock_axis_lvlh_hist(:, end);
                L = max(abs(r_hist(ix(1),:))) * 1.2;
                plot([0 d0(ix(1))*L], [0 d0(ix(2))*L], 'k--', 'LineWidth', 1.5);
                plot([0 dN(ix(1))*L], [0 dN(ix(2))*L], 'k-.', 'LineWidth', 1.5);
            end
        end

        xlabel(sprintf('%s_{%s} [m]', al{1}, tag));
        ylabel(sprintf('%s_{%s} [m]', al{2}, tag));
        title(sprintf('2D %s Projection -- %s Frame', upper(planes{pp}), frame_name));
        legend('Trajectory','Start','End','Target','Location','best');
        saveas(fig, fullfile(results_dir, sprintf('fig_2d_%s_%s.png', planes{pp}, lower(tag))));
    end
end
