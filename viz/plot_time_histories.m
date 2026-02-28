function plot_time_histories(log, p)
%PLOT_TIME_HISTORIES  Plot time histories of key variables.
%   plot_time_histories(log, p)

    Nsteps = length(log.t_hist) - 1;
    t = log.t_hist;

    % Compute derived quantities
    y_T = log.r_tb_hist(2, :);
    rad_dev = sqrt(log.r_tb_hist(1,:).^2 + log.r_tb_hist(3,:).^2);
    u_mag = sqrt(sum(log.u_hist.^2, 1));

    % ---- Figure: Approach profile ----
    figure('Name', 'Time Histories', 'Position', [150 150 900 700]);

    subplot(4,1,1);
    plot(t, y_T, 'b-', 'LineWidth', 1.5);
    ylabel('y_T [m]'); grid on;
    title('Target-Body Frame Approach Profile');

    subplot(4,1,2);
    plot(t, rad_dev, 'r-', 'LineWidth', 1.5);
    ylabel('\surd(x_T^2+z_T^2) [m]'); grid on;
    title('Radial Deviation from Docking Axis');

    subplot(4,1,3);
    hold on;
    plot(t(1:end-1), log.u_hist(1,:), 'r-', 'LineWidth', 1);
    plot(t(1:end-1), log.u_hist(2,:), 'g-', 'LineWidth', 1);
    plot(t(1:end-1), log.u_hist(3,:), 'b-', 'LineWidth', 1);
    ylabel('u_{TB} [m/s^2]'); grid on;
    legend('u_x','u_y','u_z','Location','best');
    title('Control Inputs (Target Body Frame)');

    subplot(4,1,4);
    plot(t(1:end-1), u_mag, 'k-', 'LineWidth', 1.5);
    ylabel('||u|| [m/s^2]'); xlabel('Time [s]'); grid on;
    title('Control Magnitude');

    % ---- Figure: 3D trajectory in TB frame ----
    figure('Name', '3D Approach in TB Frame', 'Position', [200 200 800 700]);
    hold on; grid on; axis equal;
    plot3(log.r_tb_hist(1,:), log.r_tb_hist(2,:), log.r_tb_hist(3,:), ...
        'k-', 'LineWidth', 2);
    plot3(log.r_tb_hist(1,1), log.r_tb_hist(2,1), log.r_tb_hist(3,1), ...
        'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot3(log.r_tb_hist(1,end), log.r_tb_hist(2,end), log.r_tb_hist(3,end), ...
        'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    plot3(0, 0, 0, 'bp', 'MarkerSize', 15, 'MarkerFaceColor', 'b');

    % Draw LOS cone wireframe
    draw_los_tetra([0;0;0], eye(3), p.cone_k, max(y_T)*1.1, ...
        p.cone_nfaces, [0.8 0.6 0]);

    xlabel('x_{TB} [m]'); ylabel('y_{TB} [m]'); zlabel('z_{TB} [m]');
    title('3D Approach Trajectory in Target Body Frame');
    legend('Trajectory','Start','End','Target','Location','best');
    view(135, 25);
end
