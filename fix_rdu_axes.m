%FIX_RDU_AXES  Regenerate fig_rdu_control_tb.png and fig_rdu_positions.png
%  with corrected z-axis limits from saved data.
%
%  Run:  >> fix_rdu_axes

clear; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(this_dir, 'results');

load(fullfile(results_dir, 'rdu_comparison.mat'), ...
     'results', 'p_all', 'rdu_vals', 'labels');

n_cases = length(rdu_vals);
colors  = {'b', 'r', [0 0.6 0], [0.8 0.4 0]};
styles  = {'-', '--', '-.', ':'};

%% ===== Figure: Control TB with z-axis +/-0.2 =====
fig4 = figure('Name','Rdu Comparison - Control TB','Position',[50 50 1100 850]);
for ax = 1:3
    subplot(3,1,ax); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        Nu = size(lg.u_hist, 2);
        plot(lg.t_hist(1:Nu), lg.u_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1);
    end
    ax_lab = {'x','y','z'};
    ylabel(sprintf('u_{%s,TB} [m/s^2]', ax_lab{ax}));
    if ax == 1
        title('Control Input -- Target Body Frame');
        legend(labels, 'Location', 'best', 'FontSize', 8);
    end
    if ax == 3
        xlabel('Time [s]');
        ylim([-0.2 0.2]);
    end
end
saveas(fig4, fullfile(results_dir, 'fig_rdu_control_tb.png'));
fprintf('Saved fig_rdu_control_tb.png with z-axis +/-0.2\n');

%% ===== Figure: Positions with z-axis +/-50 m =====
fig11 = figure('Name','Rdu Comparison - Positions','Position',[50 50 1100 900]);
labels_tb   = {'x_{TB}','y_{TB}','z_{TB}'};
labels_lvlh = {'x_{LVLH}','y_{LVLH}','z_{LVLH}'};
for ax = 1:3
    % Target body
    subplot(3,2,(ax-1)*2+1); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        plot(lg.t_hist, lg.r_tb_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.2);
    end
    ylabel([labels_tb{ax} ' [m]']);
    title([labels_tb{ax} ' -- Target Body']);
    if ax == 1, legend(labels, 'Location', 'best', 'FontSize', 7); end
    if ax == 3
        xlabel('Time [s]');
        ylim([-50 50]);
    end

    % LVLH
    subplot(3,2,(ax-1)*2+2); hold on; grid on;
    for ii = 1:n_cases
        lg = results{ii};
        plot(lg.t_hist, lg.r_lvlh_hist(ax,:), styles{ii}, 'Color', colors{ii}, 'LineWidth', 1.2);
    end
    ylabel([labels_lvlh{ax} ' [m]']);
    title([labels_lvlh{ax} ' -- LVLH']);
    if ax == 3
        xlabel('Time [s]');
        ylim([-50 50]);
    end
end
sgtitle('Position Components -- All Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig11, fullfile(results_dir, 'fig_rdu_positions.png'));
fprintf('Saved fig_rdu_positions.png with z-axis +/-50 m\n');

fprintf('\nDone. Close figures when ready.\n');
