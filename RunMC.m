%RUNMC  Entry point: Monte Carlo feasibility sweep (original resolution).
%
%   Runs the MC sweep at original resolution:
%     - 15 y-levels x 13 x-points = 195 ICs per combo
%     - omega = [1, 2, 3, 4, 5] deg/s
%     - a_max = [0.2, 0.1, 0.05, 0.02] m/s^2
%
%   Outputs saved to results/ directory.
%   Figures regenerated via plot_mc_figures.m.
%
%   Usage:  >> RunMC

this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(this_dir));

run(fullfile(this_dir, 'run_monte_carlo.m'));
