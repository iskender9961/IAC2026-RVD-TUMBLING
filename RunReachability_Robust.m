%RUNREACHABILITY_ROBUST  Entry point: Robust bounded-disturbance reachability.
%
%   Computes safe sets with worst-case disturbance tightening for all combos.
%   Guarantees feasibility for ALL disturbances in the bounded set.
%
%   Outputs:
%     - data/robust_reachability.mat
%     - figures/robust_grid.png
%
%   Usage:  >> RunReachability_Robust

this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(this_dir));

run(fullfile(this_dir, 'reachability_robust', 'run_robust_reachability.m'));
