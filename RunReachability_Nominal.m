%RUNREACHABILITY_NOMINAL  Entry point: Nominal deterministic reachability.
%
%   Computes forward and backward safe sets for all omega/a_max combos
%   using directional per-constraint erosion + sync range bound.
%
%   Outputs:
%     - data/nominal_reachability.mat
%     - figures/nominal_forward_grid.png
%     - figures/nominal_backward_grid.png
%     - figures/nominal_w*_a*.png (individual)
%
%   Usage:  >> RunReachability_Nominal

this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(this_dir));

run(fullfile(this_dir, 'reachability_nominal', 'run_nominal_reachability.m'));
