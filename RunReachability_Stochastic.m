%RUNREACHABILITY_STOCHASTIC  Entry point: Stochastic chance-constrained reachability.
%
%   Computes safe sets with Gaussian noise tightening for all combos.
%   Uses Bonferroni-corrected chance constraints (alpha = 0.05, 95% confidence).
%
%   Outputs:
%     - data/stochastic_reachability.mat
%     - figures/stochastic_grid.png
%
%   Usage:  >> RunReachability_Stochastic

this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(this_dir));

run(fullfile(this_dir, 'reachability_stochastic', 'run_stochastic_reachability.m'));
