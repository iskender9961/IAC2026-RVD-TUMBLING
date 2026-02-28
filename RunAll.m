%RUNALL  Master pipeline: execute entire analysis suite.
%
%   Runs sequentially:
%     1. Nominal reachability
%     2. Stochastic reachability
%     3. Robust reachability
%     4. Comparative analysis (overlay + tables)
%
%   Note: MC simulations are NOT included (run RunMC or RunMC_HighRes separately).
%
%   Usage:  >> RunAll

this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(this_dir));

run(fullfile(this_dir, 'run_all_reachability.m'));
