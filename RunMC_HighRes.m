%RUNMC_HIGHRES  Entry point: High-resolution Monte Carlo sweep.
%
%   Runs MC at increased grid density (SAME domain bounds):
%     - 25 y-levels x 21 x-points = 525 ICs per combo
%     - Same omega/a_max sweep as original
%     - Reproducible seed = 2026
%
%   Outputs saved to results/mc_hires_*.mat.
%
%   Usage:  >> RunMC_HighRes

this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(this_dir));

run(fullfile(this_dir, 'run_mc_hires.m'));
