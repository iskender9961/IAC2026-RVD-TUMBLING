%RUN_ALL_REACHABILITY  Master script: run all reachability analyses.
%
%   Executes in order:
%     1. Nominal reachability (forward + backward)
%     2. Stochastic reachability
%     3. Robust reachability
%     4. Comparative analysis (overlay + tables)
%
%   Does NOT run Monte Carlo (use run_monte_carlo.m separately).
%   Does NOT modify any existing files.
%
%   Run:  >> run_all_reachability

clear; close all; clc;
CALLED_FROM_MASTER = true; %#ok<NASGU>

MASTER_ROOT = fileparts(mfilename('fullpath'));
if isempty(MASTER_ROOT), MASTER_ROOT = pwd; end

% Add all reachability paths
addpath(genpath(MASTER_ROOT));

fprintf('========================================================\n');
fprintf('  REACHABILITY ANALYSIS SUITE — IAC 2026\n');
fprintf('========================================================\n\n');

t_master = tic;

%% Step 1: Nominal
fprintf('>>> Step 1/4: Running nominal reachability...\n\n');
run(fullfile(MASTER_ROOT, 'reachability_nominal', 'run_nominal_reachability.m'));

%% Step 2: Stochastic
fprintf('\n\n>>> Step 2/4: Running stochastic reachability...\n\n');
run(fullfile(MASTER_ROOT, 'reachability_stochastic', 'run_stochastic_reachability.m'));

%% Step 3: Robust
fprintf('\n\n>>> Step 3/4: Running robust reachability...\n\n');
run(fullfile(MASTER_ROOT, 'reachability_robust', 'run_robust_reachability.m'));

%% Step 4: Comparison
fprintf('\n\n>>> Step 4/4: Running comparative analysis...\n\n');
run(fullfile(MASTER_ROOT, 'reachability_compare', 'run_comparison.m'));

elapsed_master = toc(t_master);
fprintf('\n========================================================\n');
fprintf('  ALL REACHABILITY ANALYSES COMPLETE\n');
fprintf('  Total elapsed: %.1f s (%.1f min)\n', elapsed_master, elapsed_master/60);
fprintf('========================================================\n');
