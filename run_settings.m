%% General Run Settings

% choose settings for MC_analysis run

function [n_iterates,pert,poly,poly_degree,freq_save,freq_plot,...
    A_bnd, T_step_factor,A_step_factor,tfnx_step_factor,tfnx_step_nopoly] = run_settings

% Determine length of run.
n_iterates = 1000; % number of accepted iterations to complete run (suggested value 100000)

% Choose frequency of saving and plotting results. Saving frequency should
% reflect degree of autocorrelation. In other words, you want to save
% infrequently enough that neighboring saved results are not correlated.
% For a full run, freq_save = 300 meets this requirement. Frequency of
% plotting determines how often the results will be plotted on the
% in-progress data visualization. Freq_plot does not impact the saved
% results.
freq_save = 3; % choose how often to save accepted models
freq_plot = 9; % how often to plot accepted models. must be a multiple of freq_save

% Choose type of perturbation. Single-parameter perturbation only alters
% one perturbation per iteration. All-parameter perturbation alters every
% parameter in each iteration. Suggested setting is all-parameter (2)
pert = 2; % choose single-parameter (1) or all-parameter (2) perturbation

% Choose type of perturbation for thinning function. If poly=1, a polynomial
% perturbation uses a polynomial of degree poly_degree to perturb the
% thinning function input. If poly=0, independent perturbations are made to
% each parameter in the thinning function.
poly = 0; % set polynomial fit for thinning function. 0 = off, 1 = on.
poly_degree = 15; % set the order of the polynomial for tfnx perturbations

% Set bounds for accumulation rate model space as deviation from the
% initial guess
A_bnd = 0.02;

% Set stepsize factors. These values determine the base size of step used
% in each perturbation for each parameter.
if pert ==  1       % use these values when single-parameter perturbation is selected
T_step_factor = 0.5;
A_step_factor = 0.005;
tfnx_step_factor = 10^-5;
tfnx_step_nopoly = .005;
elseif pert == 2    % use these values when all-parameter perturbation is selected
T_step_factor = 0.05;
A_step_factor = 0.0005;
tfnx_step_factor = 10^-5;
tfnx_step_nopoly = .005;
end

end         % end function