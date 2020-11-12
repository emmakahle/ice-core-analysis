%% Main Code: Monte Carlo Bayesian Inference Analysis - South Pole Ice Core

% Purpose: Produce climate data sets from measurements of the South Pole
% ice core

% Code Requirements: MATLAB license, curve fitting toolbox, signal processing
% toolbox, statistics and machine learning toolbox, symbolic math toolbox

% Input measurements: high-resolution water isotopes, water isotope
% diffusion length, annual-layer thickness, depth-age scale for ice core,
% delta-age (difference between ice and gas age at given depth)

% Output data sets: temperature history, accumulation rate history,
% thinning function for ice core record

%% Start timers and counters

% Start timer to track length of run
tic
% Start iteration counter
q = 1; 

%% Import Run Settings

% Select run settings in run_settings.m function.

[n_iterates,pert,poly,poly_degree,freq_save,freq_plot,...
    A_bnd,T_step_factor,A_step_factor,tfnx_step_factor,tfnx_step_nopoly] = run_settings;

%% Import data and initiate variables for perturbation

% Import ice core data using load_data.m function.
[depth,age,Dage_obs,Dage_up,Dage_dn,sigma_obs,sigma_up,sigma_dn,...
    lambda_obs,T_init_0,T_data,A_init_0,A_data,tfnx_init_0, tfnx_CB] = load_data;

% Create vector of all data observations
d_obs = [Dage_obs; sigma_obs; lambda_obs];     % column vector of all measured observables
N_d = length(d_obs)/3;   % number of data points in grid

% If using a polynomial perturbation for thinning function, create polynomial 
% coefficients for thinning function
if poly == 1 
    tfnx_P = polyfit(age,tfnx_init_0,poly_degree)';    % generate polynomial coefficients for thinning function
end

% Initiate variables to perturb
T_init = T_init_0(:,1); 
T_init_dn = T_init_0(:,2);
T_init_up = T_init_0(:,3);
A_init = A_init_0; 
if poly ==1
    tfnx_init = tfnx_P; 
elseif poly ==0
    tfnx_init = tfnx_init_0;
end

% Create vector of entire model to perturb
M_init = [T_init; A_init; tfnx_init]; 

%% Define Data Uncertainty

% Uncertainty is calculated from the estimate of diffusion length (1
% sigma)
% Uncertainty is calculated for Dage from Epifanio et al. 2020 (2 sigma)
% Uncertainty is 3% to 10% for layer thickness (2 sigma)

[d_unc_up,d_unc_dn,std_obs1] = data_uncertainty(N_d,age,d_obs,Dage_up,...
    Dage_obs,sigma_up,sigma_obs,Dage_dn,sigma_dn);

%% Calculate and Evaluate Initial Guess
% Use the forward model (icecore_forward) to calculate the initial guess
% from our initial inputs. Calculate likelihood function to evaluate the
% initial guess.

% Run forward model
[d_cal] = icecore_forward(M_init,N_d,age,poly);

% Evaluate likelihood (i.e., how well does the current iteration output
% match the data observations?)
log_like = NaN(n_iterates,1);   % initiate vector to hold likelihood values
log_like(q) = likelihood(d_cal, d_obs, std_obs1);

%% Monte Carlo Algorithm
% Run iterations through algorithm. Evaluate each iteration and determine
% whether to accept or reject.

% Run alogrithm
[count_accept,loops,M_mat,d_cal_mat]...
    = MC_algorithm(N_d,T_step_factor,A_step_factor,tfnx_step_factor,tfnx_init,...
    tfnx_step_nopoly,age,d_obs,d_unc_up,d_unc_dn,d_cal,T_data,A_data,...
    tfnx_init_0,n_iterates,poly,pert,log_like,M_init,A_init,A_bnd,T_init_up,...
    T_init_dn,std_obs1,freq_save,freq_plot);

% Print total run time and acceptance percentage
elapsed_time = toc
acceptance_rate = count_accept/loops

%% Histogram Results
% Create and plot histograms (as gray shading) showing a posteriori result
% for each parameter

% Create back up variables for results
M_mat_full = M_mat;
d_cal_mat_backup = d_cal_mat;

% Allow user to manually determine the burn-in time based on when likelihood values
% flatten out
prompt = 'Visually, what is the burn-in index?';
burn_in = input(prompt);
% burn_in = 10000; % alternatively, set burn-in time to expected value
burn_in_in = burn_in/freq_save; % divide by save frequency to get index for M_mat and d_cal_mat

% Calculate and plot histogram results
hist_results(M_mat_full,d_cal_mat,burn_in_in,age,N_d,d_obs,d_unc_up,d_unc_dn,...
    T_data,T_init,T_init_0,A_data,A_init,A_bnd,tfnx_CB,tfnx_init_0);


%% Summary Results
% Calculate and plot mean and standard deviation of a posteriori results
% for each parameter

% Calculate and plot summary results
[mean_M,std_M,mean_D,std_D]=summary_results(M_mat_full,d_cal_mat,...
    age,N_d,d_obs,d_unc_up,d_unc_dn,T_data,T_init,T_init_0,A_data,A_init,...
    A_bnd,tfnx_CB,tfnx_init_0);

%% Save Final Climate and Ice Dynamics Histories
% Outputs are temperature, accumulation rate, and ice thinning function
% Write outputs as .txt files with columns of age, mean, and standard
% deviation

T_result = [age,mean_M(1:N_d),std_M(1:N_d)];  % temperature result
A_result = [age,mean_M(N_d+1:2*N_d),std_M(N_d+1:2*N_d)];  % accumulation result
tfnx_result = [age,mean_M(2*N_d+1:3*N_d),std_M(2*N_d+1:3*N_d)];  % thinning function result

dlmwrite('temperature.txt',T_result);
dlmwrite('accumulation_rate.txt',A_result);
dlmwrite('thinning_function.txt',tfnx_result);
