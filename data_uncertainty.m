%% Data Uncertainty

% Uncertainty is calculated from the estimate of diffusion length (1
% sigma)
% Uncertainty is calculated for Dage from Epifanio et al. 2020 (2 sigma)
% Uncertainty is 3% to 10% for layer thickness (2 sigma)


function [d_unc_up,d_unc_dn,std_obs1] = data_uncertainty(N_d,age,d_obs,Dage_up,...
    Dage_obs,sigma_up,sigma_obs,Dage_dn,sigma_dn)

% Separate uncertainty for layer thickness into 3 regimes: Holocene,
% glacial-interglacial transition, and glacial periods. sliding_percent
% defines different uncertainties in each regime.
in1 = find(age>12000, 1, 'first');
in2 = find(age>15000, 1, 'first');
sliding_percent = [linspace(0.03,0.03,length(age(1:in1)))';linspace(0.03,0.1,length(age(in1+1:in2)))';linspace(0.1,0.1,length(age(in2+1:end)))'];

% Define uncertainty bounds vector for all parameters. Diffusion length
% uncertainty is multiplied by 2 because it is imported as only 1 sigma.
d_unc_up = d_obs + [Dage_up-Dage_obs; (sigma_up-sigma_obs)*2; sliding_percent.*d_obs(2*N_d+1:end)]; 
d_unc_dn = d_obs - [Dage_obs-Dage_dn; (sigma_obs-sigma_dn)*2; sliding_percent.*d_obs(2*N_d+1:end)];

% Define vectors for standard deviation of each d_obs parameter
std_obs2 = d_unc_up-d_unc_dn;   % 2-sigma uncertainty at each data point
std_obs1 = std_obs2/2;          % 1-sigma uncertainty at each data point


end