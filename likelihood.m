%% Likelihood Function
% Calculate the log of the likelihood value, which estimates misfit between
% observations and model output

% Based on Mosegaard et al. 2002

function [log_like] = likelihood(d_cal, d_obs, std_obs1)

S_m = sum(abs(d_cal - d_obs)./std_obs1);    % Misfit function

log_like = -S_m;  % Log of the likelihood function

end

