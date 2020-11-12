%% Monte Carlo algorithm

function [count_accept,loops,M_mat,d_cal_mat]...
    = MC_algorithm(N_d,T_step_factor,A_step_factor,tfnx_step_factor,tfnx_init,...
    tfnx_step_nopoly,age,d_obs,d_unc_up,d_unc_dn,d_cal,T_data,A_data,...
    tfnx_init_0,n_iterates,poly,pert,log_like,M_init,A_init,A_bnd,T_init_up,...
    T_init_dn,std_obs1,freq_save,freq_plot)


% ****************************************************************************
% ************************** INITIATION **************************************
% ****************************************************************************

% Set initiation variables
q = 2; % main iteration variable
count_accept = 0; % count number of accepted models
l = 0; % count number of saved results
M = M_init; % set model to start from initial guess
loops = 0; % count total iterations computed to track overall acceptance rate 

% Set Step Size

% The step size for perturbations is based on a scaling of the step factors
% defined in run_settings.m
slide_T = linspace(1,2,N_d)';   % Temperature perturbations grow with time
step_T = zeros(N_d,1) + T_step_factor*slide_T ; % temperature step size
slide_A = linspace(1,2,N_d)';   % Accumulation perturbations grow with time
step_A = zeros(N_d,1) + A_step_factor*slide_A ; % accumulation step size
if poly == 1
    step_tfnx = tfnx_step_factor*tfnx_init; % step size of polynomial thinning function
elseif poly == 0
    step_tfnx = zeros(N_d,1) + tfnx_step_nopoly; % step size of discrete thinning function
end
stepsize = [step_T; step_A; step_tfnx] ;  % full step size vector



% ***************************************************************************
% *********************** In-progress Update Figures ***********************
% ***************************************************************************

% Initiate figures for visualizing in-progress updates of algorithm 

% Figure for visualizing input and output parameters
F = fig('units','inches','width',11,'height',11,'font','Helvetica','fontsize',16);
subplot(321)
plot(age/1000, d_obs(1:N_d),'r','LineWidth',1); % data observation
hold on
plot(age/1000, d_unc_up(1:N_d),'--r')   % data uncertainty bounds
plot(age/1000, d_unc_dn(1:N_d),'--r')   % data uncertainty bounds
plot(age/1000, d_cal(1:N_d),'k') % current model
ylabel('\Deltaage [yr]')
title('Data Space')
subplot(323)
plot(age/1000, d_obs(N_d+1:2*N_d),'r','LineWidth',1);   % data observation
hold on
plot(age/1000, d_unc_up(N_d+1:2*N_d),'--r')   % data uncertainty bounds
plot(age/1000, d_unc_dn(N_d+1:2*N_d),'--r')   % data uncertainty bounds
plot(age/1000, d_cal(N_d+1:2*N_d),'k') % current model
ylabel('Diffusion Length [m]')
subplot(325)
plot(age/1000, d_obs(2*N_d+1:end),'r','LineWidth',1); % data observation
hold on
plot(age/1000, d_unc_up(2*N_d+1:end),'--r')   % data uncertainty bounds
plot(age/1000, d_unc_dn(2*N_d+1:end),'--r')   % data uncertainty bounds
plot(age/1000, d_cal(2*N_d+1:end),'k') % current model
ylabel('Layer Thickness [m]')
xlabel('Age [ka]')

subplot(322)
plot(age/1000, T_data(:,1), 'g','LineWidth',1); % comparison with naive T guess
ylabel('Temperature [C]')
title('Model Space')
subplot(324)
plot(age/1000, A_data, 'g','LineWidth',1);  % comparison with naive A guess
ylabel('Accumulation [m/yr]')
subplot(326)
if poly == 1
plot(age/1000, polyval(tfnx_init,age), 'b','LineWidth',1);  % comparison with naive tfnx guess
elseif poly ==0
  plot(age/1000, tfnx_init_0, 'g','LineWidth',1); % comparison with naive tfnx guess
end
ylabel('Thinning Function')
xlabel('Age [ka]')
drawnow limitrate

% Figure for likelihood values
G = figure; 
title('Likelihood values')


% ***************************************************************************
% ******************************* RUN ITERATIONS ****************************
% ***************************************************************************

while q < n_iterates + 1
    
loops = loops + 1; % keep track of all iterations computed
M_pert = M;        % copy M to apply perturbation


% **************** GENERATE PERTURBATIONS ****************

% **** IF USING POLYNOMIAL THINNING FUNCTION ****
if poly == 1    % use polynomial fit for thinning function
if pert == 1    % single-parameter perturbation 
 
    ran = ceil(rand(1)*(length(M_pert)-(poly_degree-1))); % The index of M which will be perturbed (only allow first index of tfnx polynomial to be perturbed)
    r = normrnd(0,.5); % direction and length from normal distribution
                       % approx max and min -1 to 1, mean 0.

    M_pert(ran) = M_pert(ran) + r*stepsize(ran);    % apply perturbation

elseif pert == 2  % all-parameter perturbation
    % Note: with polynomial thinning function, must perturb T and A
    % parameters separately from the thinning function.
    
    % 2/3 chance of perturbing T, A; 1/3 chance of perturbing tfnx
    P_variable = rand;
    
if P_variable < 0.66     % perturb T and A 
    
% Create red noise for temperature perturbation
stdv = 1; % standard deviation
ar1 = 0.999; % memory
white = stdv*randn(N_d,1); % white noise
white(1) = 0; % set end point to 0
white(end) = 0; % set end point to 0
red_T = white; % create red noise from white noise
for ii = 2:N_d
    red_T(ii) = red_T(ii-1)*ar1+0.1*randn;    % red noise
end
red_smooth_T = lowpass(red_T,1/3000,1/250); % lowpass frequencies below 3000 years

% Create red noise for accumulation perturbation
red_A = white; % use same white noise as above
for ii = 2:N_d
    red_A(ii) = red_A(ii-1)*ar1+0.01*randn;    % red noise
end
red_smooth_A = lowpass(red_A,1/3000,1/250); % lowpass frequencies below 3000 years

% Apply perturbations to T and A
M_pert(1:N_d) = M_pert(1:N_d) + red_smooth_T.*stepsize(1:N_d);    % T perturbation
M_pert(N_d+1:2*N_d) = M_pert(N_d+1:2*N_d) + red_smooth_A.*stepsize(N_d+1:2*N_d);    % A perturbation

elseif P_variable >= 0.67 % perturb thinning function
    ran = ceil(poly_degree*rand(1)); % index to be perturbed
    r = normrnd(0,.5); % direction and length from normal distribution
                       % approx max and min -1 to 1, mean 0.
    M_pert(2*N_d+ran) = M_pert(2*N_d+ran) + r*stepsize(2*N_d+ran);    % apply perturbation to tfnx
end   
end % Single-parameter perturbation selection


% **** IF USING DISCRETE THINNING FUNCTION ****
elseif poly == 0
if pert == 1  % single-parameter perturbation 
    ran = ceil(rand(1)*(length(M_pert))); % The index of M which will be perturbed
    r = normrnd(0,.5); % direction and length from normal distribution
                       % approx max and min -1 to 1, mean 0.
    M_pert(ran) = M_pert(ran) + r*stepsize(ran);    % apply perturbation

elseif pert == 2  % all-parameter perturbation
  
% Create red noise for temperature perturbation
stdv = 1; % standard deviation
ar1 = 0.9; % memory
white = stdv*randn(N_d,1); % white noise
red_T = white; % create red noise from white noise
for ii = 2:N_d
    red_T(ii) = red_T(ii-1)*ar1+0.1*randn;    % red noise
end
red_smooth_T = lowpass(red_T,1/3000,1/250); % lowpass frequencies below 3000 years

% Create red noise for accumulation perturbation
red_A = white; % create red noise from white noise
for ii = 2:N_d
    red_A(ii) = red_A(ii-1)*ar1+0.1*randn;    % red noise
end
red_smooth_A = lowpass(red_A,1/3000,1/250); % lowpass frequencies below 3000 years

% Create red noise for thinning function perturbation
red_tfnx = white; % create red noise from white noise
for ii = 2:N_d
    red_tfnx(ii) = red_tfnx(ii-1)*ar1+0.1*randn;    % red noise
end
red_smooth_tfnx = lowpass(red_tfnx,1/10000,1/250); % lowpass frequencies below 10000 years
red_smooth_tfnx(1:3) = 0; % do not allow top of thinning function to vary

% Apply perturbations to each parameter
M_pert(1:N_d) = M_pert(1:N_d) + red_smooth_T.*stepsize(1:N_d);    % apply perturbation to T
M_pert(N_d+1:2*N_d) = M_pert(N_d+1:2*N_d) + red_smooth_A.*stepsize(N_d+1:2*N_d);    % apply perturbation to A
M_pert(2*N_d+1:3*N_d) = M_pert(2*N_d+1:3*N_d) + red_smooth_tfnx.*stepsize(2*N_d+1:3*N_d);    % apply perturbation to tfnx

end % type of perturbation selection
end % choice between polynomial or discrete fit for thinning function


% *************** RULES FOR PERTURBATIONS ****************    
% Require that perturbations have resulted in acceptable model inputs.
% Reject perturbations if bound test isn't satisfied.

% Accumulation must be within model bounds
A_bnd_test1 = find(M_pert(N_d+1:N_d*2) > A_init + A_bnd);
A_bnd_test2 = find(M_pert(N_d+1:N_d*2) < A_init - A_bnd);
if A_bnd_test1 > 0
        continue            % reject if out of bounds
elseif A_bnd_test2 > 0
        continue            % reject if out of bounds
end
    
% Temperature must be within model bounds
T_bnd_test1 = find(M_pert(1:N_d) > T_init_up);
T_bnd_test2 = find(M_pert(1:N_d) < T_init_dn);
if T_bnd_test1 > 0
        continue            % reject if out of bounds
elseif T_bnd_test2 > 0
        continue            % reject if out of bounds
end

    
% *************** RUN FORWARD MODEL ***************
% Calculate outputs based on the perturbed inputs

[d_cal] = icecore_forward(M_pert,N_d,age,poly);

% *************** EVALUATE MODEL OUTPUT ***************
% Store likelihood function in vector (value will be overwritten if 
% perturbation is not accepted).
log_like(q) = likelihood(d_cal, d_obs, std_obs1);


% **************** ACCEPT OR REJECT PERTURBATION **************** 
    frac = exp((log_like(q) - log_like(q-1)));  % calculate ratio of improvement on last iteration
    P = min(1,frac);    % probability of acceptance based on improvement ratio
    num = rand(1);
    if num<=P   % ACCEPT: M_pert is a better model or we met probability criterion
        M = M_pert;   % accept this model for next iteration
        count_accept = count_accept + 1; % track number of accepted models
    if mod(q,freq_save) == 0   % every freq_save models, save the model and the data output
        l=l+1;      % track number of saved models
        M_mat(:,l) = M;    % save model
        d_cal_mat(:,l) = d_cal;     % save d_cal (calculated data parameters)
        if mod(q,freq_plot) == 0 % every freq_plot models, update the plots
            figure(F)           % update input/output figure
            subplot(321)        
            hold on; plot(age/1000, d_cal(1:N_d))
            subplot(323)
            hold on; plot(age/1000, d_cal(N_d+1:2*N_d))
            subplot(325)
            hold on; plot(age/1000, d_cal(2*N_d+1:3*N_d))
            subplot(322)
            hold on; plot(age/1000, M(1:N_d))
            subplot(324)
            hold on; plot(age/1000, M(N_d+1:2*N_d))
            subplot(326)
            if poly == 1
                hold on; plot(age/1000, polyval(M(2*N_d+1:end),age))
            elseif poly ==0 
                hold on; plot(age/1000, M(2*N_d+1:3*N_d))  
            end
            drawnow limitrate
        
            figure(G)           % update likelihood value figure
            hold on
            plot(q,log_like(q),'x')
            drawnow limitrate
        end
    end
    q = q + 1; % only advance if we have accepted iteration
    end % if num<=P

% display completion percentage while computing
if mod(q,n_iterates/100) == 0  % display percentage every 1%
    disp([num2str(q./n_iterates*100) '% complete'])     
end  
end % while q < n_iterates

end