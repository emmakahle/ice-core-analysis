%% Load South Pole Ice Core Data

function [depth,age,Dage,Dage_up,Dage_dn,d18o_sigma,d18o_sigma_up,...
    d18o_sigma_dn,lambda,T_init,T_init_var,A_init,A_init_var,tfnx_init,...
    tfnx_CB] = load_data()

datalength = 212; % length of data

% *****************************************************************
% *********************** Time Scale ******************************
% *****************************************************************

% Ice Timescale
ts_data = importdata('SP19_Depth_Age.csv');
ts_depth = ts_data.data(:,1); % Depth of timescale
ts_age = ts_data.data(:,2); % age of timescale

% Gas Timescale
gas_age_data = load('SP_age_v4_prelim08192019.mat'); 
% preliminary depth, ice age, gas age
depth_age_gasage = [gas_age_data.SP_age_v4(:,1), gas_age_data.SP_age_v4(:,2),gas_age_data.SP_age_v4(:,3)]; 

% Interpolate gas ages to ice timescale
ts_gas_age = interp1(depth_age_gasage(:,1),depth_age_gasage(:,3),ts_depth);
ts = [ts_depth, ts_age, ts_gas_age];

% *****************************************************************
% *********************** Diffusion Length ************************
% *****************************************************************

% Load diffusion length record (3 isotopes combined to d18o estimate)
% Assign variables starting at 5th index to skip the modern firn where 
% diffusion is still in progress.
sigma = importdata('SPC_sigmaFINAL_CFAandSOLIDcorrected_depth_age_sigma18_up_dn.txt');
depth = sigma(5:datalength,1); 
age = sigma(5:datalength,2); 
d18o_sigma = sigma(5:datalength,3); 
d18o_sigma_up = sigma(5:datalength,4);
d18o_sigma_dn = sigma(5:datalength,5);


% *****************************************************************
% *********************** Temperature *****************************
% *****************************************************************

% **** Import water isotope data for naive estimate of temperature ****
d18O_data = importdata('SPIceCoreDatahalfcmAverage01-08-2020.txt');
depth_wi_i = d18O_data.data(:,1); % depth of water isotopes
d18O_i = d18O_data.data(:,2); % water isotope values

% **** Average temperature to match depth values defined above ****
bin_edges = zeros(length(depth)+1,1);   % create edges to define averaging bins
for n = 2:length(depth)
bin_edges(n) = (depth(n-1)+depth(n))/2;
end
bin_edges(end) = depth(end);

d18O_avg = NaN(length(depth),1); % Preallocate vector for average values

for n = 1:length(bin_edges)-1 % Take mean of each data set within bins defined by depths
d18O_avg(n) = nanmean(d18O_i(find(depth_wi_i>=bin_edges(n),1,'first'):find(depth_wi_i>=bin_edges(n+1),1,'first')));
end
d18O_avg(end) = d18O_avg(end-1); % Get rid of last value, which is NaN

% **** Correct for glacial-interglacial changes in water isotope value of
% seawater ****
SWcorr = importdata('SW_correction_kyrBP_meanoceand18O.csv');
SW_age_i = SWcorr(:,1);
SW_d18O_i = SWcorr(:,2);
% Interpolate to correct age points
SW_d18O = interp1(SW_age_i*1000, SW_d18O_i, age);
% Apply seawater correction
d18O_SWcorr = d18O_avg-SW_d18O;

% **** Create temperature estimate and initial guess from water isotopes

% Initiate temperature vector
T_init_var = NaN(length(depth),3); % variable starting guess (will turn into step function) 
% Scale water isotopes with surface slope 0.8 +/- 0.2
% Choose intercepts such that modern T=-51
T_init_var(:,1) = 1/1.2*d18O_SWcorr-8;
T_init_var(:,2) = 1/0.7*d18O_SWcorr+19; % lower bound
T_init_var(:,3) = 1/1.5*d18O_SWcorr-15; % upper bound

% Create step_function for initial temperature guess from temperature
% estimate. Calculate three periods (Holocene, transition, glacial).
T_init = NaN(length(depth),3); % initiate vector
p1a = polyfit(age(1:43),T_init_var(1:43,1),1);  % Create linear version of T_init(:,1)
T_init(1:43,1) = p1a(1)*age(1:43)+p1a(2);
p1b = polyfit(age(44:71),T_init_var(44:71,1),1);
T_init(44:71,1) = p1b(1)*age(44:71)+p1b(2);
p1c = polyfit(age(72:end),T_init_var(72:end,1),1);
T_init(72:end,1) = p1c(1)*age(72:end)+p1c(2);
p2a = polyfit(age(1:43),T_init_var(1:43,2),1);  % Create linear version of T_init(:,2)
T_init(1:43,2) = p2a(1)*age(1:43)+p2a(2);
p2b = polyfit(age(44:71),T_init_var(44:71,2),1);
T_init(44:71,2) = p2b(1)*age(44:71)+p2b(2);
p2c = polyfit(age(72:end),T_init_var(72:end,2),1);
T_init(72:end,2) = p2c(1)*age(72:end)+p2c(2);
p3a = polyfit(age(1:43),T_init_var(1:43,3),1);  % Create linear version of T_init(:,3)
T_init(1:43,3) = p3a(1)*age(1:43)+p3a(2);
p3b = polyfit(age(44:71),T_init_var(44:71,3),1);
T_init(44:71,3) = p3b(1)*age(44:71)+p3b(2);
p3c = polyfit(age(72:end),T_init_var(72:end,3),1);
T_init(72:end,3) = p3c(1)*age(72:end)+p3c(2);

% *****************************************************************
% ************************* Delta Age *****************************
% *****************************************************************

% Find Dage values halfway in between gas and ice timescales

% **** Import final Dage from Epifanio et al. 2020 ****
Dage_i = importdata('depth_Dage_Dageunc_101520.txt');
JE_depth = Dage_i(:,1);
JE_Dage = Dage_i(:,2:3); % Dage and uncertainty
JE_iceage = interp1(ts(:,1), ts(:,2), JE_depth); % interp ice age to Dage depths
JE_gasage = interp1(ts(:,1), ts(:,3), JE_depth); % interp gas age to Dage depths

% **** Create ages for Delta age that are halfway between gas and ice
% timescales ****
ts_half = (JE_iceage+JE_gasage)/2;
for nn = 1088:length(ts_half)-1
age_temp_in = find(JE_iceage<=ts_half(nn),1,'last');
Dage_depth(nn) = JE_depth(age_temp_in,1);
end
depth_Dage_i = [Dage_depth',JE_Dage(1:end-1,:)];
    
% Average values at same depths
[depth_Dage_avg,~,z] = unique(depth_Dage_i(:,1),'stable');
depth_Dage_avg(:,2) = accumarray(z,depth_Dage_i(:,2),[],@mean);
depth_Dage_avg(:,3) = accumarray(z,depth_Dage_i(:,3),[],@mean);

% Interpolate back to main grid
Dage_unc(:,1) = interp1(depth_Dage_avg(1088:end,1),depth_Dage_avg(1088:end,2),depth); % Dage on mid point timescale
Dage_unc(:,2) = interp1(depth_Dage_avg(1088:end,1),depth_Dage_avg(1088:end,3),depth); % Dage on mid point timescale

% **** Assign variables to output ****
Dage_unc(1:7,:) = repmat(Dage_unc(8,:),7,1); % Set first seven values to avoid NaN values
Dage = Dage_unc(:,1);
Dage_up = Dage_unc(:,1) + Dage_unc(:,2);
Dage_dn = Dage_unc(:,1) - Dage_unc(:,2);


% *****************************************************************
% **** Layer Thickness, Accumulation, Thinning Function ***********
% *****************************************************************

% ***** Import layer thickness data and naive estimates for thinning
% function *****

% Layer Thickness observations from ice timescale
ts_lambda = ts_data.data(:,7);

% Naive estimate of thinning function (Dansgaard-Johnsen 1969 model)
tf_data = importdata('DJ_hpt2_age_depth_tf.txt');
tf_depth = tf_data(:,2);
tf_i = tf_data(:,3);

% Naive estimate of thinning function from d15N data
d = load('SP_ties_TJ_09242019.mat');
cspline = csaps([zeros(10000,1);d.SP_ties(:,10)],[ones(10000,1);d.SP_ties(:,11)],0.00001,0:10:1750);
tfnx_CB = interp1(0:10:1750, cspline, depth); %interpolate to my depth scale

% **** Average layer thickness and thinning function estimates to defined
% depth values ****
bin_edges = zeros(length(depth)+1,1);   % create edges to define averaging bins
for n = 2:length(depth)
bin_edges(n) = (depth(n-1)+depth(n))/2;
end
bin_edges(end) = bin_edges(end-1)+3.25;

% Preallocate vectors
lambda = zeros(length(depth),1);
tfnx_DJ = zeros(length(depth),1);

% Average each data set within bins defined by depths
for n = 1:length(bin_edges)-1
    lambda(n) = mean(ts_lambda(find(ts_depth>=bin_edges(n),1,'first'):find(ts_depth>=bin_edges(n+1),1,'first')));
    tfnx_DJ(n) = mean(tf_i(find(tf_depth>=bin_edges(n),1,'first'):find(tf_depth>=bin_edges(n+1),1,'first')));
end

tfnx_init_var = tfnx_DJ; % variable starting guess
lambda(1) = (lambda(1)+lambda(2))/2; % Clean first value of layer thickness

% **** Create initial guess for accumulation from layer thickness and naive
% thinning function estimate ****
A_init_var = lambda./tfnx_DJ;

% Create step-function version of accumulation for starting guess
A_init = nan(length(lambda),1); % initiate vector
% Split step-function into three periods (Holocene, transition, glacial)
p1a = polyfit(age(1:43),A_init_var(1:43),1);  % create linear version of A_init
A_init(1:43) = p1a(1)*age(1:43)+p1a(2);
p1b = polyfit(age(44:71),A_init_var(44:71),1);
A_init(44:71) = p1b(1)*age(44:71)+p1b(2);
p1c = polyfit(age(72:end),A_init_var(72:end),1);
A_init(72:end) = p1c(1)*age(72:end)+p1c(2);

% Create polynomial starting guess for tfnx_init
p = polyfit(age,tfnx_init_var,4); 
tfnx_init = p(1)*age.^4+p(2)*age.^3+p(3)*age.^2+p(4)*age+p(5);

end
