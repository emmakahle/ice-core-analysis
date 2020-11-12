%% Herron and Langway Model
% Herron and Langway 1980 firn densification model

% INPUTS:
% Average T in degrees celcius
% Average Accumulation in meters w.e. per yr

% OUTPUTS:
% Depth profile (m)
% Density profile (kg/m^3)
% Age profile (year)

function [depth,rho,age,drho_dt] = HL_model(rho_0,T,A)

% Define parameters
rho_i = .917; % Mg/m^3
R = 8.314; % J K-1 mol-1
T_K = T+273.15; % convert to Kelvin

% Rate constants
k_0 = 11*exp(-10160/(R*T_K));
k_1 = 575*exp(-21400/(R*T_K));
h_crit = 1/(rho_i.*k_0)*(log(.550/(rho_i - .550)) - log(rho_0/(rho_i-rho_0)));
a = 1;
b = 0.5;

% Define depth grid for two zones
h_1 = 1:.1:h_crit;
h_2 = h_crit+.1:.1:2000;     % expanded range of deepest depth to avoid error with too cold of temperatures
depth = [h_1,h_2];

% Density profiles for two density zones (correct to ice equivalent)
Z_0 = exp(rho_i*k_0.*h_1 + log(rho_0/(rho_i-rho_0))); % densification parameter
rho_h1 = (rho_i*Z_0)./(1+Z_0); % zone 1 densification profile

Z_1 = exp(rho_i*k_1.*(h_2 - h_crit)./A^0.5 + log(.550/(rho_i - .550))); % densification parameter
rho_h2 = (rho_i*Z_1)./(1+Z_1); % zone 2 densification profile

% Use temperature-dependent close-off density
rho_co = (1/(rho_i*1000) + 6.95*10^(-7)*T_K - 4.3*10^(-5))^(-1)/1000; % Martinerie et al. 1992, 1994

% Calculate density profile
rho = [rho_h1, rho_h2]; % full density profile
co_depth_in = find(rho >= rho_co, 1, 'first');   % find close off depth, for rho_co = 804.3 kg/m3

% Trim density and depth profiles
rho = rho(1:co_depth_in);   % cut off density vector at close off depth
depth = depth(1:co_depth_in); % cut off depth vector at close off depth

% Age scales for two zones
t_crit = 1/(k_0*A)*log((rho_i - rho_0)/(rho_i - .550));
t_h1 = 1/(k_0*A)*log((rho_i - rho_0)./(rho_i - rho_h1));
t_h2 = 1/(k_1*A^.5)*log((rho_i - .550)./(rho_i - rho_h2)) + t_crit;

% Age profile
age = [t_h1,t_h2];
age = age(1:co_depth_in);   % cut off age vector at close off depth

% Calculate densification rate for two zones
drho_dt_h1 = k_0*A^a*(rho_i-rho_h1);
drho_dt_h2 = k_1*A^b*(rho_i-rho_h2);
drho_dt = [drho_dt_h1,drho_dt_h2];
drho_dt = drho_dt(1:co_depth_in); % cut off vector at close off depth

% Convert density and densification rate to kg/m^3
rho = rho.*1000;
drho_dt = drho_dt*1000;

end


    





