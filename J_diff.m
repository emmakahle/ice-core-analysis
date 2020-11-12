%% Johnsen Water Isotope Diffusion
% Use Johnsen et al. 2000 formulation of water isotope diffusion

function [sigma2_18] = J_diff(rho, drho_dt, T, rho_0, P)


% ****** Define diffusivity for each isotopic species ******

% Establish values needed for diffusivity calculation
T_K = T +273.15;    % convert to temp in Kelvin
m = 0.018;    % kg/mol; molar mass of water
pz = 3.454*10^12 * exp(-6133./T_K);  % Pa; saturation vapor pressure over ice
R = 8.314;  % J mol-1 kg-1 Gas constant
rho_i = 917;
rho_co = (1/(rho_i) + 6.95*10^(-7)*T_K - 4.3*10^(-5))^(-1); % Martinerie et al 1992, 1994
alpha_18_z = exp(11.839./T_K - 28.224*10^-3); % fractionation factor
Po = 1.0;   % reference pressure in atm

% Set diffusivity in air (units of m^2/s)
Da = .211 *10 ^-4 * (T_K./273.15)^1.94 * (Po/P);
Da_18 = Da ./ 1.0285;   % account for fractionation factor for 18_O (Johnsen et al. 2000)

% Calculate tortuosity
b = 1.3;      % Tortuosity parameter
invtau= 1 - b *(rho./rho_i).^2;

% Diffusivity
D_18 = m.*pz.*invtau.*Da_18.*(1./rho - 1./rho_i) ./ (R.*T_K*alpha_18_z); % Johnsen et al. 2000
D_18 = D_18*3.154*10^7;  % convert units to m^2/year


% ******* Calculate diffusion length *******

% Calculate square of diffusion length as a function of density (Gkinis et al 2014)
drho = [rho(1)-rho_0*1000 diff(rho)];    % find drho to integrate over
sigma2_18 = 1/rho_co^2*sum(2.*rho.^2.*(drho_dt).^-1.*D_18.*drho);

end