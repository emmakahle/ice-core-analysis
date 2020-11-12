%% Steady-State Ice Core Forward Model
% Forward model of ice core physics. 
% Guess inputs are temperature T, accumulation A, and thinning funciton tfnx
% Calculcate resulting delta age Dage, diffusion length sigma, and layer thickness, lambda

function [d_cal] = icecore_forward(M,N_d,age,poly)

% ************ INITIALIZE ************

rho_surf = .35*ones(N_d,1);  % specify surface density

% Specify variables from model M
T = M(1:N_d);
A = M(N_d+1:2*N_d);
if poly == 1
tfnx = polyval(M(2*N_d+1:end),age);
elseif poly ==0
    tfnx = M(2*N_d+1:end);
end

% Initialize result vectors
Dage = zeros(size(T));
sigma = zeros(size(T));
lambda = zeros(size(T));

% **************** FORWARD MODEL ****************
% For each depth point n, run forward model
for n = 1:length(T)

    % ***** Firn Densification (Herron and Langway 1980) *****
        [depth,rho,age,drho_dt] = HL_model(rho_surf(n),T(n),A(n)*.917); % convert A to water-equivalent values
        Dage_temp = age(end);   % temporary Delta age
        Ddepth_temp = depth(end);   % temporary Delta depth

    % ***** Calculate Dage *****
        % Calculate ice age at lock-in
        lock_in_rho = rho(end) - 10; % lock-in density is 10kg/m3 less than close-off density
        ice_age = age(find(rho >= lock_in_rho, 1, 'first')); % ice age at lock-in density
        % Calculate gas age at lock-in
        D_co2 = 1/0.7 * 1.81 * 10^(-2)*(273.15+T(n))^1.81; % constant from Buizert et al. 2013
        DCH_temp = depth(find(rho >= lock_in_rho, 1, 'first')) - 3; % diffusive column height = lock in depth minus a 3m convective zone
        gas_age = 1/1.367*(0.934 * DCH_temp^2/D_co2 + 4.05); % gas age at lock-in
        Dage_temp = ice_age - gas_age;
        
    % ***** Calculate Diffusion Length *****
        [sigma2_18] = J_diff(rho, drho_dt, T(n), rho_surf(n), 0.7); % Johnsen et al. 2000 diffusion formulation
        % Correct for firn densification below the depth where diffusion
        % stops
        rho_co_temp = (1/(917) + 6.95*10^(-7)*(273.15-T(n)) - 4.3*10^(-5))^(-1)/1000; % Martinerie et al 1992, 1994
        sigma2_18 = (rho_co_temp*1000/917)^2 * sigma2_18; 
        sigma_firn = sqrt(sigma2_18);

        % Adjust for thinning of ice
        sigma_temp = sigma_firn*tfnx(n);

    % ***** Calculate layer thickness ******
        lambda_temp = A(n).*tfnx(n);
    
    % Assign temporary values to full vectors
    Dage(n) = Dage_temp;
    sigma(n) = sigma_temp;
    lambda(n) = lambda_temp;

end

% Forward model output
d_cal = [Dage; sigma; lambda];

end