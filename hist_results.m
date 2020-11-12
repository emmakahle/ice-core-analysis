%% Histogram Results

function [] = hist_results(M_mat_full,d_cal_mat,burn_in_in,age,N_d,d_obs,...
    d_unc_up,d_unc_dn,T_data,T_init,T_init_0,A_data,A_init,A_bnd,tfnx_CB,tfnx_init_0)

% Trim out burn-in sampling from results
M_mat_full = M_mat_full(:,burn_in_in:end);
d_cal_mat = d_cal_mat(:,burn_in_in:end);

% set number of bins for parameter histograms
num_bins = 50;

% ************************************
% Define histograms for each parameter
% ************************************

% *** Temperature ***
edges = linspace(min(M_mat_full(1:N_d,:),[],'all'),max(M_mat_full(1:N_d,:),[],'all'),num_bins+1);
% define meshgrid for each variable across domain of age and min/max of
% variable
[X_T,Y_T] = meshgrid(age,edges(1:end-1));  
% create histograms for temperature
for n = 1:N_d
    Z_T(:,n) = histcounts(M_mat_full(n,:),edges);
end

% *** Accumulation ***
edges = linspace(min(M_mat_full(N_d+1:2*N_d,:),[],'all'),max(M_mat_full(N_d+1:2*N_d,:),[],'all'),num_bins+1);
% define meshgrid for each variable across domain of age and min/max of
% variable
[X_A,Y_A] = meshgrid(age,edges(1:end-1));  
% create histograms for accumulation
for n = 1:N_d
    Z_A(:,n) = histcounts(M_mat_full(N_d+n,:),edges);
end

% *** Thinning Function ***
edges = linspace(min(M_mat_full(2*N_d+1:3*N_d,:),[],'all'),max(M_mat_full(2*N_d+1:3*N_d,:),[],'all'),num_bins+1);
% define meshgrid for each variable across domain of age and min/max of
% variable
[X_tfnx,Y_tfnx] = meshgrid(age,edges(1:end-1));   
% create histograms for thinning function
for n = 1:N_d
    Z_tfnx(:,n) = histcounts(M_mat_full(2*N_d+n,:),edges);
end

% *** Dage ***
edges = linspace(min(d_cal_mat(1:N_d,:),[],'all'),max(d_cal_mat(1:N_d,:),[],'all'),num_bins+1);
% define meshgrid for each variable across domain of age and min/max of
% variable
[X_Dage,Y_Dage] = meshgrid(age,edges(1:end-1)); 
% create histograms for Dage
for n = 1:N_d
    Z_Dage(:,n) = histcounts(d_cal_mat(n,:),edges);
end

% *** Diffusion Length ***
edges = linspace(min(d_cal_mat(N_d+1:2*N_d,:),[],'all'),max(d_cal_mat(N_d+1:2*N_d,:),[],'all'),num_bins+1);
% define meshgrid for each variable across domain of age and min/max of
% variable
[X_sigma,Y_sigma] = meshgrid(age,edges(1:end-1));  
% create histograms for diffusion length
for n = 1:N_d
    Z_sigma(:,n) = histcounts(d_cal_mat(N_d+n,:),edges);
end

% *** Layer Thickness ***
edges = linspace(min(d_cal_mat(2*N_d+1:3*N_d,:),[],'all'),max(d_cal_mat(2*N_d+1:3*N_d,:),[],'all'),num_bins+1);
% define meshgrid for each variable across domain of age and min/max of
% variable
[X_lambda,Y_lambda] = meshgrid(age,edges(1:end-1));   
% create histograms for layer thickness
for n = 1:N_d
    Z_lambda(:,n) = histcounts(d_cal_mat(2*N_d+n,:),edges);
end

% ************************************
% ****** Plot Histogram Reults *******
% ************************************

fig('units','inches','width',11,'height',11,'font','Helvetica','fontsize',16);
% DATA SPACE
subplot(3,2,2)
s = pcolor(X_Dage/1000,Y_Dage,Z_Dage);
colormap(flipud(gray))
s.EdgeColor = 'none';
hold on;
plot(age/1000, d_obs(1:N_d),'r','LineWidth',1); % data observation
plot(age/1000, d_unc_up(1:N_d),'--r')
plot(age/1000, d_unc_dn(1:N_d),'--r')
ylabel('\Deltaage [yr]')
title('Data Parameters')
subplot(3,2,4)
s = pcolor(X_sigma/1000,Y_sigma,Z_sigma);
colormap(flipud(gray))
s.EdgeColor = 'none';
hold on;
plot(age/1000, d_obs(N_d+1:2*N_d),'r','LineWidth',1);
plot(age/1000, d_unc_up(N_d+1:2*N_d),'--r')
plot(age/1000, d_unc_dn(N_d+1:2*N_d),'--r')
ylabel('Diffusion Length [m]')
legend('Model Output','Ice Core Data','Uncertainty Bounds','Location','NorthEast')
subplot(3,2,6)
s = pcolor(X_lambda/1000,Y_lambda,Z_lambda);
colormap(flipud(gray))
s.EdgeColor = 'none';
hold on;
plot(age/1000, d_obs(2*N_d+1:end),'r','LineWidth',1);
plot(age/1000, d_unc_up(2*N_d+1:end),'--r')
plot(age/1000, d_unc_dn(2*N_d+1:end),'--r')
ylabel('Layer Thickness [m]')
xlabel('Age [ka]')

% MODEL SPACE
subplot(3,2,1)
s = pcolor(X_T/1000,Y_T,Z_T);
colormap(flipud(gray))
s.EdgeColor = 'none';
hold on;
plot(age/1000, T_data(:,1), 'g','LineWidth',1);
plot(age/1000, T_init, '--m', 'LineWidth', 1);
plot(age/1000, T_init_0(:,2), '--m', 'LineWidth', 0.5);
plot(age/1000, T_init_0(:,3), '--m', 'LineWidth', 0.5);
ylabel('Temperature [C]')
title('Model Parameters')
subplot(3,2,3)
s = pcolor(X_A/1000,Y_A,Z_A);
colormap(flipud(gray))
s.EdgeColor = 'none';
hold on;
plot(age/1000, A_data, 'g','LineWidth',1);
plot(age/1000, A_init, '--m', 'LineWidth', 1);
plot(age/1000, A_init+A_bnd, '--m', 'LineWidth', 0.5);
plot(age/1000, A_init-A_bnd, '--m', 'LineWidth', 0.5);
legend('Model Input','Independent Data Estimate','Initial Guess','Location','NorthEast') 
ylabel('Accumulation [m/yr]')
subplot(3,2,5)
s = pcolor(X_tfnx/1000,Y_tfnx,Z_tfnx);
colormap(flipud(gray))
s.EdgeColor = 'none';
hold on;
plot(age/1000, tfnx_CB, 'g','LineWidth',1);
plot(age/1000, tfnx_init_0, '--m', 'LineWidth', 1);
ylim([0.2 1])
ylabel('Thinning Function')
xlabel('Age [ka]')

end