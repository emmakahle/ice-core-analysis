%% Summary Results

function [mean_M,std_M,mean_D,std_D]=summary_results(M_mat_full,d_cal_mat,...
    age,N_d,d_obs,d_unc_up,d_unc_dn,T_data,T_init,T_init_0,A_data,A_init,...
    A_bnd,tfnx_CB,tfnx_init_0)


% find mean and standard deviation
mean_M = mean(M_mat_full,2);
std_M = 2*std(M_mat_full,[],2);
mean_D = mean(d_cal_mat,2);
std_D = 2*std(d_cal_mat,[],2);

% **** Plot results for Posteriori ****

fig('units','inches','width',11,'height',11,'font','Helvetica','fontsize',16);
% DATA SPACE
subplot(3,2,2)
shadedErrorBar(age/1000,mean_D(1:N_d),std_D(1:N_d),'lineprops','b');
hold on
plot(age/1000, d_obs(1:N_d),'r','LineWidth',1); % data observation
plot(age/1000, d_unc_up(1:N_d),'--r')
plot(age/1000, d_unc_dn(1:N_d),'--r')
ylabel('\Deltaage [yr]')
title('Data Space - Posteriori')
subplot(3,2,4)
shadedErrorBar(age/1000,mean_D(N_d+1:2*N_d),std_D(N_d+1:2*N_d),'lineprops','b');
hold on
plot(age/1000, d_obs(N_d+1:2*N_d),'r','LineWidth',1);
plot(age/1000, d_unc_up(N_d+1:2*N_d),'--r')
plot(age/1000, d_unc_dn(N_d+1:2*N_d),'--r')
legend('Mean Model, 2 Sigma STD','Ice Core Data','Uncertainty Bound','Uncertainty Bound','Location','NorthEast')
ylabel('Diffusion Length [m]')
subplot(3,2,6)
shadedErrorBar(age/1000,mean_D(2*N_d+1:end),std_D(2*N_d+1:end),'lineprops','b');
hold on
plot(age/1000, d_obs(2*N_d+1:end),'r','LineWidth',1);
plot(age/1000, d_unc_up(2*N_d+1:end),'--r')
plot(age/1000, d_unc_dn(2*N_d+1:end),'--r')
ylabel('Layer Thickness [m]')
xlabel('Age [ka]')

% MODEL SPACE
subplot(3,2,1)
shadedErrorBar(age/1000,mean_M(1:N_d),std_M(1:N_d),'lineprops','b');
hold on
plot(age/1000, T_data(:,1), 'g','LineWidth',1);
plot(age/1000, T_init, '--m', 'LineWidth', 1);
plot(age/1000, T_init_0(:,2), '--m', 'LineWidth', 0.5);
plot(age/1000, T_init_0(:,3), '--m', 'LineWidth', 0.5); 
ylabel('Temperature [C]')
ylim([-65 -45])
title('Model Space - Posteriori')
subplot(3,2,3)
shadedErrorBar(age/1000,mean_M(N_d+1:2*N_d),std_M(N_d+1:2*N_d),'lineprops','b');
hold on
plot(age/1000, A_data, 'g','LineWidth',1);
plot(age/1000, A_init, '--m', 'LineWidth', 1);
plot(age/1000, A_init+A_bnd, '--m', 'LineWidth', 0.5);
plot(age/1000, A_init-A_bnd, '--m', 'LineWidth', 0.5);
legend('Model Mean, 2 sigma STD','Data Estimate','Initial Guess','Location','NorthEast')
ylabel('Accumulation [m/yr]')
ylim([0.02 0.1])
subplot(3,2,5)
shadedErrorBar(age/1000,mean_M(2*N_d+1:end),std_M(2*N_d+1:end),'lineprops','b');
hold on
plot(age/1000, tfnx_CB, 'g','LineWidth',1);
plot(age/1000, tfnx_init_0, '--m', 'LineWidth', 1);
ylim([0.2 1])
ylabel('Thinning Function')
xlabel('Age [ka]')

end