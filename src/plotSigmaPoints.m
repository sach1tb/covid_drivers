%% 
clearvars
close all

load allTEcal_win84_detrend35.mat  %size is 24, 24*2+1
%% sigma points plot
subplot(3,3,1)
plot(phi1_tedata(2:end,:)', 'color', [0 0 0] + 0.5*(1 - [0 0 0]), 'linewidth', .5)
hold on;
plot(phi1_tedata(1,:)', 'color', [0 0 0], 'linewidth', 2)
title('$\phi_1$ (masking)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,2)
plot(phi2_tedata(2:end,:)', 'color', [0 0 0] + 0.5*(1 - [0 0 0]), 'linewidth', .5)
hold on;
plot(phi2_tedata(1,:)', 'color', [0 0 0], 'linewidth', 2)
title('$\phi_2$ (unmasking)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,3)
plot(xi1_tedata(2:end,:)', 'color', [0 0 0] + 0.5*(1 - [0 0 0]), 'linewidth', .5)
hold on;
plot(xi1_tedata(1,:)', 'color', [0 0 0], 'linewidth', 2)
title('$\xi_1$ (mobility)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,4)
plot(xi2_tedata(2:end,:)', 'color', [0 0 0] + 0.5*(1 - [0 0 0]), 'linewidth', .5)
hold on;
plot(xi2_tedata(1,:)', 'color', [0 0 0], 'linewidth', 2)
title('$\xi_2$ (isolation)','Interpreter','latex')
set(gca, 'fontsize', 18);


subplot(3,3,5)
plot(sigma_tedata(2:end,:)', 'color', [0 0 0] + 0.5*(1 - [0 0 0]), 'linewidth', .5)
hold on;
plot(sigma_tedata(1,:)', 'color', [0 0 0], 'linewidth', 2)
title('$\sigma$ (loss of immunity, post vacc.)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,6)
plot(kappa_tedata(2:end,:)', 'color', [0 0 0] + 0.5*(1 - [0 0 0]), 'linewidth', .5)
hold on;
plot(kappa_tedata(1,:)', 'color', [0 0 0], 'linewidth', 2)
title('$\kappa$ (loss of immunity, post sickness)','Interpreter','latex')
set(gca, 'fontsize', 18);


subplot(3,3,7)
plot(alpha_tedata(2:end,:)', 'color', [0 0 0] + 0.5*(1 - [0 0 0]), 'linewidth', .5)
hold on;
plot(alpha_tedata(1,:)', 'color', [0 0 0], 'linewidth', 2)
title('$\alpha$ (vaccination)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,8)
plot(Itotdot_tedata(2:end,:)', 'color', [1 0 0] + 0.5*(1 - [1 0 0]), 'linewidth', .5)
hold on;
plot(Itotdot_tedata(1,:)', 'color', [1 0 0], 'linewidth', 2)
title('$\dot{I}$','Interpreter','latex')
set(gca, 'fontsize', 18);


%% plot autocorrelation
xi_autocorr1 = xcorr(mean(xi1_tedata),mean(xi1_tedata));
xi_autocorr2 = xcorr(mean(xi2_tedata),mean(xi2_tedata));
phi_autocorr1 = xcorr(mean(phi1_tedata),mean(phi1_tedata));
phi_autocorr2 = xcorr(mean(phi1_tedata),mean(phi2_tedata));
sigma_autocorr = xcorr(mean(sigma_tedata),mean(sigma_tedata));
alpha_autocorr = xcorr(mean(alpha_tedata),mean(alpha_tedata));
[kappa_autocorr, lags] = xcorr(mean(kappa_tedata),mean(kappa_tedata));

plot(lags,xi_autocorr1,'linewidth',1.2)
hold on
plot(lags,xi_autocorr2,'linewidth',1.2)
plot(lags,phi_autocorr1,'linewidth',1.2)
plot(lags,phi_autocorr2,'linewidth',1.2)
plot(lags,sigma_autocorr,'linewidth',1.2)
plot(lags,alpha_autocorr,'linewidth',1.2)
plot(lags,kappa_autocorr,'linewidth',1.2)
grid on
legend('$\xi_1$' ...
    ,'$\xi_2$' ...
    ,'$\phi_1$' ...
    ,'$\phi_2$' ...
    ,'$\sigma$' ...
    ,'$\alpha$' ...
    ,'$\kappa$','interpreter','latex');
xlabel('lag')
ylabel('Autocorrelation')
set(gca, 'fontsize', 20);
