%% 
clearvars
close all

load allTECal_win84.mat  %size is 24, 24*2+1
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
