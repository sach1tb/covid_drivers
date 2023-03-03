%% 
clearvars
close all
clc

load ukfOutput.mat  %size is 24, 24*2+1
%% sigma points plot
subplot(3,3,1)
plot(squeeze(sigmaPointAccumulutor(18,2:end,:))', 'color', [1 0 0] + 0.5*(1 - [1 0 0]), 'linewidth', .5)
hold on;
plot(squeeze(sigmaPointAccumulutor(18,1,:))', 'color', [1 0 0], 'linewidth', 2)
title('$\phi_1$ (masking)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,2)
plot(squeeze(sigmaPointAccumulutor(19,2:end,:))', 'color', [1 0 0] + 0.5*(1 - [1 0 0]), 'linewidth', .5)
hold on;
plot(squeeze(sigmaPointAccumulutor(19,1,:))', 'color', [1 0 0], 'linewidth', 2)
title('$\phi_2$ (unmasking)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,3)
plot(squeeze(sigmaPointAccumulutor(15,2:end,:))', 'color', [1 0 0] + 0.5*(1 - [1 0 0]), 'linewidth', .5)
hold on;
plot(squeeze(sigmaPointAccumulutor(15,1,:))', 'color', [1 0 0], 'linewidth', 2)
title('$\xi_1$ (mobility)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,4)
plot(squeeze(sigmaPointAccumulutor(16,2:end,:))', 'color', [1 0 0] + 0.5*(1 - [1 0 0]), 'linewidth', .5)
hold on;
plot(squeeze(sigmaPointAccumulutor(16,1,:))', 'color', [1 0 0], 'linewidth', 2)
title('$\xi_2$ (isolation)','Interpreter','latex')
set(gca, 'fontsize', 18);


subplot(3,3,5)
plot(squeeze(sigmaPointAccumulutor(20,2:end,:))', 'color', [1 0 0] + 0.5*(1 - [1 0 0]), 'linewidth', .5)
hold on;
plot(squeeze(sigmaPointAccumulutor(20,1,:))', 'color', [1 0 0], 'linewidth', 2)
title('$\sigma$ (loss of immunity, post vacc.)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,6)
plot(squeeze(sigmaPointAccumulutor(21,2:end,:))', 'color', [1 0 0] + 0.5*(1 - [1 0 0]), 'linewidth', .5)
hold on;
plot(squeeze(sigmaPointAccumulutor(21,1,:))', 'color', [1 0 0], 'linewidth', 2)
title('$\kappa$ (loss of immunity, post sickness)','Interpreter','latex')
set(gca, 'fontsize', 18);


subplot(3,3,7)
plot(squeeze(sigmaPointAccumulutor(17,2:end,:))', 'color', [1 0 0] + 0.5*(1 - [1 0 0]), 'linewidth', .5)
hold on;
plot(squeeze(sigmaPointAccumulutor(17,1,:))', 'color', [1 0 0], 'linewidth', 2)
title('$\alpha$ (vaccination)','Interpreter','latex')
set(gca, 'fontsize', 18);

subplot(3,3,8)
sigmaPointI = squeeze(sigmaPointAccumulutor(7,2:end,:))' + squeeze(sigmaPointAccumulutor(8,2:end,:))' + squeeze(sigmaPointAccumulutor(9,2:end,:))';
meanPointI = squeeze(sigmaPointAccumulutor(7,1,:))' + squeeze(sigmaPointAccumulutor(8,1,:))' + squeeze(sigmaPointAccumulutor(9,1,:))';
plot(sigmaPointI, 'color', [1 0 0] + 0.5*(1 - [1 0 0]), 'linewidth', .5)
hold on;
plot(meanPointI, 'color', [1 0 0], 'linewidth', 2)
title('Total $I$','Interpreter','latex')
set(gca, 'fontsize', 18);
