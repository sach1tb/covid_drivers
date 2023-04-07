function plotTE()
%rise
rise = [201, 279;
        627, 694;
        775, 891];

%fall
fall = [280, 376;
        694, 761;
        892, 958 ];



startDate = datenum('02-04-2020');
endDate = datenum('11-01-2022');
dateData = linspace(startDate,endDate,1002);
load('allTEcal_win84_detrend56.mat');

%% falls
figure(1); gcf; clf;
subplot(3,3,1);
plot_TE_param2Itotdot(NetTE_phi1_Itotdot, '$\phi_1 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+fall(1,1), dateData(1)+fall(1,2)]);

subplot(3,3,2);
plot_TE_param2Itotdot(NetTE_phi1_Itotdot, '$\phi_1 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+fall(2,1), dateData(1)+fall(2,2)]);

subplot(3,3,3);
plot_TE_param2Itotdot(NetTE_phi1_Itotdot, '$\phi_1 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+fall(3,1), dateData(1)+fall(3,2)]);

subplot(3,3,4);
plot_TE_param2Itotdot(NetTE_xi2_Itotdot, '$\xi_2 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+fall(1,1), dateData(1)+fall(1,2)]);

subplot(3,3,5);
plot_TE_param2Itotdot(NetTE_xi2_Itotdot, '$\xi_2 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+fall(2,1), dateData(1)+fall(2,2)]);

subplot(3,3,6);
plot_TE_param2Itotdot(NetTE_xi2_Itotdot, '$\xi_2 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+fall(3,1), dateData(1)+fall(3,2)]);

subplot(3,3,7);
plot_TE_param2Itotdot(NetTE_alpha_Itotdot, '$\alpha_2 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+fall(1,1), dateData(1)+fall(1,2)]);

subplot(3,3,8);
plot_TE_param2Itotdot(NetTE_alpha_Itotdot, '$\alpha_2 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+fall(2,1), dateData(1)+fall(2,2)]);

subplot(3,3,9);
plot_TE_param2Itotdot(NetTE_alpha_Itotdot, '$\alpha_2 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+fall(3,1), dateData(1)+fall(3,2)]);


%% rises
figure(2); gcf; clf;
subplot(4,3,1);
plot_TE_param2Itotdot(NetTE_phi2_Itotdot, '$\phi_2 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(1,1), dateData(1)+rise(1,2)]);

subplot(4,3,2);
plot_TE_param2Itotdot(NetTE_phi2_Itotdot, '$\phi_2 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(2,1), dateData(1)+rise(2,2)]);

subplot(4,3,3);
plot_TE_param2Itotdot(NetTE_phi2_Itotdot, '$\phi_2 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(3,1), dateData(1)+rise(3,2)]);


subplot(4,3,4);
plot_TE_param2Itotdot(NetTE_xi1_Itotdot, '$\xi_1 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(1,1), dateData(1)+rise(1,2)]);

subplot(4,3,5);
plot_TE_param2Itotdot(NetTE_xi1_Itotdot, '$\xi_1 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(2,1), dateData(1)+rise(2,2)]);

subplot(4,3,6);
plot_TE_param2Itotdot(NetTE_xi1_Itotdot, '$\xi_1 \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(3,1), dateData(1)+rise(3,2)]);


subplot(4,3,7);
plot_TE_param2Itotdot(NetTE_sigma_Itotdot, '$\sigma \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(1,1), dateData(1)+rise(1,2)]);

subplot(4,3,8);
plot_TE_param2Itotdot(NetTE_sigma_Itotdot, '$\sigma \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(2,1), dateData(1)+rise(2,2)]);

subplot(4,3,9);
plot_TE_param2Itotdot(NetTE_sigma_Itotdot, '$\sigma \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(3,1), dateData(1)+rise(3,2)]);

subplot(4,3,10);
plot_TE_param2Itotdot(NetTE_kappa_Itotdot, '$\kappa \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(1,1), dateData(1)+rise(1,2)]);

subplot(4,3,11);
plot_TE_param2Itotdot(NetTE_kappa_Itotdot, '$\kappa \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(2,1), dateData(1)+rise(2,2)]);

subplot(4,3,12);
plot_TE_param2Itotdot(NetTE_kappa_Itotdot, '$\kappa \rightarrow \dot{I}$ (bits)',T);
set(gca, 'xlim', [dateData(1)+rise(3,1), dateData(1)+rise(3,2)]);

end






function plot_TE_param2Itotdot(NetTE_parameter, label,T)

infectious = csvread('data/infectiousIllinois_ci.csv');
infectious=infectious(1:1002,2);
infectious(isnan(infectious))=0;
infectious = infectious*(9.7/12.8);
addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
windowSizeDays = 12*7;

% to show dates on x-axis
startDate = datenum('02-04-2020');
endDate = datenum('11-01-2022');
dateData = linspace(startDate,endDate,T);
numberOfXTicks = 30;

gca;


TErange=(dateData(1):(dateData(end)-windowSizeDays))+windowSizeDays/2;

meanNetTE_param_Itotdot= mean(NetTE_parameter);
stdNetTE_param_Itotdot=std(NetTE_parameter);

h2 = plot(TErange,meanNetTE_param_Itotdot,'-b', 'linewidth', 2);
hold on
boundedline(TErange,meanNetTE_param_Itotdot,stdNetTE_param_Itotdot, '-b','alpha','linewidth',1.2);
set(gca, 'ytick', [-0.1 0 0.1 0.2 0.3]);
ylim([-0.1 0.3]);
ylabel(label, 'interpreter', 'latex');
yyaxis right
h1 = plot(dateData, infectious,'r--','LineWidth',2);
ylabel('$\hat{I}$', 'interpreter', 'latex');
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';

set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
set(gca, 'fontsize', 20);

grid on

end

% meanNetTE_phi1_Itot = mean(NetTE_phi1_Itotdot);
% stdNetTE_phi1_Itot = std(NetTE_phi1_Itotdot);
% h2 = plot(TErange,meanNetTE_phi1_Itot,'-r');
% hold on
% boundedline(TErange,meanNetTE_phi1_Itot,stdNetTE_phi1_Itot, '-r','alpha','linewidth',1.2);
% yyaxis right
% h1 = plot(dateData, infectious,'k--','LineWidth',2);
%
% subplot(4,3,2);
% yyaxis left
% meanNetTE_xi2_Itot = mean(NetTE_xi2_Itotdot);
% stdNetTE_xi2_Itot = std(NetTE_xi2_Itotdot);
% h3 = plot(TErange,meanNetTE_xi2_Itot,'-g');
% hold on
% boundedline(TErange,meanNetTE_xi2_Itot,stdNetTE_xi2_Itot, '-g','alpha','linewidth',1.2);
% yyaxis right
% plot(dateData, infectious,'k--','LineWidth',2);
%
% subplot(4,3,3);
% yyaxis left
% meanNetTE_alpha_Itot = mean(NetTE_alpha_Itotdot);
% stdNetTE_alpha_Itot = std(NetTE_alpha_Itotdot);
% h4 = plot(TErange,meanNetTE_alpha_Itot,'-b');
% hold on
% boundedline(TErange,meanNetTE_alpha_Itot,stdNetTE_alpha_Itot, '-b','alpha','linewidth',1.2);
% yyaxis right
% plot(dateData, infectious,'k--','LineWidth',2);
%
%
% % legend([h2,h1, h4, h3],'${\phi_1 \rightarrow \dot{I}}$','${\xi_{2} \rightarrow \dot{I}}$' ...
% %     ,'${\alpha \rightarrow \dot{I}}$','$\hat{I}$','interpreter','latex', 'location', 'northeastoutside');
% ylabel('Net TE (bits)');
% set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
% set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
% datetick('x','mmm, yy', 'keepticks')
% set(gca, 'fontsize', 20);
% ylim([-0.2 0.2]);
% grid on
%
%
% %% for falls
% figure(2)
%
% meanNetTE_phi2_Itot = mean(NetTE_phi2_Itotdot);
% stdNetTE_phi2_Itot = std(NetTE_phi2_Itotdot);
% h5 = plot(TErange,meanNetTE_phi2_Itot,'-r');
% hold on
% boundedline(TErange,meanNetTE_phi2_Itot,stdNetTE_phi2_Itot, '-r','alpha','linewidth',1.2);
%
% meanNetTE_xi1_Itot = mean(NetTE_xi1_Itotdot);
% stdNetTE_xi1_Itot = std(NetTE_xi1_Itotdot);
% h6 = plot(TErange,meanNetTE_xi1_Itot,'-g');
% hold on
% boundedline(TErange,meanNetTE_xi1_Itot,stdNetTE_xi1_Itot, '-g','alpha','linewidth',1.2);
% yyaxis right
% h7 = plot(dateData, infectious,'k--','LineWidth',2);
% yyaxis left
% meanNetTE_sigma_Itot = mean(NetTE_sigma_Itotdot);
% stdNetTE_sigma_Itot = std(NetTE_sigma_Itotdot);
% h8 = plot(TErange,meanNetTE_sigma_Itot,'-b');
% hold on
% boundedline(TErange,meanNetTE_sigma_Itot,stdNetTE_sigma_Itot, '-b','alpha','linewidth',1.2);
%
% meanNetTE_kappa_Itot = mean(NetTE_kappa_Itotdot);
% stdNetTE_kappa_Itot = std(NetTE_kappa_Itotdot);
% h9 = plot(TErange,meanNetTE_kappa_Itot,'-c');
% hold on
% boundedline(TErange,meanNetTE_kappa_Itot,stdNetTE_kappa_Itot, '-c','alpha','linewidth',1.2);
%
% % for i = 1:size(rise,1)
% %     temp = rectangle('Position',[rise(i,1),-0.6,rise(i,2)-rise(i,1),1+0.6],'FaceColor',[0.3 0.3 1.0 0.3]);
% %     temp.EdgeColor = 'none';
% % end
%
%
% legend([h5,h6,h8,h9,h7],'${\phi_{2} \rightarrow \dot{I}}$','${\xi_{1} \rightarrow \dot{I}}$', ...
%     '${\sigma \rightarrow \dot{I}}$',  '${\kappa \rightarrow \dot{I}}$','$\hat{I}$','interpreter','latex', 'location', 'northeastoutside');
% ylabel('Net TE (bits)');
% set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
% set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
% datetick('x','mmm, yy', 'keepticks')
%
% ylim([-0.2 0.2]);
% grid on
% set(gca, 'fontsize', 20);
% %NetTE_sigma_Itot(i,k)
% figure(3)
%
% h1 = plot(TErange,meanNetTE_phi1_Itot,'linewidth',1.2);
% c = get(h1,'Color');
% boundedline(TErange,meanNetTE_phi1_Itot,stdNetTE_phi1_Itot,'Color', c,'alpha','linewidth',1.2);
% hold on
%
% h2 = plot(TErange,meanNetTE_xi2_Itot,'linewidth',1.2);
% c = get(h2,'Color');
% boundedline(TErange,meanNetTE_xi2_Itot,stdNetTE_xi2_Itot,'Color', c,'alpha','linewidth',1.2);
%
% h3 = plot(TErange,meanNetTE_alpha_Itot,'linewidth',1.2);
% c = get(h3,'Color');
% boundedline(TErange,meanNetTE_alpha_Itot,stdNetTE_alpha_Itot,'Color', c,'alpha','linewidth',1.2);
%
%
% h4 = plot(TErange,meanNetTE_phi2_Itot,'linewidth',1.2);
% c = get(h4,'Color');
% boundedline(TErange,meanNetTE_phi2_Itot,stdNetTE_phi2_Itot,'Color', c,'alpha','linewidth',1.2);
%
% h5 = plot(TErange,meanNetTE_xi1_Itot,'linewidth',1.2);
% c = get(h5,'Color');
% boundedline(TErange,meanNetTE_xi1_Itot,stdNetTE_xi1_Itot,'Color', c,'alpha','linewidth',1.2);
%
% h6 = plot(TErange,meanNetTE_sigma_Itot,'linewidth',1.2);
% c = get(h6,'Color');
% boundedline(TErange,meanNetTE_sigma_Itot,stdNetTE_sigma_Itot,'Color', c,'alpha','linewidth',1.2);
%
% meanNetTE_kappa_Itot = mean(NetTE_kappa_Itotdot);
% stdNetTE_kappa_Itot = std(NetTE_kappa_Itotdot);
%
% h7 = plot(TErange,meanNetTE_kappa_Itot,'linewidth',1.2);
% c = get(h7,'Color');
% boundedline(TErange,meanNetTE_kappa_Itot,stdNetTE_kappa_Itot,'Color', c,'alpha','linewidth',1.2);
% yyaxis right
% h8 = plot(dateData, infectious,'k--','LineWidth',2);
% yyaxis left
% ylim([-0.2 0.2]);
% grid on
% legend([h1,h2,h3,h4,h5,h6,h7, h8],'${\phi_{1} \rightarrow \dot{I}}$' ...
%     ,'${\xi_{2} \rightarrow \dot{I}}$' ...
%     ,'${\alpha \rightarrow \dot{I}}$' ...
%     ,'${\phi_{2} \rightarrow \dot{I}}$' ...
%     ,'${\xi_{1} \rightarrow \dot{I}}$' ...
%     ,'${\sigma \rightarrow \dot{I}}$' ...
%     ,'${\kappa \rightarrow \dot{I}}$' ...
%     , '$\hat{I}$','interpreter','latex', 'location', 'northeastoutside');
% set(gca, 'fontsize', 20);
%
% set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
% set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
% datetick('x','mmm, yy', 'keepticks')


% end