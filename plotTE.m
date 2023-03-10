function plotTE()
%rise
rise = [201, 279;
    627, 694;
    775, 891];

%fall
fall = [280, 376;
    694, 761;
    892, 958 ];


plotTE1(rise,fall);

%% 
startDate = datenum('02-04-2020');
endDate = datenum('11-01-2022');
dateData = linspace(startDate,endDate,1002);
numberOfXTicks = 3;

figure(1); gcf;

set(gca, 'xlim', [dateData(1)+fall(1,1), dateData(1)+fall(1,2)]);
% set(gca, 'xtick', ceil(linspace(fall(1,1), fall(1,2), numberOfXTicks)));
% datetick('x','mmm, yy', 'keepticks')
title('Wave 1 (fall)');

print('-dpng', 'plots\h-fall1.png')

set(gca, 'xlim', [dateData(1)+fall(2,1), dateData(1)+fall(2,2)]);
title('Wave 2 (fall)');
print('-dpng', 'plots\h-fall2.png')


set(gca, 'xlim', [dateData(1)+fall(3,1), dateData(1)+fall(3,2)]);
title('Wave 3 (fall)');
print('-dpng', 'plots\h-fall3.png')


figure(2); gcf;
set(gca, 'xlim', [dateData(1)+rise(1,1), dateData(1)+rise(1,2)]);
title('Wave 1 (rise)');
print('-dpng', 'plots\h-rise1.png')

set(gca, 'xlim', [dateData(1)+rise(2,1), dateData(1)+rise(2,2)]);
title('Wave 2 (rise)');
print('-dpng', 'plots\h-rise2.png')

set(gca, 'xlim', [dateData(1)+rise(3,1), dateData(1)+rise(3,2)]);
title('Wave 3 (rise)');
print('-dpng', 'plots\h-rise3.png')



figure(3); gcf;

set(gca, 'xlim', [dateData(1)+rise(1,1), dateData(1)+rise(1,2)]);
print('-dpng', 'plots\a-rise1.png')

set(gca, 'xlim', [dateData(1)+rise(2,1), dateData(1)+rise(2,2)]);
print('-dpng', 'plots\a-rise2.png')

set(gca, 'xlim', [rise(3,1), rise(3,2)]);
print('-dpng', 'plots\a-rise3.png')

set(gca, 'xlim', [dateData(1)+fall(1,1), dateData(1)+fall(1,2)]);
print('-dpng', 'plots\a-fall1.png')

set(gca, 'xlim', [dateData(1)+fall(2,1), dateData(1)+fall(2,2)]);
print('-dpng', 'plots\a-fall2.png')

set(gca, 'xlim', [dateData(1)+fall(3,1), dateData(1)+fall(3,2)]);
print('-dpng', 'plots\a-fall3.png')




end






function plotTE1(rise,fall)
close all
load('allTECal_win84.mat');
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

%% for rises
figure(1)
clf;




TErange=(dateData(1):(dateData(end)-windowSizeDays))+windowSizeDays/2;


meanNetTE_phi1_Itot = mean(NetTE_phi1_Itotdot);
stdNetTE_phi1_Itot = std(NetTE_phi1_Itotdot);
h1 = plot(TErange,meanNetTE_phi1_Itot,'-r');
hold on
boundedline(TErange,meanNetTE_phi1_Itot,stdNetTE_phi1_Itot, '-r','alpha','linewidth',1.2);

meanNetTE_xi2_Itot = mean(NetTE_xi2_Itotdot);
stdNetTE_xi2_Itot = std(NetTE_xi2_Itotdot);
h2 = plot(TErange,meanNetTE_xi2_Itot,'-g');
hold on
boundedline(TErange,meanNetTE_xi2_Itot,stdNetTE_xi2_Itot, '-g','alpha','linewidth',1.2);
yyaxis right
h3 = plot(dateData, infectious,'k--','LineWidth',2);

yyaxis left
meanNetTE_alpha_Itot = mean(NetTE_alpha_Itotdot);
stdNetTE_alpha_Itot = std(NetTE_alpha_Itotdot);
h4 = plot(TErange,meanNetTE_alpha_Itot,'-b');
hold on
boundedline(TErange,meanNetTE_alpha_Itot,stdNetTE_alpha_Itot, '-b','alpha','linewidth',1.2);
% for i = 1:size(fall,1)
%     temp = rectangle('Position',[fall(i,1),-0.6,fall(i,2)-fall(i,1),1+0.6],'FaceColor',[1.0 0.3 0.3 0.3]);
%     temp.EdgeColor = 'none';
% end

legend([h1,h2, h4, h3],'${\phi_1 \rightarrow \dot{I}}$','${\xi_{2} \rightarrow \dot{I}}$' ...
    ,'${\alpha \rightarrow \dot{I}}$','$\hat{I}$','interpreter','latex', 'location', 'northeastoutside');
ylabel('Net TE (bits)');
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
set(gca, 'fontsize', 20);
ylim([-0.2 0.2]);
grid on


%% for falls
figure(2)

meanNetTE_phi2_Itot = mean(NetTE_phi2_Itotdot);
stdNetTE_phi2_Itot = std(NetTE_phi2_Itotdot);
h5 = plot(TErange,meanNetTE_phi2_Itot,'-r');
hold on
boundedline(TErange,meanNetTE_phi2_Itot,stdNetTE_phi2_Itot, '-r','alpha','linewidth',1.2);

meanNetTE_xi1_Itot = mean(NetTE_xi1_Itotdot);
stdNetTE_xi1_Itot = std(NetTE_xi1_Itotdot);
h6 = plot(TErange,meanNetTE_xi1_Itot,'-g');
hold on
boundedline(TErange,meanNetTE_xi1_Itot,stdNetTE_xi1_Itot, '-g','alpha','linewidth',1.2);
yyaxis right
h7 = plot(dateData, infectious,'k--','LineWidth',2);
yyaxis left
meanNetTE_sigma_Itot = mean(NetTE_sigma_Itotdot);
stdNetTE_sigma_Itot = std(NetTE_sigma_Itotdot);
h8 = plot(TErange,meanNetTE_sigma_Itot,'-b');
hold on
boundedline(TErange,meanNetTE_sigma_Itot,stdNetTE_sigma_Itot, '-b','alpha','linewidth',1.2);

meanNetTE_kappa_Itot = mean(NetTE_kappa_Itotdot);
stdNetTE_kappa_Itot = std(NetTE_kappa_Itotdot);
h9 = plot(TErange,meanNetTE_kappa_Itot,'-c');
hold on
boundedline(TErange,meanNetTE_kappa_Itot,stdNetTE_kappa_Itot, '-c','alpha','linewidth',1.2);

% for i = 1:size(rise,1)
%     temp = rectangle('Position',[rise(i,1),-0.6,rise(i,2)-rise(i,1),1+0.6],'FaceColor',[0.3 0.3 1.0 0.3]);
%     temp.EdgeColor = 'none';
% end


legend([h5,h6,h8,h9,h7],'${\phi_{2} \rightarrow \dot{I}}$','${\xi_{1} \rightarrow \dot{I}}$', ...
    '${\sigma \rightarrow \dot{I}}$',  '${\kappa \rightarrow \dot{I}}$','$\hat{I}$','interpreter','latex', 'location', 'northeastoutside');
ylabel('Net TE (bits)');
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')

ylim([-0.2 0.2]);
grid on
set(gca, 'fontsize', 20);
%NetTE_sigma_Itot(i,k)
figure(3)

h1 = plot(TErange,meanNetTE_phi1_Itot,'linewidth',1.2);
c = get(h1,'Color');
boundedline(TErange,meanNetTE_phi1_Itot,stdNetTE_phi1_Itot,'Color', c,'alpha','linewidth',1.2);
hold on

h2 = plot(TErange,meanNetTE_xi2_Itot,'linewidth',1.2);
c = get(h2,'Color');
boundedline(TErange,meanNetTE_xi2_Itot,stdNetTE_xi2_Itot,'Color', c,'alpha','linewidth',1.2);

h3 = plot(TErange,meanNetTE_alpha_Itot,'linewidth',1.2);
c = get(h3,'Color');
boundedline(TErange,meanNetTE_alpha_Itot,stdNetTE_alpha_Itot,'Color', c,'alpha','linewidth',1.2);


h4 = plot(TErange,meanNetTE_phi2_Itot,'linewidth',1.2);
c = get(h4,'Color');
boundedline(TErange,meanNetTE_phi2_Itot,stdNetTE_phi2_Itot,'Color', c,'alpha','linewidth',1.2);

h5 = plot(TErange,meanNetTE_xi1_Itot,'linewidth',1.2);
c = get(h5,'Color');
boundedline(TErange,meanNetTE_xi1_Itot,stdNetTE_xi1_Itot,'Color', c,'alpha','linewidth',1.2);

h6 = plot(TErange,meanNetTE_sigma_Itot,'linewidth',1.2);
c = get(h6,'Color');
boundedline(TErange,meanNetTE_sigma_Itot,stdNetTE_sigma_Itot,'Color', c,'alpha','linewidth',1.2);

meanNetTE_kappa_Itot = mean(NetTE_kappa_Itotdot);
stdNetTE_kappa_Itot = std(NetTE_kappa_Itotdot);

h7 = plot(TErange,meanNetTE_kappa_Itot,'linewidth',1.2);
c = get(h7,'Color');
boundedline(TErange,meanNetTE_kappa_Itot,stdNetTE_kappa_Itot,'Color', c,'alpha','linewidth',1.2);
yyaxis right
h8 = plot(dateData, infectious,'k--','LineWidth',2);
yyaxis left
ylim([-0.2 0.2]);
grid on
legend([h1,h2,h3,h4,h5,h6,h7, h8],'${\phi_{1} \rightarrow \dot{I}}$' ...
    ,'${\xi_{2} \rightarrow \dot{I}}$' ...
    ,'${\alpha \rightarrow \dot{I}}$' ...
    ,'${\phi_{2} \rightarrow \dot{I}}$' ...
    ,'${\xi_{1} \rightarrow \dot{I}}$' ...
    ,'${\sigma \rightarrow \dot{I}}$' ...
    ,'${\kappa \rightarrow \dot{I}}$' ...
    , '$\hat{I}$','interpreter','latex', 'location', 'northeastoutside');
set(gca, 'fontsize', 20);

set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')

% for i = 1:size(fall,1)
%     temp = rectangle('Position',[fall(i,1),-0.6,fall(i,2)-fall(i,1),1+0.6],'FaceColor',[1.0 0.3 0.3 0.3]);
%     temp.EdgeColor = 'none';
% end
% for i = 1:size(rise,1)
%     temp = rectangle('Position',[rise(i,1),-0.6,rise(i,2)-rise(i,1),1+0.6],'FaceColor',[0.3 0.3 1.0 0.3]);
%     temp.EdgeColor = 'none';
% end
figure(4)
% xi2 =
% xi1 =
% phi1
% phi2
% sigma
% alpha
% kappa
xi2 = zeros(nSigmaPoints,1001);
xi1 = zeros(nSigmaPoints,1001);
phi1 = zeros(nSigmaPoints,1001);
phi2 = zeros(nSigmaPoints,1001);
sigma= zeros(nSigmaPoints,1001);
alpha= zeros(nSigmaPoints,1001);
kappa= zeros(nSigmaPoints,1001);
for i  = 1:nSigmaPoints
    xi2(i,:) = detrend(normalize(squeeze(sigmaPoints(16,i,2:end)))');
    xi1(i,:) = detrend(normalize(squeeze(sigmaPoints(15,i,2:end)))');
    phi1(i,:) = detrend(normalize(squeeze(sigmaPoints(18,i,2:end)))');
    phi2(i,:) = detrend(normalize(squeeze(sigmaPoints(19,i,2:end)))');
    sigma(i,:) = detrend(normalize(squeeze(sigmaPoints(20,i,2:end)))');
    alpha(i,:) = detrend(normalize(squeeze(sigmaPoints(17,i,2:end)))');
    kappa(i,:) = detrend(normalize(squeeze(sigmaPoints(21,i,2:end)))');
end

xi_autocorr1 = xcorr(mean(xi1),mean(xi1));
xi_autocorr2 = xcorr(mean(xi2),mean(xi2));
phi_autocorr1 = xcorr(mean(phi1),mean(phi1));
phi_autocorr2 = xcorr(mean(phi1),mean(phi2));
sigma_autocorr = xcorr(mean(sigma),mean(sigma));
alpha_autocorr = xcorr(mean(alpha),mean(alpha));
[kappa_autocorr, lags] = xcorr(mean(kappa),mean(kappa));

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
end