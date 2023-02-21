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
figure(3); gcf;

set(gca, 'xlim', [rise(1,1), rise(1,2)]);
print('-dpng', 'manuscript\plots\rise1.png')



set(gca, 'xlim', [rise(2,1), rise(2,2)]);
print('-dpng', 'manuscript\plots\rise2.png')

set(gca, 'xlim', [rise(3,1), rise(3,2)]);
print('-dpng', 'manuscript\plots\rise3.png')

set(gca, 'xlim', [fall(1,1), fall(1,2)]);
print('-dpng', 'manuscript\plots\fall1.png')

set(gca, 'xlim', [fall(2,1), fall(2,2)]);
print('-dpng', 'manuscript\plots\fall2.png')

set(gca, 'xlim', [fall(3,1), fall(3,2)]);
print('-dpng', 'manuscript\plots\fall3.png')




end






function plotTE1(rise,fall)
load('allTECal_win84.mat');

addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
windowSizeDays = 12*7;



figure(1)
clf;
meanNetTE_phi1_Itot = mean(NetTE_phi1_Itot);
stdNetTE_phi1_Itot = std(NetTE_phi1_Itot);
h1 = plot(windowSizeDays/2+(1:numel(meanNetTE_phi1_Itot)),meanNetTE_phi1_Itot,'-r');
hold on
boundedline(windowSizeDays/2+(1:numel(meanNetTE_phi1_Itot)),meanNetTE_phi1_Itot,stdNetTE_phi1_Itot, '-r','alpha','linewidth',1.2);

meanNetTE_xi2_Itot = mean(NetTE_xi2_Itot);
stdNetTE_xi2_Itot = std(NetTE_xi2_Itot);
h2 = plot(windowSizeDays/2+(1:numel(meanNetTE_xi2_Itot)),meanNetTE_xi2_Itot,'-g');
hold on
boundedline(windowSizeDays/2+(1:numel(meanNetTE_xi2_Itot)),meanNetTE_xi2_Itot,stdNetTE_xi2_Itot, '-g','alpha','linewidth',1.2);
h3 = plot(infectious/max(infectious),'k--','LineWidth',2);


meanNetTE_alpha_Itot = mean(NetTE_alpha_Itot);
stdNetTE_alpha_Itot = std(NetTE_alpha_Itot);
h4 = plot(windowSizeDays/2+(1:numel(meanNetTE_alpha_Itot)),meanNetTE_alpha_Itot,'-b');
hold on
boundedline(windowSizeDays/2+(1:numel(meanNetTE_alpha_Itot)),meanNetTE_alpha_Itot,stdNetTE_alpha_Itot, '-b','alpha','linewidth',1.2);
for i = 1:size(fall,1)
    temp = rectangle('Position',[fall(i,1),-0.6,fall(i,2)-fall(i,1),1+0.6],'FaceColor',[1.0 0.3 0.3 0.3]);
    temp.EdgeColor = 'none';
end

legend([h1,h2, h4, h3],'${\phi_1 \rightarrow \dot{I}}$','${\xi_{2} \rightarrow \dot{I}}$' ...
    ,'${\alpha \rightarrow \dot{I}}$','$\hat{I}$','interpreter','latex', 'location', 'northeastoutside');
ylabel('Net TE (bits)');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-0.6 1]);
grid on

figure(2)

meanNetTE_phi2_Itot = mean(NetTE_phi2_Itot);
stdNetTE_phi2_Itot = std(NetTE_phi2_Itot);
h5 = plot(windowSizeDays/2+(1:numel(meanNetTE_phi2_Itot)),meanNetTE_phi2_Itot,'-r');
hold on
boundedline(windowSizeDays/2+(1:numel(meanNetTE_phi2_Itot)),meanNetTE_phi2_Itot,stdNetTE_phi2_Itot, '-r','alpha','linewidth',1.2);

meanNetTE_xi1_Itot = mean(NetTE_xi1_Itot);
stdNetTE_xi1_Itot = std(NetTE_xi1_Itot);
h6 = plot(windowSizeDays/2+(1:numel(meanNetTE_xi1_Itot)),meanNetTE_xi1_Itot,'-g');
hold on
boundedline(windowSizeDays/2+(1:numel(meanNetTE_xi1_Itot)),meanNetTE_xi1_Itot,stdNetTE_xi1_Itot, '-g','alpha','linewidth',1.2);
h7 = plot(infectious/max(infectious),'k--','LineWidth',2);

meanNetTE_sigma_Itot = mean(NetTE_sigma_Itot);
stdNetTE_sigma_Itot = std(NetTE_sigma_Itot);
h8 = plot(windowSizeDays/2+(1:numel(meanNetTE_sigma_Itot)),meanNetTE_sigma_Itot,'-b');
hold on
boundedline(windowSizeDays/2+(1:numel(meanNetTE_sigma_Itot)),meanNetTE_sigma_Itot,stdNetTE_sigma_Itot, '-b','alpha','linewidth',1.2);

for i = 1:size(rise,1)
    temp = rectangle('Position',[rise(i,1),-0.6,rise(i,2)-rise(i,1),1+0.6],'FaceColor',[0.3 0.3 1.0 0.3]);
    temp.EdgeColor = 'none';
end


legend([h5,h6,h8,h7],'${\phi_{2} \rightarrow \dot{I}}$','${\xi_{1} \rightarrow \dot{I}}$', ...
    '${\sigma \rightarrow \dot{I}}$','$\hat{I}$','interpreter','latex', 'location', 'northeastoutside');
ylabel('Net TE (bits)');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-0.6 1]);
grid on
set(gca, 'fontsize', 20);
%NetTE_sigma_Itot(i,k)
figure(3)

h1 = plot(windowSizeDays/2+(1:numel(meanNetTE_phi1_Itot)),meanNetTE_phi1_Itot,'linewidth',1.2);
c = get(h1,'Color');
boundedline(windowSizeDays/2+(1:numel(meanNetTE_phi1_Itot)),meanNetTE_phi1_Itot,stdNetTE_phi1_Itot,'Color', c,'alpha','linewidth',1.2);
hold on

h2 = plot(windowSizeDays/2+(1:numel(meanNetTE_xi2_Itot)),meanNetTE_xi2_Itot,'linewidth',1.2);
c = get(h2,'Color');
boundedline(windowSizeDays/2+(1:numel(meanNetTE_xi2_Itot)),meanNetTE_xi2_Itot,stdNetTE_xi2_Itot,'Color', c,'alpha','linewidth',1.2);

h3 = plot(windowSizeDays/2+(1:numel(meanNetTE_alpha_Itot)),meanNetTE_alpha_Itot,'linewidth',1.2);
c = get(h3,'Color');
boundedline(windowSizeDays/2+(1:numel(meanNetTE_alpha_Itot)),meanNetTE_alpha_Itot,stdNetTE_alpha_Itot,'Color', c,'alpha','linewidth',1.2);


h4 = plot(windowSizeDays/2+(1:numel(meanNetTE_phi2_Itot)),meanNetTE_phi2_Itot,'linewidth',1.2);
c = get(h4,'Color');
boundedline(windowSizeDays/2+(1:numel(meanNetTE_phi2_Itot)),meanNetTE_phi2_Itot,stdNetTE_phi2_Itot,'Color', c,'alpha','linewidth',1.2);

h5 = plot(windowSizeDays/2+(1:numel(meanNetTE_xi1_Itot)),meanNetTE_xi1_Itot,'linewidth',1.2);
c = get(h5,'Color');
boundedline(windowSizeDays/2+(1:numel(meanNetTE_xi1_Itot)),meanNetTE_xi1_Itot,stdNetTE_xi1_Itot,'Color', c,'alpha','linewidth',1.2);

h6 = plot(windowSizeDays/2+(1:numel(meanNetTE_sigma_Itot)),meanNetTE_sigma_Itot,'linewidth',1.2);
c = get(h6,'Color');
boundedline(windowSizeDays/2+(1:numel(meanNetTE_sigma_Itot)),meanNetTE_sigma_Itot,stdNetTE_sigma_Itot,'Color', c,'alpha','linewidth',1.2);

meanNetTE_kappa_Itot = mean(NetTE_kappa_Itot);
stdNetTE_kappa_Itot = std(NetTE_kappa_Itot);

h7 = plot(windowSizeDays/2+(1:numel(meanNetTE_kappa_Itot)),meanNetTE_kappa_Itot,'linewidth',1.2);
c = get(h7,'Color');
boundedline(windowSizeDays/2+(1:numel(meanNetTE_kappa_Itot)),meanNetTE_kappa_Itot,stdNetTE_kappa_Itot,'Color', c,'alpha','linewidth',1.2);

h8 = plot(infectious/max(infectious),'k--','LineWidth',2);
ylim([-0.6 1]);
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