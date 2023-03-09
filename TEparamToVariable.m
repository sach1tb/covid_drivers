clearvars
close all
addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
addpath(['cteUpdatedFiles', filesep])

load ukfOutput.mat  %size is 24, 24*2+1
sigmaPoints = sigmaPointAccumulutor;
datasetLength = size(sigmaPointAccumulutor,3);
P = covarianceMatrix;
meanValues = squeeze(sigmaPoints(:,1,:));
dxk = diff(sigmaPoints,1,3);


% windowSizeDays size for TE
windowSizeDays = 12*7;
% # of samples from the UKF
nSigmaPoints = 9;
focussedSigmaPoints = ceil(linspace(1,49,nSigmaPoints));
%CE_Stot_Itot_condE = zeros(nSigmaPoints,datasetLength-windowSizeDays); % (CE_X_Y_Z) conditional TE X->Y conditioned on Z
%CE_Stot_Itot_condEm = zeros(nSigmaPoints,datasetLength-windowSizeDays);
%CE_Stot_Itot_condEh = zeros(nSigmaPoints,datasetLength-windowSizeDays);
%TE_Stot_Itot= zeros(nSigmaPoints,datasetLength-windowSizeDays); % (TE_X_Y) TE X->Y


TE_phi1_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_phi1 = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_phi1_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_xi2_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_xi2 = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_xi2_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);


TE_phi2_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_phi2 = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_phi2_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_xi1_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_xi1 = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_xi1_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_sigma_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_sigma = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_sigma_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_alpha_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_alpha = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_alpha_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_kappa_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_kappa = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_kappa_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);


for i =1:nSigmaPoints
    j = focussedSigmaPoints(i);
    str = sprintf('Sigma Point # %d',i);
    disp(str);
    %normalizing over the whole dataset length
    S = normalize(squeeze(dxk(1,j,:)))';
    Sm = normalize(squeeze(dxk(2,j,:)))';
    Sh= normalize(squeeze(dxk(3,j,:)))';

    E = normalize(squeeze(dxk(4,j,:)))';
    Em = normalize(squeeze(dxk(5,j,:)))';
    Eh = normalize(squeeze(dxk(6,j,:)))';

    I = normalize(squeeze(dxk(7,j,:)))';
    Im = normalize(squeeze(dxk(8,j,:)))';
    Ih = normalize(squeeze(dxk(9,j,:)))';

    R = normalize(squeeze(dxk(10,j,:)))';
    D = normalize(squeeze(dxk(11,j,:)))';
    U = normalize(squeeze(dxk(12,j,:)))';
    V = normalize(squeeze(dxk(13,j,:)))';
    Stot = normalize(squeeze(dxk(1,j,:))+squeeze(dxk(2,j,:))+squeeze(dxk(3,j,:)))';
    Itot = normalize(squeeze(dxk(7,j,:))+squeeze(dxk(8,j,:))+squeeze(dxk(9,j,:)))';

    %     xi2 = detrend(normalize(squeeze(sigmaPoints(16,j,2:end)))');
    %     xi1 = detrend(normalize(squeeze(sigmaPoints(15,j,2:end)))');
    %     phi1 = detrend(normalize(squeeze(sigmaPoints(18,j,2:end)))');
    %     phi2 = detrend(normalize(squeeze(sigmaPoints(19,j,2:end)))');
    %     sigma = detrend(normalize(squeeze(sigmaPoints(20,j,2:end)))');
    %     alpha = detrend(normalize(squeeze(sigmaPoints(17,j,2:end)))');
    %     kappa = detrend(normalize(squeeze(sigmaPoints(21,j,2:end)))');

    xi2 = detrendInAWindow(normalize(squeeze(sigmaPoints(16,j,2:end)))',windowSizeDays);
    xi1 = detrendInAWindow(normalize(squeeze(sigmaPoints(15,j,2:end)))',windowSizeDays);
    phi1 = detrendInAWindow(normalize(squeeze(sigmaPoints(18,j,2:end)))',windowSizeDays);
    phi2 = detrendInAWindow(normalize(squeeze(sigmaPoints(19,j,2:end)))',windowSizeDays);
    sigma = detrendInAWindow(normalize(squeeze(sigmaPoints(20,j,2:end)))',windowSizeDays);
    alpha = detrendInAWindow(normalize(squeeze(sigmaPoints(17,j,2:end)))',windowSizeDays);
    kappa = detrendInAWindow(normalize(squeeze(sigmaPoints(21,j,2:end)))',windowSizeDays);

    

    %raw values, figure 2 add vaccination rate, figure 3 add loss of immunity
    %after vaccination.

    for k = 1:1:datasetLength-windowSizeDays

        %TE_Stot_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),Stot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);

        Hy = ent(Itot(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');

        Hx = ent(phi1(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
        TE_phi1_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),phi1(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        TE_Itot_phi1(i,k) = ete_hist(phi1(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        if Hx*Hy > eps
            NetTE_phi1_Itot(i,k) =  (TE_phi1_Itot(i,k) - TE_Itot_phi1(i,k))/sqrt(Hx*Hy);
        else
            NetTE_phi1_Itot(i,k) =  0;
        end


        Hx = ent(xi2(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
        TE_xi2_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),xi2(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        TE_Itot_xi2(i,k) = ete_hist(xi2(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);

        if Hx*Hy > eps
            NetTE_xi2_Itot(i,k) = (TE_xi2_Itot(i,k)-TE_Itot_xi2(i,k))/sqrt(Hx*Hy) ;
        else
            NetTE_xi2_Itot(i,k) = 0;
        end

        Hx = ent(phi2(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
        TE_phi2_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),phi2(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        TE_Itot_phi2(i,k) = ete_hist(phi2(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);

        if Hx*Hy > eps
            NetTE_phi2_Itot(i,k) = (TE_phi2_Itot(i,k) - TE_Itot_phi2(i,k))/sqrt(Hx*Hy);
        else
            NetTE_phi2_Itot(i,k) = 0;
        end

        Hx = ent(xi1(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
        TE_xi1_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),xi1(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        TE_Itot_xi1(i,k) = ete_hist(xi1(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        if Hx*Hy > eps
            NetTE_xi1_Itot(i,k) = (TE_xi1_Itot(i,k) - TE_Itot_xi1(i,k))/sqrt(Hx*Hy);
        else
            NetTE_xi1_Itot(i,k) = 0;
        end

        Hx = ent(sigma(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
        TE_sigma_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),sigma(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        TE_Itot_sigma(i,k) = ete_hist(sigma(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        if Hx*Hy > eps
            NetTE_sigma_Itot(i,k) = (TE_sigma_Itot(i,k) - TE_Itot_sigma(i,k))/sqrt(Hx*Hy);
        else
            NetTE_sigma_Itot(i,k) = 0;
        end

        Hx = ent(alpha(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
        TE_alpha_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),alpha(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        TE_Itot_alpha(i,k) = ete_hist(alpha(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        if Hx*Hy > eps
            NetTE_alpha_Itot(i,k) = (TE_alpha_Itot(i,k) - TE_Itot_alpha(i,k))/sqrt(Hx*Hy);
        else
            NetTE_alpha_Itot(i,k) = 0;
        end

        Hx = ent(kappa(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
        TE_kappa_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),kappa(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        TE_Itot_kappa(i,k) = ete_hist(kappa(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        if Hx*Hy > eps
            NetTE_kappa_Itot(i,k) = (TE_kappa_Itot(i,k) - TE_Itot_kappa(i,k))/sqrt(Hx*Hy);
        else
            NetTE_kappa_Itot(i,k) = 0;
        end


    end

end
save(sprintf('allTEcal_win%d.mat',windowSizeDays ));
% str = sprintf("teData_winSize%d.mat",windowSizeDays);
% save(str,"windowSizeDays",'Hy', 'TE_phi1_Itot', 'TE_Itot_phi1', 'NetTE_phi1_Itot', ...
%                           'TE_xi2_Itot', 'TE_Itot_xi2', 'NetTE_xi2_Itot', ...
%                           'TE_phi2_Itot', 'TE_Itot_phi2', 'NetTE_phi2_Itot',...
%                           'TE_xi1_Itot', 'TE_Itot_xi1', 'NetTE_xi1_Itot', ...
%                           'TE_sigma_Itot', 'TE_Itot_sigma','NetTE_sigma_Itot', ...
%                           'TE_alpha_Itot', 'TE_Itot_alpha', 'NetTE_alpha_Itot', ...
%                           'TE_kappa_Itot', 'TE_Itot_kappa', 'NetTE_kappa_Itot');

% stdSigmPoints
stdOfMean = zeros(size(sigmaPoints,1),size(sigmaPoints,3));
for i = 1:size(P,1)
    stdOfMean(i,:) = sqrt(P(i,i,:));
end

%% SEI mean plots with std
startDate = datenum('02-04-2020');
endDate = datenum('11-01-2022');
dateData = linspace(startDate,endDate,T);
numberOfXTicks = 5;

figure(1)
clf;

subplot(1,3,1)
h1 = plot(dateData,meanValues(1,:),'r','LineWidth',1.3); %S
hold on
h2 = plot(dateData,meanValues(2,:),'g','LineWidth',1.3); %Sm
h3 = plot(dateData,meanValues(3,:),'b','LineWidth',1.3); %Sh
boundedline(dateData,meanValues(1,:),stdOfMean(1,:), '-r','alpha','linewidth',1.2);
boundedline(dateData,meanValues(2,:),stdOfMean(2,:), '-g','alpha','linewidth',1.2);
boundedline(dateData,meanValues(3,:),stdOfMean(3,:), '-b','alpha','linewidth',1.2);
legend([h1 h2 h3],'$S$','$S_m$','$S_h$','Interpreter','latex');
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Susceptible')

subplot(1,3,2)
h1 = plot(dateData,meanValues(4,:),'r','LineWidth',1.3); %S
hold on
h2 = plot(dateData,meanValues(5,:),'g','LineWidth',1.3); %Sm
h3 = plot(dateData,meanValues(6,:),'b','LineWidth',1.3); %Sh
boundedline(dateData,meanValues(4,:),stdOfMean(4,:), '-r','alpha','linewidth',1.2);
boundedline(dateData,meanValues(5,:),stdOfMean(5,:), '-g','alpha','linewidth',1.2);
boundedline(dateData,meanValues(6,:),stdOfMean(6,:), '-b','alpha','linewidth',1.2);
legend([h1 h2 h3],'$E$','$E_m$','$E_h$','Interpreter','latex');
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Exposed')



subplot(1,3,3)
h1 = plot(dateData,meanValues(7,:),'r','LineWidth',1.3); %S
hold on
h2 = plot(dateData,meanValues(8,:),'g','LineWidth',1.3); %Sm
h3 = plot(dateData,meanValues(9,:),'b','LineWidth',1.3); %Sh
boundedline(dateData,meanValues(7,:),stdOfMean(7,:), '-r','alpha','linewidth',1.2);
boundedline(dateData,meanValues(8,:),stdOfMean(8,:), '-g','alpha','linewidth',1.2);
boundedline(dateData,meanValues(9,:),stdOfMean(9,:), '-b','alpha','linewidth',1.2);
legend([h1 h2 h3],'$I$','$I_m$','$I_h$','Interpreter','latex');
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Infectious')


function y = detrendInAWindow(y, windowSizeDays)
% Detrending along a window whose size is same as the TE window
% Alternatively we can decide the window size based on statioinarity
% criteria
y_smoothed = sma(y,windowSizeDays);
y = y-y_smoothed;

end

