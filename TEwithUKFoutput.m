clearvars

addpath('..\..\..\boundedline\boundedline')
addpath('..\..\..\Inpaint_nans')
addpath('..\..\..\cteUpdatedFiles\')

load ukfOutput.mat  %size is 32, 32*2+1
xkk = sigmaPointAccumulutor;
P = covarianceMatrix;
meanValues = squeeze(xkk(:,1,:));
xk = diff(xkk,1,3);
for i = 1:size(P,1)
    stdOfMean(i,:) = sqrt(P(i,i,:));
end

windowSizeDays = 56;
nPoints = 10;
TE_SmtoIm = zeros(nPoints,1001-windowSizeDays);
CE_SmtoImS = zeros(nPoints,1001-windowSizeDays);
CE_SmtoImSh = zeros(nPoints,1001-windowSizeDays);

CE_Stot_Itot_condE = zeros(nPoints,1001-windowSizeDays);
CE_Stot_Itot_condEm = zeros(nPoints,1001-windowSizeDays);
CE_Stot_Itot_condEh = zeros(nPoints,1001-windowSizeDays);
% CE_SmtoImEm = zeros(nPoints,1001-windowSizeDays);
% CE_SmtoImIm = zeros(nPoints,1001-windowSizeDays);
for i =1:nPoints

    str = sprintf('Sigma Point # %d',i);
    disp(str);
    % S = xk(1);Sm = xk(2);Sh= xk(3); E = xk(4);Em = xk(5);Eh = xk(6);I = xk(7);
    % Im = xk(8); Ih = xk(9); R = xk(10); D = xk(11); U = xk(12); V = xk(13);


    S = normalize(squeeze(xk(1,i,:)))';Sm = normalize(squeeze(xk(2,i,:)))';
    Sh= normalize(squeeze(xk(3,i,:)))';
    E = normalize(squeeze(xk(4,i,:)))';Em = normalize(squeeze(xk(5,i,:)))';
    Eh = normalize(squeeze(xk(6,i,:)))';I = normalize(squeeze(xk(7,i,:)))';
    Im = normalize(squeeze(xk(8,i,:)))';
    Ih = normalize(squeeze(xk(9,i,:)))'; R = normalize(squeeze(xk(10,i,:)))';
    D = normalize(squeeze(xk(11,i,:)))'; U = normalize(squeeze(xk(12,i,:)))';
    V = normalize(squeeze(xk(13,i,:)))';
    Stot = normalize(squeeze(xk(1,i,:))+squeeze(xk(2,i,:))+squeeze(xk(3,i,:))');
    Itot = normalize(squeeze(xk(7,i,:))+squeeze(xk(8,i,:))+squeeze(xk(9,i,:))');
    for k = 1:1:1001-windowSizeDays
        %CE_ItoSE(i,k)= cte('hist',tempS,tempE,tempI,1,ceil(sqrt(windowSizeDays)),[-1 1]);
        TE_SmtoIm(i,k) = ete_hist(Sm(k:k+windowSizeDays-1),Im(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_SmtoImS(i,k)= cte('hist',Im(k:k+windowSizeDays-1),Sm(k:k+windowSizeDays-1),S(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_SmtoImSh(i,k) = cte('hist',Im(k:k+windowSizeDays-1),Sm(k:k+windowSizeDays-1),Sh(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_Stot_Itot_condE(i,k) = cte('hist',Itot(k:k+windowSizeDays-1),Stot(k:k+windowSizeDays-1),E(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_Stot_Itot_condEm(i,k) = cte('hist',Itot(k:k+windowSizeDays-1),Stot(k:k+windowSizeDays-1),Em(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_Stot_Itot_condEh(i,k) = cte('hist',Itot(k:k+windowSizeDays-1),Stot(k:k+windowSizeDays-1),Eh(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        %CE_SmtoIEm(i,k) = cte('hist',I(k:k+windowSizeDays-1),Em(k:k+windowSizeDays-1),Sm(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        %CE_SmtoIIm(i,k) = cte('hist',I(k:k+windowSizeDays-1),Im(k:k+windowSizeDays-1),Sm(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
    end

end

%% SEI mean plots with std
figure(1)
clf;

subplot(1,3,1)
h1 = plot(meanValues(1,:),'r','LineWidth',1.3); %S
hold on
h2 = plot(meanValues(2,:),'g','LineWidth',1.3); %Sm
h3 = plot(meanValues(3,:),'b','LineWidth',1.3); %Sh
boundedline(1:numel(meanValues(1,:)),meanValues(1,:),stdOfMean(1,:), '-r','alpha','linewidth',1.2);
boundedline(1:numel(meanValues(2,:)),meanValues(2,:),stdOfMean(2,:), '-g','alpha','linewidth',1.2);
boundedline(1:numel(meanValues(3,:)),meanValues(3,:),stdOfMean(3,:), '-b','alpha','linewidth',1.2);
legend([h1 h2 h3],'$S$','$S_m$','$S_h$','Interpreter','latex');



subplot(1,3,2)
h1 = plot(meanValues(4,:),'r','LineWidth',1.3); %S
hold on
h2 = plot(meanValues(5,:),'g','LineWidth',1.3); %Sm
h3 = plot(meanValues(6,:),'b','LineWidth',1.3); %Sh
boundedline(1:numel(meanValues(4,:)),meanValues(4,:),stdOfMean(4,:), '-r','alpha','linewidth',1.2);
boundedline(1:numel(meanValues(5,:)),meanValues(5,:),stdOfMean(5,:), '-g','alpha','linewidth',1.2);
boundedline(1:numel(meanValues(6,:)),meanValues(6,:),stdOfMean(6,:), '-b','alpha','linewidth',1.2);
legend([h1 h2 h3],'$E$','$E_m$','$E_h$','Interpreter','latex');


subplot(1,3,3)
h1 = plot(meanValues(7,:),'r','LineWidth',1.3); %S
hold on
h2 = plot(meanValues(8,:),'g','LineWidth',1.3); %Sm
h3 = plot(meanValues(9,:),'b','LineWidth',1.3); %Sh
boundedline(1:numel(meanValues(7,:)),meanValues(7,:),stdOfMean(7,:), '-r','alpha','linewidth',1.2);
boundedline(1:numel(meanValues(8,:)),meanValues(8,:),stdOfMean(8,:), '-g','alpha','linewidth',1.2);
boundedline(1:numel(meanValues(9,:)),meanValues(9,:),stdOfMean(9,:), '-b','alpha','linewidth',1.2);
legend([h1 h2 h3],'$I$','$I_m$','$I_h$','Interpreter','latex');

%%
figure(2)

clf;

% TE_SmtoIm(i,k)
% CE_SmtoImS(i,k)
% CE_SmtoImSh(i,k)

meanTE_SmtoIm = mean(TE_SmtoIm);
stdTE_SmtoIm = std(TE_SmtoIm);
h1 = plot(meanTE_SmtoIm,'-r');
hold on
boundedline(1:numel(meanTE_SmtoIm),meanTE_SmtoIm,stdTE_SmtoIm, '-r','alpha','linewidth',1.2);

meanCE_SmtoImS = mean(CE_SmtoImS);
stdCE_SmtoImS = std(CE_SmtoImS);
h2 = plot(meanCE_SmtoImS,'-g');
hold on
boundedline(1:numel(meanCE_SmtoImS),meanCE_SmtoImS,stdCE_SmtoImS, '-g','alpha','linewidth',1.2);

meanCE_SmtoImSh = mean(CE_SmtoImSh);
stdCE_SmtoImSh = std(CE_SmtoImSh);
h3 = plot(meanCE_SmtoImSh,'-b');
hold on
boundedline(1:numel(meanCE_SmtoImSh),meanCE_SmtoImSh,stdCE_SmtoImSh, '-b','alpha','linewidth',1.2);


legend([h1,h2,h3],'$TE_{\dot{S}_m\rightarrow \dot{I}_m}$','$CTE_{\dot{S}_m\rightarrow \dot{I}_m | \dot{S}}$' ...
    ,'$CTE_{\dot{S}_m\rightarrow \dot{I}_m | \dot{S}_h}$','interpreter','latex');
ylabel('CTE (bits)');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on



figure(3)

clf;
% CE_Stot_Itot_condEm
% CE_Stot_Itot_condE 
% CE_Stot_Itot_condEh

meanCE_Stot_Itot_condEm = mean(CE_Stot_Itot_condEm);
stdCE_Stot_Itot_condEm = std(CE_Stot_Itot_condEm);
h1 = plot(meanCE_Stot_Itot_condEm,'-r');
hold on
boundedline(1:numel(meanCE_Stot_Itot_condEm),meanCE_Stot_Itot_condEm,stdCE_Stot_Itot_condEm, '-r','alpha','linewidth',1.2);

meanCE_Stot_Itot_condE = mean(CE_Stot_Itot_condE);
stdCE_Stot_Itot_condE = std(CE_Stot_Itot_condE);
h2 = plot(meanCE_Stot_Itot_condE,'-g');
hold on
boundedline(1:numel(meanCE_Stot_Itot_condE),meanCE_Stot_Itot_condE,stdCE_Stot_Itot_condE, '-g','alpha','linewidth',1.2);

meanCE_Stot_Itot_condEh = mean(CE_Stot_Itot_condEh);
stdCE_Stot_Itot_condEh = std(CE_Stot_Itot_condEh);
h3 = plot(meanCE_Stot_Itot_condEh,'-b');
hold on
boundedline(1:numel(meanCE_Stot_Itot_condEh),meanCE_Stot_Itot_condEh,stdCE_Stot_Itot_condEh, '-b','alpha','linewidth',1.2);


legend([h1,h2,h3],'$CTE_{\dot{S}_{tot} \rightarrow \dot{I}_{tot} | \dot{E}_m}$','$CTE_{\dot{S}_{tot} \rightarrow \dot{I}_{tot} | \dot{E}}$' ...
    ,'$CTE_{\dot{S}_{tot} \rightarrow \dot{I}_{tot} | \dot{E}_h}$','interpreter','latex');
ylabel('CTE (bits)');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on


% %%
% meanCE_SmtoIE = mean(CE_SmtoIE);
% stdCE_SmtoIE = std(CE_SmtoIE);
%
%
% subplot(2,2,2)
% h1 = plot(meanCE_SmtoIE,'-r');
% hold on
%
% boundedline(1:numel(meanCE_SmtoIE),meanCE_SmtoIE,stdCE_SmtoIE, '-r','alpha','linewidth',1.2);
%
% legend('$CTE_{\dot{S}_m\rightarrow \dot{I}|\dot{E}}$','interpreter','latex');
% ylabel('CTE (bits)');
% xlabel('Window #')
% set(gca, 'fontsize', 20);
% ylim([-inf inf]);
% grid on
%
% %%
% meanCE_SmtoIEm = mean(CE_SmtoIEm);
% stdCE_SmtoIEm = std(CE_SmtoIEm);
%
%
% subplot(2,2,3)
% h1 = plot(meanCE_SmtoIEm,'-r');
% hold on
%
% boundedline(1:numel(meanCE_SmtoIEm),meanCE_SmtoIEm,stdCE_SmtoIEm, '-r','alpha','linewidth',1.2);
%
% legend('$CTE_{\dot{S}_m\rightarrow \dot{I}|\dot{E}_m}$','interpreter','latex');
% ylabel('CTE (bits)');
% xlabel('Window #')
% set(gca, 'fontsize', 20);
% ylim([-inf inf]);
% grid on
%
% %%
% meanCE_SmtoIIm = mean(CE_SmtoIIm);
% stdCE_SmtoIIm = std(CE_SmtoIIm);
%
%
% subplot(2,2,4)
% h1 = plot(meanCE_SmtoIIm,'-r');
% hold on
%
% boundedline(1:numel(meanCE_SmtoIIm),meanCE_SmtoIIm,stdCE_SmtoIIm, '-r','alpha','linewidth',1.2);
%
% legend('$CTE_{\dot{S}_m\rightarrow \dot{I}|\dot{I}_m}$' ...
%     ,'interpreter','latex');
% ylabel('CTE (bits)');
% xlabel('Window #')
% set(gca, 'fontsize', 20);
% ylim([-inf inf]);
% grid on
%
% figure(3)
% clf;
% % CE_SmtoIS
% % CE_SmtoIE
% % CE_SmtoIEm
% % CE_SmtoIIm
%
% CE_Sm_S_E_I = meanCE_SmtoIS + meanCE_SmtoIE;
% CE_Sm_Em_E_I = meanCE_SmtoIEm + meanCE_SmtoIE;
% CE_Sm_Im_I = meanCE_SmtoIEm + meanCE_SmtoIIm;
%
%
% h1 = plot(CE_Sm_S_E_I,'-r');
% hold on
% h2 = plot(CE_Sm_Em_E_I,'-g');
% h3 = plot(CE_Sm_Im_I,'-b');
% %  (S_m, S, E, I)
% %  (S_m, E_m, E, I)
% %  (S_m, I_m, I)
% legend([h1,h2,h3],'$CTE_{\dot{S}_m \rightarrow \dot{I}$', ...
%     '$CTE_{\dot{S}_m \rightarrow \dot{I}$ ',...
%     '$CTE_{\dot{S}_m \rightarrow \dot{I}$ ',...
%     'interpreter','latex');
% ylabel('CTE (bits)');
% xlabel('Window #')
%
%
% figure(3)
% clf;
% plot(squeeze(xkk(7,1,:))');
% hold on
% plot(squeeze(xkk(8,1,:))');
% plot(squeeze(xkk(9,1,:))');
% legend('I','Im','Ih');





