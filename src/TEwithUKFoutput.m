clearvars

addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
addpath(['cteUpdatedFiles', filesep])

load ukfOutput.mat  %size is 24, 24*2+1
sigmaPoints = sigmaPointAccumulutor;
datasetLength = size(sigmaPointAccumulutor,3);
P = covarianceMatrix;
meanValues = squeeze(sigmaPoints(:,1,:));
dxk = diff(sigmaPoints,1,3);


% window size for TE
windowSizeDays = 12*7;
% # of samples from the UKF
nSigmaPoints = 10;

CE_Stot_Itot_condE = zeros(nSigmaPoints,datasetLength-windowSizeDays); % (CE_X_Y_Z) conditional TE X->Y conditioned on Z
CE_Stot_Itot_condEm = zeros(nSigmaPoints,datasetLength-windowSizeDays);
CE_Stot_Itot_condEh = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Stot_Itot= zeros(nSigmaPoints,datasetLength-windowSizeDays); % (TE_X_Y) TE X->Y

for i =1:nSigmaPoints

    str = sprintf('Sigma Point # %d',i);
    disp(str);
    %normalizing over the whole dataset length
    S = normalize(squeeze(dxk(1,i,:)))';
    Sm = normalize(squeeze(dxk(2,i,:)))';
    Sh= normalize(squeeze(dxk(3,i,:)))';

    E = normalize(squeeze(dxk(4,i,:)))';
    Em = normalize(squeeze(dxk(5,i,:)))';
    Eh = normalize(squeeze(dxk(6,i,:)))';
    
    I = normalize(squeeze(dxk(7,i,:)))';
    Im = normalize(squeeze(dxk(8,i,:)))';
    Ih = normalize(squeeze(dxk(9,i,:)))';
    
    R = normalize(squeeze(dxk(10,i,:)))';
    D = normalize(squeeze(dxk(11,i,:)))'; 
    U = normalize(squeeze(dxk(12,i,:)))';
    V = normalize(squeeze(dxk(13,i,:)))';
    Stot = normalize(squeeze(dxk(1,i,:))+squeeze(dxk(2,i,:))+squeeze(dxk(3,i,:))');
    Itot = normalize(squeeze(dxk(7,i,:))+squeeze(dxk(8,i,:))+squeeze(dxk(9,i,:))');
    for k = 1:1:datasetLength-windowSizeDays

        TE_Stot_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),Stot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_Stot_Itot_condE(i,k) = cte('hist',Itot(k:k+windowSizeDays-1),E(k:k+windowSizeDays-1),Stot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_Stot_Itot_condEm(i,k) = cte('hist',Itot(k:k+windowSizeDays-1),Em(k:k+windowSizeDays-1),Stot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_Stot_Itot_condEh(i,k) = cte('hist',Itot(k:k+windowSizeDays-1),Eh(k:k+windowSizeDays-1),Stot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);

    end

end
str = sprintf("teData_winSize%d.mat",windowSizeDays);
save(str,"windowSizeDays",'CE_Stot_Itot_condEh', ...
    "CE_Stot_Itot_condE","CE_Stot_Itot_condEm","TE_Stot_Itot");

% stdSigmPoints
stdOfMean = zeros(size(sigmaPoints,1),size(sigmaPoints,3));
for i = 1:size(P,1)
    stdOfMean(i,:) = sqrt(P(i,i,:));
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


figure(2)

clf;

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

figure(4)
clf;
%TE_Stot_Itot
meanTE_Stot_Itot = mean(TE_Stot_Itot);
stdTE_Stot_Itot = std(TE_Stot_Itot);
h1 = plot(meanTE_Stot_Itot,'-r');
hold on
boundedline(1:numel(meanTE_Stot_Itot),meanTE_Stot_Itot,stdTE_Stot_Itot, '-r','alpha','linewidth',1.2);





