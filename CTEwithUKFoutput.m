clearvars

addpath('..\..\..\boundedline\boundedline')
addpath('..\..\..\Inpaint_nans')
addpath('..\..\..\cteUpdatedFiles\')

load ukfOutput.mat  %size is 32, 32*2+1
xkk = sigmaPointAccumulutor;
xk = diff(xkk,1,3);
windowSizeDays = 56;
nPoints = 5;
CE_SmtoIS = zeros(nPoints,1001-windowSizeDays);
CE_SmtoIE = zeros(nPoints,1001-windowSizeDays);
CE_SmtoIEm = zeros(nPoints,1001-windowSizeDays);
CE_SmtoIIm = zeros(nPoints,1001-windowSizeDays);
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
    for k = 1:1:1001-windowSizeDays
        %CE_ItoSE(i,k)= cte('hist',tempS,tempE,tempI,1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_SmtoIS(i,k)= cte('hist',I(k:k+windowSizeDays-1),S(k:k+windowSizeDays-1),Sm(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_SmtoIE(i,k) = cte('hist',I(k:k+windowSizeDays-1),E(k:k+windowSizeDays-1),Sm(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_SmtoIEm(i,k) = cte('hist',I(k:k+windowSizeDays-1),Em(k:k+windowSizeDays-1),Sm(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
        CE_SmtoIIm(i,k) = cte('hist',I(k:k+windowSizeDays-1),Im(k:k+windowSizeDays-1),Sm(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
    end

end

figure(1)

clf;
%%
meanCE_SmtoIS = mean(CE_SmtoIS);
stdCE_SmtoIS = std(CE_SmtoIS);


subplot(2,2,1)
h1 = plot(meanCE_SmtoIS,'-r');
hold on

boundedline(1:numel(meanCE_SmtoIS),meanCE_SmtoIS,stdCE_SmtoIS, '-r','alpha','linewidth',1.2);

legend('$CTE_{\dot{S}_m\rightarrow \dot{I}|\dot{S}}$','interpreter','latex');
ylabel('CTE (bits)');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on


%%
meanCE_SmtoIE = mean(CE_SmtoIE);
stdCE_SmtoIE = std(CE_SmtoIE);


subplot(2,2,2)
h1 = plot(meanCE_SmtoIE,'-r');
hold on

boundedline(1:numel(meanCE_SmtoIE),meanCE_SmtoIE,stdCE_SmtoIE, '-r','alpha','linewidth',1.2);

legend('$CTE_{\dot{S}_m\rightarrow \dot{I}|\dot{E}}$','interpreter','latex');
ylabel('CTE (bits)');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on

%%
meanCE_SmtoIEm = mean(CE_SmtoIEm);
stdCE_SmtoIEm = std(CE_SmtoIEm);


subplot(2,2,3)
h1 = plot(meanCE_SmtoIEm,'-r');
hold on

boundedline(1:numel(meanCE_SmtoIEm),meanCE_SmtoIEm,stdCE_SmtoIEm, '-r','alpha','linewidth',1.2);

legend('$CTE_{\dot{S}_m\rightarrow \dot{I}|\dot{E}_m}$','interpreter','latex');
ylabel('CTE (bits)');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on

%%
meanCE_SmtoIIm = mean(CE_SmtoIIm);
stdCE_SmtoIIm = std(CE_SmtoIIm);


subplot(2,2,4)
h1 = plot(meanCE_SmtoIIm,'-r');
hold on

boundedline(1:numel(meanCE_SmtoIIm),meanCE_SmtoIIm,stdCE_SmtoIIm, '-r','alpha','linewidth',1.2);

legend('$CTE_{\dot{S}_m\rightarrow \dot{I}|\dot{I}_m}$' ...
    ,'interpreter','latex');
ylabel('CTE (bits)');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on

figure(2)
clf;
% CE_SmtoIS
% CE_SmtoIE
% CE_SmtoIEm
% CE_SmtoIIm

CE_Sm_S_E_I = meanCE_SmtoIS + meanCE_SmtoIE;
CE_Sm_Em_E_I = meanCE_SmtoIEm + meanCE_SmtoIE;
CE_Sm_Im_I = meanCE_SmtoIEm + meanCE_SmtoIIm;


h1 = plot(CE_Sm_S_E_I,'-r');
hold on
h2 = plot(CE_Sm_Em_E_I,'-g');
h3 = plot(CE_Sm_Im_I,'-b');
%  (S_m, S, E, I)
%  (S_m, E_m, E, I)
%  (S_m, I_m, I)
legend([h1,h2,h3],'$CTE_{\dot{S}_m \rightarrow \dot{I}$', ...
    '$CTE_{\dot{S}_m \rightarrow \dot{I}$ ',...
    '$CTE_{\dot{S}_m \rightarrow \dot{I}$ ',...
    'interpreter','latex');
ylabel('CTE (bits)');
xlabel('Window #')


figure(3)
clf;
plot(squeeze(xkk(7,1,:))');
hold on
plot(squeeze(xkk(8,1,:))');
plot(squeeze(xkk(9,1,:))');
legend('I','Im','Ih');





