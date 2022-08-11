clearvars

close all


addpath('..\..\..\boundedline\boundedline')
addpath('..\..\..\Inpaint_nans')
addpath('..\..\cteUpdatedFiles\')

maxT = 1000; %maximum is 1000
infectious = readmatrix('infectious.csv');
infectious(isnan(infectious))=0;

death = readmatrix('death.csv');
death(isnan(death))=0;
death = cumsum(death);

vax = readmatrix('vaccinated.csv');
vax(isnan(vax))=0;

mask= readmatrix("mask.csv");

mobility= readmatrix("mobility.csv");
mobility(isnan(mobility))=0;

tspan = 1:maxT;

windowInWeeks = 8;
windowSizeDays = windowInWeeks*7;
chunkSize = 21;

maskDot = diff(mask);
maskDot = detrend(maskDot,1);

deathDot = diff(death);
deathDot = detrend(deathDot,1);

vaxDot  = diff(vax);
vaxDot = detrend(vaxDot,1);

infectiousDot = diff(infectious);
infectiousDot = detrend(infectiousDot,1);

mobilityDot = diff(mobility);
mobilityDot = detrend(mobilityDot,1);

%% normalize

nmaskDot = normalize(maskDot');
ndeathDot = normalize(deathDot');
nvaxDot  = normalize(vaxDot');
ninfectiousDot = normalize(infectiousDot');
nmobilityDot = normalize(mobilityDot');


TE_MtoI = [];
TE_VaxtoI = [];
TE_MobtoI = [];
nullTE_MtoI = [];
nullTE_VaxtoI = [];
nullTE_MobtoI = [];
for k = 1:numel(tspan)-windowSizeDays
        tempM = nmaskDot(k:k+windowSizeDays-1);
        nullTempM = circshift(tempM,chunkSize);

        tempI= ninfectiousDot(k:k+windowSizeDays-1);
        nullTempI = circshift(tempI,chunkSize);

        tempVax = nvaxDot(k:k+windowSizeDays-1);
        nullTempVax = circshift(tempVax,chunkSize);

        tempMob = nmobilityDot(k:k+windowSizeDays-1);
        nullTempMob = circshift(tempMob,chunkSize);

        TE_MtoI(k)= ete_hist(tempI,tempM,1,ceil(sqrt(windowSizeDays))*2,[-1 1]);
        nullTE_MtoI(k)= ete_hist(tempI,nullTempM,1,ceil(sqrt(windowSizeDays)),[-1 1]);

        TE_VaxtoI(k)= ete_hist(tempI,tempVax,1,ceil(sqrt(windowSizeDays))*2,[-1 1]);
        nullTE_VaxtoI(k)= ete_hist(tempI,nullTempVax,1,ceil(sqrt(windowSizeDays)),[-1 1]);

        TE_MobtoI(k)= ete_hist(tempI,tempMob,1,ceil(sqrt(windowSizeDays))*2,[-1 1]);
        nullTE_MobtoI(k)= ete_hist(tempI,nullTempMob,1,ceil(sqrt(windowSizeDays)),[-1 1]);


end



figure(1)

clf;
%%
meanTE_MtoI = TE_MtoI;

meanNullTE_MtoI = nullTE_MtoI;

subplot(2,2,1)
h1 = plot(meanTE_MtoI,'-r');
hold on
h2 = plot(meanNullTE_MtoI,'-g');

legend([h1,h2],'$TE_{\dot{M}\rightarrow \dot{I}}$','Null $TE_{\dot{M} \rightarrow \dot{I}}$','$\lambda$S','interpreter','latex');
ylabel('TE (bits)');
xlabel('Window #');
title("Mask use");
set(gca, 'fontsize', 20);
ylim([0 0.7])
grid on

%%
meanTE_VaxtoI =TE_VaxtoI;

meanNullTE_VaxtoI = nullTE_VaxtoI;

subplot(2,2,2)
h1 = plot(meanTE_VaxtoI,'-r');
hold on
h2 = plot(meanNullTE_VaxtoI,'-g');

legend([h1,h2],'$TE_{\dot{Vax}\rightarrow \dot{I}}$','Null $TE_{\dot{Vax} \rightarrow \dot{I}}$','$\lambda$S','interpreter','latex');
ylabel('TE (bits)');
xlabel('Window #')
title("Vaccination");
set(gca, 'fontsize', 20);
ylim([0 0.7])
grid on

%%
meanTE_MobtoI = TE_MobtoI;

meanNullTE_MobtoI = nullTE_MobtoI;

subplot(2,2,3)
h1 = plot(meanTE_MobtoI,'-r');
hold on
h2 = plot(meanNullTE_MobtoI,'-g');

legend([h1,h2],'$TE_{\dot{Mob}\rightarrow \dot{I}}$','Null $TE_{\dot{Mob} \rightarrow \dot{I}}$','$\lambda$S','interpreter','latex');
ylabel('TE (bits)');
xlabel('Window #')
set(gca, 'fontsize', 20);
title("Mobility");
ylim([0 0.7])
grid on

%%

% maskDot = diff(mask);
% deathDot = diff(death);
% vaxDot  = diff(vax);
% infectiousDot = diff(infectious);
% mobilityDot = diff(mobility);



figure(2)


subplot(2,2,1)
h1 = plot(mask,'-r');

legend(h1,'Mask use','interpreter','latex');
ylabel("Mask use %");
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

subplot(2,2,2)
h1 = plot(vax,'-r');

legend(h1,'Vaccinated','interpreter','latex');
ylabel('Vaccinated');
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

subplot(2,2,3)
h1 = plot(infectious,'-r');

legend(h1,'Infections','interpreter','latex');
ylabel('Infections');
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

subplot(2,2,4)
h1 = plot(mobility/100,'-r');

legend(h1,'Mobility (Social Distancing)','interpreter','latex');
ylabel('Changing Mobility %');
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

%%
figure(3)

% maskDot = diff(mask);
% maskDot = detrend(maskDot,1);
% 
% deathDot = diff(death);
% deathDot = detrend(deathDot,1);
% 
% vaxDot  = diff(vax);
% vaxDot = detrend(vaxDot,1);
% 
% infectiousDot = diff(infectious);
% infectiousDot = detrend(infectiousDot,1);
% 
% mobilityDot = diff(mobility);
% mobilityDot = detrend(mobilityDot,1);

subplot(2,2,1)
h1 = plot(xcorr(maskDot,maskDot),'-r');

legend(h1,'Mask use','interpreter','latex');
ylabel("Autocorrelation");
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

subplot(2,2,2)
h1 = plot(xcorr(vaxDot,vaxDot),'-r');

legend(h1,'Vaccinated','interpreter','latex');
ylabel('Autocorrelation');
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

subplot(2,2,3)
h1 = plot(xcorr(infectiousDot,infectiousDot),'-r');

legend(h1,'Infections','interpreter','latex');
ylabel('Autocorrelation');
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

subplot(2,2,4)
h1 = plot(xcorr(mobilityDot,mobilityDot),'-r');

legend(h1,'Mobility (Social Distancing)','interpreter','latex');
ylabel('Autocorrelation');
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

figure(4)


subplot(2,2,1)
h1 = plot(maskDot,'-r');

legend(h1,'Mask use','interpreter','latex');
ylabel("Mask use %");
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

subplot(2,2,2)
h1 = plot(vaxDot,'-r');

legend(h1,'Vaccinated','interpreter','latex');
ylabel('Vaccinated');
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

subplot(2,2,3)
h1 = plot(infectiousDot,'-r');

legend(h1,'Infections','interpreter','latex');
ylabel('Infections');
xlabel('Day')
set(gca, 'fontsize', 20);
grid on

subplot(2,2,4)
h1 = plot(mobilityDot/100,'-r');

legend(h1,'Mobility (Social Distancing)','interpreter','latex');
ylabel('Changing Mobility %');
xlabel('Day')
set(gca, 'fontsize', 20);
grid on