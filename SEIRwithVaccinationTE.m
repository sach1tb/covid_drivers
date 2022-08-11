clearvars

addpath('..\..\..\boundedline\boundedline')
addpath('..\..\..\Inpaint_nans')
addpath('..\..\cteUpdatedFiles\')

frac_masks = 0.5;

S0 = (1-frac_masks)*(1e6-50);
E0 = 0;
I0 = (1-frac_masks)*50;
Sm0 = frac_masks*(1e6-50);
Em0 = 0;
Im0 = frac_masks*50;
eta_s = 0.35;
eta_i = 0.45;
sigma = 0.001;
beta = 0.312;
epsilon = 1/4.5;
gamma = 0.0602;
mu = 3.5/100;
alpha = 0.001;
a = [0.001 0.005];
b = [0.001 0.005];
omega = 0.005;
tau = 100;
tspan = 1:721;
kappa = 0.0*ones(1,numel(tspan));
phi_SMult = 0.01;
phi_EMult = 0;
phi_IMult = 1;
phi_SmMult = 0;
phi_EmMult = 0;
phi_ImMult = 1;


noise = 10;

runs = 10;
S = zeros(numel(tspan),runs);
E = S;
I = S;
Sm = S;
Em = S;
Im = S;
R = S;
D = S;
V = S;
phi_Sm =zeros(numel(tspan)-1,runs);
phi_S = phi_Sm;
lambda_m = phi_Sm;
lambda = phi_Sm;

Sdot = zeros(numel(tspan)-1,runs);
Edot = Sdot;
Idot = Sdot;
Smdot = Sdot;
Emdot = Sdot;
Imdot = Sdot;
Rdot = Sdot;
Ddot = Sdot;
Vdot = Sdot;
iter = 1;

while iter <= runs

    [S1,E1,I1,Sm1,Em1,Im1,R1,D1,V1,phi_S1,phi_Sm1,lambda_m1,lambda1]=seir_sim_wrapper(S0,E0,I0,Sm0, ...
        Em0,Im0,eta_s,eta_i,sigma,beta,epsilon,gamma,mu,alpha,a,b,tau,omega,kappa,noise,tspan, phi_SMult,phi_EMult,phi_IMult,phi_SmMult,phi_EmMult,phi_ImMult);
    N1 = S1+E1+I1+Sm1+Em1+R1+D1+V1;
    if (sum(isnan(S1))==0 && max(N1)<2*N1(1))
        S(:,iter) = S1;
        E(:,iter) = E1;
        I(:,iter) = I1;
        Sm(:,iter) = Sm1;
        Em(:,iter) = Em1;
        Im(:,iter) = Im1;
        R(:,iter) = R1;
        D(:,iter) = D1;
        V(:,iter) = V1;
        phi_Sm(:,iter) = phi_Sm1;
        phi_S(:,iter) = phi_S1;
        lambda_m(:,iter) = lambda_m1;
        lambda(:,iter) = lambda1;

        Sdot(:,iter) = diff(S1); %divide by dt
        Edot(:,iter) = diff(E1);
        Idot(:,iter) = diff(I1);
        Smdot(:,iter) = diff(Sm1);
        Emdot(:,iter) = diff(Em1);
        Imdot(:,iter) = diff(Im1);
        Rdot(:,iter) = diff(R1);
        Ddot(:,iter) = diff(D1);
        Vdot(:,iter) = diff(V1);


        iter = iter+1;
        fprintf('.');

        if ~mod(iter,10)
            fprintf('\n');
        end
    end
end

% [S,E,I,Sm,Em,Im,R,D,V]=seir_sim(S0,E0,I0,Sm0,Em0,Im0,eta_s,eta_i,sigma,beta,epsilon,gamma,mu,alpha,a,b,tau,omega,noise,tspan);
TE_StoE = [];
nullTE_StoE = [];

TE_EtoI = [];
nullTE_EtoI = [];

TE_SmtoS = [];
nullTE_SmtoS = [];

TE_VtoE = [];
nullTE_VtoE = [];
windowInWeeks = 8;
windowSizeDays = windowInWeeks*7;
chunkSize = 21;
for i = 1:runs
    SSdot = normalize(Sdot(:,i)');
    EEdot = normalize(Edot(:,i)');
    IIdot = normalize(Idot(:,i)');
    SSmdot = normalize(Smdot(:,i)');
    EEmdot = normalize(Emdot(:,i)');
    IImdot = normalize(Imdot(:,i)');
    RRdot = normalize(Rdot(:,i)');
    DDdot = normalize(Ddot(:,i)');
    VVdot = normalize(Vdot(:,i)');
    for k = 1:1:numel(tspan)-windowSizeDays
        tempS = SSdot(k:k+windowSizeDays-1);
        nullTempS = circshift(tempS,chunkSize);

        tempSm = SSmdot(k:k+windowSizeDays-1);
        nullTempSm = circshift(tempS,chunkSize);

        tempE = EEdot(k:k+windowSizeDays-1);
        nullTempE = circshift(tempE,chunkSize);

        tempI = IIdot(k:k+windowSizeDays-1);
        nullTempI = circshift(tempI,chunkSize);

        tempV = VVdot(k:k+windowSizeDays-1);
        nullTempV = circshift(tempV,chunkSize);

        tempR = RRdot(k:k+windowSizeDays-1);
        nullTempR = circshift(tempR,chunkSize);

        TE_StoE(i,k)= ete_hist(tempE,tempS,1,ceil(sqrt(windowSizeDays))*2,[-1 1]);
        nullTE_StoE(i,k)= ete_hist(tempE,nullTempS,1,ceil(sqrt(windowSizeDays)),[-1 1]);

        TE_SmtoS(i,k)= ete_hist(tempS,tempSm,1,ceil(sqrt(windowSizeDays)),[-1 1]);
        nullTE_SmtoE(i,k)= ete_hist(tempS,nullTempSm,1,ceil(sqrt(windowSizeDays)),[-1 1]);

        TE_EtoI(i,k)= ete_hist(tempI,tempE,1,ceil(sqrt(windowSizeDays)),[-1 1]);
        nullTE_EtoI(i,k)= ete_hist(tempI,nullTempE,1,ceil(sqrt(windowSizeDays)),[-1 1]);

        TE_VtoE(i,k)= ete_hist(tempE,tempV,1,ceil(sqrt(windowSizeDays)),[-1 1]);
        nullTE_VtoE(i,k)= ete_hist(tempE,nullTempV,1,ceil(sqrt(windowSizeDays)),[-1 1]);

        TE_RtoS(i,k)= ete_hist(tempS,tempR,1,ceil(sqrt(windowSizeDays)),[-1 1]);
        nullTE_RtoS(i,k)= ete_hist(tempS,nullTempR,1,ceil(sqrt(windowSizeDays)),[-1 1]);

    end

end


figure(1)

clf;
%%
meanTE_StoE = mean(TE_StoE);
stdTE_StoE = std(TE_StoE);

meanNullTE_StoE = mean(nullTE_StoE);
stdNullTE_StoE = std(nullTE_StoE);

coupling =lambda.*S(1:end-1,:);
meanCoupling = mean(coupling');
stdCoupling = std(coupling');

subplot(3,2,1)
h1 = plot(meanTE_StoE,'-r');
hold on
h2 = plot(meanNullTE_StoE,'-g');

boundedline(1:numel(meanTE_StoE),meanTE_StoE,stdTE_StoE, '-r','alpha','linewidth',1.2);
boundedline(1:numel(meanNullTE_StoE),meanNullTE_StoE,stdNullTE_StoE, '-g','alpha','linewidth',1.2);
ylabel('TE (bits)');
ylim([0 1]);
yyaxis right
h3 = plot(meanCoupling,'-b');
legend([h1,h2,h3],'$TE_{\dot{S}\rightarrow \dot{E}}$','Null $TE_{\dot{S} \rightarrow \dot{E}}$','$\lambda$S','interpreter','latex');
ylabel('Coupling Strength');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on

%%
meanTE_SmtoS = mean(TE_SmtoS);
stdTE_SmtoS = std(TE_SmtoS);

meanNullTE_SmtoE = mean(nullTE_SmtoE);
stdNullTE_SmtoE = std(nullTE_SmtoE);

coupling =phi_S.*Sm(1:end-1,:);
meanCoupling = mean(coupling');
stdCoupling = std(coupling');

subplot(3,2,2)
h1 = plot(meanTE_SmtoS,'-r');
hold on
h2 = plot(meanNullTE_SmtoE,'-g');

boundedline(1:numel(meanTE_SmtoS),meanTE_SmtoS,stdTE_SmtoS, '-r','alpha','linewidth',1.2);
boundedline(1:numel(meanNullTE_SmtoE),meanNullTE_SmtoE,stdNullTE_SmtoE, '-g','alpha','linewidth',1.2);
ylabel('TE (bits)');
ylim([0 1]);
yyaxis right
h3 = plot(meanCoupling,'-b');
legend([h1,h2,h3],'$TE_{\dot{Sm}\rightarrow \dot{S}}$','Null $TE_{\dot{Sm} \rightarrow \dot{S}}$','$\phi_{S} S_m$','interpreter','latex');
ylabel('Coupling Strength');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on

%%
meanTE_EtoI = mean(TE_EtoI);
stdTE_EtoI = std(TE_EtoI);

meanNullTE_EtoI= mean(nullTE_EtoI);
stdNullTE_EtoI = std(nullTE_EtoI);

coupling =epsilon.*E(1:end-1,:);
meanCoupling = mean(coupling');
stdCoupling = std(coupling');

subplot(3,2,3)
h1 = plot(meanTE_EtoI,'-r');
hold on
h2 = plot(meanNullTE_EtoI,'-g');

boundedline(1:numel(meanTE_EtoI),meanTE_EtoI,stdTE_EtoI, '-r','alpha','linewidth',1.2);
boundedline(1:numel(meanNullTE_EtoI),meanNullTE_EtoI,stdNullTE_EtoI, '-g','alpha','linewidth',1.2);
ylabel('TE (bits)');
ylim([0 1]);
yyaxis right
h3 = plot(meanCoupling,'-b');
legend([h1,h2,h3],'$TE_{\dot{E}\rightarrow \dot{I}}$','Null $TE_{\dot{E} \rightarrow \dot{I}}$','$\epsilon$E','interpreter','latex');
ylabel('Coupling Strength');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on

%%


meanTE_VtoE = mean(TE_VtoE);
stdTE_VtoE = std(TE_VtoE);

meanNullTE_VtoE= mean(nullTE_VtoE);
stdNullTE_VtoE = std(nullTE_VtoE);

coupling =(1-omega)*sigma.*V(1:end-1,:);
meanCoupling = mean(coupling');
stdCoupling = std(coupling');

subplot(3,2,4)
h1 = plot(meanTE_VtoE,'-r');
hold on
h2 = plot(meanNullTE_VtoE,'-g');

boundedline(1:numel(meanTE_VtoE),meanTE_VtoE,stdTE_VtoE, '-r','alpha','linewidth',1.2);
boundedline(1:numel(meanNullTE_VtoE),meanNullTE_VtoE,stdNullTE_VtoE, '-g','alpha','linewidth',1.2);
ylabel('TE (bits)');
ylim([0 1]);
yyaxis right
h3 = plot(meanCoupling,'-b');
legend([h1,h2,h3],'$TE_{\dot{V}\rightarrow \dot{E}}$','Null $TE_{\dot{V} \rightarrow \dot{E}}$','$(1-\omega)\sigma V$','interpreter','latex');
ylabel('Coupling Strength');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on

%%

meanTE_RtoS = mean(TE_RtoS);
stdTE_RtoS = std(TE_RtoS);

meanNullTE_RtoS= mean(nullTE_RtoS);
stdNullTE_RtoS = std(nullTE_RtoS);
temp=repmat(kappa',1,runs);
coupling =(1-omega)*(temp(1:end-1,:).*R(1:end-1,:));
meanCoupling = mean(coupling');
stdCoupling = std(coupling');

subplot(3,2,5)
h1 = plot(meanTE_RtoS,'-r');
hold on
h2 = plot(meanNullTE_RtoS,'-g');

boundedline(1:numel(meanTE_RtoS),meanTE_RtoS,stdTE_RtoS, '-r','alpha','linewidth',1.2);
boundedline(1:numel(meanNullTE_RtoS),meanNullTE_RtoS,stdNullTE_RtoS, '-g','alpha','linewidth',1.2);
ylabel('TE (bits)');
ylim([0 1]);
yyaxis right
h3 = plot(meanCoupling,'-b');
legend([h1,h2,h3],'$TE_{\dot{R}\rightarrow \dot{S}}$','Null $TE_{\dot{R} \rightarrow \dot{S}}$','$(1-\omega_R) \kappa R$','interpreter','latex');
ylabel('Coupling Strength');
xlabel('Window #')
set(gca, 'fontsize', 20);
ylim([-inf inf]);
grid on


figure(2)
clf;
subplot(2,2,1)


plot(S)
hold on
title('Susceptible')


subplot(2,2,2)
plot(E)
hold on
title('Exposed')

subplot(2,2,3)
plot(I)
hold on
title('Infectious')


subplot(2,2,4)
plot(R)
hold on
title('Recovered')
