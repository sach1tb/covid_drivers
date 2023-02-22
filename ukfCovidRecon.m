clearvars
% close all
addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
load fminconOptimisedParameters.mat

nc=13; np=11;
n=nc+np;%number of state
m = 6; %number of measurements
dt = 1;
ddt= dt; % smaller timestep for stable dynamics

% covariance of measurement
%[infectious,death,vax,mask,mobility,Total Population]
% multipliers: infectious is std,
% for population Census.gov says 90% confidence level
% For a 90% confidence level, the critical factor or z-value is 1.645
% MOE = Z*std/sqrt(n)
% https://www2.census.gov/programs-surveys/acs/tech_docs/accuracy/2019_ACS_Accuracy_Document_Worked_Examples.pdf
Rp=[0.00,0.1^2,0.01^2,0.1^2,2^2,(0.001/1.645)^2];

%-- fmnincon optimal parameters (taking the first element of each parameter vector)
% beta0 = modelParams.beta;
% mu0 = modelParams.mu;
% epsilon0 = modelParams.epsilon; %1.0 - optimised
% gamma0 = modelParams.gamma;
% eta_Ih0 = modelParams.eta_Ih;
% eta_Im0 = modelParams.eta_Im;
% eta_Sh0 = modelParams.eta_Sh;
% eta_Sm0 = modelParams.eta_Sm;
% kappa0=modelParams.kappa0;
% sigma0=modelParams.sigma0;
% alpha0 = modelParams.alpha;
%
% phi10 = modelParams.phi1;
% phi20 = modelParams.phi2;
% xi20 = modelParams.xi2;
% xi10 = modelParams.xi1;


beta0 = 0.5;  %transmissibility
mu0 = 2e-4; %Mortality
epsilon0 = 1.0; %1.0 - optimised, rate of progression to infectious class
gamma0 = 0.0110; % Rate of recovery
eta_Ih0 = 0.99; % used in lambda calculation
eta_Im0 = 0.69; % used in lambda calculation
eta_Sh0 = eta_Ih0;  % used in lambda calculation
eta_Sm0 = eta_Im0/2;  % used in lambda calculation
kappa0 = 0.1; % because all kappas are same but may depend on compartment sizes
sigma0 = 0.0001; % because all sigmas are same
alpha0 = 0.001;  % Vaccination rate [Maged2022]

phi10 = 0.01; % S -> S_m
phi20 = 0.01;  % S_m -> S
xi20 = 0.001; % S -> S_h
xi10 = 0.001;  % S_h -> S


% initial estimate

parameterInit = [beta0,xi20,xi10 ...
    alpha0,phi10,phi20,sigma0,kappa0, ...
    mu0, gamma0, epsilon0];


% sigma limits for constrained ukf
sigmaLimitsMax = [1e7*ones(1,nc), ones(1,np)];
sigmaLimitsMin = [0*ones(1,nc), 0*ones(1,np)];

Q=diag([1e-12*ones(1,nc), 0.1*parameterInit] + eps); % covariance of process
f=@(x) seirDynamics(x,eta_Ih0,eta_Im0,eta_Sm0,eta_Sh0,dt);  % nonlinear state equations
h=@(x) seirObservation(x);                               % measurement equation
s=[zeros(1,nc),parameterInit]';  %
s(1) = 12830632*9.7/12.8-24; %9.7e6-24; % initial state susceptible
s(7) = 24; % initial infections
x=s; %initial state          % initial state with noise
P = Q;                        % initial state covraiance
% total dynamic steps


% data/observations

% infectious = csvread('data/infectiousIllinois.csv');
infectious = csvread('data/infectiousIllinois_ci.csv');
infectious_95p_upper_bound=infectious(:,3);
infectious_95p_lower_bound=infectious(:,4);
infectious=infectious(1:1002,2);
% 95% bounds are 3.92 standard deviations apart
infectious_std=(infectious_95p_upper_bound-infectious_95p_lower_bound)/3.92;
infectious(isnan(infectious))=0;
infectious = infectious*(9.7/12.8);
T=size(infectious, 1);

infectious = interp1(1:dt:T, infectious, 1:ddt:T);

days = 1:1:T;

death = csvread('data/deathIllinois.csv');
death(isnan(death))=0;
death = cumsum(death);
death = death*(9.7/12.8);
death = interp1(1:dt:T, death, 1:ddt:T);

vax = csvread('data/vaccinatedIllinois.csv');
vax(isnan(vax))=0;
vax = vax*(9.7/12.8);
vax = interp1(1:dt:T, vax, 1:ddt:T);


mask= csvread('data/maskIllinois.csv');
mask = interp1(1:dt:T, mask, 1:ddt:T);


mobility = csvread('data/mobilityIllinois.csv');
mobility(isnan(mobility))=0;
mobility = interp1(1:dt:T, mobility, 1:ddt:T);


dayStops = [1 332 697 1002];
% data source https://www.statista.com/statistics/815172/chicago-metro-area-population/
% popChicagoMetro = [9684738 9601605 9509934 9433330];

% census.gov
popChicagoMetro = [12830632, 12812508, 12671469, 12582032]*9.7/12.8;


popDays = interp1(dayStops,popChicagoMetro,days);
popDays = interp1(1:dt:T, popDays, 1:ddt:T);


T=size(infectious, 2);

xV = zeros(n,T);          %estmate        % allocate memory
sV = zeros(n,T);          %actual
zV = zeros(m,T);

z = [infectious;death;vax;mask;mobility;popDays]; % measurements

pmat = zeros(n,n,T);
Xprev = zeros((np+nc),2*(np+nc)+1);
sigmaPointAccumulutor = zeros(size(Xprev,1),size(Xprev,2),n);
covarianceMatrix = zeros(n,n,T);
for k=1:T
    zk=z(:,k);                            % save actual state
    zV(:,k)  = zk;                         % save measurment


    % z(1) = xk(7)+xk(8)+xk(9); %Infectious
    % z(2) = xk(11); %death
    % z(3) = xk(13); % Vax
    % z(4) = (xk(2)+xk(5)+xk(8))/totPop; % Mask
    % z(5) = -100*(xk(3)+xk(6)+xk(9))/totPop; % Mobility
    % z(6) = totPop; % population


    R = diag([infectious_std(k)^2, zk(2)*Rp(2), zk(3)*Rp(3),Rp(4), Rp(5), zk(6)*Rp(6)]);
    [x, P, Xprev] = ukfConstrained(f,x,P,h,zk,Q,R,sigmaLimitsMin,sigmaLimitsMax);            % ekf
    sigmaPointAccumulutor(:,:,k) = Xprev;
    covarianceMatrix(:,:,k) = P;
    pmat(:,:,k) = P;
    xV(:,k) = x;                            % save estimate
    fprintf('.');
    if ~mod(k,28)
        fprintf('%d\n', k);
    end
end

% Post calculation filters (smoothing and outlier removal)
smooth_yes=1;
if smooth_yes
    for jj = 1:n
        %xV(jj,:) = filloutliers(xV(jj,:),'nearest',2);
        % smooth smoothens dynamocs as well
        %     xV(jj,:)= smooth(xV(jj,:),28,'lowess');
        xV(jj,:)= sma(xV(jj,:),7);
    end
end

save('ukfOutput.mat','sigmaPointAccumulutor','covarianceMatrix','xV', ...
                    'infectious', 'death', 'vax', 'mask', 'mobility', ...
                    'popDays', 'T', 'nc', 'np');


%% plot results
load('ukfOutput.mat');
startDate = datenum('02-04-2020');
endDate = datenum('11-01-2022');
dateData = linspace(startDate,endDate,T);
numberOfXTicks = 5;

figure(1); gcf; clf;

subplot(3,3,1)

%plot(1:size(xV,2),(xV(1,:)+xV(2,:)+xV(3,:)),'r--','LineWidth',2)
plot(dateData,(xV(1,:)+xV(2,:)+xV(3,:)),'r--','LineWidth',2)
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')

%xlabel("Day");
ylabel("Population")
title('Susceptible')

subplot(3,3,2)

plot(dateData,(xV(4,:)+xV(5,:)+xV(6,:)),'r--','LineWidth',2)
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
%xlabel("Day");
ylabel("Population")
title('Exposed')

subplot(3,3,3)

%waves
waves = [201, 376;
    627, 761;
    775, 958];

plot(dateData,infectious,'k-','LineWidth',2);
grid on;
hold on
plot(dateData,xV(7,:)+xV(8,:)+xV(9,:),'r--','LineWidth',2)
for i = 1:size(waves,1)
    temp = rectangle('Position',[waves(i,1)+dateData(1),0,waves(i,2)-waves(i,1),max(infectious)],'FaceColor',[0.5 0.5 0.5 0.3]);
    temp.EdgeColor = 'none';
end

set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Infectious')


subplot(3,3,4)
plot(dateData,xV(10,:),'r--','LineWidth',2)
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Recovered')

subplot(3,3,5)

plot(dateData,death,'k-','LineWidth',2);
grid on;

hold on
plot(dateData,xV(11,:),'r--','LineWidth',2)
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Deaths (Cumulative)')


subplot(3,3,6)
plot(dateData,vax,'k-','LineWidth',2);
grid on;
hold on
plot(dateData,xV(13,:),'r--','LineWidth',2)
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Vaccinated (Cumulative)')


subplot(3,3,7)
plot(dateData,mask,'k-','LineWidth',2);
grid on;
hold on

plot(dateData,(xV(2,:)+xV(5,:)+xV(8,:))./sum(xV([1:10,12],:),1) ,'r--','LineWidth',2)
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Masked')

subplot(3,3,8)
plot(dateData,mobility,'k-','LineWidth',2);
grid on;
hold on

plot(dateData,-100*(xV(3,:)+xV(6,:)+xV(9,:))./sum(xV([1:10,12],:),1),'r--','LineWidth',2)
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Mobility')

subplot(3,3,9)
plot(dateData,popDays,'k-','LineWidth',2);
grid on;
hold on;
plot(dateData,sum(xV([1:10, 12],:),1),'r--','LineWidth',2)
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Total population')



figure(2); gcf; clf;


legendStr={"\beta", "\xi_1 (mobility)"  ...
    , "\xi_2 (isolation)", "\alpha (vaccination)", "\phi_1 (masking)", ...
    "\phi_2 (unmasking)", "\sigma (loss of imm, post vacc.)",  ...
    "\kappa (loss of imm, post sick.)", "\mu (mortality)", "\gamma (recovery)", "\epsilon (incubation)"};

for ii = 1:np

    subplot(3,4,ii);

    plot(dateData, xV(ii+nc,:),'k-','LineWidth',2) ;
    hold on
    grid on;
    set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
    set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
    datetick('x','mmm, yy', 'keepticks');
    ylabel("Value")
    title(legendStr{ii});
end



figure(3); gcf; clf;

subplot(1,3,1)

plot(dateData,xV(1,:),'LineWidth',2);
hold on
plot(dateData,xV(2,:),'LineWidth',2);
plot(dateData,xV(3,:),'LineWidth',2);
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks');
title('Susceptible')
legend("S","S_m","S_h", 'location','northeastoutside');


subplot(1,3,2)

plot(dateData,xV(4,:),'LineWidth',2);
hold on
plot(dateData,xV(5,:),'LineWidth',2);
plot(dateData,xV(6,:),'LineWidth',2);
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks');
title('Exposed')
legend("E","E_m","E_h", 'location','northeastoutside');


subplot(1,3,3)

plot(dateData,xV(7,:),'LineWidth',2);
hold on
plot(dateData,xV(8,:),'LineWidth',2);
plot(dateData,xV(9,:),'LineWidth',2);
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks');
title('Infectious')
legend("I","I_m","I_h", 'location','northeastoutside');
%%


function z = seirObservation(xk)

totPop = sum(xk([1:10,12]));
z(1) = xk(7)+xk(8)+xk(9); %Infectious
z(2) = xk(11); %death
z(3) = xk(13); % Vax
z(4) = (xk(2)+xk(5)+xk(8))/totPop; % Mask
z(5) = -100*(xk(3)+xk(6)+xk(9))/totPop; % Mobility
z(6) = totPop; % population

end