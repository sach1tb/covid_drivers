clearvars
% close all
addpath(['..', filesep, 'boundedline', filesep, 'boundedline'])
addpath(['..', filesep, 'boundedline', filesep, 'Inpaint_nans'])
load fminconOptimisedParameters.mat

nc=13; np=11;
n=nc+np;%number of state
m = 6; %number of measurements
dt = 1;
ddt= dt; % smaller timestep for stable dynamics

% covariance of measurement
%[infectious,death,vax,mask,mobility,Total Population]

Rp=[0.005,0.01,0.01,0.1,5,0.05]; 

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
eta_Ih0 = 0.1; % used in lambda calculation
eta_Im0 = 0.1; % used in lambda calculation
eta_Sh0 = 0.1;  % used in lambda calculation
eta_Sm0 = 0.1;  % used in lambda calculation
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

Q=diag([1e-6*ones(1,nc), 0.1*parameterInit] + eps); % covariance of process
f=@(x) seirDynamics(x,eta_Ih0,eta_Im0,eta_Sm0,eta_Sh0,dt);  % nonlinear state equations
h=@(x) seirObservation(x);                               % measurement equation
s=[zeros(1,nc),parameterInit]';  %
s(1) = 9684738-24; %9.7e6-24; % initial state susceptible
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
popChicagoMetro = [9684738 9601605 9509934 9433330]; 

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
    disp(k);
end



save('ukfOutput.mat','sigmaPointAccumulutor','covarianceMatrix','xV');

% Post calculation filters (smoothing and outlier removal)
for jj = 1:n
    %xV(jj,:) = filloutliers(xV(jj,:),'nearest',2); 
    % smooth smoothens dynamocs as well
%     xV(jj,:)= smooth(xV(jj,:),28,'lowess');
end

%% plot results
figure(1); gcf; clf;

subplot(3,3,1)


plot(1:size(xV,2),(xV(1,:)+xV(2,:)+xV(3,:)),'r--','LineWidth',2)

title('Susceptible')

subplot(3,3,2)

plot((xV(4,:)+xV(5,:)+xV(6,:)),'r--','LineWidth',2)
title('Exposed')

subplot(3,3,3)
plot(infectious,'k-','LineWidth',2);
hold on
plot(xV(7,:)+xV(8,:)+xV(9,:),'r--','LineWidth',2)

title('Infectious')


subplot(3,3,4)
plot(xV(10,:),'r--','LineWidth',2)
title('Recovered')

subplot(3,3,5)

plot(death,'-');

hold on
plot(xV(11,:),'r--','LineWidth',2)

title('Deaths')


subplot(3,3,6)
plot(vax,'k-','LineWidth',2);
hold on
plot(xV(13,:),'r--','LineWidth',2)

title('Vaccinated')


subplot(3,3,7)
plot(mask,'k-','LineWidth',2);
hold on

plot((xV(2,:)+xV(5,:)+xV(8,:))./sum(xV([1:10,12],:),1) ,'r--','LineWidth',2)
title('Masked')

subplot(3,3,8)
plot(mobility,'k-','LineWidth',2);
hold on

plot(-100*(xV(3,:)+xV(6,:)+xV(9,:))./sum(xV([1:10,12],:),1),'r--','LineWidth',2)
title('Mobility')

subplot(3,3,9)
plot(popDays,'k-','LineWidth',2);
hold on;
plot(sum(xV([1:10, 12],:),1),'r--','LineWidth',2)
title('Total population')



figure(2); gcf; clf;


legendStr={"\beta", "\xi_2 (isolation)"  ...
    , "\xi_1 (mobility)", "\alpha (Vax)", "\phi_1 (Masking)", ...
    "\phi_2 (unmasking)", "\sigma (Vax lost)",  ...
    "\kappa (Imm lost)", "\mu (mortality)", "\gamma (Rec rate)", "\epsilon (E to I)"};

for ii = 1:np

subplot(3,4,ii);

plot(xV(ii+nc,:)) ;
hold on

title(legendStr{ii});
end



figure(3); gcf; clf;

subplot(1,3,1)

plot(1:size(xV,2),xV(1,:));
hold on
plot(1:size(xV,2),xV(2,:));
plot(1:size(xV,2),xV(3,:));
title('Susceptible')
legend("S","S_m","S_h");


subplot(1,3,2)

plot(1:size(xV,2),xV(4,:));
hold on
plot(1:size(xV,2),xV(5,:));
plot(1:size(xV,2),xV(6,:));
title('Exposed')
legend("E","E_m","E_h");


subplot(1,3,3)

plot(1:size(xV,2),xV(7,:));
hold on
plot(1:size(xV,2),xV(8,:));
plot(1:size(xV,2),xV(9,:));
title('Infectious')
legend("I","I_m","I_h");
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