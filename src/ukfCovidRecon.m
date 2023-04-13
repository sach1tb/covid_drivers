clearvars
addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
load fminconOptimisedParameters.mat


% some important parameters
nc=13; np=11;
n=nc+np;%number of state
m = 6; %number of measurements
dt = 1;
ddt= dt; % smaller timestep for stable dynamics


% data/observations
% infectious = csvread('../data/infectiousIllinois.csv');
infectious = csvread('../data/infectiousIllinois_ci.csv');
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

death = csvread('../data/deathIllinois.csv');
death(isnan(death))=0;
death = cumsum(death);
death = death*(9.7/12.8);
death = interp1(1:dt:T, death, 1:ddt:T);

vax = csvread('../data/vaccinatedIllinois.csv');
vax(isnan(vax))=0;
% find when vaccination started v
vax_start_day=find(vax==0);
vax_start_day=vax_start_day(end);
vax = vax*(9.7/12.8);
vax = interp1(1:dt:T, vax, 1:ddt:T);



% data source https://www.statista.com/statistics/815172/chicago-metro-area-population/
% popChicagoMetro = [9684738 9601605 9509934 9433330];

% census.gov
popChicagoMetro = [12830632, 12812508, 12671469, 12582032]*9.7/12.8;
dayStops = [1 332 697 1002];

popDays = interp1(dayStops,popChicagoMetro,days);
popDays = interp1(1:dt:T, popDays, 1:ddt:T);

mask= csvread('../data/maskIllinois.csv');
mask = interp1(1:dt:T, mask, 1:ddt:T);

mobility = csvread('../data/mobilityIllinois.csv');
mobility(isnan(mobility))=0;
mobility = interp1(1:dt:T, mobility, 1:ddt:T);


% "Between 2019 and 2021, the number of people primarily working from home
% tripled from 5.7% (roughly 9 million people) to 17.9% (27.6 million people), 
% according to new 2021 American Community Survey (ACS) 1-year estimates 
% released today by the U.S. Census Bureau."
% https://www.census.gov/newsroom/press-releases/2022/people-working-from-home.html
% 
% This means that we should set 5.7% of population to work from home as
% default. A positive mobility, which we see initially, could then mean that
% people who worked mainly from home came out or simply that people
% started buying things in anticipation of lockdowns. 
% 
% look at the update in the measurement model where we use
normal_wfh=popDays(1)*5.7/100;

% initialize the filter

beta0 = 0.312;  % transmissibility [Maged2022, table 2]
mu0 = 0.035; % Mortality [Maged2022, table 2]
epsilon0 = 1.0/4.5; % [Maged2022, pp 5]
gamma0 = 0.0602; % Rate of recovery, [Maged2022, pp 5]
kappa0 = 1/240; % loss of immunity, initial should be low
sigma0 = 1/240; % loss of immunity, initial should be low
alpha0 = 0.00001;  % Vaccination rate [Maged2022 has 0.001 but that doesn't make sense!]

% next 4 are simply small to begin with represent 10% rate of movement on
% either side
phi10 = 0.1; % S -> S_m % maksing
phi20 = 0.1;  % S_m -> S % unmasking
xi10 = 0.1; % S -> S_h % isolation
xi20 = 0.1;  % S_h -> S % mobility

% constant
eta_Ih0 = 0.99; % decrease in infectivity because of isolation
eta_Im0 = 0.69; % decrease in infectivity because of masks [Brooks and Butler 2021]
eta_Sh0 = eta_Ih0;  % decrease in susceptibility because of isolation
eta_Sm0 = eta_Im0/2;  % decrease in susceptibility because of masks

% initial estimate

parameterInit = [beta0,xi10,xi20 ...
                alpha0,phi10,phi20, ...
                sigma0,kappa0, mu0,...
                gamma0, epsilon0];


s=[zeros(1,nc),parameterInit]';  %
% to initialize we set the initial value of mask and mobility so that it
% matches the data
s(1) = popDays(1)-24-mask(1)*popDays(1)-(5.7+mobility(1))/100*popDays(1); 
s(2)= mask(1)*popDays(1);
s(3)= (5.7+mobility(1))/100*popDays(1);
s(7) = 24; % initial infections
                       % initial state covraiance


% measurement model
h=@(x) seirObservation(x,normal_wfh); 
z = [infectious;death;vax;mask;mobility;popDays]; % measurements
% covariance of measurement
%[infectious,death,vax,mask,mobility,Total Population]
% multipliers: 
% infectious is std based on upper and lower bounds
% death is 10% because it is large
% vaccination is 1%
% 
% for population Census.gov says 90% confidence level
% For a 90% confidence level, the critical factor or z-value is 1.645
% MOE = Z*std/sqrt(n)
% https://www2.census.gov/programs-surveys/acs/tech_docs/accuracy/2019_ACS_Accuracy_Document_Worked_Examples.pdf
Rp=[0.00,0.1^2,0.01^2,0.1^2,2^2,(0.001/1.645)^2];%

% process model
f=@(x) seirDynamics(x,eta_Ih0,eta_Im0,eta_Sm0,eta_Sh0,dt);  % nonlinear state equations

% process noise 
Q=diag([1e-6*ones(1,nc), 0.1^2*beta0, 0.1^2*xi10, 0.1^2*xi20 ...
    0.1^2*alpha0,0.1^2*phi10,0.1^2*phi20,0.1^2*sigma0,0.1^2*kappa0, ...
    0.1^2*mu0, 0.1^2*gamma0, 0.1^2*epsilon0]); % covariance of process
                              % measurement equation
% total dynamic steps

% sigma limits for constrained ukf
sigmaLimitsMax = [popDays(1)*ones(1,nc), 1 1 1 1 1 1 1 1 0.2 1 0.25];
sigmaLimitsMin = [0*ones(1,nc), 0.3 0*ones(1,np-1)]+1e-12;

% initialize the arrays, allocate memory
xV = zeros(n,T);          %estmate        
sV = zeros(n,T);          %actual
zV = zeros(m,T);
pmat = zeros(n,n,T);
Xsigma = zeros((np+nc),2*(np+nc)+1);
sigmaPointAccumulutor = zeros(size(Xsigma,1),size(Xsigma,2),n);
covarianceMatrix = zeros(n,n,T);


x=s; %initial state          
P = Q; 

% run the filter
for k=1:T
    zk=z(:,k);                            % save actual state
    zV(:,k)  = zk;                         % save measurment


    % z(1) = xk(7)+xk(8)+xk(9); %Infectious
    % z(2) = xk(11); %death
    % z(3) = xk(13); % Vax
    % z(4) = (xk(2)+xk(5)+xk(8))/totPop; % Mask
    % z(5) = -100*(xk(3)+xk(6)+xk(9))/totPop; % Mobility
    % z(6) = totPop; % population
    % only let sigma increase if 240 days have passed since vaccination started
    % 240 days from Truszkowska, revax, supplementary, (254-14, or 284-44)
    % Centers for Disease Control and Prevention. Joint statement from HHS public health and medical 
    % experts on COVID-19 booster shots; 2021. Available from: 
    % https://www.cdc.gov/media/releases/2021/s0818-covid-19-booster-shots.html.
    
    % don't let people get vaccinated before it actually starts
    % ukf doesn't like this one!!
%     if k < vax_start_day
%         sigmaLimitsMax(4+nc)=eps;
%     else
%         sigmaLimitsMax(4+nc)=1;
%     end
    
    % don't let people recover from vaccination before 8 months have passed
    % since the first vaccination happened
    if k>240
        if vax(k-240)<1 
            sigmaLimitsMax(7+nc)=0;
        else
            sigmaLimitsMax(7+nc)=1;
        end
    else
        sigmaLimitsMax(7+nc)=0;
    end
    
    % only let kappa grow after 240 days assuming first person infected
    % loses immunity
    if k<240
        sigmaLimitsMax(8+nc)=0;
    else
        sigmaLimitsMax(8+nc)=1;
    end
    
    Rt = diag([infectious_std(k)^2, zk(2)*Rp(2), zk(3)*Rp(3),Rp(4), Rp(5), zk(6)*Rp(6)]);
    [x, P, Xsigma] = ukfConstrained(f,x,P,h,zk,Q,Rt,sigmaLimitsMin,sigmaLimitsMax);            % ekf
    sigmaPointAccumulutor(:,:,k) = Xsigma;
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
    swin=14; % 2 weeks
    % smooth the estimate
    for jj = 1:n
        xV(jj,:)= sma(xV(jj,:),swin);
    end
    
    % smooth all sigma points
    for ii=1:n
        for jj=1:(2*n+1)
            sigmaPointAccumulutor(ii,jj,:)=...
                sma(sigmaPointAccumulutor(ii,jj,:),swin);
        end
    end
end

save('ukfOutput.mat','sigmaPointAccumulutor','covarianceMatrix','xV', ...
                    'infectious', 'death', 'vax', 'mask', 'mobility', ...
                    'popDays', 'T', 'pmat', 'nc', 'np','s');


%% plot results
load('ukfOutput.mat');
startDate = datenum('02-04-2020');
endDate = datenum('11-01-2022');
dateData = linspace(startDate,endDate,T);
numberOfXTicks = 7;

% waves
waves = [201, 376;
    627, 761;
    775, 958];

figure(1); gcf; clf;

subplot(3,3,1)

%plot(1:size(xV,2),(xV(1,:)+xV(2,:)+xV(3,:)),'r--','LineWidth',2)
plot(dateData,(xV(1,:)+xV(2,:)+xV(3,:)),'r--','LineWidth',2)
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')

set(gca, 'ylim', [sum(s(1:7))/3, sum(s(1:7))]);
ylabel("Population")
title('Susceptible')

subplot(3,3,2)

plot(dateData,(xV(4,:)+xV(5,:)+xV(6,:)),'r--','LineWidth',2)
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
set(gca, 'ylim', [0, 1e6]);
ylabel("Population")
title('Exposed')

subplot(3,3,3)

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
set(gca, 'ylim', [0, 15e4]);
ylabel("Population")
title('Infectious')


subplot(3,3,4)
plot(dateData,xV(10,:),'r--','LineWidth',2)
grid on;
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
datetick('x','mmm, yy', 'keepticks')
set(gca, 'ylim', [0, 10e5]);
ylabel("Population")
title('Recovered')

subplot(3,3,5)

plot(dateData,death,'k-','LineWidth',2);
grid on;

hold on
plot(dateData,xV(11,:),'r--','LineWidth',2)
set(gca, 'xtick', ceil(linspace(dateData(1), dateData(T), numberOfXTicks)));
set(gca, 'XLimSpec', 'Tight', 'fontsize', 16);
set(gca, 'ylim', [0, 4e4]);
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
set(gca, 'ylim', [0, 10e6]);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Vaccinated (Cumulative)')


subplot(3,3,7)
plot(dateData,mask,'k-','LineWidth',2);
set(gca, 'ylim', [0, 1]);
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
set(gca, 'ylim', [-100, 0]);
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
set(gca, 'ylim', [9e6, 10e6]);
datetick('x','mmm, yy', 'keepticks')
ylabel("Population")
title('Total population')



figure(2); gcf; clf;


legendStr={"\beta", "\xi_1 (isolation)"  ...
    , "\xi_2 (mobility)", "\alpha (vaccination)", "\phi_1 (masking)", ...
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
set(gca, 'ylim', [0, sum(s(1:7))]);
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
set(gca, 'ylim', [0, 8e5]);
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
set(gca, 'ylim', [0, 15e4]);
datetick('x','mmm, yy', 'keepticks');
title('Infectious')
legend("I","I_m","I_h", 'location','northeastoutside');


% plot all sigma points
figure(4); gcf; clf;
n=size(sigmaPointAccumulutor,1);
for ii=1:n
    subplot(4,6,ii);
    plot(squeeze(sigmaPointAccumulutor(ii,:,:))', 'color', ones(1,3)*0.75, ...
            'linewidth', 1);
    hold on;
    plot(squeeze(sigmaPointAccumulutor(ii,1,:))', 'k',  ...
            'linewidth', 2); 
    if ii> nc
        title(legendStr{ii-nc});
    end
end

%%


function z = seirObservation(xk, normal_wfh)

totPop = sum(xk([1:10,12]));
z(1) = xk(7)+xk(8)+xk(9); %Infectious
z(2) = xk(11); %death
z(3) = xk(13); % Vax
z(4) = (xk(2)+xk(5)+xk(8))/totPop; % Mask
z(5) = -100*(xk(3)+xk(6)+xk(9)-normal_wfh)/totPop; % Mobility
z(6) = totPop; % population

end