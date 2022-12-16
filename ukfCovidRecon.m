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

R=diag([1000,100,500,0.3,10,1000]); % covariance of measurement,%[infectious,death,vax,mask,mobility,Total Population]

%-- fmnincon optimal parameters (taking the first element of each parameter vector)
beta0 = modelParams.betaVector(1)*0.01;
mu0 = modelParams.muVector(1)*0.01;
epsilon0 = modelParams.epsilonVector(1)*0.01; %1.0 - optimised
gamma0 = modelParams.gammaVector(1)*0.01;
eta_Ih0 = modelParams.eta_Ih(1)*0.01;
eta_Im0 = modelParams.eta_Im(1)*0.01;
eta_Sh0 = modelParams.eta_Sh(1)*0.001;
eta_Sm0 = modelParams.eta_Sm(1)*0.01;
kappa_R0 = modelParams.kappa_RVector(1)*0.001;
kappa_Rm0 = modelParams.kappa_RmVector(1)*0.01;
kappa_Rh0 = modelParams.kappa_RhVector(1)*0.01;
kappa0=kappa_R0; % because all kappas are same but may depend on compartment sizes
sigma_S0 = modelParams.sigma_SVector(1)*0.01;
sigma_Sm0 = modelParams.sigma_SmVector(1)*0.01;
sigma_Sh0 = modelParams.sigma_ShVector(1)*0.01;
sigma0=sigma_S0; % because all sigmas are same
alpha0 = 0.0246;

phi_Sm0 = modelParams.phi1Vector(1)*0.01;
phi_S0 = modelParams.phi2Vector(1)*0.01;
xi_Sh0 = modelParams.xi2Vector(1)*0.01;
xi_S0 = modelParams.xi1Vector(1)*0.01;

% initial estimate

parameterInit = [beta0,xi_Sh0,xi_S0 ...
    alpha0,phi_Sm0,phi_S0,sigma0,kappa0, ...
    mu0, gamma0, epsilon0];


% sigma limits for constrained ukf
sigmaLimitsMax = [1e7*ones(1,nc), ones(1,np)];
sigmaLimitsMin = [0*ones(1,nc), 0*ones(1,np)];

Q=diag([1e-6*ones(1,nc), 0.1*parameterInit]); % covariance of process
f=@(x) seirDynamics(x,eta_Ih0,eta_Im0,eta_Sm0,eta_Sh0,dt);  % nonlinear state equations
h=@(x) seirObservation(x);                               % measurement equation
s=[zeros(1,nc),parameterInit]';  %
s(1) = 9684738-24; %9.7e6-24; % initial state susceptible
s(7) = 24; % initial infections
x=s; %initial state          % initial state with noise
P = Q;                        % initial state covraiance
                                    % total dynamic steps

 
% data/observations

infectious = csvread('data/infectiousIllinois.csv');
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
    % xV(jj,:) = filloutliers(xV(jj,:),'nearest',2); 
    % smooth smoothens dynamocs as well
    xV(jj,:)= smooth(xV(jj,:),14,'lowess');
end

%% plot results
figure(1); gcf; clf;

subplot(3,3,1)


plot(1:size(xV,2),(xV(1,:)+xV(2,:)+xV(3,:)),'--')

title('Susceptible')

subplot(3,3,2)

plot((xV(4,:)+xV(5,:)+xV(6,:)),'--')
title('Exposed')

subplot(3,3,3)
plot(infectious,'-','Color','g');
hold on
plot(xV(7,:)+xV(8,:)+xV(9,:),'--','Color','r','LineWidth',1.4)

title('Infectious')


subplot(3,3,4)
plot(xV(10,:),'--')
title('Recovered')

subplot(3,3,5)

plot(death,'-');

hold on
plot(xV(11,:),'--','Color','r')

title('Deaths')


subplot(3,3,6)
plot(vax,'-');
hold on
plot(xV(13,:),'--')

title('Vaccinated')


subplot(3,3,7)
plot(mask,'-');
hold on

plot((xV(2,:)+xV(5,:)+xV(8,:))./sum(xV([1:10,12],:),1) ,'--')
title('Masked')

subplot(3,3,8)
plot(mobility,'-');
hold on

plot(-100*(xV(3,:)+xV(6,:)+xV(9,:))./sum(xV([1:10,12],:),1),'--')
title('Mobility')

subplot(3,3,9)
plot(popDays,'-');
hold on;
plot(sum(xV([1:10, 12],:),1),'--')
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

function [x_kp1] = seirDynamics(xk,eta_Ih,eta_Im,eta_Sm,eta_Sh,dt)
x_kp1=xk;
% state variables
S = xk(1);
Sm = xk(2);
Sh= xk(3);
E = xk(4);
Em = xk(5);
Eh = xk(6);
I = xk(7);
Im = xk(8);
Ih = xk(9);
R = xk(10);
D = xk(11);
U = xk(12);
V = xk(13);

beta = xk(14);

xi_Sh = xk(15);
xi_Eh = xk(15);
xi_Ih = xk(15);

xi_S = xk(16);
xi_E = xk(16);
xi_I = xk(16);

alpha = xk(17);

phi_Sm = xk(18);
phi_Em = xk(18);
phi_Im = xk(18);

phi_S = xk(19);
phi_E = xk(19);
phi_I = xk(19);
% sigma is a single value and is scaled with the proportion of the
% susceptible population
sigma_S = xk(20)*S/(S+Sh+Sm);
sigma_Sm = xk(20)*Sm/(S+Sh+Sm);
sigma_Sh = xk(20)*Sh/(S+Sh+Sm);

kappa_R= xk(21)*S/(S+Sh+Sm);
kappa_Rm= xk(21)*Sm/(S+Sh+Sm);
kappa_Rh= xk(21)*Sh/(S+Sh+Sm);

mu = xk(22);
gamma = xk(23);
epsilon = xk(24);
n = 1000;
Dt = dt/n;


for i = 1:n
    T = S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U;
    lambda = (beta/T)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);
    lambda_m = (beta/T)*(1-eta_Sm)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);
    lambda_h = (beta/T)*(1-eta_Sh)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);

    dS = (-(xi_Sh + alpha + phi_Sm + lambda)*S + kappa_R * R + xi_S *Sh + sigma_S * U + phi_S*Sm)*Dt;
    dSm =(-(lambda_m + phi_S + alpha)*Sm + sigma_Sm * U + kappa_Rm * R + phi_Sm * S)*Dt;
    dSh =(-(alpha + xi_S + lambda_h)*Sh + sigma_Sh * U + xi_Sh * S + kappa_Rh * R)*Dt;

    dE =(-(epsilon + phi_Em + xi_Eh)*E + xi_E * Eh + phi_E * Em + lambda*S )*Dt;
    dEm =(-(epsilon + phi_E) * Em + phi_Em*E + lambda_m* Sm)*Dt;
    dEh =(-(epsilon + xi_E)*Eh + lambda_h * Sh + xi_Eh * E)*Dt;

    dI =(-(mu + gamma + xi_Ih + phi_Im)*I + epsilon*E + phi_I* Im + xi_I *Ih)*Dt;
    dIm =(-(phi_I + gamma + mu)*Im + phi_Im*I + epsilon*Em)*Dt;
    dIh = (-(xi_I + gamma + mu)*Ih + epsilon*Eh + xi_Ih*I)*Dt;

    dR =  (-(kappa_R + kappa_Rm + kappa_Rh)*R+(I + Im + Ih)*gamma)*Dt;
    dD = ((I + Im + Ih)*mu)*Dt;
    dU = (-(sigma_S + sigma_Sm + sigma_Sh)*U+(S + Sh + Sm)*alpha)*Dt;

    totGrads = dS + dSm + dSh + dE + dEm + dEh + dI + dIm + dIh + dR + dD + dU;

    % check if the sum of gradients is below a threshold
    if abs(totGrads) > 1
        disp(totGrads);
        keyboard
    end

    dV =  alpha*(S + Sm+ Sh)*Dt;

    S = S+dS;
    Sm = Sm + dSm;
    Sh = Sh + dSh;

    E = E + dE;
    Em = Em + dEm;
    Eh = Eh + dEh;

    I = I + dI;
    Im = Im + dIm;
    Ih = Ih + dIh;

    R = R + dR;
    D = D + dD;
    U = U + dU;

    V = V + dV;

end


x_kp1(1) = S;
x_kp1(2)=Sm;
x_kp1(3)=Sh;
x_kp1(4)=E ;
x_kp1(5)=Em;
x_kp1(6)=Eh;
x_kp1(7)=I ;
x_kp1(8)=Im;
x_kp1(9)=Ih;
x_kp1(10)=R ;
x_kp1(11)=D ;
x_kp1(12)=U ;
x_kp1(13)=V ;



end

function z = seirObservation(xk)

totPop = sum(xk([1:10,12]));
z(1) = xk(7)+xk(8)+xk(9); %Infectious
z(2) = xk(11); %death
z(3) = xk(13); % Vax
z(4) = (xk(2)+xk(5)+xk(8))/totPop; % Mask
z(5) = -100*(xk(3)+xk(6)+xk(9))/totPop; % Mobility
z(6) = totPop; % population

end