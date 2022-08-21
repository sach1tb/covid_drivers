clear all
%close all
%clc
global infectious death vax mask maxT numParts mobility magedExposed weights initialInfection deltaVariant omicronVariant
global S0 Sm0 Sh0 E0 Em0 Eh0
global I0 Im0 Ih0 R0 D0 U0 Vp0 N0
global objective counter

objective =[];
counter = 1;
maxT = 1000; %maximum is 1000
infectious = readmatrix('infectiousIllinois.csv');
infectious(isnan(infectious))=0;
infectious = infectious*(9.7/12.8);

death = readmatrix('deathIllinois.csv');
death(isnan(death))=0;
death = cumsum(death);
death = death*(9.7/12.8);


vax = readmatrix('vaccinatedIllinois.csv');
vax(isnan(vax))=0;
vax = vax*(9.7/12.8);
mask= readmatrix("maskIllinois.csv");

mobility = readmatrix("mobilityIllinois.csv");
mobility(isnan(mobility))=0;

magedExposed = readmatrix("authorExposed.csv");
magedTime = magedExposed(:,1);
magedExposed = magedExposed(:,2);
magedExposed = interp1(magedTime,magedExposed,1:700);

death = death(1:maxT);
vax = vax(1:maxT);
infectious =  infectious(1:maxT);
mask =  mask(1:maxT);
mobility =  mobility(1:maxT);


initialInfection = 1:332+151;
deltaVariant = 332+151+1:332+151+153;
omicronVariant = 332+151+153+1:1000;

% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1, 𝑤_𝑀=0.1 
weights = [0.25, 0.3, 0.25, 0.1, 0.1];
% weights = [0.35, 0.1, 0.35, 0.1, 00.1];
numParts = 3;
xmin=[0.0001*ones(1,8) 0.0001 0.0001 0.0001  0.1 0.1*ones(1,7)];
xmax = [1.0*ones(1,8) 0.006 0.1 1 0.5 1*ones(1,7)];

% phi_1
% phi_2
% xi_1
% xi_2
% sigma_S
% sigma_Sm
% sigma_Sh
% alpha ...8
% mu
% gamma
% epsilon
% beta
% eta_Ih
% eta_Im
% eta_Sh
% eta_Sm
% kappa_R
% kappa_Rm
% kappa_Rh ...19

x0 = 0.01*(xmin(1:numel(xmin))+xmax(1:numel(xmin)));
options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',40000,'Algorithm','sqp');

%% Initialize variables #1

S0 = 9.7e6-24;
Sm0 = 0;
Sh0 = 0;
E0 = 0;
Em0 = 0;
Eh0 = 0;
I0 = 24;
Im0 =0 ;
Ih0 = 0;
R0 = 0;
D0 = 0;
U0 = 0;
Vp0 = 0;
N0 = 9.7e6;

xopt1 = fmincon(@objectiveFirstInfection,x0,[],[],[],[],xmin,xmax,[],options);
disp('Objective 1')
beep
beep
pause(5)
[S0,Sm0,Sh0,E0,Em0,Eh0,I0,Im0,Ih0,R0,D0,U0,Vp0,N0] ...
    = getSEIRIsolationEndValues( ...
    S0,Sm0,Sh0,E0,Em0,Eh0,I0,Im0,Ih0,R0,D0,U0,Vp0,N0,xopt1,initialInfection);

xopt2 = fmincon(@objectiveDeltaVariant,x0,[],[],[],[],xmin,xmax,[],options);
disp('Objective 2')
beep
beep
pause(5)
[S0,Sm0,Sh0,E0,Em0,Eh0,I0,Im0,Ih0,R0,D0,U0,Vp0,N0] ...
    = getSEIRIsolationEndValues( ...
    S0,Sm0,Sh0,E0,Em0,Eh0,I0,Im0,Ih0,R0,D0,U0,Vp0,N0,xopt2,deltaVariant);

xopt3 = fmincon(@objectiveOmicronVariant,x0,[],[],[],[],xmin,xmax,[],options);
disp('Objective 3')
beep
beep
pause(5)
plotResultsConcatenated(xopt1,xopt2,xopt3,weights,initialInfection,deltaVariant, omicronVariant,objective)

%% Initial Infection
function obj = objectiveFirstInfection(params)

global S0 Sm0 Sh0 E0 Em0 Eh0
global I0 Im0 Ih0 R0 D0 U0 Vp0 N0
global objective counter

tspan = 0:1:maxT-1;
dt = 1;


x = [];
x(1,1) = S0;
x(1,2) = Sm0;
x(1,3) = Sh0 ;
x(1,4) = E0;
x(1,5) = Em0;
x(1,6) = Eh0;
x(1,7) = I0 ;
x(1,8) = Im0 ;
x(1,9) = Ih0 ;
x(1,10) = R0;
x(1,11) = D0;
x(1,12) = U0 ;
x(1,13) = Vp0;
x(1,14) = N0;

phi_1 = params(1);
phi_2= params(2);
xi_1= params(3);
xi_2= params(4);
sigma_S= params(5);
sigma_Sm= params(6);
sigma_Sh= params(7);
alpha= params(8);
mu= params(9);
gamma= params(10);
epsilon= params(11);
beta= params(12);
eta_Ih= params(13);
eta_Im= params(14);
eta_Sh= params(15);
eta_Sm= params(16);
kappa_R = params(17);
kappa_Rm = params(18);
kappa_Rh = params(19);

for k=1:numel(initialInfection)-1

    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)]...
    = seir_simIsolation( ...
    x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14),phi_1,phi_2,xi_1,xi_2, ...
    sigma_S,sigma_Sm,sigma_Sh,alpha,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh);


end

S= x(:,1);
Sm= x(:,2);
Sh= x(:,3);
E= x(:,4);
Em= x(:,5);
Eh= x(:,6);
I= x(:,7);
Im= x(:,8);
Ih= x(:,9);
R= x(:,10);
D= x(:,11);
U= x(:,12);
Vp= x(:,13);
N = 9.7e6;

obj1 = sum((D/N-ddeath/N).^2)^0.5; %w_D
obj3 = sum((Vp/N-vvax/N).^2)^0.5;  %w_U
obj2 = sum(((I+Im+Ih)/N-iinfectious/N).^2)^0.5;  %w_I
maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
maskedReference = mmask;
obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
obj4 = sum(((estimatedMobility-mmobility/100)).^2)^0.5; % w_H
% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1, 𝑤_𝑀=0.1 
objective(counter) = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
counter = counter +1;
obj = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
%obj = sum((E+Em+Eh - magedExposed').^2)^0.5;
%obj = obj1*obj3;
end

function obj = objectiveDeltaVariant(params)

global infectious death mask vax maxT numParts mobility  weights deltaVariant
global S0 Sm0 Sh0 E0 Em0 Eh0
global I0 Im0 Ih0 R0 D0 U0 Vp0 N0
global objective counter

iinfectious = infectious(deltaVariant);
ddeath = death(deltaVariant);
mmask = mask(deltaVariant);
vvax = vax(deltaVariant);
mmobility = mobility(deltaVariant);

tspan = 0:1:maxT-1;
dt = 1;


x = [];
x(1,1) = S0;
x(1,2) = Sm0;
x(1,3) = Sh0 ;
x(1,4) = E0;
x(1,5) = Em0;
x(1,6) = Eh0;
x(1,7) = I0 ;
x(1,8) = Im0 ;
x(1,9) = Ih0 ;
x(1,10) = R0;
x(1,11) = D0;
x(1,12) = U0 ;
x(1,13) = Vp0;
x(1,14) = N0;

phi_1 = params(1);
phi_2= params(2);
xi_1= params(3);
xi_2= params(4);
sigma_S= params(5);
sigma_Sm= params(6);
sigma_Sh= params(7);
alpha= params(8);
mu= params(9);
gamma= params(10);
epsilon= params(11);
beta= params(12);
eta_Ih= params(13);
eta_Im= params(14);
eta_Sh= params(15);
eta_Sm= params(16);
kappa_R = params(17);
kappa_Rm = params(18);
kappa_Rh = params(19);

for k=1:numel(deltaVariant)-1

    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)]...
    = seir_simIsolation( ...
    x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14),phi_1,phi_2,xi_1,xi_2, ...
    sigma_S,sigma_Sm,sigma_Sh,alpha,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh);

end

S= x(:,1);
Sm= x(:,2);
Sh= x(:,3);
E= x(:,4);
Em= x(:,5);
Eh= x(:,6);
I= x(:,7);
Im= x(:,8);
Ih= x(:,9);
R= x(:,10);
D= x(:,11);
U= x(:,12);
Vp= x(:,13);
N = 9.7e6;

obj = errorCalc(S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,UVp,N);

obj1 = sum((D/N-ddeath/N).^2)^0.5; %w_D
obj3 = sum((Vp/N-vvax/N).^2)^0.5;  %w_U
obj2 = sum(((I+Im+Ih)/N-iinfectious/N).^2)^0.5;  %w_I
maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
maskedReference = mmask;
obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
obj4 = sum(((estimatedMobility-mmobility/100)).^2)^0.5; % w_H
% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1, 𝑤_𝑀=0.1 
objective(counter) = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
counter = counter +1;
obj = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
%obj = sum((E+Em+Eh - magedExposed').^2)^0.5;
%obj = obj1*obj3;
end

function obj = objectiveOmicronVariant(params)

global infectious death mask vax maxT numParts mobility  weights omicronVariant
global S0 Sm0 Sh0 E0 Em0 Eh0
global I0 Im0 Ih0 R0 D0 U0 Vp0 N0
global objective counter

iinfectious = infectious(omicronVariant);
ddeath = death(omicronVariant);
mmask = mask(omicronVariant);
vvax = vax(omicronVariant);
mmobility = mobility(omicronVariant);

tspan = 0:1:maxT-1;
dt = 1;

x = [];
x(1,1) = S0;
x(1,2) = Sm0;
x(1,3) = Sh0 ;
x(1,4) = E0;
x(1,5) = Em0;
x(1,6) = Eh0;
x(1,7) = I0 ;
x(1,8) = Im0 ;
x(1,9) = Ih0 ;
x(1,10) = R0;
x(1,11) = D0;
x(1,12) = U0 ;
x(1,13) = Vp0;
x(1,14) = N0;

phi_1 = params(1);
phi_2= params(2);
xi_1= params(3);
xi_2= params(4);
sigma_S= params(5);
sigma_Sm= params(6);
sigma_Sh= params(7);
alpha= params(8);
mu= params(9);
gamma= params(10);
epsilon= params(11);
beta= params(12);
eta_Ih= params(13);
eta_Im= params(14);
eta_Sh= params(15);
eta_Sm= params(16);
kappa_R = params(17);
kappa_Rm = params(18);
kappa_Rh = params(19);

for k=1:numel(omicronVariant)-1

    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)]...
    = seir_simIsolation( ...
    x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14),phi_1,phi_2,xi_1,xi_2, ...
    sigma_S,sigma_Sm,sigma_Sh,alpha,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh);



% 
%     [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
%         x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)] = ...
%         ...
%         seir_simIsolation     ( x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
%         x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14) ...
%         ,a1,a2,b1,b2,partT, sigma_Smax,sigma_Sm_max,sigma_Sh_max, c, kFactor ...
%         ,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,t0,k,numParts );


end

S= x(:,1);
Sm= x(:,2);
Sh= x(:,3);
E= x(:,4);
Em= x(:,5);
Eh= x(:,6);
I= x(:,7);
Im= x(:,8);
Ih= x(:,9);
R= x(:,10);
D= x(:,11);
U= x(:,12);
Vp= x(:,13);
N = 9.7e6;

obj1 = sum((D/N-ddeath/N).^2)^0.5; %w_D
obj3 = sum((Vp/N-vvax/N).^2)^0.5;  %w_U
obj2 = sum(((I+Im+Ih)/N-iinfectious/N).^2)^0.5;  %w_I
maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
maskedReference = mmask;
obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
obj4 = sum(((estimatedMobility-mmobility/100)).^2)^0.5; % w_H
% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1, 𝑤_𝑀=0.1 
objective(counter) = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
counter = counter +1;
obj = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
%obj = sum((E+Em+Eh - magedExposed').^2)^0.5;
%obj = obj1*obj3;
end
