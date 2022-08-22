clear all
%close all
%clc
global infectious death vax mask maxT numParts mobility magedExposed weights initialInfection deltaVariant omicronVariant
global S0 Sm0 Sh0 E0 Em0 Eh0
global I0 Im0 Ih0 R0 D0 U0 Vp0 N0
global objective counter dayCounter

dayCounter = 1;

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
weights = [0.30, 0.2, 0.3, 0.1, 0.1];
% weights = [0.35, 0.1, 0.35, 0.1, 00.1];
numParts = 3;
xmin=[0.0001*ones(1,8) 0.0001 0.0001 0.0001  0.001 0.001*ones(1,7)];
xmax = [1.0*ones(1,8) 0.01 0.1 1 1.0 1*ones(1,7)];

% phi_1
% phi_2
% xi_1
% xi_2
% sigma_S
% sigma_Sm
% sigma_Sh
% c ...8
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

% x0 = 0.01*(xmin(1:numel(xmin))+xmax(1:numel(xmin)));
x0 =[0.0334093116883648	0.000100000000000000	0.883632333755389 ...
    0.638394824749586	0.204445550120683	0.00398370549320305 ...
    0.000100000000000000	0.00276238063762399 ...
    0.000165558383109743	0.100000000000000	0.311245339544135 ...
    0.932802081022819	0.0680480998745700	0.963223381107490	...
    0.839156884468624	0.401332688593568	0.799977066838738 ...	
    0.00100000000000000	1];
x1 = [1	0.705961434080331	0.294259587563625	1	0.000100000000000000	...
    1	0.000100000000000000	0.000100000000000000	0.000100000000000000 ...
    0.100000000000000	0.000100000000000000	1	0.00100000000000000 ...
    0.00100000000000000	0.00100000000000000	0.00100000000000000 ...
    0.00100000000000000	0.00100000000000000	0.00946664522173071];
x2 = [0.000100000000000000	0.00379504767884131	0.000100000000000000	...
    0.00797967618357940	0.000100000000000000	0.000100000000000000 ...
    1	0.00227633460741759	0.000100000000000000	0.00412435537665440 ...
    0.000100000000000000	1	0.00100000000000000	0.00100000000000000 ...
    0.00100000000000000	1	0.00100000000000000	0.104472268472839 ...
    0.00100000000000000];

options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',40000,'Algorithm','interior-point');

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
    S0,Sm0,Sh0,E0,Em0,Eh0,I0,Im0,Ih0,R0,D0,U0,Vp0,N0,xopt1,initialInfection,400);

xopt2 = fmincon(@objectiveDeltaVariant,x1,[],[],[],[],xmin,xmax,[],options);
disp('Objective 2')
beep
beep
pause(5)
[S0,Sm0,Sh0,E0,Em0,Eh0,I0,Im0,Ih0,R0,D0,U0,Vp0,N0] ...
    = getSEIRIsolationEndValues( ...
    S0,Sm0,Sh0,E0,Em0,Eh0,I0,Im0,Ih0,R0,D0,U0,Vp0,N0,xopt2,deltaVariant,400);

xopt3 = fmincon(@objectiveOmicronVariant,x2,[],[],[],[],xmin,xmax,[],options);
disp('Objective 3')
beep
beep
pause(5)

plotResultsConcatenated(xopt1,xopt2,xopt3,weights,initialInfection,deltaVariant, omicronVariant,objective)

%% Initial Infection
function obj = objectiveFirstInfection(params)

global initialInfection
global S0 Sm0 Sh0 E0 Em0 Eh0
global I0 Im0 Ih0 R0 D0 U0 Vp0 N0
global objective counter


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
c= params(8);
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
    sigma_S,sigma_Sm,sigma_Sh,c,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh, k);

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


objective(counter) = errorCalc(S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,initialInfection);
counter = counter +1;
obj = errorCalc(S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,initialInfection);
end

function obj = objectiveDeltaVariant(params)

global S0 Sm0 Sh0 E0 Em0 Eh0
global I0 Im0 Ih0 R0 D0 U0 Vp0 N0
global objective counter deltaVariant


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
c= params(8);
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
    sigma_S,sigma_Sm,sigma_Sh,c,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh, 400);


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

objective(counter) = errorCalc(S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,deltaVariant);
counter = counter +1;
obj = errorCalc(S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,deltaVariant);

end

function obj = objectiveOmicronVariant(params)

global  omicronVariant
global S0 Sm0 Sh0 E0 Em0 Eh0
global I0 Im0 Ih0 R0 D0 U0 Vp0 N0
global objective counter

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
c= params(8);
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
    sigma_S,sigma_Sm,sigma_Sh,c,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh,400);
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


objective(counter) = errorCalc(S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,omicronVariant);
counter = counter +1;
obj = errorCalc(S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,omicronVariant);

end


function obj = errorCalc(S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,duration)

global infectious death mask vax mobility  weights

iinfectious = infectious(duration);
ddeath = death(duration);
mmask = mask(duration);
vvax = vax(duration);
mmobility = mobility(duration);


obj1 = sum((D/max(ddeath)-ddeath/max(ddeath)).^2)^0.5; %w_D
obj3 = sum((Vp/max(vvax)-vvax/max(vvax)).^2)^0.5;  %w_U
obj2 = sum(((I+Im+Ih)/max(iinfectious)-iinfectious/max(iinfectious)).^2)^0.5;  %w_I
maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
maskedReference = mmask;
obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
obj4 = sum(((estimatedMobility-mmobility/100)).^2)^0.5; % w_H
obj = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
% obj = max([obj1,obj2,obj3,obj4,obj5]);

% obj1 = max(abs(D/N-ddeath/N)); %w_D
% obj3 = max(abs(Vp/N-vvax/N));  %w_U
% obj2 = max(abs((I+Im+Ih)/N-iinfectious/N));  %w_I
% maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
% maskedReference = mmask;
% obj5 = max(abs(maskedEstimated-maskedReference)); % w_M
% estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
% obj4 = max(abs((estimatedMobility-mmobility/100))); % w_H
% 
% obj = obj1;

end
