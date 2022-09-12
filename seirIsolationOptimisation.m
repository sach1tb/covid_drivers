clear all
%close all
%clc
global infectious death vax mask maxT mobility magedExposed weights
global error iterCounter
iterCounter = 1;
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


% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1, 𝑤_𝑀=0.1 
% weights = [0.9, 0.0, 0.1, 0.0, 0.0];
weights = [0.3, 0.2, 0.3, 0.1, 0.1];

% weights = [0.5, 0.2, 0.5, 0.0, 0.0];

%%
% sigma_Sh, sigma_S, sigma_Sm : t1 t2
% xi1, xi2, phi1, phi2 : t3 t4 t5
% kappa_Rh, kappa_R, kappa_Rm : t6  ...10
% gamma, mu : t1 t2
% alpha, epsilon : t1 t2
% beta : t1 t2  ...15
% eta_Ih, eta_Im, eta_Sh,eta_Sm
%%
xmin=[0.001*ones(1,46),  0.001*ones(1,4),  0.001,0.6360,0.001,0.001,0.001,0.001 ];
xmax = [1*ones(1,46),  ones(1,4),  0.482, 0.8,1,1,1,1];

x0 = [0.01*(xmin(1:numel(xmin)-6)+xmax(1:numel(xmin)-6)), 0.2,0.6530,0.3,0.3,0.3,0.3];
options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations', 30000, 'Algorithm', 'interior-point');
xopt = fmincon(@objective,x0,[],[],[],[],xmin,xmax,[],options);
% options = optimset('Display','iter-detailed','MaxFunEvals', 30000);
% xopt = fminsearch(@objective,x0,options);

% options = optimoptions('patternsearch','Display','iter','MaxFunctionEvaluations', 100000);
% xopt = patternsearch(@objective,x0,[],[],[],[],xmin,xmax,[],options);

% options = optimoptions('gamultiobj','Display','iter');
% xopt = gamultiobj(@objective,56,[],[],[],[],xmin,xmax,options);


xopt
objective(xopt)
plotResults(xopt, weights,error)

function obj = objective(params)
global infectious death mask vax maxT numParts mobility  weights
global error iterCounter
tspan = 0:1:maxT-1;
dt = 1;

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



%%
% sigma_Sh, sigma_S, sigma_Sm : t1 t2
% xi1, xi2, phi1, phi2 : t3 t4 t5
% kappa_Rh, kappa_R, kappa_Rm : t6  ...10
% gamma, mu : t1 t2
% alpha, epsilon : t1 t2
% beta : t1 t2  ...15
% eta_Ih, eta_Im, eta_Sh,eta_Sm
%%
sigma_ShVector = params(1:3);
sigma_SVector= params(4:6);
sigma_SmVector= params(7:9);
xi1Vector= params(10:13); 
xi2Vector= params(14:17); 
phi1Vector = params(18:21); 
phi2Vector = params(22:25);
gammaVector = params(26:28);
muVector = params(29:31);
kappa_RhVector= 0*params(32:33);
kappa_RVector= 0*params(34:35);
kappa_RmVector= 0*params(36:37);
alphaVector= params(38:40);
epsilonVector= params(41:43);
betaVector= params(44:46);
eta_Ih = params(47);
eta_Im = params(48);
eta_Sh = params(49);
eta_Sm = params(50);


t1 = 1000*params(51);
t2 = 1000*params(52);
t3 = 1000*params(53);
t4 = 1000*params(54);
t5 = 1000*params(55);
t6 = 1000*params(56);


for k=1:numel(tspan)-1

%%
% sigma_Sh, sigma_S, sigma_Sm : t1 t2
% xi1, xi2, phi1, phi2 : t3 t4 t5
% kappa_Rh, kappa_R, kappa_Rm : t6  ...10
% gamma, mu : t1 t2
% alpha, epsilon : t1 t2
% beta : t1 t2  ...15
% eta_Ih, eta_Im, eta_Sh,eta_Sm
%%


    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)] = ...
        seir_simIsolation( x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14), ...
        sigma_ShVector, sigma_SVector, sigma_SmVector, xi1Vector, xi2Vector, phi1Vector, phi2Vector ...
        ,kappa_RhVector, kappa_RVector, kappa_RmVector, ...
        gammaVector, muVector, alphaVector ...
        , epsilonVector, betaVector,eta_Ih,eta_Im,eta_Sh,eta_Sm,t1, t2, t3, t4, t5, t6,k);
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

% obj1 = sum((D/N-death/N).^2)^0.5; %w_D
% obj3 = sum((Vp/N-vax/N).^2)^0.5;  %w_U
% obj2 = sum(((I+Im+Ih)/N-infectious/N).^2)^0.5;  %w_I
% maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
% maskedReference = mask;
% obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
% estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
% obj4 = sum(((estimatedMobility-mobility/100)).^2)^0.5; % w_H

obj1 = sum((D/max(death)-death/max(death)).^2)^0.5; %w_D
obj3 = sum((Vp/max(vax)-vax/max(vax)).^2)^0.5;  %w_U
obj2 = sum(((I+Im+Ih)/max(infectious)-infectious/max(infectious)).^2)^0.5;  %w_I
maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
maskedReference = mask;
obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
obj4 = sum(((estimatedMobility-mobility/100)).^2)^0.5; % w_H


% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1, 𝑤_𝑀=0.1 
error(iterCounter) = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
iterCounter = iterCounter +1;
obj = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
% obj = max([weights(1)*obj1,weights(2)*obj2,weights(3)*obj3,weights(4)*obj4,weights(5)*obj5]);
%obj = sum((E+Em+Eh - magedExposed').^2)^0.5;
%obj = obj1*obj3;
end