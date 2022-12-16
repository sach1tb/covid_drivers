clearvars
close all

% This script operates on Illinois data scaled to Chicago population and
% uses a linear optimizer to estimate parameters in 2, 3, 4 partitions. It also tries to find optimal partitions. This
% script needs external functions: seir_simIsolation, piecewiseContinuous.m and
% plotResults.m

% Define global parameters 
global infectious death vax mask maxT mobility weights
global error iterCounter

iterCounter = 1;  % Keeps track of the number of days passed
infectious = readmatrix('data/infectiousIllinois.csv');
infectious(isnan(infectious))=0; % replace nan data with 0
infectious = infectious*(9.7/12.8); % Scaling to Chicago population
maxT =size(infectious, 1); % length of dataset

death = readmatrix('data/deathIllinois.csv');
death(isnan(death))=0;
death = cumsum(death);
death = death*(9.7/12.8);


vax = readmatrix('data/vaccinatedIllinois.csv');
vax(isnan(vax))=0;
vax = vax*(9.7/12.8); % Scaling to Chicago population
mask= readmatrix("data/maskIllinois.csv");

mobility = readmatrix("data/mobilityIllinois.csv");
mobility(isnan(mobility))=0;


death = death(1:maxT);
vax = vax(1:maxT);
infectious =  infectious(1:maxT);
mask =  mask(1:maxT);
mobility =  mobility(1:maxT);

% define weights for death, vaccination, infection, mask and mobility
weights = [0.3, 0.3, 0.2, 0.1, 0.1]; % Death and vaccinations are more 
% reliable than infection, mask and mobility data


%%
% sigma_Sh, sigma_S, sigma_Sm : t1 t2 % 3 partitions
% xi1, xi2, phi1, phi2 : t3 t4 t5 % 4 partitions
% gamma, mu : t1 t2 % three partitions
% kappa_Rh, kappa_R, kappa_Rm : t6    % 2 partition
% alpha, epsilon : t1 t2 % three partitions
% beta : t1 t2 % three partitions
% eta_Ih, eta_Im, eta_Sh,eta_Sm % no partitions 
% t1 thru t6 : partitions in time
%%
% 1-46 => sigma_Sh - beta, 47-50 => eta_Ih - eta_Sm, 51-56 => t1-t6 
xmin=[0.001*ones(1,46),  0.001*ones(1,4),  0.001,0.6360,0.001,0.001,0.001,0.001 ];
xmax = [1*ones(1,46),  ones(1,4),  0.482, 0.8,1,1,1,1];

x0 = [0.01*(xmin(1:numel(xmin)-6)+xmax(1:numel(xmin)-6)), 0.2,0.6530,0.3,0.3,0.3,0.3];
options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations', 30000, 'Algorithm', 'interior-point');
xopt = fmincon(@objective,x0,[],[],[],[],xmin,xmax,[],options);



xopt = xopt';
% Create structure to export optimised parameters
modelParams.sigma_ShVector = xopt(1:3);
modelParams.sigma_SVector= xopt(4:6);
modelParams.sigma_SmVector= xopt(7:9);
modelParams.xi1Vector= xopt(10:13); 
modelParams.xi2Vector= xopt(14:17); 
modelParams.phi1Vector = xopt(18:21); 
modelParams.phi2Vector = xopt(22:25);
modelParams.gammaVector = xopt(26:28);
modelParams.muVector = xopt(29:31);
modelParams.kappa_RhVector= xopt(32:33);
modelParams.kappa_RVector= xopt(34:35);
modelParams.kappa_RmVector= xopt(36:37);
modelParams.alphaVector= xopt(38:40);
modelParams.epsilonVector= xopt(41:43);
modelParams.betaVector= xopt(44:46);
modelParams.eta_Ih = xopt(47);
modelParams.eta_Im = xopt(48);
modelParams.eta_Sh = xopt(49);
modelParams.eta_Sm = xopt(50);



save("fminconOptimisedParameters.mat","modelParams",'-mat')
objective(xopt)
plotResults(xopt, weights,error)

function obj = objective(params)
global infectious death mask vax maxT mobility  weights
global error iterCounter
tspan = 0:1:maxT-1;
dt = 1;

% state variable initialisation
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
% sigma_Sh, sigma_S, sigma_Sm : t1 t2 % 3 partitions
% xi1, xi2, phi1, phi2 : t3 t4 t5 % 4 partitions
% gamma, mu : t1 t2 % three partitions
% kappa_Rh, kappa_R, kappa_Rm : t6   % 2 partition
% alpha, epsilon : t1 t2 % three partitions
% beta : t1 t2  % three partitions
% eta_Ih, eta_Im, eta_Sh,eta_Sm % no partitions since exposed->infected
% rate is constant
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
kappa_RhVector= params(32:33);
kappa_RVector= params(34:35);
kappa_RmVector= params(36:37);
alphaVector= params(38:40);
epsilonVector= params(41:43);
betaVector= params(44:46);
eta_Ih = params(47);
eta_Im = params(48);
eta_Sh = params(49);
eta_Sm = params(50);

% estimate the partitions
t1 = 1000*params(51);
t2 = 1000*params(52);
t3 = 1000*params(53);
t4 = 1000*params(54);
t5 = 1000*params(55);
t6 = 1000*params(56);


for k=1:numel(tspan)-1

%%
% sigma_Sh, sigma_S, sigma_Sm : t1 t2 % 3 partitions
% xi1, xi2, phi1, phi2 : t3 t4 t5 % 4 partitions
% gamma, mu : t1 t2 % three partitions
% kappa_Rh, kappa_R, kappa_Rm : t6   % 2 partition
% alpha, epsilon : t1 t2 % three partitions
% beta : t1 t2  % three partitions
% eta_Ih, eta_Im, eta_Sh,eta_Sm % 
% 
%%

    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)] = ...
        seir_simIsolation( x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14), ...
        sigma_ShVector, sigma_SVector, sigma_SmVector, xi1Vector, xi2Vector, phi1Vector, phi2Vector, ...
        gammaVector, muVector, ...
        kappa_RhVector, kappa_RVector, kappa_RmVector, ...
        alphaVector, ...
        epsilonVector, betaVector,eta_Ih,eta_Im,eta_Sh,eta_Sm, ...
        t1, t2, t3, t4, t5, t6, k);

end
% state variables
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

% sub objectives
obj1 = sum((D/max(death)-death/max(death)).^2)^0.5; %w_D
obj3 = sum((Vp/max(vax)-vax/max(vax)).^2)^0.5;  %w_U
obj2 = sum(((I+Im+Ih)/max(infectious)-infectious/max(infectious)).^2)^0.5;  %w_I
maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
maskedReference = mask;
obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
obj4 = sum(((estimatedMobility-mobility/100)).^2)^0.5; % w_H


error(iterCounter) = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5; % calculate the objective function value and store it in an array
iterCounter = iterCounter +1; %advance day
obj = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5; % Calculate the objective function value for optimiser

end