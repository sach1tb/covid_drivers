clearvars
close all
clc

global infectious death vax mask maxT mobility weights
global error 

infectious = csvread('data/infectiousIllinois.csv');
infectious(isnan(infectious))=0; % replace nan data with 0
infectious = infectious*(9.7/12.8); % Scaling to Chicago population
maxT =335; % length of dataset

death = csvread('data/deathIllinois.csv');
death(isnan(death))=0;
death = cumsum(death);
death = death*(9.7/12.8);


vax = csvread('data/vaccinatedIllinois.csv');
vax(isnan(vax))=0;
vax = vax*(9.7/12.8); % Scaling to Chicago population
mask= csvread("data/maskIllinois.csv");

mobility = csvread("data/mobilityIllinois.csv");
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
% sigma_Sh, sigma_S, sigma_Sm
% xi1, xi2, phi1, phi2 
% gamma, mu
% kappa_Rh, kappa_R, kappa_Rm
% alpha, epsilon 
% beta
% eta_Ih, eta_Im, eta_Sh,eta_Sm 
%%

xmin = zeros(1,19);
xmax = 1*ones(1,19);
xmax(5)=0.1; 
xmax(1:3) = 0.00; % all sigmas for vaccination
xmax(13) = 0.00; % alpha for vaccination
xmax(14) = 0.5; % epsilon
xmax(15) = 0.5; % beta
x0 = 0.5*(xmin+xmax); %initial guess values
options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',...
                3000, 'Algorithm', 'interior-point');
xopt = fmincon(@objective,x0,[],[],[],[],xmin,xmax,[],options);



xopt = xopt';
% Create structure to export optimised parameters
modelParams.sigma0 = xopt(1);
modelParams.xi1= xopt(2); 
modelParams.xi2= xopt(3); 
modelParams.phi1 = xopt(4); 
modelParams.phi2 = xopt(5);
modelParams.gamma = xopt(6);
modelParams.mu = xopt(7);
modelParams.kappa0= xopt(8);
modelParams.alpha= xopt(9);
modelParams.epsilon= xopt(10);
modelParams.beta= xopt(11);
modelParams.eta_Ih = xopt(12);
modelParams.eta_Im = xopt(13);
modelParams.eta_Sh = xopt(14);
modelParams.eta_Sm = xopt(15);


save("fminconOptimisedParameters.mat","modelParams",'-mat')
objective(xopt)
plotResultsOptimization(xopt, weights)

function obj = objective(params)
global infectious death mask vax maxT mobility  weights
global error
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
% sigma_Sh, sigma_S, sigma_Sm
% xi1, xi2, phi1, phi2 
% gamma, mu
% kappa_Rh, kappa_R, kappa_Rm
% alpha, epsilon 
% beta
% eta_Ih, eta_Im, eta_Sh,eta_Sm 
%%

sigma0 = params(1);
xi1= params(2); 
xi2= params(3); 
phi1 = params(4); 
phi2 = params(5);
gamma = params(6);
mu = params(7);
kappa0= params(8);
alpha= params(9);
epsilon= params(10);
beta= params(11);
eta_Ih = params(12);
eta_Im = params(13);
eta_Sh = params(14);
eta_Sm = params(15);


for k=1:numel(tspan)-1

%%
% sigma_Sh, sigma_S, sigma_Sm
% xi1, xi2, phi1, phi2 
% gamma, mu
% kappa_Rh, kappa_R, kappa_Rm
% alpha, epsilon 
% beta
% eta_Ih, eta_Im, eta_Sh,eta_Sm 
%%


    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)] = ...
        seirDynamicsforOptimization( x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14), ...
        sigma0, xi1, xi2, phi1, phi2, ...
        gamma, mu, ...
        kappa0, ...
        alpha, ...
        epsilon, beta,eta_Ih,eta_Im,eta_Sh,eta_Sm);

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
% obj1 = sum((D/max(death)-death/max(death)).^2)^0.5; %w_D
% obj3 = sum((Vp/max(vax)-vax/max(vax)).^2)^0.5;  %w_U
% obj2 = sum(((I+Im+Ih)/max(infectious)-infectious/max(infectious)).^2)^0.5;  %w_I
% maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U);
% maskedReference = mask;
% obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
% estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U);
% obj4 = sum(((estimatedMobility-mobility/100)).^2)^0.5; % w_H
% error = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5; % calculate the objective function value and store it in an array
% obj = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5; % Calculate the objective function value for optimiser


obj = sum((D-death).^2)^0.5; %w_D
end