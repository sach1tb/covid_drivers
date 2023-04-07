clearvars
close all

global infectious death vax mask maxT mobility weights
global error 

infectious = csvread('data/infectiousIllinois_ci.csv');
infectious=infectious(:,2);
infectious(isnan(infectious))=0; % replace nan data with 0
infectious = infectious*(9.7/12.8); % Scaling to Chicago population
maxT =335; % length of dataset

death = csvread('data/deathIllinois.csv');
% death = death(:,2);
death(isnan(death))=0;
death = cumsum(death);
death = death*(9.7/12.8);


vax = csvread('data/vaccinatedIllinois.csv');
% vax = vax(:,2);
vax(isnan(vax))=0;
vax = vax*(9.7/12.8); % Scaling to Chicago population


mask= csvread("data/maskIllinois.csv");
% mask= mask(:,2);

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
% beta (1), xi2 (2), xi1 (3), alpha (4), phi2 (5), 
% phi1 (6), sigma0 (7), kappa0 (8), mu (9), gamma (10), epsilon (11)
%%

xmin = zeros(1,15);
% sigma0=0, alpha=0, for vaccination
% mu = 0.01 for mortality
% epsilon = 0.5;
xmax = [0.5, .1, 1, 0, 1, ...
        1, 0, 1, 0.01, 0.1, ...
        0.5, 1, 1, 1, 1];


x0 = 0.5*(xmin+xmax); %initial guess values
options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',...
                2000);
xopt = fmincon(@objective_fun,x0,[],[],[],[],xmin,xmax,[],options);



% Create structure to export optimised parameters

[  beta, xi1, xi2, alpha, ...
            phi1, phi2, sigma0, kappa0, ...
            mu, gamma, epsilon]=vec2params(xopt(1:11));

modelParams.beta=beta;
modelParams.xi1=xi1;
modelParams.xi2=xi2;
modelParams.alpha=alpha;
modelParams.phi1=phi1;
modelParams.phi2=phi2;
modelParams.sigma0=sigma0;
modelParams.kappa0=kappa0;
modelParams.mu=mu;
modelParams.gamma=gamma; 
modelParams.epsilon=epsilon;
modelParams.eta_Ih = xopt(12);
modelParams.eta_Im = xopt(13);
modelParams.eta_Sh = xopt(14);
modelParams.eta_Sm = xopt(15);


save("fminconOptimisedParameters.mat","modelParams",'-mat')
modelParams
[val, x]=objective_fun(xopt);
plotResultsOptimization(x', maxT)

function [obj, x] = objective_fun(params)
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

x=x';




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


%     [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
%         x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)] = ...
%         seirDynamicsforOptimization( x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
%         x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14), ...
%         sigma0, xi1, xi2, phi1, phi2, ...
%         gamma, mu, ...
%         kappa0, ...
%         alpha, ...
%         epsilon, beta,eta_Ih,eta_Im,eta_Sh,eta_Sm);

    xkp1=seirDynamics([x(:,k); params(1:11)'], params(12), params(13), params(14), params(15),1);
    x(:,k+1)=xkp1(1:13,1);

end
% state variables
% S= x(:,1);
% Sm= x(:,2);
% Sh= x(:,3);
% E= x(:,4);
% Em= x(:,5);
% Eh= x(:,6);
% I= x(:,7);
% Im= x(:,8);
% Ih= x(:,9);
% R= x(:,10);
% D= x(:,11);
% U= x(:,12);
% Vp= x(:,13);
% N = 9.7e6;

% sub objective_funs
% obj1 = sum((D/max(death)-death/max(death)).^2)^0.5; %w_D
% obj3 = sum((Vp/max(vax)-vax/max(vax)).^2)^0.5;  %w_U
% obj2 = sum(((I+Im+Ih)/max(infectious)-infectious/max(infectious)).^2)^0.5;  %w_I
% maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U);
% maskedReference = mask;
% obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
% estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U);
% obj4 = sum(((estimatedMobility-mobility/100)).^2)^0.5; % w_H
% error = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5; % calculate the objective_fun function value and store it in an array
% obj = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5; % Calculate the objective_fun function value for optimiser


obj1 = sqrt(sum((x(11,:)'-death).^2)); %w_D
% obj2 = sqrt(sum((sum(x(7:9,:))'-infectious).^2));
obj=obj1;
end