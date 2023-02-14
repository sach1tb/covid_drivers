function [  beta, xi1, xi2, alpha, ...
            phi1, phi2, sigma0, kappa0, ...
            mu, gamma, epsilon]=vec2params(vec)

beta = vec(1);
xi1  = vec(2);
xi2  = vec(3);


% xi_Sh = vec(2);
% xi_Eh = vec(2);
% xi_Ih = vec(2);
% 
% xi_S = vec(3);
% xi_E = vec(3);
% xi_I = vec(3);

alpha = vec(4);

phi1 = vec(5);
phi2 = vec(6);

% phi_Sm = vec(5);
% phi_Em = vec(5);
% phi_Im = vec(5);
% 
% phi_S = vec(6);
% phi_E = vec(6);
% phi_I = vec(6);

% sigma is a single value and is scaled with the proportion of the
% susceptible population
sigma0 = vec(7);
% sigma_S = vec(7)*S/(S+Sh+Sm);
% sigma_Sm = vec(7)*Sm/(S+Sh+Sm);
% sigma_Sh = vec(7)*Sh/(S+Sh+Sm);

kappa0 = vec(8);
% kappa_R= vec(8)*S/(S+Sh+Sm);
% kappa_Rm= vec(8)*Sm/(S+Sh+Sm);
% kappa_Rh= vec(8)*Sh/(S+Sh+Sm);

mu = vec(9);
gamma = vec(10);
epsilon = vec(11);