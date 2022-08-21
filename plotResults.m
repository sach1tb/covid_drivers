function plotResults(params,numTerms,maxT, weights,part,figNo)

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



death = death(part);
vax = vax(part);
infectious =  infectious(part);
mask =  mask(part);
mobility =  mobility(part);

tspan = part;
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

phi1 = [];
phi2 = [];
eta1 = [];
eta2 = [];
lambda = [];
lambda_m = [];
lambda_h = [];

for k=1:numel(part)-1

    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)...
        ,phi1(k),phi2(k),xi1(k),xi2(k),lambda(k),lambda_m(k),lambda_h(k)]...
    = seir_simIsolation( ...
    x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14),phi_1,phi_2,xi_1,xi_2, ...
    sigma_S,sigma_Sm,sigma_Sh,alpha,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh );



% 
%     [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
%         x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)] = ...
%         ...
%         seir_simIsolation     ( x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
%         x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14) ...
%         ,a1,a2,b1,b2,partT, sigma_Smax,sigma_Sm_max,sigma_Sh_max, c, kFactor ...
%         ,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,t0,k,numParts );
end



% for k=1:numel(tspan)-1
% 
% 
%     [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
%         x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)...
%         ,phi1(k),phi2(k),xi1(k),xi2(k),lambda(k),lambda_m(k),lambda_h(k)] = ...
%         seir_simIsolation(x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
%         x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14) ...
%         ,a1,a2,b1,b2,Omega1,Omega2,Theta1,Theta2,zeta1,zeta2,psi1,psi2, ...
%         sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,kFactor,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,t0,k,numTerms);
% end
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
figure(figNo);
clf;
subplot(3,3,1)
plot(tspan,S+Sm+Sh)
hold on
maskedEstimated = (Sm+Em+Im)./(S+E+I+Sm+Em+Im+Sh+Eh+Ih+R+D+U);
maskedReference = mask;
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R);

title('Susceptible')


subplot(3,3,2)
plot(D)
hold on
plot(death)
legend('Estimated','Real Data')
title('Deaths')


subplot(3,3,3)
plot(I+Im+Ih)
hold on
plot(infectious)
title('Infectious')


subplot(3,3,4)
plot(Vp)
hold on
plot(vax)
title('Vaccinated Perpetual')

subplot(3,3,5)
plot(E+Em+Eh)
title('Exposed')

subplot(3,3,6)
plot(maskedEstimated)
hold on
plot(maskedReference)
title('Masked Estimated')

subplot(3,3,7)
plot(estimatedMobility)
hold on
plot(mobility/100)
title('Mobility')

subplot(3,3,8)
plot(R)
hold on

title('Recovered')


subplot(3,3,9)
plot(phi1);
hold on
plot(phi2);
plot(xi1);
plot(xi2);

legend('\phi_1','\phi_2','\xi_1','\xi_2')
% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1,𝑤_𝑀=0.1 
str = sprintf("w_D = %0.3f, w_I = %0.3f,w_U = %0.3f, w_H = %0.3f, w_M = %0.3f", ...
    weights(1),weights(2),weights(3),weights(4),weights(5));
str2 =  sprintf("\n \\sigma_{S} = %0.4f, \\sigma_{Sm} = %0.4f, \\sigma_{Sh} = %0.4f, \\alpha= %0.4f," + ...
    "\n" + "\\beta = %0.4f, \\mu = %0.4f, \\epsilon = %0.4f, \\gamma = %0.4f, \\eta_{I_h}=%0.4f," + ...
    "\\eta_{I_m} = %0.4f,\\eta_{S_h} = %0.4f,\\eta_{S_m} = %0.4f, \\kappa_R= %0.4f, \\kappa_R_m= %0.4f, \\kappa_R_h= %0.4f",sigma_S,sigma_Sm ...
    ,sigma_Sh,alpha,beta,mu ,epsilon,gamma,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh);

sgtitle(str+str2);



end