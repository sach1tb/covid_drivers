function plotResults(params,numTerms,maxT, weights)

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



death = death(1:maxT);
vax = vax(1:maxT);
infectious =  infectious(1:maxT);
mask =  mask(1:maxT);
mobility =  mobility(1:maxT);

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

a1 = params(1+0*numTerms:1*numTerms);
a2 = params(1+1*numTerms:2*numTerms);
b1 = params(1+2*numTerms:3*numTerms);
b2 = params(1+3*numTerms:4*numTerms);


static = 4*numTerms+1;
nTerms = numTerms-1;
Omega1 = params(static:static+nTerms-1);
Omega2 = params(static+nTerms:static+2*nTerms-1);
Theta1 = params(static+2*nTerms:static+3*nTerms-1);
Theta2 = params(static+3*nTerms:static+4*nTerms-1);
zeta1 = params(static+4*nTerms:static+5*nTerms-1);
zeta2 = params(static+5*nTerms:static+6*nTerms-1);
psi1 = params(static+6*nTerms:static+7*nTerms-1);
psi2 = params(static+7*nTerms:static+8*nTerms-1);

static = static+8*nTerms-1;
sigma_Smax  = params(static+1);
sigma_Sm_max = params(static+2);
sigma_Sh_max= params(static+3);
c = params(static+4);
t0  = params(static+5);
beta  = params(static+6);
mu  = params(static+7); 
epsilon  = params(static+8);
gamma  = params(static+9);
eta_Ih  = params(static+10);
eta_Im  = params(static+11);
eta_Sh  = params(static+12);
eta_Sm  = params(static+13);

kFactor = 5;

a1
a2
b1
b2
% % Omega1
% % Omega2
% % Theta1
% % Theta2
% % zeta1
% % zeta2
% % psi1 
% % psi2 
sigma_Smax
sigma_Sm_max
sigma_Sh_max
c
t0 
beta
mu 
epsilon
gamma
eta_Ih
eta_Im
eta_Sh
eta_Sm


phi1 = [];
phi2 = [];
eta1 = [];
eta2 = [];
lambda = [];
lambda_m = [];
lambda_h = [];
for k=1:numel(tspan)-1


    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)...
        ,phi1(k),phi2(k),xi1(k),xi2(k),lambda(k),lambda_m(k),lambda_h(k)] = ...
        seir_simIsolation(x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14) ...
        ,a1,a2,b1,b2,Omega1,Omega2,Theta1,Theta2,zeta1,zeta2,psi1,psi2, ...
        sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,kFactor,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,t0,k,numTerms);
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
figure(1);
clf;
subplot(3,3,1)
plot(tspan,S+Sm+Sh)
hold on
maskedEstimated = (Sm+Em+Im)./(S+E+I+Sm+Em+Im+Sh+Eh+Ih+R+D+U);
maskedReference = mask;
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R);

title('Susceptible')


subplot(3,3,2)
plot(tspan,D)
hold on
plot(tspan,death)
legend('Estimated','Real Data')
title('Deaths')


subplot(3,3,3)
plot(tspan,I+Im+Ih)
hold on
plot(infectious)
title('Infectious')


subplot(3,3,4)
plot(tspan,Vp)
hold on
plot(vax)
title('Vaccinated Perpetual')

subplot(3,3,5)
plot(tspan,E+Em+Eh)
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
str2 =  sprintf("\n \\sigma_{S_{max}} = %0.3f, \\sigma_{Sm_{max}} = %0.3f, \\sigma_{Sh_{max}} = %0.3f, c= %0.3f, t_0 = %0.3f," + ...
    "\n" + "\\beta = %0.3f, \\mu = %0.3f, \\epsilon = %0.3f, \\gamma = %0.3f, \\eta_{I_h}=%0.3f," + ...
    "\\eta_{I_m} = %0.3f,\\eta_{S_h} = %0.3f,\\eta_{S_m} = %0.3f",sigma_Smax,sigma_Sm_max ...
    ,sigma_Sh_max,c,t0,beta,mu ,epsilon,gamma,eta_Ih,eta_Im,eta_Sh,eta_Sm);

sgtitle(str+str2);



end