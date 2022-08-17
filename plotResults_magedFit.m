function plotResults_magedFit(params,numTerms,maxT)

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

tspan = 0:1:maxT-1;
dt = 1;

S0 = 1e6-50;
Sm0 = 0;
Sh0 = 0;
E0 = 0;
Em0 = 0;
Eh0 = 0;
I0 = 50;
Im0 =0 ;
Ih0 = 0;
R0 = 0;
D0 = 0;
U0 = 0;
Vp0 = 0;
N0 = 1e6;

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

% a1 = params(1+0*numTerms:1*numTerms)*0;
% a2 = params(1+1*numTerms:2*numTerms)*0;
% b1 = params(1+2*numTerms:3*numTerms)*0;
% b2 = params(1+3*numTerms:4*numTerms)*0;
% 
% Omega1 = params(1+4*numTerms:5*numTerms)*0;
% Omega2 = params(1+5*numTerms:6*numTerms)*0;
% Theta1 = params(1+6*numTerms:7*numTerms)*0;
% Theta2 = params(1+7*numTerms:8*numTerms)*0;
% zeta1 = params(1+8*numTerms:9*numTerms)*0;
% zeta2 = params(1+9*numTerms:10*numTerms)*0;
% psi1 = params(1+10*numTerms:11*numTerms)*0;
% psi2 = params(1+11*numTerms:12*numTerms)*0;
% sigma_Smax = params(12*numTerms+1)*0;
% sigma_Sm_max= params(12*numTerms+2)*0;
% sigma_Sh_max= params(12*numTerms+3)*0;
% c= params(12*numTerms+4)*0;
mu = 3.5/100;
gamma = 0.0602;
epsilon = 1/4.5;
beta = 0.312;
eta_i= params(12*numTerms+5);
eta_s = 0.35;
eta_h = 0.45;
t0 = 300;
%kappa = 0*ones(1,numel(tspan));

a1 = 0;
a2 = 0;
b1 = 0;
b2 = 0;

Omega1 = 0;
Omega2 = 0;
Theta1 = 0;
Theta2 = 0;
zeta1 = 0;
zeta2 = 0;
psi1 = 0;
psi2 = 0;
sigma_Smax = 0;
sigma_Sm_max= 0;
sigma_Sh_max= 0;
c= 0;
eta_i= 0;




phi1 = [];
phi2 = [];
eta1 = [];
eta2 = [];
labmda = [];
lambda_m = [];
lambda_h = [];
for k=1:numel(tspan)-1


    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)...
        ,phi1(k),phi2(k),eta1(k),eta2(k),lambda(k),lambda_m(k),lambda_h(k)] = ...
        seir_simIsolation_magedFit(x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14) ...
        ,a1,a2,b1,b2,Omega1,Omega2,Theta1,Theta2,zeta1,zeta2,psi1,psi2, ...
        sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,mu,gamma,epsilon,beta,eta_i,eta_s,eta_h,t0,k,numTerms);
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
V= x(:,12);
Vp= x(:,13);
N = 9.7e6;
figure(1);
clf;
subplot(3,3,1)
plot(tspan,S+Sm+Sh)
hold on
maskedEstimated = (Sm+Em+Im)./(S+E+I+Sm+Em+Im+Sh+Eh+Ih+R+D+V);
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
hold on
plot(magedExposed)
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
plot(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+V+D)
hold on

title('Total Population')


subplot(3,3,9)
plot(phi1);
hold on
plot(phi2);
plot(eta1);
plot(eta2);
plot(lambda);
plot(lambda_m);
plot(lambda_h);
legend('phi1','phi2','eta1','eta2','lambda','lambda_m','lambda_h')



end