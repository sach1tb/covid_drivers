function plotResultsOptimization(params, weights)
part = 1:335;
infectious = csvread('data/infectiousIllinois.csv');
infectious(isnan(infectious))=0;
infectious = infectious*(9.7/12.8);

death = csvread('data/deathIllinois.csv');
death(isnan(death))=0;
death = cumsum(death);
death = death*(9.7/12.8);



vax = csvread('data/vaccinatedIllinois.csv');
vax(isnan(vax))=0;
vax = vax*(9.7/12.8);
mask= csvread("data/maskIllinois.csv");

mobility = csvread("data/mobilityIllinois.csv");
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


sigma_ShVector = params(1);
sigma_SVector= params(2);
sigma_SmVector= params(3);
xi1Vector= params(4);
xi2Vector= params(5);
phi1Vector = params(6);
phi2Vector = params(7);
gammaVector = params(8);
muVector = params(9);
kappa_RhVector= params(10);
kappa_RVector= params(11);
kappa_RmVector= params(12);
alphaVector= params(13);
epsilonVector= params(14);
betaVector= params(15);

eta_Ih = params(16);
eta_Im = params(17);
eta_Sh = params(18);
eta_Sm = params(19);


for k=1:numel(part)-1

    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)...
        ,phi1(k),phi2(k),xi1(k),xi2(k),lambda(k),lambda_m(k),lambda_h(k),alpha1(k),...
            beta1(k),eta_Ih1(k),eta_Im1(k),eta_Sh1(k) ...
        ,eta_Sm1(k),kappa0(k)]...
        = seirDynamicsforOptimization(x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14), ...
        sigma0, xi1Vector, xi2Vector ...
        , phi1Vector, phi2Vector, ...
        gammaVector, muVector,kappa0, ...
        alphaVector ...
        , epsilonVector, betaVector,eta_Ih,eta_Im,eta_Sh,eta_Sm );

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

close all

figure(1)

subplot(4,3,1)
plot(S+Sm+Sh)
hold on
maskedEstimated = (Sm+Em+Im)./(S+E+I+Sm+Em+Im+Sh+Eh+Ih+R+U);
maskedReference = mask;
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U);

title('Susceptible')


subplot(4,3,2)
plot(D)
hold on
plot(death)

legend('Estimated','Real Data')
title('Deaths')


subplot(4,3,3)
plot(I+Im+Ih)
hold on
plot(infectious)
title('Infectious')


subplot(4,3,4)
plot(Vp)
hold on
plot(vax)
title('Vaccinated Perpetual')

subplot(4,3,5)
plot(E+Em+Eh)
title('Exposed')

subplot(4,3,6)
plot(maskedEstimated)
hold on
plot(maskedReference)
title('Masked Estimated')

subplot(4,3,7)
plot(estimatedMobility)
hold on
plot(mobility/100)
title('Mobility')

subplot(4,3,8)
plot(R)
hold on

title('Recovered')


subplot(4,3,9)
plot(phi1);
hold on
plot(phi2);
plot(xi1);
plot(xi2);

legend('\phi_1','\phi_2','\xi_1','\xi_2')

subplot(4,3,10)
plot(alpha1);
hold on
plot(beta1);
legend('\alpha','\beta')

subplot(4,3,11)
plot(eta_Ih1);
hold on
plot(eta_Im1);
plot(eta_Sh1);
plot(eta_Sm1);
legend('\eta_{I_h}','\eta_{I_m}','\eta_{S_h}','\eta_{S_m}')

subplot(4,3,12)
plot(kappa_R1);
hold on
plot(kappa_Rm1);
plot(kappa_Rh1);
legend('\kappa_{R}','\kappa_{R_m}','\kappa_{R_h}')


% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1,𝑤_𝑀=0.1
str = sprintf("w_D = %0.3f, w_I = %0.3f,w_U = %0.3f, w_H = %0.3f, w_M = %0.3f", ...
    weights(1),weights(2),weights(3),weights(4),weights(5));

sgtitle(str);

end