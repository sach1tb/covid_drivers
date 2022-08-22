function plotResultsConcatenated(params1, params2, params3, weights,part1,part2,part3,objective)
close all
param = [params1;params2;params3];
parts=[part1(1), part1(end);part2(1), part2(end); part3(1), part3(end)];
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

phi1 = [];
phi2 = [];
eta1 = [];
eta2 = [];
lambda = [];
lambda_m = [];
lambda_h = [];
for run = 1:3
params = param(run,:);
part = parts(run,1):parts(run,2);


ddeath = death(part);
vvax = vax(part);
iinfectious =  infectious(part);
mmask =  mask(part);
mmobility =  mobility(part);


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



for k=parts(run,1):parts(run,2)

    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)...
        ,phi1(k),phi2(k),xi1(k),xi2(k),lambda(k),lambda_m(k),lambda_h(k),alpha1(k),beta1(k),eta_Ih1(k),eta_Im1(k),eta_Sh1(k) ...
        ,eta_Sm1(k),kappa_R1(k),kappa_Rm1(k),kappa_Rh1(k)]...
    = seir_simIsolation( ...
    x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14),phi_1,phi_2,xi_1,xi_2, ...
    sigma_S,sigma_Sm,sigma_Sh,c,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh,k );

end

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

subplot(4,3,1)
plot(S+Sm+Sh)
hold on
plot([parts(1,1) parts(1,1)],[0 5e6],'k-.')
plot([parts(1,2) parts(1,2)],[0 5e6],'k-.')
plot([parts(2,1) parts(2,1)],[0 5e6],'k-.')
plot([parts(2,2) parts(2,2)],[0 5e6],'k-.')
plot([parts(3,1) parts(3,1)],[0 5e6],'k-.')
plot([parts(3,2) parts(3,2)],[0 5e6],'k-.')
maskedEstimated = (Sm+Em+Im)./(S+E+I+Sm+Em+Im+Sh+Eh+Ih+R+D+U);
maskedReference = mask;
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R);

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
% str2 =  sprintf("\n \\sigma_{S} = %0.4f, \\sigma_{Sm} = %0.4f, \\sigma_{Sh} = %0.4f, \\alpha= %0.4f," + ...
%     "\n" + "\\beta = %0.4f, \\mu = %0.4f, \\epsilon = %0.4f, \\gamma = %0.4f, \\eta_{I_h}=%0.4f," + ...
%     "\\eta_{I_m} = %0.4f,\\eta_{S_h} = %0.4f,\\eta_{S_m} = %0.4f, \\kappa_R= %0.4f, \\kappa_R_m= %0.4f, \\kappa_R_h= %0.4f",sigma_S,sigma_Sm ...
%     ,sigma_Sh,alpha,beta,mu ,epsilon,gamma,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh);

sgtitle(str);

figure(2)
plot(objective,'LineWidth',1.3)

end