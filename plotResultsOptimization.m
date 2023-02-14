function plotResultsOptimization(x, maxT)
part = 1:maxT;
infectious = csvread('data/infectiousIllinois_ci.csv');
infectious=infectious(:,2);
infectious(isnan(infectious))=0;
infectious = infectious*(9.7/12.8);

death = csvread('data/deathIllinois.csv');
% death = death(:,2);
death(isnan(death))=0;
death = cumsum(death);
death = death*(9.7/12.8);



vax = csvread('data/vaccinatedIllinois.csv');
% vax = vax(:,2);
vax(isnan(vax))=0;
vax = vax*(9.7/12.8);

mask= csvread("data/maskIllinois.csv");
% mask= mask(:,2);

mobility = csvread("data/mobilityIllinois.csv");
mobility(isnan(mobility))=0;



death = death(part);
vax = vax(part);
infectious =  infectious(part);
mask =  mask(part);
mobility =  mobility(part);



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


close all

figure(1)

subplot(3,3,1)
plot(S+Sm+Sh)
hold on
maskedEstimated = (Sm+Em+Im)./(S+E+I+Sm+Em+Im+Sh+Eh+Ih+R+U);
maskedReference = mask;
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U);

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
plot(E+Em+Eh)
title('Exposed')

subplot(3,3,5)
plot(maskedEstimated)
hold on
plot(maskedReference)
title('Masked Estimated')

subplot(3,3,6)
plot(estimatedMobility)
hold on
plot(mobility/100)
title('Mobility')

subplot(3,3,7)
plot(R)
hold on

title('Recovered')


% subplot(4,3,9)
% plot(phi1);
% hold on
% plot(phi2);
% plot(xi1);
% plot(xi2);
% 
% legend('\phi_1','\phi_2','\xi_1','\xi_2')
% 
% subplot(4,3,10)
% plot(alpha1);
% hold on
% plot(beta1);
% legend('\alpha','\beta')
% 
% subplot(4,3,11)
% plot(eta_Ih1);
% hold on
% plot(eta_Im1);
% plot(eta_Sh1);
% plot(eta_Sm1);
% legend('\eta_{I_h}','\eta_{I_m}','\eta_{S_h}','\eta_{S_m}')
% 
% subplot(4,3,12)
% plot(kappa_R1);
% hold on
% plot(kappa_Rm1);
% plot(kappa_Rh1);
% legend('\kappa_{R}','\kappa_{R_m}','\kappa_{R_h}')

% 
% str = sprintf("w_D = %0.3f, w_I = %0.3f,w_U = %0.3f, w_H = %0.3f, w_M = %0.3f", ...
%     weights(1),weights(2),weights(3),weights(4),weights(5));
% 
% sgtitle(str);

end