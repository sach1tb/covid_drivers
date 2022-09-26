close all
clearvars
%fmincon optimised parameters
beta0 = 0.5;
mu0 = 2e-4;
epsilon0 = 0.8; %1.0 - optimised
gamma0 = 0.0110;
eta_Ih0 = 0.1;
eta_Im0 = 0.1;
eta_Sh0 = 0.1;
eta_Sm0 = 0.1;
kappa_R0 = 0.1;
kappa_Rm0 = 0.1;
kappa_Rh0 = 0.1;
sigma_S0 = 0.0001;
sigma_Sm0 = 0.0001;
sigma_Sh0 = 0.0001;
alpha0 = 0.0246;

phi_Sm0 = 0.01;
phi_S0 = 0.05;
xi_Sh0 = 0.001;
xi_S0 = 0.001;

parameterInit = [beta0,eta_Ih0,eta_Im0,eta_Sm0,eta_Sh0,xi_Sh0,xi_S0 ...
    alpha0,phi_Sm0,phi_S0,sigma_S0,sigma_Sm0,sigma_Sh0,kappa_R0,kappa_Rm0 ...
    kappa_Rh0, mu0, gamma0, epsilon0];
variableObserved = zeros(1,13);
variableObserved(1) = 9684738-24; %9.7e6-24; % initial state susceptible
variableObserved(7) = 24;
xV = [variableObserved,parameterInit]';
for i = 1:365

    [xV(:,i+1)] = seirDynamics(xV(:,i),1);

end

figure(1); gcf; clf;

subplot(3,3,1)


plot(1:size(xV,2),(xV(1,:)+xV(2,:)+xV(3,:)),'--')

title('Susceptible')

subplot(3,3,2)

plot((xV(4,:)+xV(5,:)+xV(6,:)),'--')
title('Exposed')

subplot(3,3,3)

plot(xV(7,:)+xV(8,:)+xV(9,:),'--','Color','r','LineWidth',1.4)
title('Infectious')


subplot(3,3,4)
plot(xV(10,:),'--')
title('Recovered')

subplot(3,3,5)
plot(xV(11,:),'--','Color','r')
title('Deaths')



subplot(3,3,6)
plot(xV(13,:),'--')

title('Vaccinated')


subplot(3,3,7)
plot((xV(2,:)+xV(5,:)+xV(8,:))./sum(xV(1:12,:),1) ,'--')
title('Masked')

subplot(3,3,8)

plot(-100*(xV(3,:)+xV(6,:)+xV(9,:))./sum(xV(1:12,:),1),'--')
title('Mobility')

subplot(3,3,9)
plot(totPopulation,'-');
hold on;
plot(sum(xV(1:12,:),1),'--')
title('Total population')




function [x_kp1] = seirDynamics(xk,dt)
x_kp1=xk;

S = xk(1);
Sm = xk(2);
Sh= xk(3);
E = xk(4);
Em = xk(5);
Eh = xk(6);
I = xk(7);
Im = xk(8);
Ih = xk(9);
R = xk(10);
D = xk(11);
U = xk(12);
V = xk(13);

beta = xk(14);
eta_Ih = xk(15);
eta_Im = xk(16);
eta_Sm = xk(17);
eta_Sh = xk(18);

xi_Sh = xk(19);
xi_Eh = xk(19);
xi_Ih = xk(19);

xi_S = xk(20);
xi_E = xk(20);
xi_I = xk(20);

alpha = xk(21);

phi_Sm = xk(22);
phi_Em = xk(22);
phi_Im = xk(22);

phi_S = xk(23);
phi_E = xk(23);
phi_I = xk(23);

sigma_S = xk(24);
sigma_Sm = xk(25);
sigma_Sh = xk(26);

kappa_R= xk(27);
kappa_Rm= xk(28);
kappa_Rh= xk(29);

mu = xk(30);
gamma = xk(31);
epsilon = xk(32);
n = 1000;
Dt = dt/n;


for i = 1:n
    N = S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+D+U;
    lambda = (beta/N)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);
    lambda_m = (beta/N)*(1-eta_Sm)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);
    lambda_h = (beta/N)*(1-eta_Sh)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);

    dS = (-(xi_Sh + alpha + phi_Sm + lambda)*S + kappa_R * R + xi_S *Sh + sigma_S * U + phi_S*Sm)*Dt;
    dSm =(-(lambda_m + phi_S + alpha)*Sm + sigma_Sm * U + kappa_Rm * R + phi_Sm * S)*Dt;
    dSh =(-(alpha + xi_S + lambda_h)*Sh + sigma_Sh * U + xi_Sh * S + kappa_Rh * R)*Dt;

    dE =(-(epsilon + phi_Em + xi_Eh)*E + xi_E * Eh + phi_E * Em + lambda*S )*Dt;
    dEm =(-(epsilon + phi_E) * Em + phi_Em*E + lambda_m* Sm)*Dt;
    dEh =(-(epsilon + xi_E)*Eh + lambda_h * Sh + xi_Eh * E)*Dt;

    dI =(-(mu + gamma + xi_Ih + phi_Im)*I + epsilon*E + phi_I* Im + xi_I *Ih)*Dt;
    dIm =(-(phi_I + gamma + mu)*Im + phi_Im*I + epsilon*Em)*Dt;
    dIh = (-(xi_I + gamma + mu)*Ih + epsilon*Eh + xi_Ih*I)*Dt;

    dR =  (-(kappa_R + kappa_Rm + kappa_Rh)*R+(I + Im + Ih)*gamma)*Dt;
    dD = ((I + Im + Ih)*mu)*Dt;
    dU = (-(sigma_S + sigma_Sm + sigma_Sh)*U+(S + Sh + Sm)*alpha)*Dt;

    totGrads = dS + dSm + dSh + dE + dEm + dEh + dI + dIm + dIh + dR + dD + dU;
    if abs(totGrads) > 1
        disp(totGrads);
        keyboard
    end

    dV =  alpha*(S + Sm+ Sh)*Dt;

    S = S+dS;
    Sm = Sm + dSm;
    Sh = Sh + dSh;

    E = E + dE;
    Em = Em + dEm;
    Eh = Eh + dEh;

    I = I + dI;
    Im = I + dIm;
    Ih = Ih + dIh;

    R = R + dR;
    D = D + dD;
    U = U + dU;

    V = V + dV;

end
end