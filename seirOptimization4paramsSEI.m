clear all
%close all
%clc

global authorY1 authorY2 authorY3;
authorData = readmatrix('authorExposed.csv');
authorData2 = authorData(:,2);
authorData2(authorData2<0)=0;

authorData1 = authorData(:,1);
authorY1 = interp1(authorData1,authorData2,0:1:700);

authorData = readmatrix('authorSusceptible.csv');
authorData2 = authorData(:,2);


authorData1 = authorData(:,1);
authorY2 = interp1(authorData1,authorData2,0:1:700);



authorData = readmatrix('authorInfectious.csv');
authorData2 = authorData(:,2);

authorData2(authorData2<0)=0;
authorData1 = authorData(:,1);
authorY3 = interp1(authorData1,authorData2,0:1:700);



x0 = zeros(1,5);
options = optimoptions('fmincon','Display','final');
xopt = fmincon(@objective,x0,[],[],[],[],zeros(1,5),0.1*ones(1,5),[],options);
% ,[],[],[],[],[],[],[],[]
xopt
objective(xopt)
plotResults(xopt)


function plotResults(xx)
global authorY1 authorY2 authorY3;
tspan = 0:1:700;
dt = 1;

S0 = 1*(1e6-50);
E0 = 0;
I0 = 50;
Sm0 = 0.0*(1e6-50);
Em0 = 0;
Im0 = 0;
R0 = 0;
D0 = 0;
V0 = 0 ;

x = [];
x(1,1) = S0;
x(1,2) = E0;
x(1,3) =I0 ;
x(1,4) =Sm0;
x(1,5) =Em0;
x(1,6) =Im0;
x(1,7) =R0 ;
x(1,8) =D0 ;
x(1,9) =V0 ;
temp = [];
temp(1,1) = E0+Em0;
for k=1:numel(tspan)-1
    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), x(k+1,8), x(k+1,9),~,~,~,~] = ...
        seir_sim(x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), x(k,7), x(k,8), x(k,9), dt, xx(end), 100, xx(1:2), xx(3:4));
end

S = x(:,1);
E = x(:,2);
I = x(:,3);
Sm= x(:,4);
Em = x(:,5);
Im= x(:,6);
R = x(:,7);
D = x(:,8);
v = x(:,9);

figure(1)
clf;
subplot(2,2,1)
plot(tspan,S)
hold on
plot(tspan,authorY2,'r-.')
title('Susceptible')

subplot(2,2,2)
plot(tspan,E)
hold on
plot(tspan,authorY1,'r-.')
title('Exposed')

subplot(2,2,3)
plot(tspan,I)
hold on
plot(tspan,authorY3,'r-.')
title('Infectious')

subplot(2,2,4)
plot(tspan,R)
title('Recovered')




end


function obj = objective(xx)
global authorY1 authorY2 authorY3

tspan = 0:1:700;
dt = 1;

S0 = 1e6-50;
E0 = 0;
I0 = 50;
Sm0 = 0;
Em0 = 0;
Im0 = 0;
R0 = 0;
D0 = 0;
V0 = 0 ;

x = [];
x(1,1) = S0;
x(1,2) = E0;
x(1,3) =I0 ;
x(1,4) =Sm0;
x(1,5) =Em0;
x(1,6) =Im0;
x(1,7) =R0 ;
x(1,8) =D0 ;
x(1,9) =V0 ;
Expos = [];
Expos(1,1) = E0+Em0;

Sus = [];
Sus(1,1) = S0;

Infect = [];
Infect(1,1) = I0;
for k=1:numel(tspan)-1
    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), x(k+1,8), x(k+1,9),~,Sus(k+1),Expos(k+1),Infect(k+1)] = ...
        seir_sim(x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), x(k,7), x(k,8), x(k,9), dt, xx(end), 100, xx(1:2), xx(3:4));
end
obj = sum((Expos-authorY1).^2)^0.5; %Exposed
%obj = sum((Sus-authorY2).^2)^0.5; %Sus
% obj = sum((Infect-authorY3).^2)^0.5; %Infect
end



function [S,E,I,Sm,Em,Im,R,D,V,N,SS,EE,II] = seir_sim(S,E,I,Sm,Em,Im,R,D,V,dt,omega,tau,a,b)
eta_s = 0.35;
eta_i = 0.45;
sigma = 0.001;
beta = 0.312;
epsilon = 1/4.5;
gamma = 0.0602;
mu = 3.5/100; %death rate
alpha = 0.001*1;
% a = 0.3;
% b = 0.1;
%tau = 100000;

phi_S = (b(1)-a(1))/(1+exp(-10*(I+Im-tau)))+a(1);
phi_E = (b(1)-a(1))/(1+exp(-10*(I+Im-tau)))+a(1);
phi_I = (b(1)-a(1))/(1+exp(-10*(I+Im-tau)))+a(1);

phi_Sm = (b(2)-a(2))/(1+exp(-10*(I+Im-tau)))+a(2);

phi_Em = (b(2)-a(2))/(1+exp(-10*(I+Im-tau)))+a(2);
phi_Im = (b(2)-a(2))/(1+exp(-10*(I+Im-tau)))+a(2);



n = 100;
Dt = dt/n;
for i = 1:n

    N = S+E+I+Sm+Em+Im+R+D+V;
    lambda = (beta/N)*(I+(1-eta_i)*Im);
    lambda_m = (beta/N)*(1-eta_s)*(I+(1-eta_i)*Im);

    dS  = (-(phi_Sm+lambda)*S + phi_S*Sm - alpha*S)*Dt;
    dE  = (-(phi_Em + epsilon)*E + phi_E*Em + lambda*S + omega*sigma*V)*Dt;
    dI  = (-(phi_Im+gamma+mu)*I + phi_I*Im + epsilon*E)*Dt;
    dSm = (-(phi_S + lambda_m)*Sm + phi_Sm*S - alpha*Sm)*Dt;
    dEm = (-(phi_E+epsilon)*Em + phi_Em*E+lambda_m*Sm + omega*sigma*V)*Dt;
    dIm = (-(phi_I + gamma + mu)*Im + phi_Im*I + epsilon*Em)*Dt;
    dR  = (gamma*(I+Im))*Dt;
    dD  = (mu*(I+Im))*Dt;
    dV  = (alpha*(S+Sm)-sigma*V)*Dt;

    S = S + dS  ;
    E = E + dE  ;
    I = I + dI  ;
    Sm = Sm + dSm;
    Em = Em + dEm ;
    Im = Im + dIm;
    R = R + dR  ;
    D = D + dD  ;
    V = V + dV  ;
    SS = S+Sm;
    EE = E+Em;
    II = I+Im;
end


end

