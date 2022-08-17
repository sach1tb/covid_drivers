clear all
%close all
%clc
global infectious death vax mask maxT numTerms mobility magedExposed weights
maxT = 1000; %maximum is 1000
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


% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1, 𝑤_𝑀=0.1 
% weights = [0.9, 0.0, 0.1, 0.0, 0.0];
weights = [0.0, 0.0, 0.0, 0.0, 1.0];
numTerms = 1;
xmin=[0.00001*ones(1,4*numTerms+8*(numTerms-1)) 0.001 0.001 0.001 0.001  300  0.2  0.001  0.1    0.01   0.01   0.01  0.01  0.01];
xmax = [1.0*ones(1,4*numTerms+8*(numTerms-1))   1     1     1     1      1000  0.4  0.01  0.3    0.3    1.0   1.0  1.0  1.0];
% a1
% a2
% b1
% b2
% Omega1
% Omega2
% Theta1
% Theta2
% zeta1
% zeta2
% psi1 
% psi2 
% % sigma_Smax
% % sigma_Sm_max
% % sigma_Sh_max
% % c
% t0 
% beta
% mu 
% epsilon
% gamma
% eta_Ih
% eta_Im
% eta_Sh
% eta_Sm
x0 = 0.01*(xmin(1:numel(xmin))+xmax(1:numel(xmin)));
options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',50000,'Algorithm','interior-point');
xopt = fmincon(@objective,x0,[],[],[],[],xmin,xmax,[],options);

xopt
objective(xopt)
plotResults(xopt,numTerms,maxT, weights)

function obj = objective(params)
global infectious death mask vax maxT numTerms mobility  weights

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
%%
% sigma_Smax = params(static+1);
% sigma_Sm_max = params(static+2);
% sigma_Sh_max = params(static+3);
% c = params(static+4);
% eta_i= params(static+5);
% 
% t0 = params(static+6);
% beta =  params(static+7);
% mu = params(static+8);
% epsilon = params(static+9);
% kFactor = 5;
% gamma = params(static+10);
% eta_s = params(static+11); %0.35;
% eta_h = params(static+12); %0.45;
% %kappa = 0*ones(1,numel(tspan)););

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
for k=1:numel(tspan)-1

    % disp(k)
    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)] = ...
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

obj1 = sum((D/N-death/N).^2)^0.5; %w_D
obj3 = sum((Vp/N-vax/N).^2)^0.5;  %w_U
obj2 = sum(((I+Im+Ih)/N-infectious/N).^2)^0.5;  %w_I
maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
maskedReference = mask;
obj5 = sum((maskedEstimated-maskedReference).^2)^0.5; % w_M
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
obj4 = sum(((estimatedMobility-mobility/100)).^2)^0.5; % w_H
% 𝑤_𝐷=0.35, 𝑤_𝐼=0.1, 𝑤_𝑈=0.35, 𝑤_𝐻=0.1, 𝑤_𝑀=0.1 
obj = weights(1)*obj1+weights(2)*obj2+weights(3)*obj3+weights(4)*obj4+weights(5)*obj5;
%obj = sum((E+Em+Eh - magedExposed').^2)^0.5;
%obj = obj1*obj3;
end

function [c, ceq] = constraint(params)
global infectious death mask vax maxT numTerms mobility exposed
tspan = 0:1:maxT-1;
dt = 1;


S0 = 9.7e6-240;
Sm0 = 0;
Sh0 = 0;
E0 = 0;
Em0 = 0;
Eh0 = 0;
I0 = 240;
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
sigma_Smax = params(static+1);
sigma_Sm_max= params(static+2);
sigma_Sh_max= params(static+3);
c= params(static+4);


gamma = 0.0602;


eta_i= params(static+5);
eta_s = 0.35;
eta_h = 0.45;
t0 = params(static+6);
beta =  params(static+7);
mu = params(static+8);
epsilon = params(static+9);
kFactor =  params(static+10);
for k=1:numel(tspan)-1

    %     [S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,V,Vp,N] ... = seir_simIsolation( ...
    %    S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,V,Vp,N,a1,a2,b1,b2, ...
    %    Omega1,Omega2,Theta1,Theta2,zeta1,zeta2,psi1,psi2, ...
    %    sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,mu,gamma,epsilon,beta,eta_i,eta_s,eta_h,t0
    %    )
    % disp(k)
    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14), ...
        phi1(k),phi2(k),eta1(k),eta2(k),lambda(k),lambda_m(k),lambda_h(k)] = ...
        seir_simIsolation(x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14) ...
        ,a1,a2,b1,b2,Omega1,Omega2,Theta1,Theta2,zeta1,zeta2,psi1,psi2, ...
        sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,kFactor,mu,gamma,epsilon,beta,eta_i,eta_s,eta_h,t0,k,numTerms);
end

c(1)=abs(min(phi1));
c(2)=abs(min(phi2));
c(3)=abs(min(eta1));
c(4)=abs(min(eta2));


c(5) = abs(max(phi1)) -1;
c(6) = abs(max(phi2)) - 1;
c(7) = abs(max(eta1)) - 1;
c(8) = abs(max(eta2)) - 1;
% c(9) = mean(lambda) - 1;
% c(10) = mean(lambda_m)-1;
% c(11) = mean(lambda_h)-1;

ceq = [];

end