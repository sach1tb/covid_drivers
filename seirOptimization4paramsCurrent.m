clear all
%close all
%clc
global infectious death vax mask maxT numTerms mobility magedExposed
maxT = 700; %maximum is 1000
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



numTerms = 3;
xmin=0.01*ones(1,12*numTerms+5);
xmax = [10.0*ones(1,12*numTerms) 1.0*ones(1,5)];
x0 = 0.01*(xmin+xmax);

options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',5000,'Algorithm','interior-point');
xopt = fmincon(@objective,x0,[],[],[],[],xmin,xmax,@constraint,options);
% ,[],[],[],[],[],[],[],[]
xopt
objective(xopt)
plotResults(xopt,numTerms,maxT)

function obj = objective(params)
global infectious death mask vax maxT numTerms mobility magedExposed

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

Omega1 = params(1+4*numTerms:5*numTerms);
Omega2 = params(1+5*numTerms:6*numTerms);
Theta1 = params(1+6*numTerms:7*numTerms);
Theta2 = params(1+7*numTerms:8*numTerms);
zeta1 = params(1+8*numTerms:9*numTerms);
zeta2 = params(1+9*numTerms:10*numTerms);
psi1 = params(1+10*numTerms:11*numTerms);
psi2 = params(1+11*numTerms:12*numTerms);
sigma_Smax = params(12*numTerms+1)*0;
sigma_Sm_max= params(12*numTerms+2)*0;
sigma_Sh_max= params(12*numTerms+3)*0;
c= params(12*numTerms+4);
mu = 3.5/100;
gamma = 0.0602;
epsilon = 1/4.5;
beta = 0.312;
eta_i= params(12*numTerms+5);
eta_s = 0.35;
eta_h = 0.45;
t0 = 24;
%kappa = 0*ones(1,numel(tspan));
for k=1:numel(tspan)-1

    %     [S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,V,Vp,N] ... = seir_simIsolation( ...
    %    S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,V,Vp,N,a1,a2,b1,b2, ...
    %    Omega1,Omega2,Theta1,Theta2,zeta1,zeta2,psi1,psi2, ...
    %    sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,mu,gamma,epsilon,beta,eta_i,eta_s,eta_h,t0
    %    )
    % disp(k)
    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)] = ...
        seir_simIsolation(x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
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
U= x(:,12);
Vp= x(:,13);
N = 9.7e6;

obj1 = sum((D-death).^2)^0.5;
obj2 = sum((Vp-vax).^2)^0.5;
obj3 = sum((I+Im+Ih-infectious).^2)^0.5;
maskedEstimated = (Sm+Em+Im)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
maskedReference = mask;
obj4 = sum((maskedEstimated*N-maskedReference*N).^2)^0.5;
estimatedMobility = -(Sh+Eh+Ih)./(S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U+D);
obj5 = sum((N*(estimatedMobility-mobility/100)).^2)^0.5;

obj = 0.35*obj1+0.35*obj2+0.1*obj3+0.1*obj4+0.1*obj5;
obj = sum((E+Em+Eh - magedExposed').^2)^0.5;
%obj = obj1*obj3;
end

function [c, ceq] = constraint(params)
global infectious death mask vax maxT numTerms mobility exposed
tspan = 0:1:maxT-1;
dt = 1;

S0 = 9.7e6-24;
Sm0 = 0;
Sh0 = 0;
E0 = 0;
Em0 = 0;
Eh0 = 0;
I0 = 24;
Im0 = 0;
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
Omega1 = params(1+4*numTerms:5*numTerms);
Omega2 = params(1+5*numTerms:6*numTerms);
Theta1 = params(1+6*numTerms:7*numTerms);
Theta2 = params(1+7*numTerms:8*numTerms);
zeta1 = params(1+8*numTerms:9*numTerms);
zeta2 = params(1+9*numTerms:10*numTerms);
psi1 = params(1+10*numTerms:11*numTerms);
psi2 = params(1+11*numTerms:12*numTerms);
sigma_Smax = params(12*numTerms+1);
sigma_Sm_max= params(12*numTerms+2);
sigma_Sh_max= params(12*numTerms+3);
c= params(12*numTerms+4);
mu = 3.5/100;
gamma = 0.0602;
epsilon = 1/4.5;
beta = 0.312;
eta_i= params(12*numTerms+5);
eta_s = 0.35;
eta_h = 0.45;
t0 = 24;
%kappa = 0*ones(1,numel(tspan));
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
        sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,mu,gamma,epsilon,beta,eta_i,eta_s,eta_h,t0,k,numTerms);
end

c(1)=abs(min(phi1));
c(2)=abs(min(phi2));
c(3)=abs(min(eta1));
c(4)=abs(min(eta2));


c(5) = abs(max(phi1)) -1;
c(6) = abs(max(phi2)) - 1;
c(7) = abs(max(eta1)) - 1;
c(8) = abs(max(eta2)) - 1;
c(9) = mean(lambda)-1;
c(10) = mean(lambda_m)-1;
c(11) = mean(lambda_h)-1;

ceq = [];



end

