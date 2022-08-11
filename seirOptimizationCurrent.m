clear all
close all
%clc
global infectious death vax mask maxT
maxT = 1000; %maximum is 1000
infectious = readmatrix('infectious.csv');
infectious(isnan(infectious))=0;

death = readmatrix('death.csv');
death(isnan(death))=0;
death = cumsum(death);

vax = readmatrix('vaccinated.csv');
vax(isnan(vax))=0;

mask= readmatrix("mask.csv");

death = death(1:maxT);
vax = death(1:maxT);
infectious =  infectious(1:maxT);
mask =  mask(1:maxT);


xmin=[0 0 0 0 0 0 0.001*ones(1,6) 0.0001*ones(1,6) 0 0 0 0.000001 0 0];
xmax = [0.3*ones(1,6) 0.005*ones(1,6) 0.01*ones(1,6) 0.2 0.2 0.2 0.5 0.5 0.5];
x0 = 0.1*(xmin+xmax);

options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',50000);
xopt = fmincon(@objective,x0,[],[],[],[],xmin,xmax,[],options);
% ,[],[],[],[],[],[],[],[]
xopt
objective(xopt)
plotResults(xopt)


function plotResults(xx)
global infectious death vax mask maxT
frac = [0];

figure(1)
clf;
for ii = 1:numel(frac)

    tspan = 0:1:maxT-1;
    dt = 1;


    S0 = 1*(330e6-100);
    E0 = 0;
    I0 = 1*100;
    Sm0 = 0*(330e6-100);
    Em0 = 0;
    Im0 = 0*50;

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
    x(1,10)=V0 ;
    x(1,11)=I0+Im0;
    temp = [];
    temp(1,1) = E0+Em0;

    phi_SMult = 1;
    phi_EMult = 1;
    phi_IMult = 1;
    phi_SmMult = 1;
    phi_EmMult = 1;
    phi_ImMult = 1;

    eta_s = 0.35;
    eta_i = 0.45;

    beta = 0.312;
    % epsilon = 1/4.5;
    epsilon = xx(17);

    a = xx(1:6);
f = xx(7:12);
freq_i = xx(13:18);
sigma_Vm_max = xx(19);
sigma_V_max = xx(20);
c = xx(21);
mu = xx(22);
epsilon = xx(23);
gamma = xx(24);
    tau = 24;
    %kappa = 0*ones(1,numel(tspan));
    for k=1:numel(tspan)-1

        noise = 0;
        timeDay = k;
        [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), x(k+1,8), x(k+1,9),x(k+1,10),x(k+1,11),~,~,~,~] ...
            = seir_sim(eta_s,eta_i, sigma_Vm_max,sigma_V_max,c,beta,epsilon,gamma, mu,a,f,freq_i,timeDay,tau, ...
            x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), x(k,7), x(k,8), x(k,9),x(k,10),x(k,11), ...
            dt,noise, phi_SMult,phi_EMult,phi_IMult,phi_SmMult,phi_EmMult,phi_ImMult);
    end
    S= x(:,1);
    E= x(:,2);
    I= x(:,3);
    Sm= x(:,4);
    Em= x(:,5);
    Im= x(:,6);
    R= x(:,7);
    D = x(:,8);
    V= x(:,9);
    Vp= x(:,10);
    Ip= x(:,11);
    maskedEstimated = (Sm+Em+Im)./(S+E+I+Sm+Em+Im);
    maskedReference = mask;
    Vcumulative = cumsum(V);
    subplot(3,2,1)


    plot(tspan,S+Sm)
    hold on
    %plot(tspan,S)
    %plot(tspan,Sm)
    title('Susceptible')


    subplot(3,2,2)
    plot(tspan,D)
    hold on
    plot(tspan,death)
    title('Deaths')


    subplot(3,2,3)
    plot(tspan,I+Im)
    hold on
    plot(infectious)
    title('Infectious')


    subplot(3,2,4)
    plot(tspan,Vp)
    hold on
    plot(vax)
    title('Vaccinated Perpetual')

    subplot(3,2,5)
    plot(tspan,E+Em)
    hold on
    title('Exposed')

    subplot(3,2,6)
    plot(maskedEstimated)
    hold on
    plot(maskedReference)
    title('Masked Estimated')

end





end


function obj = objective(xx)
global infectious death mask vax maxT;

tspan = 0:1:maxT-1;
dt = 1;

S0 = 1*(330e6-100);
E0 = 0;
I0 = 1*100;
Sm0 = 0*(330e6-100);
Em0 = 0;
Im0 = 0*50;

R0 = 0;
D0 = 0;
V0 = 0 ;

x = [];
x(1,1) = S0;
x(1,2) = E0;
x(1,3) = I0 ;
x(1,4) = Sm0;
x(1,5) = Em0;
x(1,6) = Im0;
x(1,7) = R0 ;
x(1,8) = D0 ;
x(1,9) = V0 ;
x(1,10) =V0 ;
x(1,11) =I0+Im0 ;

phi_SMult = 1;
phi_EMult = 1;
phi_IMult = 1;
phi_SmMult = 1;
phi_EmMult = 1;
phi_ImMult = 1;

eta_s = 0.35;
eta_i = 0.45;
beta = 0.312;
%epsilon = 1/4.5;

% gamma = 0.0602;

 %death rateS
a = xx(1:6);
f = xx(7:12);
freq_i = xx(13:18);
sigma_Vm_max = xx(19);
sigma_V_max = xx(20);
c = xx(21);
mu = xx(22);
epsilon = xx(23);
gamma = xx(24);
tau = 24;
%kappa = 0*ones(1,numel(tspan));
for k=1:numel(tspan)-1


    % [S,E,I,Sm,Em,Im,R,D,V,N,SS,EE,II,phi_S,phi_Sm,lambda,lambda_m] = seir_sim(eta_s,eta_i, ...
    %     sigma_Vm_max,sigma_V_max,c,beta,epsilon,gamma, ...
    %     mu,a,f,timeDay,tau, ...
    %     S,E,I,Sm,Em,Im,R,D,V,dt,noise, phi_SMult,phi_EMult,phi_IMult,phi_SmMult,phi_EmMult,phi_ImMult )

    noise = 0;
    timeDay = k;
    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), x(k+1,8), x(k+1,9),x(k+1,10),x(k+1,11),~,~,~,~] ...
        = seir_sim(eta_s,eta_i, sigma_Vm_max,sigma_V_max,c,beta,epsilon,gamma, mu,a,f,freq_i,timeDay,tau, ...
        x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), x(k,7), x(k,8), x(k,9), x(k,10), x(k,11), ...
        dt,noise, phi_SMult,phi_EMult,phi_IMult,phi_SmMult,phi_EmMult,phi_ImMult);
end
S= x(:,1);
E= x(:,2);
I= x(:,3);
Sm= x(:,4);
Em= x(:,5);
Im= x(:,6);
R= x(:,7);
D = x(:,8);
V= x(:,9);
Vp= x(:,10);
Ip= x(:,11);
N = 330e6;
Vcumulative = cumsum(V);
% obj = 0.5*sum(((D./N)-(death./N)).^2)^0.5 + 0.4*sum(((V./N)-(vax./N)).^2)^0.5 + 0.1*sum((((I+Im)./N)-(infectious./N)).^2)^0.5;
obj1 = sum((D-death).^2)^0.5;
obj2 = sum((Vp-vax).^2)^0.5;
obj3 = sum((I+Im-infectious).^2)^0.5;
maskedEstimated = (Sm+Em+Im)./(S+E+I+Sm+Em+Im);
maskedReference = mask;
obj4 = sum((maskedEstimated*N-maskedReference*N).^2)^0.5;
%obj = 0.4*obj1+0*obj2+0.35*obj3+0.25*obj4;
obj = 0.5*obj1+0.1*obj2+0.3*obj3+0.1*obj4;
% obj = 0.4*obj1+0.4*obj2+0.1*obj3+0.1*obj4;

% obj = 0.5*obj1+0.2*obj2+0.20*obj3+0.1*obj4;
end



