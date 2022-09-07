
clearvars

    
nc=13; np=19;
n=nc+np;%number of state
m = 5; %number of measurements
dt = 1;
R=diag([1000,100,1000,1,10]);        % covariance of measurement,%[infectious,death,vax,mask,mobility]


%-- fmnincon optimal parameters
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

Q=diag([1e-6*ones(1,nc), 0.1*parameterInit]); % covariance of process
f=@(x) seirDynamics(x,dt);  % nonlinear state equations
h=@(x) seirObservation(x);                               % measurement equation
s=[zeros(1,nc),parameterInit]';  %
s(1) = 9.7e6-24;% initial state
s(7) = 24;
x=s; %initial state          % initial state with noise
P = Q;                        % initial state covraiance
N=1000;  %?\ check                                   % total dynamic steps
xV = zeros(n,N);          %estmate        % allocate memory
sV = zeros(n,N);          %actual
zV = zeros(m,N);

% readData
maxT = N; %maximum is 1000
infectious = csvread('infectiousIllinois.csv');
infectious(isnan(infectious))=0;
infectious = infectious*(9.7/12.8);

death = csvread('deathIllinois.csv');
death(isnan(death))=0;
death = cumsum(death);
death = death*(9.7/12.8);


vax = csvread('vaccinatedIllinois.csv');
vax(isnan(vax))=0;
vax = vax*(9.7/12.8);
mask= csvread('maskIllinois.csv');

mobility = csvread('mobilityIllinois.csv');
mobility(isnan(mobility))=0;

death = death(1:maxT);
vax = vax(1:maxT);
infectious =  infectious(1:maxT);
mask =  mask(1:maxT);
mobility =  mobility(1:maxT);

z = [infectious';death';vax';mask';mobility']; % measurments
for k=1:N
    zk=z(:,k);                            % save actual state
    zV(:,k)  = zk;                             % save measurment
    [x, P] = ukf(f,x,P,h,zk,Q,R);            % ekf
    xV(:,k) = x;                            % save estimate
end

% remove outliers
for jj=1:n 
    mu=mean(xV(jj,:));
    st=std(xV(jj,:));
    
    idx=xV(jj,:)>mu+2*st|xV(jj,:)<mu-2*st;
    xV(jj,idx)=nan;
end

 % plot results
figure(1); gcf; clf;

subplot(3,3,1)

% plot(-(xV(3,:)+xV(6,:)+xV(9,:))./sum(xV(1:12,:),1),'--')
plot((xV(1,:)+xV(2,:)+xV(3,:)),'--')
title('Susceptible')

subplot(3,3,2)

% plot(-(xV(3,:)+xV(6,:)+xV(9,:))./sum(xV(1:12,:),1),'--')
plot((xV(4,:)+xV(5,:)+xV(6,:)),'--')
title('Exposed')

subplot(3,3,3)
plot(infectious,'-');
hold on
plot(xV(7,:)+xV(8,:)+xV(9,:),'--')
title('Infectious')


subplot(3,3,4)
plot(xV(10,:),'--')
title('Recovered')

subplot(3,3,5)
% plot(zV(2,:),'-');
plot(death,'-');

hold on
plot(xV(11,:),'--')
title('Deaths')



subplot(3,3,6)
plot(vax,'-');
hold on
plot(xV(13,:),'--')
title('Vaccinated')


subplot(3,3,7)
plot(mask,'-');
hold on
% plot((xV(2,:)+xV(5,:)+xV(8,:))./sum(xV(1:12,:),1) ,'--')
plot((xV(2,:)+xV(5,:)+xV(8,:))./sum(xV(1:12,:),1) ,'--')
title('Masked')

subplot(3,3,8)
plot(mobility,'-');
hold on
% plot(-(xV(3,:)+xV(6,:)+xV(9,:))./sum(xV(1:12,:),1),'--')
plot(-(xV(3,:)+xV(6,:)+xV(9,:))./sum(xV(1:12,:),1),'--')
title('Mobility')

%use fmincon optimal parameters
%incorporate euler maryama inside the dynamics
%




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
    
    S = S + (-(xi_Sh + alpha + phi_Sm + lambda)*S + kappa_R * R + xi_S *Sh + sigma_S * U + phi_S*Sm)*Dt;
    Sm =Sm+ (-(lambda_m + phi_S + alpha)*Sm + sigma_Sm * U + kappa_Rm * R + phi_Sm * S)*Dt;
    Sh =Sh+ (-(alpha + xi_S + lambda_h)*Sh + sigma_Sh * U + xi_Sh * S + kappa_Rh * R)*Dt;
    
    E =E+ (-(epsilon + phi_Em + xi_Eh)*E + xi_E * Eh + phi_E * Em + lambda*S )*Dt;
    Em = Em+(-(epsilon + phi_E) * Em + phi_Em*E + lambda_m* Sm)*Dt;
    Eh =Eh+ (-(epsilon + xi_E)*Eh + lambda_h * Sh + xi_Eh * E)*Dt;
    
    I =I+ (-(mu + gamma + xi_Ih + phi_Im)*I + epsilon*E + phi_I* Im + xi_I *Ih)*Dt;
    Im =Im+ (-(phi_I + gamma + mu)*Im + phi_Im*I + epsilon*Em)*Dt;
    Ih =Ih+ (-(xi_I + gamma + mu)*Ih + epsilon*Eh + xi_Ih*I)*Dt;
    
    R = R+ (-(kappa_R + kappa_Rm + kappa_Rh)*R+(I + Im + Ih)*gamma)*Dt;
    D = D+((I + Im + Ih)*mu)*Dt;
    U = U+(-(sigma_S + sigma_Sm + sigma_Sh)*U+(S + Sh + Sm)*alpha)*Dt;
end

V = V+ alpha*(S + Sm+ Sh);

x_kp1(1) = S;
x_kp1(2)=Sm;
x_kp1(3)=Sh;
x_kp1(4)=E ;
x_kp1(5)=Em;
x_kp1(6)=Eh;
x_kp1(7)=I ;
x_kp1(8)=Im;
x_kp1(9)=Ih;
x_kp1(10)=R ;
x_kp1(11)=D ;
x_kp1(12)=U ;
x_kp1(13)=V ;

maxvals=[1e7, 1e7/2, 1e7/2, ...
         1e7/2, 1e7/4, 1e7/4, ...
         1e7/2, 1e7/4, 1e7/4, ...
         1e7, 5e4, 1e5, 1e5];

for jj=1:13
    x_kp1(jj,x_kp1(jj,:)<0)=0;
    x_kp1(jj,x_kp1(jj,:)>maxvals(jj))=maxvals(jj);
end


for jj=14:32
    x_kp1(jj,x_kp1(jj,:)<0)=0;
    x_kp1(jj,x_kp1(jj,:)>1)=1;
end

end

function z = seirObservation(xk)
%[infectious,death,vax,mask,mobility]
% x_kp1(1) =S;
% x_kp1(2)=Sm;
% x_kp1(3)=Sh;
% x_kp1(4)=E ;
% x_kp1(5)=Em;
% x_kp1(6)=Eh;
% x_kp1(7)=I ;
% x_kp1(8)=Im;
% x_kp1(9)=Ih;
% x_kp1(10)=R ;
% x_kp1(11)=D ;
% x_kp1(12)=U ;

z(1) = xk(7)+xk(8)+xk(9); %Infectious
z(2) = xk(11); %death
z(3) = xk(13); % Vax
z(4) = (xk(2)+xk(5)+xk(8))/sum(sum(xk(1:12))); % Mask
z(5) = -(xk(3)+xk(6)+xk(9))./(sum(xk(1:12))); % Mobility

end