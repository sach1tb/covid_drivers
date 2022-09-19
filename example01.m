function example01
% this example compares the performance of different filtering algorithms
% on a dataset obtained by simulating a lorenz attractor. A high sampling
% rate is assumed and therefore a constant velocity motion model is assumed
% The measurement model is also linear 

% NOTE
% >> in the script means a big change is needed for you to adapt it to your
% setup. e.g. a new likelihood function or measurement model
% ^^ means there are parameters that can be tuned to see if it changes the
% performance


addpath /Users/sachit/Dropbox/EASeL/common/MYLIB/filtering


infectious = csvread('infectiousIllinois.csv');
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

Z=[infectious, death, vax, mask, mobility]';


% Z=Z(:,1:300);

% ready to filter
% initialize
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

X0=zeros(31,1);
X0(1)=9.7e6-24;% initial state
X0(7)=24;
X0(13:end)=parameterInit;


% uncomment one of these below to see how it performs, see the function
% itself to see how it works


n=size(X0,1);
T=size(Z,2);



% ^^ number of particles (change this to see how it affects the
% performance)
N=1000; 

% particles n x N x T
p=zeros(n,N,T);

% ^^ pepper the initial estimate with some noise 
for jj=1:N
    p(:,jj,1)=X0 + randn(n,1).*X0/10; 
end

% >> the measurement model is now replaced by a likelihood function
% ^^ noise values these should be changed depending on the measurement
% model above
% Z=[infectious, death, vax, mask, mobility]';
eta=[100,10,10000,10,10];

hfun=@(p) seirObservation(p);
glfn=@(Z,p) glfn1(Z,p, hfun, eta);

% we still need these to show stuff
Xh=zeros(n,T);
P=zeros(n,n,T);


% >> motion model
dt=1;
gmot=@(x) seirDynamics(x, dt);

% ^^disturbance noise
w=diag([ones(1,12), parameterInit/10]);

k0=1;
kF=T;

% particle filter
wts=ones(N,T);
for k=k0:kF

    % update
    [p(:,:,k), wts(:,k)] = pf_update(p(:,:,k), wts(:,k), Z(:,k), glfn);

    % this function pulls out the estimate from the distribution
    % ^^ the flag can be 1,2, or 3 and can give different estimates
    flag=2;
    Xh(:,k)=postest(p(:,:,k), wts(:,k), flag);
    P(:,:,k)=cov(p(:,:,k)');
    
    % predict
    p(:,:,k+1) = pf_predict(p(:,:,k), gmot, w);

    fprintf('.');
    if ~mod(k,20), fprintf('\n'); end
end

figure(1); gcf; clf;

subplot(3,3,1)
plot((Xh(1,:)+Xh(2,:)+Xh(3,:)),'--')
title('Susceptible')

subplot(3,3,2)
plot((Xh(4,:)+Xh(5,:)+Xh(6,:)),'--')
title('Exposed')

subplot(3,3,3)
plot(infectious,'-');
hold on
plot(Xh(7,:)+Xh(8,:)+Xh(9,:),'--')
title('Infectious')


subplot(3,3,4)
plot(Xh(10,:),'--')
title('Recovered')

subplot(3,3,5)
% plot(zV(2,:),'-');
plot(death,'-');

hold on
plot(Xh(11,:),'--')
title('Deaths')



subplot(3,3,6)
plot(vax,'-');
hold on
plot(Xh(12,:),'--')
title('Vaccinated')


subplot(3,3,7)
plot(mask,'-');
hold on
% plot((Xh(2,:)+Xh(5,:)+Xh(8,:))./sum(Xh(1:12,:),1) ,'--')
plot((Xh(2,:)+Xh(5,:)+Xh(8,:))./sum(Xh(1:12,:),1) ,'--')
title('Masked')

subplot(3,3,8)
plot(mobility,'-');
hold on
% plot(-(Xh(3,:)+Xh(6,:)+Xh(9,:))./sum(Xh(1:12,:),1),'--')
plot(-(Xh(3,:)+Xh(6,:)+Xh(9,:))./sum(Xh(1:12,:),1),'--')
title('Mobility')

subplot(3,3,9);
plot(Xh(13,:));
title('\beta');


function wts=glfn1(Z, p, hfun, eta)

% >> convert from actual value to Z

% >> this line for example would instead consist of the full nonlinear
% measurment model like the epipolar model or the camera model
Zh=hfun(p);


% eta is actually a diagonal matrix
eta=diag(eta);
wts=normpdf(Zh(1), Z(1), eta(1,1)).*...
    normpdf(Zh(2), Z(2), eta(2,2)).*...
    normpdf(Zh(3), Z(3), eta(3,3)).*...
    normpdf(Zh(4), Z(4), eta(4,4)).*...
    normpdf(Zh(5), Z(5), eta(5,5)).*...
    unifpdf(p(1), 0, 1e12).*...
    unifpdf(p(2), 0, 1e12).*...
    unifpdf(p(3), 0, 1e12).*...
    unifpdf(p(4), 0, 1e12).*...
    unifpdf(p(5), 0, 1e12).*...
    unifpdf(p(6), 0, 1e12).*...
    unifpdf(p(7), 0, 1e12).*...
    unifpdf(p(8), 0, 1e12).*...
    unifpdf(p(9), 0, 1e12).*...
    unifpdf(p(10), 0, 1e12).*...
    unifpdf(p(11), 0, 1e12).*...
    unifpdf(p(12), 0, 1e12).*...
    unifpdf(p(13), 0.2, 1).*...
    unifpdf(p(14), 0, 1).*...
    unifpdf(p(15), 0, 1).*...
    unifpdf(p(16), 0, 1).*...
    unifpdf(p(17), 0, 1).*...
    unifpdf(p(18), 0, 1).*...
    unifpdf(p(19), 0, 1).*...
    unifpdf(p(20), 0, 1).*...
    unifpdf(p(21), 0, 1).*...
    unifpdf(p(22), 0, 1).*...
    unifpdf(p(23), 0, 1).*...
    unifpdf(p(24), 0, 1).*...
    unifpdf(p(25), 0, 1).*...
    unifpdf(p(26), 0, 1).*...
    unifpdf(p(27), 0, 1).*...
    unifpdf(p(28), 0, 1).*...
    unifpdf(p(29), 0, 1).*...
    unifpdf(p(30), 0, 1).*...
    unifpdf(p(31), 0, 1);



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
z(3) = xk(12); % Vax
z(4) = (xk(2)+xk(5)+xk(8))/sum(sum(xk(1:12))); % Mask
z(5) = -(xk(3)+xk(6)+xk(9))./(sum(xk(1:12))); % Mobility

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

beta = xk(13);
eta_Ih = xk(14);
eta_Im = xk(15);
eta_Sm = xk(16);
eta_Sh = xk(17);

xi_Sh = xk(18);
xi_Eh = xk(18);
xi_Ih = xk(18);

xi_S = xk(19);
xi_E = xk(19);
xi_I = xk(19);

alpha = xk(20);

phi_Sm = xk(21);
phi_Em = xk(21);
phi_Im = xk(21);

phi_S = xk(22);
phi_E = xk(22);
phi_I = xk(22);

sigma_S = xk(23);
sigma_Sm = xk(24);
sigma_Sh = xk(25);
kappa_R= xk(26);
kappa_Rm= xk(27);
kappa_Rh= xk(28);
mu = xk(29);
gamma = xk(30);
epsilon = xk(31);
n = 1;
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
