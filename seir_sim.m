
function [S,E,I,Sm,Em,Im,R,D,V,Vp,N,SS,EE,II,phi_S,phi_Sm,lambda,lambda_m] ...
= seir_simIsolation( ...
   S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,V,Vp,a1,a2,b1,b2,Omega1,Omega2,Theta1,Theta2,zeta1,zeta2,psi1,psi2,sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,mu,gamma,epsilon,beta,eta_i,eta_s,eta_h,t0 )

phi_1 = a1(1)+a1(2)*cos(Omega1(1)*timeDay-Theta1(1));
phi_2 = a2(1)+a2(2)*cos(Omega2(1)*timeDay-Theta2(1));
phi_Sm = phi_1;
phi_Em = phi_1;
phi_Im = phi_1;
phi_S = phi_2;
phi_E = phi_2;
phi_I = phi_2;

eta_1 = b1(1)+b1(2)*cos(zeta1(1)*timeDay-psi1(1));
eta_2 = b2(1)+b2(2)*cos(zeta2(1)*timeDay-psi2(1));
eta_Sh = eta_1;
eta_Eh = eta_1;
eta_Ih = eta_1;
eta_Sm = eta_2;
eta_Em = eta_2;
eta_Im = eta_2;

k = 1e5;
alpha = c/(exp(-k*(t-t0)));

%% sigma_S

if(t<=t0)
    sigma_S = 0;
end
if(t>t0 && t<=t0+180)
    sigma_S = sigma_Smax*((t-t0)/180);
end
if(t>t0+180)
    sigma_S = sigma_Smax;
end

%% sigma_Sm

if(t<=t0)
    sigma_Sm = 0;
end
if(t>t0 && t<=t0+180)
    sigma_Sm = sigma_Sm_max*((t-t0)/180);
end
if(t>t0+180)
    sigma_Sm = sigma_Sm_max;
end

%% sigma_Sh

if(t<=t0)
    sigma_Sh = 0;
end
if(t>t0 && t<=t0+180)
    sigma_Sh = sigma_Sh_max*((t-t0)/180);
end
if(t>t0+180)
    sigma_Sh = sigma_Sh_max;
end

kappa_R = sigma_S;
kappa_Rm = sigma_Sm;
kappa_Rh = sigma_Sh;




n = 100;
Dt = dt/n;

for i = 1:n

    N = S+E+I+Sm+Em+Im+Sh+Eh+Ih+R+D+V;
    lambda = (beta/N)*(I+Ih+(1-eta_i)*Im);
    lambda_m = (beta/N)*(1-eta_s)*(I+Ih+(1-eta_i)*Im);
    lambda_h = (beta/N)*(1-eta_h)*(I+Ih+(1-eta_i)*Im);

    dS = (-(phi_Sm + lambda + alpha)*S + sigma_S*V + kappa_R * R + phi_S*S_m)*Dt;
    dSm = (-(phi_S + lambda_m + eta_Sh +  alpha)*S_m + sigma_Sm*V + eta_S*Sh+ kappa_Rm * R + phi_Sm*S)*Dt;
    dSh = (-(eta_Sm+lambda_h+alpha)*S_h + eta_Sh*S_m+sigma_Sh*V + kappa_Rh*R)*Dt;

    dE = -((epsilon+phi_Em)*E+lambda*S+phi_E*S)*Dt;
    dEm = -((eta_Eh+epsilon+phi_E)*E_m + eta_Em * E_h + lambda_m*S_m + phi_Em *E)*Dt;
    dEh = -((epsilon + eta_Em)*Eh + lambda_h*Sh+eta_Eh*Em)*Dt;
    
    dI = -((phi_Im+mu+gamma)*I+phi_I*Im+epsilon*E)*Dt;
    dIm = -((phi_I+gamma+mu+eta_Ih)*Im+epsilon*Em+phi_Im*I+eta_Im*Ih)*Dt;
    dIh = -((mu + gamma + eta_Im)*Ih + epsilon*Eh+eta_Ih*Im)*Dt;

    dR = -((kappa_R+kappa_Rm+kappa_Rh)*R+(I + Im + Ih)*gamma)*Dt;
    dD = ((I+Im+Ih)*mu)*Dt;
    dV = -((sigma_S+sigma_Sm+sigma_Sh)*V+(S+Sh+Sm)*alpha)*Dt;

    dVp = ((S+Sh+Sm)*alpha)*Dt;

    S = S + dS  ;
    Sm = Sm + dSm;
    Sh = Sh + dSh;

    E = E + dE  ;
    Em = Em + dEm  ;
    Eh = Eh + dEh  ;

    I = I + dI  ;
    Im = Im + dIm  ;
    Ih = Ih+ dIh  ;


    
    R = R + dR  ;
    D = D + dD  ;
    V = V + dV  ;

    Vp = Vp + dVp  ;


    SS = S+Sm+Sh;
    EE = E+Em+Eh;
    II = I+Im+Ih;
end
end