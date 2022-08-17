
function [S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,phi_1,phi_2,eta_1,eta_2,lambda,lambda_m,lambda_h,alpha] ...
    = seir_simIsolation( ...
    S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,a1,a2,b1,b2, ...
    Omega1,Omega2,Theta1,Theta2,zeta1,zeta2,psi1,psi2, ...
    sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,mu,gamma,epsilon,beta,eta_i,eta_s,eta_h,t0,timeDay,numTerms )
% a1 = a1./(sum(a1));
% a2 = a2./(sum(a2));
% b1 = b1./(sum(b1));
% b2 = b2./(sum(b2));
phi_1 = a1(1);
phi_2 = a2(1);
eta_1 = b1(1);
eta_2 = b2(1);
if(numTerms>1)

    for ii = 1:numTerms-1
        phi_1 = phi_1 +a1(ii+1)*cos(Omega1(ii)*timeDay-Theta1(ii));
        phi_2 = phi_2 +a2(ii+1)*cos(Omega2(ii)*timeDay-Theta2(ii));

        eta_1 = eta_1 +b1(ii+1)*cos(zeta1(ii)*timeDay-psi1(ii));
        eta_2 = eta_2 +b2(ii+1)*cos(zeta2(ii)*timeDay-psi2(ii));
    end
end


phi_Sm = phi_1;
phi_Em = phi_1;
phi_Im = phi_1;
phi_S = phi_2;
phi_E = phi_2;
phi_I = phi_2;

phi_S = 8.8864e-10;
phi_E = 8.8864e-10;
phi_I = 8.8864e-10;
phi_Sm = 0.0050;
phi_Em = 0.0050;
phi_Im = 0.0050;
lambda = 6.0381e-06;
lambda_m = 3.9248e-06;


eta_Sh = eta_1;
eta_Eh = eta_1;
eta_Ih = eta_1;
eta_S = eta_2;
eta_E = eta_2;
eta_I = eta_2;

k = 1e5;
alpha = c/(1+exp(-k*(timeDay-t0)));

%% sigma_S

if(timeDay<=t0)
    sigma_S = 0;
end
if(timeDay>t0 && timeDay<=t0+180)
    sigma_S = sigma_Smax*((timeDay-t0)/180);
end
if(timeDay>t0+180)
    sigma_S = sigma_Smax;
end

%% sigma_Sm

if(timeDay<=t0)
    sigma_Sm = 0;
end
if(timeDay>t0 && timeDay<=t0+180)
    sigma_Sm = sigma_Sm_max*((timeDay-t0)/180);
end
if(timeDay>t0+180)
    sigma_Sm = sigma_Sm_max;
end

%% sigma_Sh

if(timeDay<=t0)
    sigma_Sh = 0;
end
if(timeDay>t0 && timeDay<=t0+180)
    sigma_Sh = sigma_Sh_max*((timeDay-t0)/180);
end
if(timeDay>t0+180)
    sigma_Sh = sigma_Sh_max;
end

kappa_R = sigma_S;
kappa_Rm = sigma_Sm;
kappa_Rh = sigma_Sh;




n = 100;
Dt = 1/n;
% lambda = 6.0381e-06;
% lambda_m = 3.9248e-06;
lambda_h = 0;
for i = 1:n

    N = S+E+I+Sm+Em+Im+Sh+Eh+Ih+R+D+U;
    %lambda = (beta/N)*(I+Ih+(1-eta_i)*Im);
    lambda_m = (beta/N)*(1-eta_s)*(I+Ih+(1-eta_i)*Im);
%     lambda_h = (beta/N)*(1-eta_h)*(I+Ih+(1-eta_i)*Im);

    lambda = (beta/N)*(I+Ih+(1-eta_i)*Im);
    %lambda_m = (beta/N)*(1-eta_s)*(I+Ih+(1-eta_i)*Im);
    lambda_h = (beta/N)*(1-eta_h)*(I+Ih+(1-eta_i)*Im)*0;

    
    dS = (-(eta_Sh + alpha + phi_Sm + lambda)*S + kappa_R * R + eta_S *Sh + sigma_S * U + phi_S*Sm)*Dt;
    dSm = (-(lambda_m + phi_S + alpha)*Sm + sigma_Sm * U + kappa_Rm * R + phi_Sm * S)*Dt;
    dSh = (-(alpha + eta_S + lambda_h)*Sh + sigma_Sh * U + eta_Sh * S + kappa_Rh * R)*Dt;

    dE = (-(epsilon + phi_Em + eta_Eh)*E + eta_E * Eh + phi_E * Em + lambda*S )*Dt;
    dEm = (-(epsilon + phi_E) * Em + phi_Em*E + lambda_m* Sm)*Dt;
    dEh = (-(epsilon + eta_E)*Eh + lambda_h * Sh + eta_Eh * E)*Dt;

    dI = (-(mu + gamma + eta_Ih + phi_Im)*I + epsilon*E + phi_I* Im + eta_I *Ih)*Dt;
    dIm = (-(phi_I + gamma + mu)*Im + phi_Im*I + epsilon*Em)*Dt;
    dIh = (-(eta_I + gamma + mu)*Ih + epsilon*Eh + eta_Ih*I)*Dt;

    dR = (-(kappa_R+kappa_Rm+kappa_Rh)*R+(I + Im + Ih)*gamma)*Dt;
    dD = ((I+Im+Ih)*mu)*Dt;
    dU = (-(sigma_S+sigma_Sm+sigma_Sh)*U+(S+Sh+Sm)*alpha)*Dt;

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
    U = U + dU  ;

    Vp = Vp + dVp  ;

end
end