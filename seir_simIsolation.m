function [S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,phi_Sm,phi_S,xi_Sh,xi_S,lambda,lambda_m,lambda_h,alpha] ...
    = seir_simIsolation( ...
    S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,a1,a2,b1,b2,partT, ...
    sigma_Smax,sigma_Sm_max,sigma_Sh_max,c,k,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,t0,timeDay,numParts )


a11 = a1(1);
a12 = a1(2);
a13 = a1(3);
a21 = a2(1);
a22 = a2(2);
a23 = a2(3);

b11 = b1(1);
b12 = b1(2);
b13 = b1(3);
b21 = b2(1);
b22 = b2(2);
b23 = b2(3);

t1 = partT(1);
t2 = partT(2);
t3 = partT(3);
syms t;
phi_1 = piecewise((0<t)&(t<=t1),a11,(t1<t)&(t<=t2),a12,(t2<t)&(t<=t3),a13);
phi_2 = piecewise((0<t)&(t<=t1),a21,(t1<t)&(t<=t2),a22,(t2<t)&(t<=t3),a23);
xi_1 = piecewise((0<t)&(t<=t1),b11,(t1<t)&(t<=t2),b12,(t2<t)&(t<=t3),b13);
xi_2 = piecewise((0<t)&(t<=t1),b21,(t1<t)&(t<=t2),b22,(t2<t)&(t<=t3),b23);


phi_Sm = double(subs(phi_1,timeDay));
phi_Em = double(subs(phi_1,timeDay));
phi_Im = double(subs(phi_1,timeDay));
phi_S = double(subs(phi_2,timeDay));
phi_E = double(subs(phi_2,timeDay));
phi_I = double(subs(phi_2,timeDay));


xi_Sh = double(subs(xi_1,timeDay));
xi_Eh = double(subs(xi_1,timeDay));
xi_Ih = double(subs(xi_1,timeDay));
xi_S = double(subs(xi_2,timeDay));
xi_E = double(subs(xi_2,timeDay));
xi_I = double(subs(xi_2,timeDay));

% ,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,

beta = piecewise((0<t)&(t<=t1),beta(1),(t1<t)&(t<=t2),beta(2),(t2<t)&(t<=t3),beta(3));
beta = double(subs(beta,timeDay));

eta_Ih = piecewise((0<t)&(t<=t1),eta_Ih(1),(t1<t)&(t<=t2),eta_Ih(2),(t2<t)&(t<=t3),eta_Ih(3));
eta_Ih = double(subs(eta_Ih,timeDay));

eta_Im = piecewise((0<t)&(t<=t1),eta_Im(1),(t1<t)&(t<=t2),eta_Im(2),(t2<t)&(t<=t3),eta_Im(3));
eta_Im = double(subs(eta_Im,timeDay));

eta_Sh = piecewise((0<t)&(t<=t1),eta_Sh(1),(t1<t)&(t<=t2),eta_Sh(2),(t2<t)&(t<=t3),eta_Sh(3));
eta_Sh = double(subs(eta_Sh,timeDay));

eta_Sm = piecewise((0<t)&(t<=t1),eta_Sm(1),(t1<t)&(t<=t2),eta_Sm(2),(t2<t)&(t<=t3),eta_Sm(3));
eta_Sm = double(subs(eta_Sm,timeDay));

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



n = 1000;
Dt = 1/n;
lambda = 0;
lambda_m = 0;
lambda_h = 0;
for i = 1:n

    N = S+E+I+Sm+Em+Im+Sh+Eh+Ih+R+D+U;
%     lambda = (beta/N)*(I+Ih+(1-eta_i)*Im);
%     lambda_m = (beta/N)*(1-eta_s)*(I+Ih+(1-eta_i)*Im);
%     lambda_h = (beta/N)*(1-eta_h)*(I+Ih+(1-eta_i)*Im);

    lambda = (beta/N)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);
    lambda_m = (beta/N)*(1-eta_Sm)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);
    lambda_h = (beta/N)*(1-eta_Sh)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);

    dS = (-(xi_Sh + alpha + phi_Sm + lambda)*S + kappa_R * R + xi_S *Sh + sigma_S * U + phi_S*Sm)*Dt;
    dSm = (-(lambda_m + phi_S + alpha)*Sm + sigma_Sm * U + kappa_Rm * R + phi_Sm * S)*Dt;
    dSh = (-(alpha + xi_S + lambda_h)*Sh + sigma_Sh * U + xi_Sh * S + kappa_Rh * R)*Dt;

    dE = (-(epsilon + phi_Em + xi_Eh)*E + xi_E * Eh + phi_E * Em + lambda*S )*Dt;
    dEm = (-(epsilon + phi_E) * Em + phi_Em*E + lambda_m* Sm)*Dt;
    dEh = (-(epsilon + xi_E)*Eh + lambda_h * Sh + xi_Eh * E)*Dt;

    dI = (-(mu + gamma + xi_Ih + phi_Im)*I + epsilon*E + phi_I* Im + xi_I *Ih)*Dt;
    dIm = (-(phi_I + gamma + mu)*Im + phi_Im*I + epsilon*Em)*Dt;
    dIh = (-(xi_I + gamma + mu)*Ih + epsilon*Eh + xi_Ih*I)*Dt;

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