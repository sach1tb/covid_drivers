function [S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,phi_1,phi_2,xi_1,xi_2,lambda,lambda_m,lambda_h,alpha,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm ...
    ,kappa_R,kappa_Rm,kappa_Rh] ...
    = seir_simIsolation( ...
    S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N,sigma_ShVector, sigma_SVector, sigma_SmVector, xi1Vector, xi2Vector, phi1Vector, phi2Vector ...
        ,kappa_RhVector, kappa_RVector, kappa_RmVector, ...
        gammaVector, muVector, alphaVector ...
        , epsilonVector, betaVector,eta_Ih,eta_Im,eta_Sh,eta_Sm,t1, t2, t3, t4, t5, t6, currentDay )



%%
% sigma_Sh, sigma_S, sigma_Sm : t1 t2
% xi1, xi2, phi1, phi2 : t3 t4 t5
% kappa_Rh, kappa_R, kappa_Rm : t6  ...10
% gamma, mu : t1 t2
% alpha, epsilon : t1 t2
% beta : t1 t2  ...15
% eta_Ih, eta_Im, eta_Sh,eta_Sm
%%



sigma_Sh = piecewiseContinuous(sigma_ShVector,[t1, t2],currentDay);
sigma_S = piecewiseContinuous(sigma_SVector,[t1, t2],currentDay);
sigma_Sm = piecewiseContinuous(sigma_SmVector,[t1, t2],currentDay);
xi_1 = piecewiseContinuous(xi1Vector,[t3,t4,t5],currentDay);
xi_2 = piecewiseContinuous(xi2Vector,[t3,t4,t5],currentDay);
phi_1 = piecewiseContinuous(phi1Vector,[t3,t4,t5],currentDay);
phi_2 = piecewiseContinuous(phi2Vector,[t3,t4,t5],currentDay);
kappa_Rh = piecewiseContinuous(kappa_RhVector,[t6],currentDay);
kappa_R = piecewiseContinuous(kappa_RVector,[t6],currentDay);
kappa_Rm = piecewiseContinuous(kappa_RmVector,[t6],currentDay);
gamma = piecewiseContinuous(gammaVector,[t1, t2],currentDay);
mu = piecewiseContinuous(muVector,[t1, t2],currentDay);
alpha = piecewiseContinuous(alphaVector,[t1, t2],currentDay);
epsilon = piecewiseContinuous(epsilonVector,[t1, t2],currentDay);
beta = piecewiseContinuous(betaVector,[t1, t2],currentDay);

phi_Sm = phi_1;
phi_Em = phi_1;
phi_Im = phi_1;
phi_S = phi_2;
phi_E = phi_2;
phi_I = phi_2;


xi_Sh = xi_1;
xi_Eh = xi_1;
xi_Ih = xi_1;
xi_S = xi_2;
xi_E = xi_2;
xi_I = xi_2;

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