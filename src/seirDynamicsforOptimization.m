function [S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,V,N,phi_1,phi_2,xi_1,xi_2,lambda,...
    lambda_m,lambda_h,alpha,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa0] ...
    = seirDynamicsforOptimization(S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,V,N,...
                    sigma0, xi_1, xi_2, phi_1, phi_2, ...
                    gamma, mu,kappa0, alpha, ...
                    epsilon, beta,eta_Ih,eta_Im,eta_Sh,eta_Sm)



%%
% sigma_Sh, sigma_S, sigma_Sm
% xi_1, xi_2, phi_1, phi_2
% gamma, mu
% kappa_Rh, kappa_R, kappa_Rm
% alpha, epsilon
% beta
% eta_Ih, eta_Im, eta_Sh,eta_Sm
%%

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



sigma_S = sigma0*S/(S+Sh+Sm);
sigma_Sm = sigma0*Sm/(S+Sh+Sm);
sigma_Sh = sigma0*Sh/(S+Sh+Sm);

kappa_R= kappa0*S/(S+Sh+Sm);
kappa_Rm= kappa0*Sm/(S+Sh+Sm);
kappa_Rh= kappa0*Sh/(S+Sh+Sm);


n = 1000;
Dt = 1/n;
lambda = 0;
lambda_m = 0;
lambda_h = 0;


for i = 1:n
    T = S+Sm+Sh+E+Em+Eh+I+Im+Ih+R+U;
    lambda = (beta/T)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);
    lambda_m = (beta/T)*(1-eta_Sm)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);
    lambda_h = (beta/T)*(1-eta_Sh)*(I+(1-eta_Ih)*Ih+(1-eta_Im)*Im);

    dS = (-(xi_Sh + alpha + phi_Sm + lambda)*S + kappa_R * R + xi_S *Sh + sigma_S * U + phi_S*Sm)*Dt;
    dSm =(-(lambda_m + phi_S + alpha)*Sm + sigma_Sm * U + kappa_Rm * R + phi_Sm * S)*Dt;
    dSh =(-(alpha + xi_S + lambda_h)*Sh + sigma_Sh * U + xi_Sh * S + kappa_Rh * R)*Dt;

    dE =(-(epsilon + phi_Em + xi_Eh)*E + xi_E * Eh + phi_E * Em + lambda*S )*Dt;
    dEm =(-(epsilon + phi_E) * Em + phi_Em*E + lambda_m* Sm)*Dt;
    dEh =(-(epsilon + xi_E)*Eh + lambda_h * Sh + xi_Eh * E)*Dt;

    dI =(-(mu + gamma + xi_Ih + phi_Im)*I + epsilon*E + phi_I* Im + xi_I *Ih)*Dt;
    dIm =(-(phi_I + gamma + mu)*Im + phi_Im*I + epsilon*Em)*Dt;
    dIh = (-(xi_I + gamma + mu)*Ih + epsilon*Eh + xi_Ih*I)*Dt;

    dR =  (-(kappa_R + kappa_Rm + kappa_Rh)*R+(I + Im + Ih)*gamma)*Dt;
    dD = ((I + Im + Ih)*mu)*Dt;
    dU = (-(sigma_S + sigma_Sm + sigma_Sh)*U+(S + Sh + Sm)*alpha)*Dt;

    totGrads = dS + dSm + dSh + dE + dEm + dEh + dI + dIm + dIh + dR + dD + dU;

    % check if the sum of gradients is below a threshold
    if abs(totGrads) > 1
        disp(totGrads);
        keyboard
    end

    dV =  alpha*(S + Sm+ Sh)*Dt;

    S = S+dS;
    Sm = Sm + dSm;
    Sh = Sh + dSh;

    E = E + dE;
    Em = Em + dEm;
    Eh = Eh + dEh;

    I = I + dI;
    Im = Im + dIm;
    Ih = Ih + dIh;

    R = R + dR;
    D = D + dD;
    U = U + dU;

    V = V + dV;

end


end