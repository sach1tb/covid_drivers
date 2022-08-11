function [S,E,I,Sm,Em,Im,R,D,V,phi_S,phi_Sm,lambda,lambda_m]=seir_sim_wrapper(S0,E0,I0,Sm0,Em0,Im0,eta_s,eta_i,sigma,beta ...
    ,epsilon,gamma,mu,alpha,a,b,tau,omega,kappa,noise,tspan, phi_SMult,phi_EMult,phi_IMult,phi_SmMult,phi_EmMult,phi_ImMult)


dt = 1;
R0 = 0;
D0 = 0;
V0 = 0;
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

phi_S = zeros(1,numel(tspan)-1);
phi_Sm = phi_S;
lambda = phi_S;
lambda_m = phi_S;
for k=1:numel(tspan)-1
    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), x(k+1,8), x(k+1,9),~,~,~,~,phi_S(k),phi_Sm(k),lambda(k),lambda_m(k)] = ...
        seir_sim(eta_s,eta_i,sigma,beta,epsilon,gamma,mu,alpha,a,b,tau,omega,kappa(k),x(k,1) ...
        , x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), x(k,7), x(k,8), x(k,9), dt,noise, phi_SMult,phi_EMult,phi_IMult,phi_SmMult,phi_EmMult,phi_ImMult);
end

S = x(:,1);
E = x(:,2);
I = x(:,3);
Sm= x(:,4);
Em = x(:,5);
Im= x(:,6);
R = x(:,7);
D = x(:,8);
V = x(:,9);
end
