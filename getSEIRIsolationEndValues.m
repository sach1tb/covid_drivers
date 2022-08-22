function [S,Sm,Sh,E,Em,Eh,I,Im,Ih,R,D,U,Vp,N] ...
    = getSEIRIsolationEndValues( ...
    S0,Sm0,Sh0,E0,Em0,Eh0,I0,Im0,Ih0,R0,D0,U0,Vp0,N0,params,infectionDuration,dayCounter)


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

phi_1 = params(1);
phi_2= params(2);
xi_1= params(3);
xi_2= params(4);
sigma_S= params(5);
sigma_Sm= params(6);
sigma_Sh= params(7);
c= params(8);
mu= params(9);
gamma= params(10);
epsilon= params(11);
beta= params(12);
eta_Ih= params(13);
eta_Im= params(14);
eta_Sh= params(15);
eta_Sm= params(16);
kappa_R = params(17);
kappa_Rm = params(18);
kappa_Rh = params(19);


for k=1:numel(infectionDuration)

    [x(k+1,1), x(k+1,2), x(k+1,3), x(k+1,4), x(k+1,5), x(k+1,6), x(k+1,7), ...
        x(k+1,8), x(k+1,9),x(k+1,10), x(k+1,11), x(k+1,12), x(k+1,13), x(k+1,14)]...
    = seir_simIsolation( ...
    x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), ...
        x(k,7), x(k,8), x(k,9),x(k,10), x(k,11), x(k,12), x(k,13), x(k,14),phi_1,phi_2,xi_1,xi_2, ...
    sigma_S,sigma_Sm,sigma_Sh,c,mu,gamma,epsilon,beta,eta_Ih,eta_Im,eta_Sh,eta_Sm,kappa_R,kappa_Rm, kappa_Rh,dayCounter );

end

S= x(end,1);
Sm= x(end,2);
Sh= x(end,3);
E= x(end,4);
Em= x(end,5);
Eh= x(end,6);
I= x(end,7);
Im= x(end,8);
Ih= x(end,9);
R= x(end,10);
D= x(end,11);
U= x(end,12);
Vp= x(end,13);
N =  x(end,14);



end