function [x_kp1] = seirDynamics(xk,eta_Ih,eta_Im,eta_Sm,eta_Sh,dt)
x_kp1=xk;
% state variables
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

xi_Sh = xk(15);
xi_Eh = xk(15);
xi_Ih = xk(15);

xi_S = xk(16);
xi_E = xk(16);
xi_I = xk(16);

alpha = xk(17);

phi_Sm = xk(18);
phi_Em = xk(18);
phi_Im = xk(18);

phi_S = xk(19);
phi_E = xk(19);
phi_I = xk(19);
% sigma is a single value and is scaled with the proportion of the
% susceptible population
sigma_S = xk(20)*S/(S+Sh+Sm);
sigma_Sm = xk(20)*Sm/(S+Sh+Sm);
sigma_Sh = xk(20)*Sh/(S+Sh+Sm);

kappa_R= xk(21)*S/(S+Sh+Sm);
kappa_Rm= xk(21)*Sm/(S+Sh+Sm);
kappa_Rh= xk(21)*Sh/(S+Sh+Sm);

mu = xk(22);
gamma = xk(23);
epsilon = xk(24);
n = 1000;
Dt = dt/n;


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



end
