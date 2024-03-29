function [x,P,Xsigma]=ukfConstrained(fstate,x,P,hmeas,z,Q,R,sigmaLimitsMin,sigmaLimitsMax)
% UKF   Unscented Kalman Filter for nonlinear dynamic systems
% [x, P] = ukf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P
% for nonlinear dynamic system (for simplicity, noises are assumed as additive):
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           h: function handle for h(x)
%           z: current measurement
%           Q: process noise covariance
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance
%
% Example:
%{
n=3;      %number of state
q=0.1;    %std of process
r=0.1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
h=@(x)x(1);                               % measurement equation
s=[0;0;1];                                % initial state
x=s+q*randn(3,1); %initial state          % initial state with noise
P = eye(n);                               % initial state covraiance
N=20;                                     % total dynamic steps
xV = zeros(n,N);          %estmate        % allocate memory
sV = zeros(n,N);          %actual
zV = zeros(1,N);
for k=1:N
  z = h(s) + r*randn;                     % measurments
  sV(:,k)= s;                             % save actual state
  zV(k)  = z;                             % save measurment
  [x, P] = ukf(f,x,P,h,z,Q,R);            % ekf
  xV(:,k) = x;                            % save estimate
  s = f(s) + q*randn(3,1);                % update process
end
for k=1:3                                 % plot results
  subplot(3,1,k)
  plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
end
%}
% Reference: Julier, SJ. and Uhlmann, J.K., Unscented Filtering and
% Nonlinear Estimation, Proceedings of the IEEE, Vol. 92, No. 3,
% pp.401-422, 2004.
%
% By Yi Cao at Cranfield University, 04/01/2008
%
L=numel(x);                                 %number of states
m=numel(z);                                 %number of measurements
alpha=1;                                 %default, tunable
ki=0;                                       %default, tunable
beta=2;                                     %default, tunable
lambda=alpha^2*(L+ki)-L;                    %scaling factor
c=L+lambda;                                 %scaling factor
Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
c=sqrt(c);

try
    chol(P);
catch
    P=nearestSPD(P); % ideally this should go
end
X=sigmas(x,P,c);                            %sigma points around x
X = constrainSigma(X,sigmaLimitsMin,sigmaLimitsMax);          % constrain
[x1,X1,P1,X2]=ut(fstate,X,Wm,Wc,L,Q);          %unscented transformation of process


try
    chol(P1);
catch
    P1=nearestSPD(P1); % ideally this should go
end
%
X1=sigmas(x1,P1,c);                         %sigma points around x1
X1 = constrainSigma(X1,sigmaLimitsMin,sigmaLimitsMax);          % constrainZ
% X2=X1-x1(:,ones(1,size(X1,2)));             %deviation of X1
[z1,Z1,P2,Z2]=ut(hmeas,X1,Wm,Wc,m,R);       %unscented transformation of measurements
P12=X2*diag(Wc)*Z2';                        %transformed cross-covariance
K=P12/P2; %*inv(P2);
x=x1+K*(z-z1);                              %state update
P=P1-K*P12';                                %covariance update

% apply the constraints to the ukf estimate as well
% Kandepu2008
% "If UKF estimate violates the constraints,
% the same projection technique can be used."
x=constrainSigma(x,sigmaLimitsMin,sigmaLimitsMax);

try
    chol(P);
catch
    P=nearestSPD(P); % ideally this should go
end
Xsigma=sigmas(x,P,c);
Xsigma=constrainSigma(Xsigma,sigmaLimitsMin,sigmaLimitsMax);
% recalculate P
% X3=Xsigma-x(:,ones(1,L));
% P=X3*diag(Wc)*X3'+0;

function [y,Y,P,Y1]=ut(f,X,Wm,Wc,n,R)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations

L=size(X,2);
y=zeros(n,1);
Y=zeros(n,L);
for k=1:L
    Y(:,k)=f(X(:,k));
    y=y+Wm(k)*Y(:,k);
end
Y1=Y-y(:,ones(1,L));
P=Y1*diag(Wc)*Y1'+R;


function X=sigmas(x,P,c)
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points

A = c*chol(P)';
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A];



function X = constrainSigma(X,sigmaLimitsMin,sigmaLimitsMax)
[~, npoints] = size(X);

% vectorize for speed
for k=1:npoints
    idx1=X(:,k)<sigmaLimitsMin';
    X(idx1,k)=sigmaLimitsMin(idx1);
    
    idx2=X(:,k)>sigmaLimitsMax';
    X(idx2,k)=sigmaLimitsMax(idx2);
end

