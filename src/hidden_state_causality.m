clear variables

addpath cteUpdatedFiles/

rng(3);
T=120;
% ulam map;
fx=@(x) 2-x.^2;
X=zeros(1,T); Y=X;
eta=10*eps;
epsilon=0.3;
% X -> Y
for k=1:T
    X(k+1)=fx(X(k))+randn*eta;
    Y(k+1)=fx(epsilon*X(k)+(1-epsilon)*Y(k));
end


numberOfBins=8; support=[-2,2];
binarize=1;
nsymbols=0;

nshuffle=1000;

% nonlinear measurement
ZX=(X.^2);

% estimate using UKF
n=1;      %number of state
q=eps;    %std of process 
r=0.1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement  
f=@(x) 2-x.^2;  % nonlinear state equations, ulam map
h=@(x)x.^2;     % measurement equation, nonlinear
s=0;            % initial state
x=s+q*randn(1,1); % initial state with noise
P = eye(n);                               % initial state covraiance
xV = zeros(n,T);          %estmate        % allocate memory
sV = zeros(n,T);          %actual
zV = zeros(1,T);
for k=1:T
  z = h(s) + r*randn;                % measurments
  sV(:,k)= s;                        % save actual state
  zV(k)  = z;                        % save measurment
  [x, P] = ukf(f,x,P,h,z,Q,R);       % ekf 
  xV(:,k) = x;                       % save estimate
  s = f(s) + q*randn;                % update process 
end
zV=[zV, zV(end)];
xV=[xV, xV(end)];

for ii=1:

[IXYhat]=emi_with_shuffle(xV,Y, nshuffle, numberOfBins, ...
    support, binarize, nsymbols);

[IZY]=emi_with_shuffle(zV,Y, nshuffle, numberOfBins, ...
    support, binarize, nsymbols);

[IXY, pIshuffle, Isup, IXYshuffle]=emi_with_shuffle(X,Y, nshuffle, numberOfBins, ...
    support, binarize, nsymbols);

figure(1); gcf; clf;
plot(X, 'r', 'linewidth', 2); hold on;
plot(xV, 'r:', 'linewidth', 2);
plot(zV, 'k*');
set(gca, 'fontsize', 16);
xlabel('time');
legend('X(t)', '\hat{X}(t)', 'Z_X(t)');
% one-tailed test visual with alpha
alpha=0.01;
cpI=cumsum(pIshuffle);
idx=find(cpI>1-alpha);

figure(2); gcf; clf;
plot(Isup,pIshuffle, 'k', 'linewidth', 2);
hold on;
plot([IXY,IXY], [0,1], 'r--', 'linewidth', 2);
plot([Isup(idx(1)),Isup(idx(1))], [0,1], 'k:', 'linewidth', 2);
plot([IZY,IZY], [0,1], 'b--', 'linewidth', 2);
plot([IXYhat,IXYhat], [0,1], 'r:', 'linewidth', 2);

legend('$I(X,Y_{sh})$', '$I(X,Y)$', '$p<0.01$',...
        '$I(Z_X;Y)$', '$I(\hat{X};Y)$');
xlabel('$I(bits)$', 'interpreter', 'latex');
ylabel('$p(I)$', 'interpreter', 'latex');
set(gca, 'fontsize', 16);