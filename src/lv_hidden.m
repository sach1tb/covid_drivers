clear variables

addpath cteUpdatedFiles/


% data=readtable('../data/hudson_bay_data.xlsx');
T=240;
dt=0.1;
% lv model
X=zeros(3,T); 

b=.5; % prey growth rate (exponential)
s=.1; % predation rate (more predators mean more predation)
d=1; % predator death rate
e=.1; % predator growth rate, (nutritional based on the number of prey)

% b=0.5*ones(1,T);

% model
X(:,1)=[d/(e*s)+10,b(1)/s, .3]';
for k=1:T
    X(:,k+1)=lvsim(X(:,k), s,e,d,k,dt);
end




numberOfBins=8; support=[-2,2];
binarize=1;
nsymbols=0;

nshuffle=1000;

% estimate using UKF
n=3;      %number of state
m=2;      %number measurements
q=[1 .1 .1]';    %std of process 
r=[1 .1];    %std of measurement
Q=diag(q.^2); % covariance of process
R=r.^2;        % covariance of measurement  
h=@(x) x([1,3]);     % measurement equation, nonlinear
x0=[10,2.5,.1]';            % initial state
x=x0+q*randn(1,1); % initial state with noise

Xsigma = zeros(n,n*2+1);
sigmaPointAccumulutor = zeros(size(Xsigma,1),size(Xsigma,2),n);
P = eye(n);                               % initial state covraiance
xV = zeros(n,T);          %estmate        % allocate memory
sV = zeros(n,T);          %actual
zV = zeros(m,T);

Z=X([1, 3],:)+diag(r)*randn(2,size(X,2));

for k=1:T
  zV(:,k)  = Z(:,k);                        % save measurment
  f=@(x) lvsim(x, s, e, d, k, dt);  % nonlinear state equations, ulam map

  [x, P] = ukf(f,x,P,h,zV(:,k),Q,R);       % ekf 
  xV(:,k) = x;                       % save estimate
  Xsigma=sigmas(x,P,0.001);
  sigmaPointAccumulutor(:,:,k) = Xsigma;
end


figure(1); gcf; clf;
labels={'prey', 'predatory', 'growth rate'};
for ii=1:n
    subplot(1,n,ii);
    plot(1:T, X(ii,1:end-1), 'k', 'linewidth', 2); hold on;
    plot(1:T, xV(ii,:), 'r:', 'linewidth', 2);
    ylabel(labels{ii});
    set(gca, 'fontsize', 24);
end

figure(1); gcf; 
zidx=[1 3];

for jj=1:m
    subplot(1,n,zidx(jj));
    plot(1:T, zV(jj,:), 'k*', 'markersize', 4);
end

[IZ1Z2, pZZIshuffle, IZZsup]=emi_with_shuffle(xV(1,:),xV(3,:), ...
        50000, [], [], 1, 0);
[IX1X2, pXXIshuffle, IXXsup]=emi_with_shuffle(xV(1,:),xV(2,:), ...
        50000, [], [], 1, 0);

alpha=0.01;
figure(2); gcf; clf;
subplot(1,2,1);
area(IZZsup, pZZIshuffle, 'facecolor', 'k', 'facealpha', 0.1);
hold on;
% one-tailed test visual with alpha
cpI=cumsum(pZZIshuffle);
idx=find(cpI>1-alpha);

plot([IZ1Z2,IZ1Z2], [0,1], 'r-', 'linewidth', 2);
plot([IZZsup(idx(1)),IZZsup(idx(1))], [0,1], 'k:', 'linewidth', 2);

subplot(1,2,2);
area(IXXsup, pXXIshuffle, 'facecolor', 'k', 'facealpha', 0.1);
hold on;
% one-tailed test visual with alpha
cpI=cumsum(pXXIshuffle);
idx=find(cpI>1-alpha);

plot([IX1X2,IX1X2], [0,1], 'r-', 'linewidth', 2);
plot([IXXsup(idx(1)),IXXsup(idx(1))], [0,1], 'k:', 'linewidth', 2);
% legend('$I(H,b_{sh})$', '$I(H,b)$',...
%         '$I(H;P)$');
xlabel('$I(bits)$', 'interpreter', 'latex');
ylabel('$p(I)$', 'interpreter', 'latex');
set(gca, 'fontsize', 16);


function xkp1=lvsim(xk, s,e,d, k, dt)


    hk=xk(1,1);
    pk=xk(2,1);
    bk=xk(3,1);

    hkp1=hk+(bk*hk - s*hk*pk)*dt;
    pkp1=pk+(e*s*hk*pk - d*pk)*dt;
    bkp1=0.5*sin(2*pi/12*k*dt-pi/12)+0.75;
%     b=0.5*sin(2*pi/12*(0:dt:T)-pi/2)+0.75;

    xkp1(1,1)=hkp1;
    xkp1(2,1)=pkp1;
    xkp1(3,1)=bkp1;
end

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
end