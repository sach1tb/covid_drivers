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

Xsigma = zeros(1,1*2+1);
sigmaPointAccumulutor = zeros(size(Xsigma,1),size(Xsigma,2),n);
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
  Xsigma=sigmas(x,P,0.001);
  sigmaPointAccumulutor(:,:,k) = Xsigma;
  s = f(s) + q*randn;                % update process 
end

zV=[zV, zV(end)];
xV=[xV, xV(end)];

% calculate I(Idot(sigma); param(sigma))
IXYsig_dist=zeros(3,1);
for aa=1:3
    Xhat=squeeze(sigmaPointAccumulutor(1,aa,:))';
    IXYsig_dist(aa)=emi_with_shuffle(Xhat, Y(1:end-1), 0, [], [], 1, 0);
end

% calculate I(Idot(shuffle);param(shuffle)) four times each shuffle
nsh=3000;



IXYsig_dist_shuffle=zeros(3,nsh);
for ss=1:nsh
    for aa=1:3
        Xhat=squeeze(sigmaPointAccumulutor(1,aa,:))';
        IXYsig_dist_shuffle(aa,ss)=emi_with_shuffle(Xhat, Y(randperm(numel(Y)-1)), 0, [], [], 1, 0);
    end
end

[Ixysig_freq, ixysig_sup]=histcounts(IXYsig_dist_shuffle(:), 0:.01:1);
ixysig_sup=ixysig_sup(1:end-1);
pxysig=Ixysig_freq/sum(Ixysig_freq);


[IXhatY]=emi_with_shuffle(xV,Y, 0, numberOfBins, ...
    support, binarize, nsymbols);

[IZY]=emi_with_shuffle(zV,Y, 0, numberOfBins, ...
    support, binarize, nsymbols);
% 
% [IXY, pIshuffle, Isup, IXYshuffle]=emi_with_shuffle(X,Y, nshuffle, numberOfBins, ...
%     support, binarize, nsymbols);

figure(1); gcf; clf;
plot(X, 'r', 'linewidth', 2); hold on;
plot(xV, 'r:', 'linewidth', 2);
plot(squeeze(sigmaPointAccumulutor)', 'color', [1,.5, .5], 'linewidth', 1);
plot(zV, 'k*');
set(gca, 'fontsize', 16);
xlabel('time');
legend('X(t)', '\hat{X}(t)', 'Z_X(t)');
% one-tailed test 
figure(2); gcf; clf;
alpha=0.01;
area(ixysig_sup, pxysig, 'facecolor', 'k', 'facealpha', 0.1);
hold on;
% one-tailed test visual with alpha
cpI=cumsum(pxysig);
idx=find(cpI>1-alpha);

[val, p_idx]=min(sqrt((ixysig_sup-mean(IXYsig_dist(:))).^2));
if val < 0.01
    pval=1-cpI(p_idx);
else
    error('fix support resolution');
end
%             pval=ranksum(mean(Ip2sig_dist(:)),Ip2sig_dist_shuffle(:),'tail', 'right', 'alpha', 0.05);
h1=plot([mean(IXYsig_dist(:)),mean(IXYsig_dist(:))], [0,1], 'r-', 'linewidth', 2);
plot([ixysig_sup(idx(1)),ixysig_sup(idx(1))], [0,1], 'k:', 'linewidth', 2);

plot([IZY,IZY], [0,1], 'b--', 'linewidth', 2);
plot([IXhatY,IXhatY], [0,1], 'r:', 'linewidth', 2);


% cpI=cumsum(pIshuffle);
% idx=find(cpI>1-alpha);
% 
% figure(2); gcf; clf;
% plot(Isup,pIshuffle, 'k', 'linewidth', 2);
% hold on;
% plot([IXY,IXY], [0,1], 'r--', 'linewidth', 2);
% plot([Isup(idx(1)),Isup(idx(1))], [0,1], 'k:', 'linewidth', 2);

legend('$I(X,Y_{sh})$', '$I(X,Y)$', '$p<0.01$',...
        '$I(Z_X;Y)$', '$I(\hat{X};Y)$');
xlabel('$I(bits)$', 'interpreter', 'latex');
ylabel('$p(I)$', 'interpreter', 'latex');
set(gca, 'fontsize', 16);


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