function [te, X, Y]=ete_kde(X,Y,bw, support)
% function te=ete1(X,Y,bw)
% 
% X and Y are d x n signals
% bw is the bandwidth or the resolution
%
% te is the transfer entropy from Y -> X
%

% get the size of each variable (it may be multi-dimensional)
[d, nx]=size(X);

% normalize data
if isempty(support)
    support=[-1 1];
    X=normalize(X);
    Y=normalize(Y);
end

% get the chains X(n), Y(n), and X(n+1)
Xn=X(:,1:end-1);
Yn=Y(:,1:end-1);
Xnp1=X(:,2:end);

% x=[Xnp1; Xn; Yn];

% we set n=10 because with three variables it becomes 1000 and that is a
% lot
% irrespective of the dimension, just pick the support to be between min
% and max for all variables
n=10;
s=ones(3*d,1)*linspace(support(1), support(2), n);


% [x1, x2, x3]=ndgrid(s(1,:), s(2,:), s(3,:));
% this part assigns the indices of the support to create the new 3D version
idx=zeros(n,n,n,3);
[idx(:,:,:,1), idx(:,:,:,2), idx(:,:,:,3)]=ndgrid(1:n);
s1=zeros(3*d,n^3);

% assign to each support
for j=1:3
    k=d*(j-1)+1:d*j;
    s1(k,:)=s(k,idx(:,:,:,j));
end
s=s1;

% compute individual components of the entropy % x=[Xnp1; Xn; Yn];
pXXY=kde([Xnp1; Xn; Yn], bw, s);
pXX=kde([Xnp1; Xn], bw, s(1:d*2,:));
pX=kde(Xn, bw, s(d+1:d*2,:));
pXY=kde([Xn; Yn], bw, s(d+1:d*3,:));

% nz=(pXXY~=0 & pXX~=0 & pX~=0 & pXY~=0);
% instead make them all more than epsilon because small values can exist
nz=(pXXY>eps & pXX>eps & pX>eps & pXY>eps);
te=sum(pXXY(nz).*log2(pXXY(nz)./((pXX(nz)./pX(nz)).*pXY(nz))));

% te1-te

% [bw p]=kde3(X,10/bw);
% h=-sum(p(p~=0).*log2(p(p~=0)));


function X=normalize(X)

n=size(X,2);
X=X-min(X,[],2)*ones(1,n);

X=X./(max(X,[],2)*ones(1,n));
X=X*2-1;
