function X=normalize(X)

% normalize
% X=X-min(X(:)); X=X/max(X(:)); X=X*2-1;
n=size(X,2);
X=X-min(X,[],2)*ones(1,n);

X=X./(max(X,[],2)*ones(1,n));
X=X*2-1;

% pXnp1XnYn=kde([Xnp1; Xn; Yn], bw);
% pXnYn=kde([Xn; Yn], bw);
% pXnp1Xn=kde([Xnp1; Xn], bw);
% pXn=kde(Xn, bw);
% 
% te=sum(pXnp1XnYn.*log2(pXnp1XnYn.*pXn/(pXnYn.*pXnp1Xn)));
