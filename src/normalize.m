function X=normalize(X)
% function X=normalize(X)
% X is n x T where n is the # of signals and T is the length

% normalize w.r.t. all signals 
% use if all signals are of the same type, e.g. turn rate of spp
X=X-min(X(:)); X=X/max(X(:)); X=X*2-1;

% normalize each signal separately
% use if signals are different e.g. heart rate and breath rate
% n=size(X,2);
% X=X-min(X,[],2)*ones(1,n);
% 
% X=X./(max(X,[],2)*ones(1,n));
% X=X*2-1;
