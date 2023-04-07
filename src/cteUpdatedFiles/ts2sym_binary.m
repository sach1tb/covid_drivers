function [sX, numberOfBins] = ts2sym_binary(X, m)
% X is the 1D data
% m is the symbol length
%
%
% Bandt, C., & Pompe, B. (2002). Permutation entropy: a natural complexity measure for time series. 
% Physical review letters, 88(17), 174102.
% 
% instead of creating symbols based on how they can be sorted, we simply
% permute the binary symbols with repetition allowed so that the symbols
% are
% permn([0 1], 3)
% 
% ans =
% 
%      0     0     0
%      0     0     1
%      0     1     0
%      0     1     1
%      1     0     0
%      1     0     1
%      1     1     0
%      1     1     1

[d, n]=size(X);
if n < d
    error('input data should be of the form d x n where d is the dimension');
end

if sum(X~=0 & X~=1)
   error('only binary sequences allowed');
end 
syms=permn([0 1], m);
numberOfBins=size(syms,1);

% create symbols
kk=1;
for jj=1:n-m
    % note: symbol is simply the location of that sequence in the time
    % series
    sX(kk)=find(ismember(syms,X(jj:jj+m-1),'rows'));
    kk=kk+1;
end