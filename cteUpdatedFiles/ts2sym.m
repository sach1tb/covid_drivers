function [sX, numberOfBins] = ts2sym(X, m)
% X is the 1D data
% m is the symbol length
%
%
% Bandt, C., & Pompe, B. (2002). Permutation entropy: a natural complexity measure for time series. 
% Physical review letters, 88(17), 174102.
% 
% Staniek, M., & Lehnertz, K. (2008). Symbolic transfer entropy. 
% Physical Review Letters, 100(15), 158101.
%
% Let us take a series with seven values: x = (4,7,9,10,6,11,3)
% We compare three consecutive values. (4,7,9) and (7,9,10) represent the 
% permutation 012 since they are in increasing order. 
% (9,10,6) and (6,11,3) correspond to the permutation 201 since x(t+2) < x(t) < x(t+1)
% while (10,6,11) has the permutation type 102 with x(t+1) < x(t) < x(t+2).
%
%
% % instead of creating symbols based on how they can be sorted, we simply
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

sdata=(unique(X));
isbin=0;

% binary time series should be set as above 
if numel(sdata)==2 && ismember(0, sdata) && ismember(1, sdata)
    syms=permn([0 1], m);
    isbin=1;
else
    syms=perms(1:m);
end
numberOfBins=size(syms,1);

% create symbols
kk=1;
for jj=1:n-m
    if ~isbin
        % note: symbol is simply the sort of the time series are sorted 
        [~, idx]=sort(X(jj:jj+m-1));
        sX(kk)=find(ismember(syms,idx,'rows'));
    else
        sX(kk)=find(ismember(syms,X(jj:jj+m-1),'rows'));
    end
    kk=kk+1;
end