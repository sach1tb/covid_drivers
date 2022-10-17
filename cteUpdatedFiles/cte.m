function ztoxy=cte(type, X,Y,Z, timeDownSample, numberOfBins, support)
% function ytox=ete(type, X,Y, timeDownSample, numberOfBins, bw, support)
%
% X, Y are 1 x n time series
%
% type is 'hist', 'kde'
% timeDownSample is how much you downsample the series 1:timeDownSample:end
% numberOfBins sets the number of bins e.g. 4, 8, 12
% bw is bandwidth
% support is a 1x2 vector with minimum and maximum values of the support
% space

% 

switch type
    case 'hist'
        ztoxy=cte_hist(X,Y, Z, timeDownSample, numberOfBins, support);
    case 'sym'
        % support is m
        ztoxy=cte_sym(X,Y, Z, timeDownSample, support);
    otherwise
        ytox=0;
end



function ztoxy=cte_hist(X,Y, Z,timeDownSample, numberOfBins, support)
% function [ytox, h]=get_te_hist(X,Y, method, timeDownSample, numberOfBins)
%
% X, Y are 1 x n time series
%
% if method is 'hist'
% timeDownSample is how much you downsample the series 1:timeDownSample:end
% numberOfBins sets the number of bins e.g. 4, 8, 12
% support is a 1x2 vector with minimum and maximum values of the support
% space

% 

if isempty(timeDownSample), timeDownSample=1; end

[d, n]=size(X);
if n < d
    error('input data should be of the form d x n where d is the dimension');
end


% normalize data
if isempty(support)
    support=[-1 1];
    X=normalize(X);
    Y=normalize(Y);
end


ztoxy=ent([X(1,2:end); X(1,1:end-1); Y(1,1:end-1)], ones(3,1)*numberOfBins, ones(3,1)*support, 'x;y;z') - ...
     ent([X(1,1:end-1); Y(1,1:end-1)], ones(2,1)*numberOfBins, ones(2,1)*support, 'x;y') - ...
     ent([X(1,2:end); X(1,1:end-1); Y(1,1:end-1); Z(1,1:end-1)], ones(4,1)*numberOfBins, ones(4,1)*support, 'w;x;y;z') + ...
     ent([X(1,1:end-1); Y(1,1:end-1); Z(1,1:end-1)], ones(3,1)*numberOfBins, ones(3,1)*support, 'x;y;z');

function X=normalize(X)

% normalize
% X=X-min(X(:)); X=X/max(X(:)); X=X*2-1;
[d, n]=size(X);
if n < d
    error('input data should be of the form d x n where d is the dimension');
end

X=X-min(X,[],2)*ones(1,n);

X=X./(max(X,[],2)*ones(1,n));
X=X*2-1;


function ztoxy=cte_sym(X,Y, Z, timeDownSample, m)
% function ytox=ete_sym(X,Y, timeDownSample, m)
%
% X, Y are 1 x n time series (assumed discretized)
%
% if method is 'hist'
% timeDownSample is how much you downsample the series 1:timeDownSample:end
% this is 'l' in the Staniek paper
% numberOfBins sets the number of bins e.g. 4, 8, 12
% m is the length of the symbol vector (-> the number of bins is perms(1:m))

% Staniek, M., & Lehnertz, K. (2008). Symbolic transfer entropy. 
% Physical Review Letters, 100(15), 158101.

if isempty(timeDownSample), timeDownSample=1; end

[d, n]=size(X);
if n < d
    error('input data should be of the form d x n where d is the dimension');
end

if m>5
    error('symbol length very high?!');
end

syms=perms(1:m);
numberOfBins=size(syms,1);

% create symbols
kk=1;
for jj=1:n-m
    [~, idx]=sort(X(jj:jj+m-1));
    sX(kk)=find(ismember(syms,idx,'rows'));

    [~, idx]=sort(Y(jj:jj+m-1));
    sY(kk)=find(ismember(syms,idx,'rows'));
    
    [~, idx]=sort(Z(jj:jj+m-1));
    sZ(kk)=find(ismember(syms,idx,'rows'));
    kk=kk+1;
end


% markovOrder=1;
support=[1, numberOfBins];
ztoxy=ent([sX(1,2:end); sX(1,1:end-1); sY(1,1:end-1)], ones(3,1)*numberOfBins, ones(3,1)*support, 'x;y;z') - ...
     ent([sX(1,1:end-1); sY(1,1:end-1)], ones(2,1)*numberOfBins, ones(2,1)*support, 'x;y') - ...
     ent([sX(1,2:end); sX(1,1:end-1); sY(1,1:end-1); sZ(1,1:end-1)], ones(4,1)*numberOfBins, ones(4,1)*support, 'w;x;y;z') + ...
     ent([sX(1,1:end-1); sY(1,1:end-1); sZ(1,1:end-1)], ones(3,1)*numberOfBins, ones(3,1)*support, 'x;y;z');
        

% 
% ytox=cte_hist(sX,sY, sZ, markovOrder,timeDownSample, numberOfBins, [1 numberOfBins]);
        



