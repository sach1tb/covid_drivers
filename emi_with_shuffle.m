function [IXY, pIshuffle, Isup, IXYshuffle]=emi_with_shuffle(X,Y, nshuffle, numberOfBins, ...
    support, binarize, nsymbols)
%
% X,Y are 1xn, 1-D time series
% nshuffle is a large number 
% 
% binarize: if 1, then binarize, and the number of bins and support will 
%           change automatically
% nsymbols: should be more than 1 to convert time series to symbolic form,
%           otherwise set to 0
%
% NOTE: 
% 'addpath cteUpdatedFiles/' before running this script
% 
% if running without any argument, the script will test with ulam map;
% may need to run several times since ulam map can diverge

if nargin < 1 
    rng(3);
    T=90;
    % ulam map;
    fx=@(x) 2-x.^2;
    X=zeros(1,T); Y=X;
    eta=1e-10;
    epsilon=0.3;
    % X -> Y
    for k=1:T
        X(k+1)=fx(X(k))+randn*eta;
        Y(k+1)=fx(epsilon*X(k)+(1-epsilon)*Y(k));
    end
    
    % diff Y so that we are looking at causality I(X,Ydot)?
%     Y=[diff(Y), Y(end)];
    
    numberOfBins=8; support=[-2,2];
    binarize=1; 
    nsymbols=0;
    
    nshuffle=1000;
end

if binarize
    X=binarize_wrt_median(X);
    Y=binarize_wrt_median(Y);
    numberOfBins=2; support=[0,1];
end

if nsymbols
    [X, ~]=ts2sym(X,nsymbols);
    [Y, nb]=ts2sym(Y,nsymbols);
    numberOfBins=nb; support=[1 nb];
end

IXYshuffle=zeros(1,nshuffle);
IXY=emi(X, Y, numberOfBins, support);

for ii=1:nshuffle
    % shuffling does not change H(X), H(Y), but it 
    % changes the joint distribution of X[k],Y[k] and therefore H(X,Y)
    % note that we shuffle here instead of the raw time series... it's ok
    % since the effect should be the same before and after binarization but
    % not necessarily when we symbolize
    IXYshuffle(1,ii)=emi(X, Y(randperm(numel(Y))), numberOfBins, support);
end

[freq,Isup]=hist(IXYshuffle,0:.01:.1);
pIshuffle=freq/sum(freq);

if nargin<1 % check if it works
    figure(1); gcf; clf;
    plot(X, 'r'); hold on;
    plot(Y, 'k');
    
    figure(2); gcf; clf;
    plot(Isup,pIshuffle, 'k', 'linewidth', 2);
    hold on;
    plot([IXY,IXY], [0,1], 'r--', 'linewidth', 2);
    legend('I(X,shuffled)', 'I(X,Y)');
    xlabel('I(bits)');
    ylabel('p(I)');
    set(gca, 'fontsize', 16);
end


