function ytox=ete(type, X,Y, timeDownSample, numberOfBins, bw, support)
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
        ytox=ete_hist(X,Y, timeDownSample, numberOfBins, support);
    case 'kde'
        ytox=ete_kde(X,Y,bw, support);
    case 'sym'
        % support is m
        ytox=ete_sym(X,Y,timeDownSample, support);
    otherwise
        ytox=0;
end