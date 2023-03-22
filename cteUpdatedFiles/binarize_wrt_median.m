function Xb = binarize_wrt_median(X)
% X is the 1D data
% Xb is the binarized data based on whether a value is more or less than
% the median
%

if nargin < 1
    % generate a random process
    X=randn(1,100);
    X=sin(0:.1:10);
end

[d, n]=size(X);
if n < d
    error('input data should be of the form d x n where d is the dimension');
end

% calculate the median
% X should be row vector
Xbar=median(X,2);

Xb=double(X>Xbar);

% check whether it works
if nargin <1
    figure(1); gcf; clf;
    plot(X, 'k:');
    hold on;
    plot(X*0+Xbar, 'g');
    plot(Xb, 'r--');
end

