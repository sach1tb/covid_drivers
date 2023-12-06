function [IXYcZ, pIshuffle, Isup, IXYcZshuffle]=ecmi_with_shuffle(X,Y,Z, nshuffle, ...
                    numberOfBins, support, binarize, nsymbols)
%
% X,Y,Z are 1xn, 1-D time series
% nshuffle is a large number 
% 
% binarize: if 1, then binarize, and the number of bins and support will 
%           change automatically
% nsymbols: should be more than 1 to convert time series to symbolic form,
%           otherwise set to 0
%
% 
% if running without any argument, the script will test with ulam map;
% may need to run several times since ulam map can diverge
if nargin < 1 
    rng(5);
    T=120;
    % ulam map;
    % we use CMI like transfer entropy so that we should pick up
    % TE(Y->X) = H(X_k+1|X_k)-H(X_k+1|X_k,Y_k)
    % by treating X_k+1 as one series and X_k as another
    % TE/cmi or I(X_k+1; Y_k|X_k) > I(X_k+1;Y_k(shuffled)|X_k)
    % because
    % TE(Y->X) = I(X_k+1; Y_k|X_k) c.f. I(X;Y|Z) so that X is X_k+1, Y is
    % Y_k, and Z is X_k
    % 
    % Shahsavari et al. (2020). Estimating conditional transfer entropy... 
    % Entropy, 22(10), 1124.
    fx=@(x) 2-x.^2;
    X1=zeros(1,T); Y1=X1; 
    eta=1e-10;
    
    epsilon=0.1;
    % X -> Y
    for k=1:T+1
        X1(k+1)=fx(X1(k))+randn*eta;
        Y1(k+1)=fx(epsilon*X1(k)+(1-epsilon)*Y1(k));
    end
    X=Y1(2:end); % x_k+1
    Y=X1(1:end-1); % y_k
    Z=Y1(1:end-1); % x_k
    
    numberOfBins=8; support=[-2,2];
    binarize=1; 
    nsymbols=0;
    
    nshuffle=1000;
end

X=detrend(X);
Y=detrend(Y);
Z=detrend(Z);

if binarize
    X=binarize_wrt_median(X);
    Y=binarize_wrt_median(Y);
    Z=binarize_wrt_median(Z);
    numberOfBins=2; support=[0,1];
end

if nsymbols
    X=ts2sym(X,nsymbols);
    Y=ts2sym(Y,nsymbols);
    [Z, nb]=ts2sym(Z,nsymbols);
    numberOfBins=nb; support=[1 nb];
end

IXYcZshuffle=zeros(1,nshuffle);
IXYcZ=ecmi(X, Y, Z, numberOfBins, support);
% IXYcZ=ete_hist(Z,Y,1,numberOfBins, support);

for ii=1:nshuffle
    % shuffling does not change H(X), H(Y), but it 
    % changes the joint distribution of X[k],Y[k] and therefore H(X,Y)
    % note that we shuffle here instead of the raw time series... it's ok
    % since the effect should be the same before and after binarization but
    % not necessarily when we symbolize
    IXYcZshuffle(1,ii)=ecmi(X, Y, Z(randperm(numel(Z))), numberOfBins, support);
%     IXYcZshuffle(1,ii)=ete_hist(Z, Y(randperm(numel(Y))), 1,numberOfBins, support);

end

[freq,Isup]=hist(IXYcZshuffle);
pIshuffle=freq/sum(freq);


if nargin<1 % check if it works
    figure(1); gcf; clf;
    plot(X, 'r'); hold on;
    plot(Y, 'k');
    plot(Z, 'g');
    
    figure(2); gcf; clf;
    plot(Isup,pIshuffle, 'k', 'linewidth', 2);
    hold on;
    plot([IXYcZ,IXYcZ], [0,1], 'r--', 'linewidth', 2);
    legend('I(X,Y|shuffled)', 'I(X,Y|Z)');
    xlabel('I(bits)');
    ylabel('p(I)');
    set(gca, 'fontsize', 16);
end


