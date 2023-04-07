function ytox=ete_sym(X,Y, timeDownSample, m)
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

% syms=perms(1:m);
% numberOfBins=size(syms,1);
% 
% % create symbols
% kk=1;
% for jj=1:n-m
%     [~, idx]=sort(X(jj:jj+m-1));
%     sX(kk)=find(ismember(syms,idx,'rows'));
% 
%     [~, idx]=sort(Y(jj:jj+m-1));
%     sY(kk)=find(ismember(syms,idx,'rows'));
%     kk=kk+1;
% end

[sX, numberOfBinsX]= ts2sym(X, m);
[sY, numberOfBinsY]= ts2sym(Y, m);


% ytox=te_hist(sX,sY,markovOrder,timeDownSample, numberOfBins, [1 numberOfBins]);
supportX=1:numberOfBinsX;
supportY=1:numberOfBinsY;
ytox=ent([sX(1,2:end); sX(1,1:end-1)], [numberOfBinsX; numberOfBinsX], [supportX; supportX], 'x;y') - ...
     ent(sX(1,1:end-1), numberOfBinsX, supportX, 'x') - ...
     ent([sX(1,2:end); sX(1,1:end-1); sY(1,1:end-1)], [numberOfBinsX; numberOfBinsX; numberOfBinsY], ...
                        [supportX; supportX; supportY], 'x;y;z') + ...
     ent([sX(1,1:end-1); sY(1,1:end-1)], [numberOfBinsX; numberOfBinsY], [supportX; supportY], 'x;y');
        

function ytox=te_hist(X,Y,markovOrder,timeDownSample, numberOfBins, support)

% Schreiber, T., 2000. Measuring information transfer. Physical Review Letters, 85(2), pp.461?464.
% ytox=p(X[n], X[n-1], ..., X[n-markovOrder], Y[n-1], Y[n-2],
% Y[n-l])log((p(X[n]|X[n-1], ..., X[n-markovOrder], Y[n-1], Y[n-2], ...,
% Y[n-l])/p(X[n]|X[n-1], ..., X[n-markovOrder]))
%
% markovOrder=1, l=1
%
% ytox=p(X[n], X[n-1], Y[n-1])log((p(X[n]|X[n-1], Y[n-1])/p(X[n]|X[n-1]))
%
%


% downsample
X=X(1:timeDownSample:end);
Y=Y(1:timeDownSample:end);

% binning
% Z=[X(:), Y(:)]; % use both datasets
% Z=X(:);
% bins=linspace(min(Z(:)), max(Z(:)), numberOfBins);
bins=linspace(support(1), support(2), numberOfBins); % because the incoming data is normalized anyway
% fprintf('binends=[%.1f, %.1f]\n', bins(1), bins(end));
binwidth=bins(2)-bins(1);
bb=bins;
bins=[-inf bins(1)-binwidth/2 bins+binwidth/2];
bins(end)=bb(end);
nbins=numel(bins);

% initialize the data
pXXY=zeros(nbins-1, nbins-1, nbins-1);
N=numel(X);

% create the markov chain
Xchain=zeros(markovOrder+1, N-markovOrder); 
for jj=markovOrder+1:-1:1
    Xchain(markovOrder+2-jj,:)=X(jj:jj+N-markovOrder-1);
end
Y=Y(1:N-markovOrder);


% compute 3D histogram
% Y=sin(2*pi*(1:numel(Y))/4); % this matches
% Y=sin(2*pi*(1:numel(Y))/4); % this is lower freq
for ii=1:nbins-1
    for jj=1:nbins-1
        for kk=1:nbins-1
            % the same markovOrder is used since we are using a logical operator
            pXXY(ii,jj,kk)=sum(Xchain(1,:)>=bins(ii) & Xchain(1,:)<=bins(ii+1) & ...
                               Xchain(2,:)>=bins(jj) & Xchain(2,:)<=bins(jj+1) & ...
                               Y>=bins(kk) & Y<=bins(kk+1));
        end
    end
end

pXXY=pXXY(2:end, 2:end, 2:end);
pXXY=pXXY/sum(pXXY(:));

% p(X[markovOrder-1], Y[markovOrder-1])
pXY=zeros(nbins-1, nbins-1);
for jj=1:nbins-1
    for kk=1:nbins-1
        pXY(jj,kk)=sum(Xchain(2,:)>=bins(jj) & Xchain(2,:)<=bins(jj+1) & ...
                               Y>=bins(kk) & Y<=bins(kk+1));
    end
end
pXY=pXY(2:end,2:end);
pXY=pXY/sum(pXY(:));


% p(X[markovOrder], X[markovOrder-1])
pXX=zeros(nbins-1,nbins-1);
for ii=1:nbins-1
    for jj=1:nbins-1
        pXX(ii,jj)=sum(Xchain(1,:)>=bins(ii) & Xchain(1,:)<=bins(ii+1) & ...
                               Xchain(2,:)>=bins(jj) & Xchain(2,:)<=bins(jj+1));
    end
end
pXX=pXX(2:end,2:end);
pXX=pXX/sum(pXX(:));


% p(X[markovOrder-1])
pXm1=zeros(nbins-1,1);
for jj=1:nbins-1
    pXm1(jj)=sum(Xchain(2,:)>=bins(jj) & Xchain(2,:)<=bins(jj+1));
end
pXm1=pXm1(2:end);
pXm1=pXm1/sum(pXm1(:));


pXXY(isnan(pXXY))=0;
pXX(isnan(pXX))=0;
pXY(isnan(pXY))=0;
pXm1(isnan(pXm1))=0;

% ytox=p(X[n], X[n-1], Y[n-1])log((p(X[n]|X[n-1], Y[n-1])/p(X[n]|X[n-1]))
ytox=0;
for ii=1:nbins-2
    for jj=1:nbins-2
        for kk=1:nbins-2
            if pXXY(ii,jj,kk) && pXX(ii,jj) && pXY(jj,kk) && pXm1(jj)
                ytox=ytox+pXXY(ii,jj,kk)*log2(pXXY(ii,jj,kk)/((pXX(ii,jj)/pXm1(jj))*pXY(jj,kk)));
            end
        end
    end
end

