function ixy=emi(X,Y, numberOfBins, support)
% mutual information

% I(X;Y) = H(X) + H(Y) - H(X;Y)
%
% 
% t=0:.1:10; delay=t.^1.1;
% x=sin(1*t); y=sin(1*t+delay); plot(t, x,t,y); emi(x, y, 10, [-1 1]);

if numel(numberOfBins) ==1
    numberOfBins=numberOfBins*ones(1,2);
end

if size(support,1)==1
    support=ones(2,1)*support;
end

ixy=    ent(X, numberOfBins(1), support(1,:), 'x') + ...
        ent(Y, numberOfBins(2), support(2,:), 'x') - ...
        ent([X; Y], numberOfBins, support, 'x;y');