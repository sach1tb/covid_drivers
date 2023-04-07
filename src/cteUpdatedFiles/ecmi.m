function ixycz=ecmi(X,Y, Z, numberOfBins, support)
%
% calculates I(X;Y|Z) 
% numberOfBins sets the number of bins e.g. 4, 8, 12
% support is a 1x2 vector with minimum and maximum values of the support
% space

% 

if numel(numberOfBins) ==1
    numberOfBins=numberOfBins*ones(1,3);
end

if size(support,1)==1
    support=ones(3,1)*support;
end


ixycz=  ent([X; Z], numberOfBins([1 3]), support([1 3],:), 'x;y') + ...
        ent([Y; Z], numberOfBins([2 3]), support([2 3],:), 'x;y') - ...
        ent([X; Y; Z], numberOfBins, support, 'x;y;z') - ...
        ent(Z, numberOfBins(3), support(3,:), 'x');
        