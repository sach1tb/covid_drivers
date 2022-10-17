function X = constrainSigma(X,sigmaLimitsMin,sigmaLimitsMax)
[len, npoints] = size(X);


% Apply constraints
for i = 1:len
    for k = 1:npoints
        if X(i,k) < sigmaLimitsMin(i)
            X(i,k) = sigmaLimitsMin(i);
        elseif X(i,k) > sigmaLimitsMax(i)
            X(i,k)  = sigmaLimitsMax(i);
        end
    end
end
end