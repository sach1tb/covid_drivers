function X = constrainSigma(X,sigmaLimitsMin,sigmaLimitsMax,P)
[len, npoints] = size(X);


% Apply constraints
for i = 1:len
    for k = 1:npoints
        if X(i,k) < sigmaLimitsMin(i)
            X(i,k) = sigmaLimitsMin(i) + 0.01*sqrt(P(i,i));
        elseif X(i,k) > sigmaLimitsMax(i)
            X(i,k)  = sigmaLimitsMax(i) - 0.001*sqrt(P(i,i)) ;
        end
    end
end
end