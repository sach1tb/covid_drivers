function X = constrainSigma(X,sigmaLimitsMin,sigmaLimitsMax)
[~, npoints] = size(X);


% Apply constraints
% for i = 1:len
%     for k = 1:npoints
%         if X(i,k) < sigmaLimitsMin(i)
%             X(i,k) = sigmaLimitsMin(i);
%         elseif X(i,k) > sigmaLimitsMax(i)
%             X(i,k)  = sigmaLimitsMax(i);
%         end
%     end
% end
% end

% vectorize for speed
for k=1:npoints
    idx1=X(:,k)<sigmaLimitsMin';
    X(idx1,k)=sigmaLimitsMin(idx1);
    
    idx2=X(:,k)>sigmaLimitsMax';
    X(idx2,k)=sigmaLimitsMax(idx2);
end