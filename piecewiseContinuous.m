function y = piecewiseContinuous(xVector,tVector,currentTime)
y = -1;
ifelseLength = numel(tVector);

if currentTime < tVector(1)
    y = xVector(1);
end
if y ==-1
    for i = 2:ifelseLength
        if (currentTime>= tVector(i-1) && currentTime < tVector(i))
            y = xVector(i);
        end
    end
end

if (y == -1 && currentTime >= tVector(end))
    y = xVector(end);
end

end