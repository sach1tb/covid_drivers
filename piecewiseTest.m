clear all
close all
clc

tspan = 1:1000;

tVector = [200 500 600];
xVector = [1 2 1 4];
for i = 1:numel(tspan)
y(i) = piecewiseContinuous(xVector,tVector,tspan(i));
end
plot(y)