clear variables

addpath cteUpdatedFiles/

% does the entropy change when a random variable is mapped through a
% nonlinear function?

n=5000;
X=randn(n)*2;

[freq, sup]=histcounts(X, -5:.1:5);

p=freq/sum(freq);

Hx=ent(normalize(X), 10,[0 1], 'x')

Y=2*(X);
Hy=ent(normalize(Y), 10, [0 1], 'x')

Y=(X)/3;
Hy=ent(normalize(Y), 10, [0 1], 'x')

fprintf('nonlinear...\n');

Y=sin(X);
Hy=ent(normalize(Y), 10, [0 1], 'x')


Y=(X).^2;
Hy=ent(normalize(Y), 10, [0 1], 'x')