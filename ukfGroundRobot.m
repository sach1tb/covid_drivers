clearvars
% close all
addpath('..\boundedline\boundedline')
addpath('..\boundedline\Inpaint_nans')

% generating 2D data from a mobile robot with unicycle kinematics
% set the parameters
% note the timestep size.
dt=1/10;
tspan = 0:dt:20; % time for simulation

% ode with initial conditions 
i1=[2 2 0];
v=5; 
% this is just a simulation so we keep the unicycle model because it makes
% nice tracks.
[~,gt] = ode45(@unicycle_ode,tspan,i1,[],v);

% measurement model is for a radar
eta=1;
z=radar1(gt(:,1:2)')+randn(2, size(gt,1))*eta;

T=size(z,2);


% ready to filter
% initialize
x(1:2,1)=gt(1,1:2); 
x(3:4,1)=1; % arbitrary velocity

n=4;

Fk= [ 1, 0, dt,  0
    0, 1, 0,  dt
    0, 0, 1,  0
    0, 0, 0,  1];

f=@(x) Fk*x;  % nonlinear state equations

Q= eye(4)*10;

h=@(x) radar1(x);                               % measurement equation
P=Q; %initial covariance
R=diag([2 2]*3);
pmat = zeros(n,n,T);
covarianceMatrix = zeros(n,n,T);
for k=1:T
    zk=z(:,k);                            % save actual state
    zV(:,k)  = zk;                         % save measurment
    [x, P] = ukf(f,x,P,h,zk,Q,R);
        pmat(:,:,k) = P;

    xV(:,k) = x;                            % save estimate
    disp(k);
end

% plot results
% show the results
figure(1); gcf; clf;

ylabels={'x_1', 'x_2', '\theta'};
gt=gt';
for ii=1:3
    subplot(2,2,ii+1); gca;
    plot(gt(ii,:), 'k', 'linewidth', 2);
    hold on;
%     plot(Z(ii,:), 'k*');
    plot(xV(ii,:), 'r', 'linewidth', 2);
    boundedline(1:T, xV(ii,:), squeeze(P(ii,ii,:))', '-g', 'alpha');
    set(gca, 'fontname', 'times', 'fontsize', 24);
    grid on;
    ylabel(ylabels{ii});
end
xlabel('time');

subplot(3,2,1);
plot(gt(1,:), gt(2,:), 'k', 'linewidth', 2);
hold on;
plot(xV(1,:), xV(2,:), 'r', 'linewidth', 2);
axis image;
set(gca, 'fontname', 'times', 'fontsize', 24);
legend('ground truth','estimate', 'location', 'northeastoutside');




function Xdot = unicycle_ode(t,X,v)
Xdot(1,1) = v*cos(X(3));
Xdot(2,1) = v*sin(X(3));
% set omega as a function of time to create interesting trajectory
Xdot(3,1) = 1*sin(.5*t);
end

function z = radar1(X)

z(1,:)=sqrt(X(1,:).^2+X(2,:).^2);
z(2,:)=atan2(X(2,:), X(1,:));
end