function [x1, x2, v1, v2, lidx]=vicsekHierarchyCoupling(N, ts, eta, etaf, ...
                        sp, ...
                        lidx, show,hierarchyFlag,cp)
% function X=vicsek2d(N, ts, eta, etaf, ...
%                         top, topval, sp, ...
%                         lidx, show)
%
% eta is the the noise value for dtheta
% ts number of timesteps
% N is the number of birds
% top is the topology ('metric', 'nn')
% topval is the value for metric, radius, for nn number of nearest
% neighbors
% show is 1 if you want to plot. set to 0 for analysis
%
% Vicsek, T., Czirak, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995). 
% Novel Type of Phase Transition in a System of Self-Driven Particles. 
% Physical Review Letters, 75, 1226-1229. doi:10.1103/PhysRevLett.75.1226

% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s);


% default params if you run it as is without arguments
if nargin ==0
    N=3; ts=800; eta=2; 
    show=1;
    boundary=0; % 0=periodic, 1=reflective
    etaf=eta;
    hierarchyFlag =1; % 0 = no coupling, 1 = hierarchical coupling
        
    % speed
    sp=0.5; % units/s
    cp = 0.7; %0 is no influence, 1 is full influence
    % leader identity, random at the random times
    lidx=ones(1,ts);

end

% the units can be anything to match whatever
% is realistic
b=[-.100 .100 -.100 .100]; % boundary (m)
dt=0.2; % seconds?

% initialize
x1=zeros(N,ts); x2=x1; v1=x1; v2=x1;
% use component-wise to speed up code
x1(:,1)=b(1)+rand(N,1)*(b(2)-b(1)); % x_1
x2(:,1)=b(3)+rand(N,1)*(b(4)-b(3)); % x_2


th=rand(N,1)*2*pi;
th_prev=th;
v1(:,1)=cos(th); % xdot_1 / v_1 (direction only)
v2(:,1)=sin(th); % xdot_2 / v_2 (direction only)



for k=1:ts
    
    % topology
    D=sqrt((x1(:,k)*ones(1,N)-ones(N,1)*x1(:,k)').^2 + ...
            (x2(:,k)*ones(1,N)-ones(N,1)*x2(:,k)').^2);
   
    % motion model 
    x1(:,k+1)=x1(:,k)+sp*v1(:,k)*dt;
    x2(:,k+1)=x2(:,k)+sp*v2(:,k)*dt;
    
    for ii=1:N
        
        if hierarchyFlag == 1
            if ii>1
                nidx = [ii ii-1];
            else 
                nidx = 1;
            end
        end
         if hierarchyFlag == 0
                nidx = [ii];
         end
        
        % average angle for particle hierarchy 1->2->3...
        sin_th=sum(v2(nidx,k))/numel(nidx);
        cos_th=sum(v1(nidx,k))/numel(nidx);
%         if ii = 2
%         end
        % follower
        if numel(nidx)>1
            th_hat_ri=atan2(sin_th, cos_th);
            % direction update 
            th=(1-cp)*th_prev(ii) + (cp)*th_hat_ri-etaf/2 + rand*etaf;
            th_prev(ii) = th;
        else % leader doesn't interact
            th_hat_ri=atan2(v2(ii,k), v1(ii,k));
            th=th_hat_ri-eta/2 + rand*eta;
            th_prev(ii) = th;
        end
        
        % update the direction for the next time-step
        v1(ii,k+1)=cos(th);
        v2(ii,k+1)=sin(th);
        
    end
    
    % periodic boundary condition
    %x1(x1(:,k+1)<b(1),k+1)=x1(x1(:,k+1)<b(1),k+1)+(b(2)-b(1));
    %x1(x1(:,k+1)>b(2),k+1)=x1(x1(:,k+1)>b(2),k+1)-(b(2)-b(1));

    %x2(x2(:,k+1)<b(3),k+1)=x2(x2(:,k+1)<b(3),k+1)+(b(4)-b(3));
    %x2(x2(:,k+1)>b(4),k+1)=x2(x2(:,k+1)>b(4),k+1)-(b(4)-b(3));
    
    % display; alternatively use the output for analysis
    if show
        figure(1); gcf; 
        clf
        plot(x1(:,k), x2(:,k), 'b.');
        hold on;
        
        plot(x1(lidx(k),k), x2(lidx(k),k), 'ro');
        b = [min(x1(:))-5 max(x1(:))+5 min(x2(:))-5 max(x2(:))+5];
        axis(b);
        axis square;

        quiver(x1(:,k), x2(:,k), v1(:,k)*.4, v2(:,k)*.4, 0, 'k');

        xlabel(sprintf('[%d (%d)]', k, lidx(:,k)));
        
        if boundary==2
            th=linspace(-pi, pi, 50);
            plot((b(2)-b(1))/2*cos(th), (b(2)-b(1))/2*sin(th), 'b');
        end
        
        drawnow;
    end
    
end