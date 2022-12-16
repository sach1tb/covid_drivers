clear all
close all
clc

addpath('..\..\..\..\boundedline\boundedline')
addpath('..\..\..\..\Inpaint_nans')
addpath('cteUpdatedFiles\')
N=3; ts=400; eta=2;
show=0; etaf =eta; sp = 0.05;hierarchyFlag = 1;
boundary=0; % 0=periodic, 1=reflective
lidx = 1*ones(1,ts);
x = [];
y = [];
z = [];
c = 1;
cp = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
for cap = 1:numel(cp)
    for run = 1:10
        [xPos, yPos, xVel, yVel, lidx]=vicsekHierarchyCoupling(N, ts, eta, etaf, sp,lidx, show, hierarchyFlag,cp(cap));
        heading1 = (atan2(yVel(1,:),xVel(1,:)));
        heading2 = (atan2(yVel(2,:),xVel(2,:)));
        heading3 = (atan2(yVel(3,:),xVel(3,:)));
        x(run,:)=downsample(heading1,1);
        y(run,:)=downsample(heading2,1);
        z(run,:)=downsample(heading3,1);
    end
    window = ts;

    CTExToyz = [];
    CTEyTozx = [];
    CTExTozy = [];
    CTEzToxy = [];
    CTEzToyx = [];
    for i = 1:size(y,1)
        c = 1;
        for j = 1:1:size(y,2)-window
            CTExToyz(i,c)= 0;
            CTEyTozx(i,c)= 0;
            CTExTozy(i,c)= 0;
            CTEzToxy(i,c) = 0;
            c = c+1;
        end
    end


    maxLength = size(y,2);
    maxHeight = size(y,1);
    j = 1:1:maxLength-window;
    innerLoopLength = numel(j);
    for i = 1:maxHeight
        xTemp = normalize(x(i,:));
        yTemp = normalize(y(i,:));
        zTemp = normalize(z(i,:));
        for k = 1
            CTExToyz(i)= cte('hist',yTemp(k:k+window-1),zTemp(k:k+window-1),xTemp(k:k+window-1),1,ceil(sqrt(window)),[-1 1]);
            CTEyTozx(i)= cte('hist',zTemp(k:k+window-1),xTemp(k:k+window-1),yTemp(k:k+window-1),1,ceil(sqrt(window)),[-1 1]);
            CTExTozy(i)= cte('hist',zTemp(k:k+window-1),yTemp(k:k+window-1),xTemp(k:k+window-1),1,ceil(sqrt(window)),[-1 1]);
            CTEzToxy(i)= cte('hist',xTemp(k:k+window-1),yTemp(k:k+window-1),zTemp(k:k+window-1),1,ceil(sqrt(window)),[-1 1]);
            CTEzToyx(i)= cte('hist',yTemp(k:k+window-1),xTemp(k:k+window-1),zTemp(k:k+window-1),1,ceil(sqrt(window)),[-1 1]);
        end
    end

    meanCTExToyz(cap) = mean(CTExToyz);
    stdCTExToyz(cap) =  std(CTExToyz);

    meanCTExTozy(cap) = mean(CTExTozy);
    stdCTExTozy(cap) =  std(CTExTozy);

    meanCTEyTozx(cap) = mean(CTEyTozx);
    stdCTEyTozx(cap) =  std(CTEyTozx);

    meanCTEzToxy(cap) = mean(CTEzToxy);
    stdCTEzToxy (cap)=  std(CTEzToxy);

    meanCTEzToyx(cap) = mean(CTEzToyx);
    stdCTEzToyx(cap) =  std(CTEzToyx);
end




figure(1)
clf;
h1 = plot(cp,meanCTExToyz,'-r');
hold on
h2 = plot(cp,meanCTEyTozx,'-b');
h3 = plot(cp,meanCTExTozy,'-g');
h4 = plot(cp,meanCTEzToxy,'-c');
h5 = plot(cp,meanCTEzToyx,'-y');
boundedline(cp,meanCTExToyz,stdCTExToyz, '-r','alpha','linewidth',1.2);
boundedline(cp,meanCTEyTozx,stdCTEyTozx, '-b','alpha','linewidth',1.2);
boundedline(cp,meanCTExTozy,stdCTExTozy, '-g','alpha','linewidth',1.2);
boundedline(cp,meanCTEzToxy,stdCTEzToxy, '-c','alpha','linewidth',1.2);
boundedline(cp,meanCTEzToyx,stdCTEzToyx, '-y','alpha','linewidth',1.2);
set(gca, 'fontsize', 24);
ylabel('CTE (bits)');
xlabel('Coupling Parameter')
ylim([0, inf])
%legend([h1,h2, h3,h4, h5],'$CE_{sin(\theta_{X})\rightarrow sin(\theta_{Y})|sin(\theta_{Z})}$','$CE_{sin(\theta_{Y})\rightarrow sin(\theta_{Z})|sin(\theta_{X})}$','$CE_{sin(\theta_{X})\rightarrow sin(\theta_{Z})|sin(\theta_{Y})}$','$CE_{sin(\theta_{Z})\rightarrow sin(\theta_{X})|sin(\theta_{Y})}$','$CE_{sin(\theta_{Z})\rightarrow sin(\theta_{Y})|sin(\theta_{X})}$','interpreter','latex');
legend([h1,h2, h3,h4, h5],'$CE_{\theta_{X}\rightarrow \theta_{Y}|\theta_{Z}}$','$CE_{\theta_{Y}\rightarrow \theta_{Z}|\theta_{X}}$','$CE_{\theta_{X}\rightarrow \theta_{Z}|\theta_{Y}}$','$CE_{\theta_{Z}\rightarrow \theta_{X}|\theta_{Y}}$','$CE_{\theta_{Z}\rightarrow \theta_{Y}|\theta_{X}}$','interpreter','latex');
fig = gcf;
fig.Position = [681 389 810 590];
