clearvars

addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
addpath(['cteUpdatedFiles', filesep])

load ukfOutput.mat  %size is 24, 24*2+1

% debug set to 1 for testing
debug=1;

%rise
rise = [201, 279;
        627, 694;
        775, 891];
rise_params_labels={'\xi_2', '\phi_2', '\sigma', '\kappa'};
rise_params=[16,19,20,21];     

%fall
fall = [280, 376;
        694, 761;
        892, 958 ];
fall_params_labels={'\xi_1', '\alpha','\phi_1'};
fall_params=[15,17,18];      

nshuffle=10000;

alpha=0.01; % for one-tailed test

Itot = sum(xV(7:9,:),1);
Itotdot=diff(Itot);


% smooth all sigma points
sigmaPointAccumulutor1=sigmaPointAccumulutor*0;
for ii=1:24
    for jj=1:49
        sigmaPointAccumulutor1(ii,jj,:)=sma(sigmaPointAccumulutor(ii,jj,:),7);
    end
end

Itot_sigma=squeeze(sum(sigmaPointAccumulutor1(7:9,:,:),1));
Itot_sigma_dot=diff(Itot_sigma,1,2);

% [  beta, xi1, xi2, alpha, ...
%             phi1, phi2, sigma0, kappa0, ...
%             mu, gamma, epsilon]=vec2params(vec)

% rises
if debug
    figure(1); gcf; clf;
    figure(2); gcf; clf; jj=0;
end
for ii=1:numel(rise_params)
    param=xV(rise_params(ii),1:end-1);
    param_sigma=squeeze(sigmaPointAccumulutor1(rise_params(ii),:,1:end-1));
    param_label=rise_params_labels{ii};
    % plot the parameter to be doubly sure
    if debug
        figure(1); gcf; 
        subplot(1,numel(rise_params),ii);
        plot(param, 'linewidth', 2);
        ylabel(param_label);
    end
    for rr=1:size(rise,1)
        win=rise(rr,1):rise(rr,2);
        
        % calculate I(Idot;param)
        [Ip2Itotdot, pIshuffle, Isup, Ip2shuffle]=emi_with_shuffle(Itotdot(win), param(win),nshuffle, [], ...
                        [], 1, 0);
                    
        % calculate I(Idot(sigma); param(sigma))            
        Ip2sig_dist=zeros(49,49);
        for aa=1:49
            for bb=1:49
                Itot_X=Itot_sigma_dot(aa,win);
                param_Y=param_sigma(bb,win);
                Ip2sig_dist(aa,bb)=emi_with_shuffle(Itot_X, param_Y, 1, [], [], 1, 0);
            end
        end
        [Ip2sig_freq, ip2sig_sup]=hist(Ip2sig_dist(:));
        pIp2sig=Ip2sig_freq/sum(Ip2sig_freq);
        
        if debug
            jj=jj+1;
            figure(2); gcf;
            subplot(4,3,jj);
            plot(Isup,pIshuffle, 'k', 'linewidth', 2);
            hold on;
            plot(ip2sig_sup, pIp2sig, 'r:', 'linewidth', 2);
            
            % one-tailed test visual with alpha
            cpI=cumsum(pIshuffle);
            idx=find(cpI>1-alpha);
            
            
            plot([Ip2Itotdot,Ip2Itotdot], [0,1], 'r-', 'linewidth', 2);
            % check, then comment
%             plot([Ip2sig_dist(1,1), Ip2sig_dist(1,1)], [0,1],'b:', 'linewidth', 2);
            plot([Isup(idx(1)),Isup(idx(1))], [0,1], 'k:', 'linewidth', 2);
            xlabel(['I($\dot{I};' param_label '$) (bits)'], 'interpreter', 'latex');
            ylabel('p');
            if jj<4
            title(sprintf('rise %d', rr));
            end
            set(gca, 'fontsize', 16, 'xlim', [0 0.35]);
        end            
    end
end


% falls
if debug
    figure(3); gcf; clf;
    figure(4); gcf; clf; jj=0;
end
for ii=1:numel(fall_params)
    param=xV(fall_params(ii),1:end-1);
    param_label=fall_params_labels{ii};
    % plot the parameter to be doubly sure
    if debug
        figure(3); gcf; 
        subplot(1,numel(fall_params),ii);
        plot(param, 'linewidth', 2);
        ylabel(param_label);
    end
    for ff=1:size(fall,1)
        win=fall(ff,1):fall(ff,2);
        [Ip2Itotdot, pIshuffle, Isup, Ip2shuffle]=emi_with_shuffle(Itotdot(win), param(win),nshuffle, [], ...
                        [], 1, 0);

        
        % calculate I(Idot(sigma); param(sigma))            
        Ip2sig_dist=zeros(49,49);
        for aa=1:49
            for bb=1:49
                Itot_X=sma(Itot_sigma_dot(aa,win),7);
                param_Y=sma(param_sigma(bb,win),7);
                Ip2sig_dist(aa,bb)=emi_with_shuffle(Itot_X, param_Y, 1, [], [], 1, 0);
            end
        end
        [Ip2sig_freq, ip2sig_sup]=hist(Ip2sig_dist(:));
        pIp2sig=Ip2sig_freq/sum(Ip2sig_freq);            
                    
        if debug
            jj=jj+1;
            figure(4); gcf;
            subplot(numel(fall_params),3,jj);
            plot(Isup,pIshuffle, 'k', 'linewidth', 2);
            hold on;
            plot(ip2sig_sup, pIp2sig, 'g:', 'linewidth', 2);
            % one-tailed test with alpha
            cpI=cumsum(pIshuffle);
            idx=find(cpI>1-alpha);
            
            plot([Ip2Itotdot,Ip2Itotdot], [0,1], 'g-', 'linewidth', 2);
            % check then comment
%             plot([Ip2sig_dist(1,1), Ip2sig_dist(1,1)], [0,1],'b:', 'linewidth', 2);
            plot([Isup(idx(1)),Isup(idx(1))], [0,1], 'k:', 'linewidth', 2);
            xlabel(['I($\dot{I};' param_label '$) (bits)'], 'interpreter', 'latex');
            ylabel('p');
            if jj<4
            title(sprintf('fall %d', ff));
            end
            set(gca, 'fontsize', 16, 'xlim', [0 0.6]);
        end            
    end
end





