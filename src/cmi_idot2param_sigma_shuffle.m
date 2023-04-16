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
fall_params_labels={'\xi_1', '\phi_1', '\alpha'};
fall_params=[15,18,17];

nshuffle=10000;

alpha=0.01; % for one-tailed test

Itot_sigma=squeeze(sum(sigmaPointAccumulutor(7:9,:,:),1));
Itot_sigma_dot=diff(Itot_sigma,1,2);

% [  beta, xi1, xi2, alpha, ...
%             phi1, phi2, sigma0, kappa0, ...
%             mu, gamma, epsilon]=vec2params(vec)

% rises
for ii=1:numel(rise_params)
    if debug
        figure(ii); gcf; clf;
        jj=0;
    end
    param=xV(rise_params(ii),1:end-1);
    param_sigma=squeeze(sigmaPointAccumulutor(rise_params(ii),:,1:end-1));
    param_label=rise_params_labels{ii};
    for kk=1:numel(rise_params)
        if ii ~=kk
            param2=xV(rise_params(kk),1:end-1);
            param2_sigma=squeeze(sigmaPointAccumulutor(rise_params(kk),:,1:end-1));
            param2_label=rise_params_labels{kk};
            
            for rr=1:size(rise,1)
                win=rise(rr,1):rise(rr,2);
                
                
                 % calculate I(Idot(sigma); param(sigma))
                Ip2sig_dist=zeros(49,49,49);
                for aa=1:49
                    for bb=1:49
                        for cc=1:49
                            Itot_X=Itot_sigma_dot(aa,win);
                            param_Y=param_sigma(bb,win);
                            param_Z=param2_sigma(cc,win);
                            Ip2sig_dist(aa,bb,cc)=ecmi_with_shuffle(Itot_X, ...
                                    param_Y,param_Z, 0, [], ...
                                    [], 1, 0);
                        end
                    end
                end

                % calculate I(Idot;param(shuffle)) four times each shuffle
                nsh=1;
                Ip2sig_dist_shuffle=zeros(49,49,49,nsh);
                for ss=1:nsh
                    for aa=1:49
                        for bb=1:49
                            for cc=1:49
                            Itot_X=Itot_sigma_dot(aa,win);
                            param_Y=param_sigma(bb,win);
                            param_Z=param2_sigma(cc,win);
                            Ip2sig_dist_shuffle(aa,bb,cc,ss)=ecmi_with_shuffle(Itot_X, ...
                                    param_Y,param_Z(randperm(numel(win))), 0, [], ...
                                    [], 1, 0);
                            end
                        end
                    end
                end
                
                [Ip2sig_shuffle_freq, ip2sig_shuffle_sup]=histcounts(Ip2sig_dist_shuffle(:), 0:.01:1);
                ip2sig_shuffle_sup=ip2sig_shuffle_sup(1:end-1);
                pIp2sig_shuffle=Ip2sig_shuffle_freq/sum(Ip2sig_shuffle_freq);
        
                
                if debug
                    jj=jj+1;
                    figure(ii); gcf;
                    subplot(numel(rise_params)-1,3,jj);
                    area(ip2sig_shuffle_sup,pIp2sig_shuffle,  'facecolor', 'k', 'facealpha', 0.1);
                    
                    % two-tailed test visual with alpha
                    cpI=cumsum(pIp2sig_shuffle);
                    idxu=find(cpI>1-alpha/2);
                    idxl=find(cpI<alpha/2);
                    
                    hold on;
                    h1=plot([mean(Ip2sig_dist(:)),mean(Ip2sig_dist(:))], [0,1], 'r-', 'linewidth', 2);
                    plot([ip2sig_shuffle_sup(idxu(1)),ip2sig_shuffle_sup(idxu(1))], [0,1], 'k:', 'linewidth', 2);
                    if ~isempty(idxl)
                        plot([ip2sig_shuffle_sup(idxl(end)),ip2sig_shuffle_sup(idxl(end))], [0,1], 'k:', 'linewidth', 2);
                    end
                    xlabel(['I($', '\dot{I}',';', param_label , '|' , param2_label, '$) (bits)'], 'interpreter', 'latex');
                    ylabel('p');
                    if jj<4
                    title(sprintf('rise %d', rr));
                    end
                    set(gca, 'fontsize', 16);
                    set(gca,'Xscale','log')
                    set(gca, 'fontsize', 16, 'ylim', [0 0.5]);
                    %                     set(gca, 'xlim', [0, 0.6]);
                end
            end
        end
    end
end


% close all
%% falls
for ii=1:numel(fall_params)
    if debug
        figure(ii+4); gcf; clf;
        jj=0;
    end
    param=xV(fall_params(ii),1:end-1);
    param_sigma=squeeze(sigmaPointAccumulutor(fall_params(ii),:,1:end-1));
    param_label=fall_params_labels{ii};
    for kk=1:numel(fall_params)
        if ii ~=kk
            param2=xV(fall_params(kk),1:end-1);
            param2_sigma=squeeze(sigmaPointAccumulutor(fall_params(kk),:,1:end-1));
            param2_label=fall_params_labels{kk};
            
            for ff=1:size(fall,1)
                win=fall(ff,1):fall(ff,2);
                
                % calculate I(Idot(sigma); param(sigma))
                Ip2sig_dist=zeros(49,49,49);
                for aa=1:49
                    for bb=1:49
                        for cc=1:49
                            Itot_X=Itot_sigma_dot(aa,win);
                            param_Y=param_sigma(bb,win);
                            param_Z=param2_sigma(cc,win);
                            Ip2sig_dist(aa,bb,cc)=ecmi_with_shuffle(Itot_X, ...
                                    param_Y,param_Z, 0, [], ...
                                    [], 1, 0);
                        end
                    end
                end
 
                % calculate I(Idot;param(shuffle)) four times each shuffle
                nsh=1;
                Ip2sig_dist_shuffle=zeros(49,49,49,nsh);
                for ss=1:nsh
                    for aa=1:49
                        for bb=1:49
                            for cc=1:49
                            Itot_X=Itot_sigma_dot(aa,win);
                            param_Y=param_sigma(bb,win);
                            param_Z=param2_sigma(cc,win);
                            Ip2sig_dist_shuffle(aa,bb,cc,ss)=ecmi_with_shuffle(Itot_X, ...
                                    param_Y,param_Z(randperm(numel(win))), 0, [], ...
                                    [], 1, 0);
                            end
                        end
                    end
                end
                
                [Ip2sig_shuffle_freq, ip2sig_shuffle_sup]=histcounts(Ip2sig_dist_shuffle(:), 0:.01:1);
                ip2sig_shuffle_sup=ip2sig_shuffle_sup(1:end-1);
                pIp2sig_shuffle=Ip2sig_shuffle_freq/sum(Ip2sig_shuffle_freq);
                
                if debug
                    jj=jj+1;
                    figure(ii+4); gcf;
                    subplot(numel(fall_params)-1,3,jj);
                    area(ip2sig_shuffle_sup,pIp2sig_shuffle,  'facecolor', 'k', 'facealpha', 0.1);
                    
                    % one-tailed test visual with alpha
                    cpI=cumsum(pIp2sig_shuffle);
                    idxu=find(cpI>1-alpha/2);
                    idxl=find(cpI<alpha/2);
                    
                    hold on;
                    h1=plot([mean(Ip2sig_dist(:)),mean(Ip2sig_dist(:))], [0,1], 'r-', 'linewidth', 2);
                    plot([ip2sig_shuffle_sup(idxu(1)),ip2sig_shuffle_sup(idxu(1))], [0,1], 'k:', 'linewidth', 2);

                    if ~isempty(idxl)
                        plot([ip2sig_shuffle_sup(idxl(end)),ip2sig_shuffle_sup(idxl(end))], [0,1], 'k:', 'linewidth', 2);
                    end
                    xlabel(['I($', '\dot{I}',';', param_label , '|' , param2_label, '$) (bits)'], 'interpreter', 'latex');
                    ylabel('p');
                    if jj<4
                    title(sprintf('fall %d', ff));
                    end
                    set(gca, 'fontsize', 16);
                    set(gca,'Xscale','log')
                    set(gca, 'fontsize', 16, 'ylim', [0 0.5]);
                end
            end
        end
    end
end





