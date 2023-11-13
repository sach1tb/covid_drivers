clearvars

addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
addpath(['cteUpdatedFiles', filesep])

load ukfOutput.mat  %size is 24, 24*2+1

% debug set to 1 for testing
debug=1;

%rise
rise = [203, 272;
        616, 691;
        772, 887];
rise_params_labels={'\xi_2', '\phi_2', '\sigma', '\kappa'};
rise_params=[16,19,20,21];

%fall
fall = [272, 383;
        691, 771;
        887, 958 ];
fall_params_labels={'\xi_1', '\phi_1', '\alpha',};
fall_params=[15,18,17];

alpha=0.01; % for one-tailed test

Itot_sigma=squeeze(sum(sigmaPointAccumulutor(7:9,:,:),1));
Itot_sigma_dot=diff(Itot_sigma,1,2);

% [  beta, xi1, xi2, alpha, ...
%             phi1, phi2, sigma0, kappa0, ...
%             mu, gamma, epsilon]=vec2params(vec)

%% rises
for ii=1:numel(rise_params) % param1
    param=xV(rise_params(ii),1:end-1);
    param_sigma=squeeze(sigmaPointAccumulutor(rise_params(ii),:,1:end-1));
    
    for kk=1:numel(rise_params) % param2
        if ii ~=kk
            param2=xV(rise_params(kk),1:end-1);
            param2_sigma=squeeze(sigmaPointAccumulutor(rise_params(kk),:,1:end-1));
            
            
            for rr=1:size(rise,1) % window
                win=rise(rr,1):rise(rr,2);
                
                
                % calculate I(Idot(sigma); param(sigma))
                CMIp2sig_dist=zeros(49,49,49);
                for aa=1:49
                    for bb=1:49
                        for cc=1:49
                            Itot_X=Itot_sigma_dot(aa,win);
                            param_Y=param_sigma(bb,win);
                            param_Z=param2_sigma(cc,win);
                            CMIp2sig_dist(aa,bb,cc)=ecmi_with_shuffle(Itot_X, ...
                                param_Y,param_Z, 0, [], ...
                                [], 1, 0);
                        end
                    end
                end
                
                % calculate I(Idot;param(shuffle)) four times each shuffle
                nsh=1;
                CMIp2sig_dist_shuffle=zeros(49,49,49,nsh);
                for ss=1:nsh
                    for aa=1:49
                        for bb=1:49
                            for cc=1:49
                                Itot_X=Itot_sigma_dot(aa,win);
                                param_Y=param_sigma(bb,win);
                                param_Z=param2_sigma(cc,win);
                                CMIp2sig_dist_shuffle(aa,bb,cc,ss)=ecmi_with_shuffle(Itot_X, ...
                                    param_Y(randperm(numel(win))),param_Z, 0, [], ...
                                    [], 1, 0);
                            end
                        end
                    end
                end
                
                [CMIp2sig_shuffle_freq, cmi_p2sig_shuffle_sup]=histcounts(CMIp2sig_dist_shuffle(:), 0:.01:1);
                cmi_data(ii,kk,rr).cmi_p2sig_shuffle_sup=cmi_p2sig_shuffle_sup(1:end-1);
                cmi_data(ii,kk,rr).pCMIp2sig_shuffle=CMIp2sig_shuffle_freq/sum(CMIp2sig_shuffle_freq);
                cmi_data(ii,kk,rr).CMIp2sig_dist_mean=mean(CMIp2sig_dist(:));
            end
        end
    end
end
save('cmi_rise.mat', 'cmi_data');

%% plot
load('cmi_rise.mat');
for ii=1:numel(rise_params)
    figure(ii); gcf; clf;
    jj=0;
    param_label=rise_params_labels{ii};
    for kk=1:numel(rise_params)
        param2_label=rise_params_labels{kk};
        if kk~=ii
        for rr=1:size(rise,1)
            jj=jj+1;
            figure(ii); gcf;
            subplot(numel(rise_params)-1,3,jj);
            area(cmi_data(ii,kk,rr).cmi_p2sig_shuffle_sup,cmi_data(ii,kk,rr).pCMIp2sig_shuffle,...
                'facecolor', 'k', 'facealpha', 0.1);
            
            % two-tailed test visual with alpha
            cpI=cumsum(cmi_data(ii,kk,rr).pCMIp2sig_shuffle);
            idxu=find(cpI>1-alpha);
            %                     idxu=find(cpI>1-alpha/2);
            %                     idxl=find(cpI<alpha/2);
            
            hold on;
            h1=plot([cmi_data(ii,kk,rr).CMIp2sig_dist_mean,cmi_data(ii,kk,rr).CMIp2sig_dist_mean], [0,1], 'r-', 'linewidth', 2);
            plot([cmi_data(ii,kk,rr).cmi_p2sig_shuffle_sup(idxu(1)),cmi_data(ii,kk,rr).cmi_p2sig_shuffle_sup(idxu(1))], [0,1], 'k:', 'linewidth', 2);
            %                     if ~isempty(idxl)
            %                         plot([cmi_p2sig_shuffle_sup(idxl(end)),cmi_p2sig_shuffle_sup(idxl(end))], [0,1], 'k:', 'linewidth', 2);
            %                     end
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
% close all
%% falls
for ii=1:numel(fall_params)
    param=xV(fall_params(ii),1:end-1);
    param_sigma=squeeze(sigmaPointAccumulutor(fall_params(ii),:,1:end-1));
    
    for kk=1:numel(fall_params)
        if ii ~=kk
            param2=xV(fall_params(kk),1:end-1);
            param2_sigma=squeeze(sigmaPointAccumulutor(fall_params(kk),:,1:end-1));
            
            
            for ff=1:size(fall,1)
                win=fall(ff,1):fall(ff,2);
                
                % calculate I(Idot(sigma); param(sigma))
                CMIp2sig_dist=zeros(49,49,49);
                for aa=1:49
                    for bb=1:49
                        for cc=1:49
                            Itot_X=Itot_sigma_dot(aa,win);
                            param_Y=param_sigma(bb,win);
                            param_Z=param2_sigma(cc,win);
                            CMIp2sig_dist(aa,bb,cc)=ecmi_with_shuffle(Itot_X, ...
                                param_Y,param_Z, 0, [], ...
                                [], 1, 0);
                        end
                    end
                end
                
                % calculate I(Idot;param(shuffle)) four times each shuffle
                nsh=1;
                CMIp2sig_dist_shuffle=zeros(49,49,49,nsh);
                for ss=1:nsh
                    for aa=1:49
                        for bb=1:49
                            for cc=1:49
                                Itot_X=Itot_sigma_dot(aa,win);
                                param_Y=param_sigma(bb,win);
                                param_Z=param2_sigma(cc,win);
                                CMIp2sig_dist_shuffle(aa,bb,cc,ss)=ecmi_with_shuffle(Itot_X, ...
                                    param_Y(randperm(numel(win))),param_Z, 0, [], ...
                                    [], 1, 0);
                            end
                        end
                    end
                end
                
                [CMIp2sig_shuffle_freq, cmi_p2sig_shuffle_sup]=histcounts(CMIp2sig_dist_shuffle(:), 0:.01:1);
                cmi_data_fall(ii,kk,ff).cmi_p2sig_shuffle_sup=cmi_p2sig_shuffle_sup(1:end-1);
                cmi_data_fall(ii,kk,ff).pCMIp2sig_shuffle=CMIp2sig_shuffle_freq/sum(CMIp2sig_shuffle_freq);
                cmi_data_fall(ii,kk,ff).CMIp2sig_dist_mean=mean(CMIp2sig_dist(:));
                
            end
        end
    end
end
save('cmi_fall.mat', 'cmi_data_fall');

%% plot
load('cmi_fall.mat');
for ii=1:numel(fall_params)
    figure(ii+4); gcf; clf;
    jj=0;
    param_label=fall_params_labels{ii};
    for kk=1:numel(fall_params)
        param2_label=fall_params_labels{kk};
        if kk~=ii
        for ff=1:size(fall,1)
            jj=jj+1;
            figure(ii+4); gcf;
            subplot(numel(fall_params)-1,3,jj);
            area(cmi_data_fall(ii,kk,ff).cmi_p2sig_shuffle_sup,...
                cmi_data_fall(ii,kk,ff).pCMIp2sig_shuffle,...
                'facecolor', 'k', 'facealpha', 0.1);
            
            % one-tailed test visual with alpha
            cpI=cumsum(cmi_data_fall(ii,kk,ff).pCMIp2sig_shuffle);
            idxu=find(cpI>1-alpha);
            %                     idxu=find(cpI>1-alpha/2);
            %                     idxl=find(cpI<alpha/2);
            
            hold on;
            h1=plot([cmi_data_fall(ii,kk,ff).CMIp2sig_dist_mean,...
                        cmi_data_fall(ii,kk,ff).CMIp2sig_dist_mean], [0,1], 'g-', 'linewidth', 2);
            plot([cmi_data_fall(ii,kk,ff).cmi_p2sig_shuffle_sup(idxu(1)),...
                cmi_data_fall(ii,kk,ff).cmi_p2sig_shuffle_sup(idxu(1))], [0,1], 'k:', 'linewidth', 2);
            
            xlabel(['I($', '\dot{I}',';', param_label , '|' , param2_label, '$) (bits)'], 'interpreter', 'latex');
            ylabel('p');
            if jj<  4
                title(sprintf('fall %d', ff));
            end
            set(gca, 'fontsize', 16);
            set(gca,'Xscale','log')
            set(gca, 'fontsize', 16, 'ylim', [0 0.5]);
        end
        end
    end
end



