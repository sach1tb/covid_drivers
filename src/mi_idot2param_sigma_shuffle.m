clearvars

addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
addpath(['cteUpdatedFiles', filesep])

load ukfOutput.mat  %size is 24, 24*2+1

% debug set to 1 for testing
debug=1;

% converting to dates
startDate = datenum('02-04-2020');
endDate = datenum('11-01-2022');
dateData = linspace(startDate,endDate,T);

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
fall_params_labels={'\xi_1', '\phi_1', '\alpha',};
fall_params=[15,18,17];

nsh=20;

alpha=0.01; % for one-tailed test

Itot = sum(xV(7:9,:),1);
Itotdot=diff(Itot);


Itot_sigma=squeeze(sum(sigmaPointAccumulutor(7:9,:,:),1));
Itot_sigma_dot=diff(Itot_sigma,1,2);

% [  beta, xi1, xi2, alpha, ...
%             phi1, phi2, sigma0, kappa0, ...
%             mu, gamma, epsilon]=vec2params(vec)

%% rises
figure(1); gcf; clf;
for ii=1:numel(rise_params)
    param=xV(rise_params(ii),1:end-1);
    param_sigma=squeeze(sigmaPointAccumulutor(rise_params(ii),:,1:end-1));
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
        
        
        % calculate I(Idot(sigma); param(sigma))
        Ip2sig_dist=zeros(49,49);
        for aa=1:49
            for bb=1:49
                Itot_X=Itot_sigma_dot(aa,win);
                param_Y=param_sigma(bb,win);
                Ip2sig_dist(aa,bb)=emi_with_shuffle(Itot_X, param_Y, 0, [], [], 1, 0);
            end
        end
        
        % calculate I(Idot;param(shuffle)) four times each shuffle
        
        Ip2sig_dist_shuffle=zeros(49,49,nsh);
        for ss=1:nsh
            for aa=1:49
                for bb=1:49
                    Itot_X=Itot_sigma_dot(aa,win);
                    param_Y=param_sigma(bb,win);
                    Ip2sig_dist_shuffle(aa,bb,ss)=emi_with_shuffle(Itot_X, param_Y(randperm(numel(win))), 0, [], [], 1, 0);
                end
            end
        end
        
        [Ip2sig_shuffle_freq, ip2sig_shuffle_sup]=histcounts(Ip2sig_dist_shuffle(:), 0:.01:1);
        mi_data_rise(ii,rr).ip2sig_shuffle_sup=ip2sig_shuffle_sup(1:end-1);
        mi_data_rise(ii,rr).pIp2sig_shuffle=Ip2sig_shuffle_freq/sum(Ip2sig_shuffle_freq);
        mi_data_rise(ii,rr).Ip2sig_dist_mean=mean(Ip2sig_dist(:));
        
    end
end
save('mi_rise.mat', 'mi_data_rise');

%% plot
load('mi_rise.mat');

figure(2); gcf; clf; jj=0;

for ii=1:numel(rise_params)
    param_label=rise_params_labels{ii};
    for rr=1:size(rise,1)
        jj=jj+1;
        figure(2); gcf;
        subplot(4,3,jj);
        area(mi_data_rise(ii,rr).ip2sig_shuffle_sup, ...
            mi_data_rise(ii,rr).pIp2sig_shuffle, 'facecolor', 'k', 'facealpha', 0.1);
        hold on;
        % one-tailed test visual with alpha
        cpI=cumsum(mi_data_rise(ii,rr).pIp2sig_shuffle);
        idx=find(cpI>1-alpha);
        
        [val, p_idx]=min(sqrt((mi_data_rise(ii,rr).ip2sig_shuffle_sup-mi_data_rise(ii,rr).Ip2sig_dist_mean).^2));
        if val < 0.01
            pval=1-cpI(p_idx);
        else
            error('fix support resolution');
        end
        %             pval=ranksum(mean(Ip2sig_dist(:)),Ip2sig_dist_shuffle(:),'tail', 'right', 'alpha', 0.05);
        h1=plot([mi_data_rise(ii,rr).Ip2sig_dist_mean,...
            mi_data_rise(ii,rr).Ip2sig_dist_mean], [0,1], 'r-', 'linewidth', 2);
        plot([mi_data_rise(ii,rr).ip2sig_shuffle_sup(idx(1)),...
            mi_data_rise(ii,rr).ip2sig_shuffle_sup(idx(1))], [0,1], 'k:', 'linewidth', 2);
        xlabel(['$\log(\mathcal{I}(\dot{I};' param_label '))$ (bits)'], 'interpreter', 'latex');
        ylabel('p');
        if jj<4
            title(sprintf('%s to %s', ...
                datestr(dateData(rise(rr,1))), datestr(dateData(rise(rr,2)))), ...
                'fontweight', 'normal');
        end
        fprintf('Rise (%s to %s), ${I}(\\dot{I};%s)$,%.3f(%.3f), p=%.4f\n',...
            datestr(dateData(rise(rr,1))), datestr(dateData(rise(rr,2))), ...
            param_label,mi_data_rise(ii,rr).Ip2sig_dist_mean, ...
            mi_data_rise(ii,rr).ip2sig_shuffle_sup(idx(1)), pval);
        %             title(sprintf('p=%.5f', pval)); % comment after checking
        set(gca, 'fontsize', 20, 'xlim', [0 1]);
        set(gca,'Xscale','log')
        set(gca, 'fontsize', 20, 'ylim', [0 0.5]);
        if pval < 0.01
            legend(h1,'*');
        end
    end
end


%% falls
figure(3); gcf; clf;
for ii=1:numel(fall_params)
    param=xV(fall_params(ii),1:end-1);
    param_sigma=squeeze(sigmaPointAccumulutor(fall_params(ii),:,1:end-1));
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
        % calculate I(Idot(sigma); param(sigma))
        Ip2sig_dist=zeros(49,49);
        for aa=1:49
            for bb=1:49
                Itot_X=Itot_sigma_dot(aa,win);
                param_Y=param_sigma(bb,win);
                Ip2sig_dist(aa,bb)=emi_with_shuffle(Itot_X, param_Y, 0, [], [], 1, 0);
            end
        end
        
        % calculate I(Idot(shuffle);param(shuffle)) four times each shuffle
        Ip2sig_dist_shuffle=zeros(49,49,nsh);
        for ss=1:nsh
            for aa=1:49
                for bb=1:49
                    Itot_X=Itot_sigma_dot(aa,win);
                    param_Y=param_sigma(bb,win);
                    Ip2sig_dist_shuffle(aa,bb,ss)=emi_with_shuffle(Itot_X, param_Y(randperm(numel(win))), 0, [], [], 1, 0);
                end
            end
        end
        
        [Ip2sig_shuffle_freq, ip2sig_shuffle_sup]=histcounts(Ip2sig_dist_shuffle(:), 0:.01:1);
        mi_data_fall(ii,ff).ip2sig_shuffle_sup=ip2sig_shuffle_sup(1:end-1);
        mi_data_fall(ii,ff).pIp2sig_shuffle=Ip2sig_shuffle_freq/sum(Ip2sig_shuffle_freq);
        mi_data_fall(ii,ff).Ip2sig_dist_mean=mean(Ip2sig_dist(:));
        
    end
end

save('mi_fall.mat', 'mi_data_fall');

%% plot
load('mi_fall.mat');


figure(4); gcf; clf; jj=0;


for ii=1:numel(fall_params)
    param_label=fall_params_labels{ii};
    for ff=1:size(fall,1)
        jj=jj+1;
        figure(4); gcf;
        subplot(3,3,jj);
        area(mi_data_fall(ii,ff).ip2sig_shuffle_sup, ...
            mi_data_fall(ii,ff).pIp2sig_shuffle, 'facecolor', 'k', 'facealpha', 0.1);
        hold on;
        % one-tailed test visual with alpha
        cpI=cumsum(mi_data_fall(ii,ff).pIp2sig_shuffle);
        idx=find(cpI>1-alpha);
        
        [val, p_idx]=min(sqrt((ip2sig_shuffle_sup-mi_data_fall(ii,ff).Ip2sig_dist_mean).^2));
        if val < 0.01
            pval=1-cpI(p_idx);
        else
            error('fix support resolution');
        end
        %             pval=ranksum(mean(Ip2sig_dist(:)),Ip2sig_dist_shuffle(:),'tail', 'right', 'alpha', 0.05);
        h2=plot([mi_data_fall(ii,ff).Ip2sig_dist_mean,...
            mi_data_fall(ii,ff).Ip2sig_dist_mean], [0,1], 'g-', 'linewidth', 2);
        plot([mi_data_fall(ii,ff).ip2sig_shuffle_sup(idx(1)),...
            mi_data_fall(ii,ff).ip2sig_shuffle_sup(idx(1))], [0,1], 'k:', 'linewidth', 2);
        xlabel(['$\log(\mathcal{I}(\dot{I};' param_label '))$ (bits)'], 'interpreter', 'latex');
        ylabel('p');
        if jj<4
            title(sprintf('%s to %s', ...
                datestr(dateData(fall(ff,1))), datestr(dateData(fall(ff,2)))), ...
                'fontweight', 'normal');
        end
        fprintf('Fall (%s to %s), ${I}(\\dot{I};%s)$,%.3f(%.3f), p=%.4f\n',...
            datestr(dateData(fall(ff,1))), datestr(dateData(fall(ff,2))), ...
            param_label,mi_data_fall(ii,ff).Ip2sig_dist_mean, ...
            mi_data_fall(ii,ff).ip2sig_shuffle_sup(idx(1)), pval);
        %             title(sprintf('p=%.5f', pval)); % comment after checking
        set(gca, 'fontsize', 20, 'xlim', [0 1]);
        set(gca, 'fontsize', 20, 'ylim', [0 0.5]);
        set(gca,'Xscale','log');
        if pval < 0.01
            legend(h2,'*');
        end
    end
    
end


