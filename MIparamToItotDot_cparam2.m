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
rise_params_labels={'\xi_1', '\phi_2', '\sigma', '\kappa'};
rise_params=[15,19,20,21];

%fall
fall = [280, 376;
    694, 761;
    892, 958 ];
fall_params_labels={'\xi_2', '\alpha','\phi_1'};
fall_params=[16,17,18];

nshuffle=10000;

alpha=0.01; % for one-tailed test

Itot = sum(xV(7:9,:),1);
Itotdot=diff(Itot);

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
    param_label=rise_params_labels{ii};
    for kk=1:numel(rise_params)
        if ii ~=kk
            param2=xV(rise_params(kk),1:end-1);
            param2_label=rise_params_labels{kk};
            
            for rr=1:size(rise,1)
                win=rise(rr,1):rise(rr,2);
                [Ip2Itotdotcparam, pIshuffle, Isup, Ip2shufflecparam]=ecmi_with_shuffle(param(win),Itotdot(win), param2(win), nshuffle, [], ...
                    [], 1, 0);
                
                if debug
                    jj=jj+1;
                    figure(ii); gcf;
                    subplot(numel(rise_params)-1,3,jj);
                    plot(Isup,pIshuffle, 'k', 'linewidth', 2);
                    
                    % one-tailed test visual with alpha
                    cpI=cumsum(pIshuffle);
                    idx=find(cpI>1-alpha);
                    
                    hold on;
                    plot([Ip2Itotdotcparam,Ip2Itotdotcparam], [0,1], 'r-', 'linewidth', 2);
                    plot([Isup(idx(1)),Isup(idx(1))], [0,1], 'k:', 'linewidth', 2);
                    xlabel(['I($' param_label ';\dot{I}|', param2_label, '$) (bits)'], 'interpreter', 'latex');
                    ylabel('p');
                    if jj<4
                    title(sprintf('rise %d', rr));
                    end
                    set(gca, 'fontsize', 16);
                    %                     set(gca, 'xlim', [0, 0.6]);
                end
            end
        end
    end
end


%% falls
for ii=1:numel(fall_params)
    if debug
        figure(ii); gcf; clf;
        jj=0;
    end
    param=xV(fall_params(ii),1:end-1);
    param_label=fall_params_labels{ii};
    for kk=1:numel(fall_params)
        if ii ~=kk
            param2=xV(fall_params(kk),1:end-1);
            param2_label=fall_params_labels{kk};
            
            for ff=1:size(fall,1)
                win=fall(ff,1):fall(ff,2);
                [Ip2Itotdotcparam, pIshuffle, Isup, Ip2shufflecparam]=ecmi_with_shuffle(param(win),Itotdot(win), param2(win), nshuffle, [], ...
                    [], 1, 0);
                
                if debug
                    jj=jj+1;
                    figure(ii); gcf;
                    subplot(numel(fall_params)-1,3,jj);
                    plot(Isup,pIshuffle, 'k', 'linewidth', 2);
                    
                    % one-tailed test visual with alpha
                    cpI=cumsum(pIshuffle);
                    idx=find(cpI>1-alpha);
                    
                    hold on;
                    plot([Ip2Itotdotcparam,Ip2Itotdotcparam], [0,1], 'g-', 'linewidth', 2);
                    plot([Isup(idx(1)),Isup(idx(1))], [0,1], 'k:', 'linewidth', 2);
                    xlabel(['I($' param_label ';\dot{I}|', param2_label, '$) (bits)'], 'interpreter', 'latex');
                    ylabel('p');
                    if jj<4
                    title(sprintf('fall %d', ff));
                    end
                    set(gca, 'fontsize', 16);
                    %                     set(gca, 'xlim', [0, 0.6]);
                end
            end
        end
    end
end





