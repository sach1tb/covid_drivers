clearvars

load ukfOutput.mat  %size is 24, 24*2+1

% debug set to 1 for testing
debug=1;

% converting to dates
startDate = datenum('02-04-2020');
endDate = datenum('11-01-2022');
dateData = linspace(startDate,endDate,T);

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

nsh=20;

alpha=0.01; % for one-tailed test

Itot = sum(xV(7:9,:),1);
Itotdot=diff(Itot);

figure(1); gcf;clf;
subplot(1,2,1);
plot(Itot, 'r');
hold on;
for ii=1:size(rise,1)
    plot([rise(ii,2), rise(ii,2)], [min(Itot), max(Itot)], 'k');
end
grid on;
ylabel('I');

subplot(1,2,2);
plot(Itotdot, 'r');
hold on;
for ii=1:size(rise,1)
    plot([rise(ii,2), rise(ii,2)], [min(Itotdot), max(Itotdot)], 'k');
end
grid on;
ylabel('Idot');




Itot_sigma=squeeze(sum(sigmaPointAccumulutor(7:9,:,:),1));
Itot_sigma_dot=diff(Itot_sigma,1,2);


z = [infectious;death;vax;mask;mobility;popDays];


figure(2); gcf; clf; jj=0;
% for rr=1:size(rise,1)
    for jj=1:size(z,1)
        for kk=jj:size(z,1)
%             win=rise(rr,1):fall(rr,2);
            win=1:size(z,2);
            figure(2); gcf;
            subplot(4,3,jj);
            [rcorr, pcorr]=corrcoef(z(jj,win), z(kk,win));
%             corrmat(jj,kk, rr)=rcorr(1,2);
%             pcorrmat(jj,kk, rr)=pcorr(1,2);
            corrmat(jj,kk)=rcorr(1,2);
            pcorrmat(jj,kk)=pcorr(1,2);

        end
    end
% end

for jj=1:size(z,1)
    for kk=1:size(z,1)
        fprintf('%.3f (%.3f),', corrmat(jj,kk), pcorrmat(jj,kk));
    end
    fprintf('\n');
end
