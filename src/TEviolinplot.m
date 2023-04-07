clearvars
close all
load('allTECal_win84.mat');

figure(1)

%rise
rise = [201, 279;
    627, 694;
    775, 891]-round(windowSizeDays/2);

%fall
fall = [280, 376;
    694, 761;
    892, 958 ]-round(windowSizeDays/2);

% rise wave 1
figure(1)
temp = mean(NetTE_phi2_Itot)';
aggregate = [];
aggregate = temp(rise(1,1):rise(1,2));
name = repmat(["\phi_2"],[numel(rise(1,1):rise(1,2)),1]);

temp = mean(NetTE_xi1_Itot)';
aggregate = [aggregate; temp(rise(1,1):rise(1,2))];
name = [name; repmat(["\xi_1"],[numel(rise(1,1):rise(1,2)),1]) ];

temp = mean(NetTE_sigma_Itot)';
aggregate = [aggregate; temp(rise(1,1):rise(1,2))];
name = [name; repmat(["\sigma"],[numel(rise(1,1):rise(1,2)),1]) ];
subplot(2,2,1)
vs = violinplot(aggregate, name);
plotdetails(1,1);

%% rise wave 2
temp = mean(NetTE_phi2_Itot)';
aggregate = [];
name = [];
aggregate = temp(rise(2,1):rise(2,2));
name = repmat(["\phi_2"],[numel(rise(2,1):rise(2,2)),1]);

temp = mean(NetTE_xi1_Itot)';
aggregate = [aggregate; temp(rise(2,1):rise(2,2))];
name = [name; repmat(["\xi_1"],[numel(rise(2,1):rise(2,2)),1]) ];

temp = mean(NetTE_sigma_Itot)';
aggregate = [aggregate; temp(rise(2,1):rise(2,2))];
name = [name; repmat(["\sigma"],[numel(rise(2,1):rise(2,2)),1]) ];
subplot(2,2,2)
vs = violinplot(aggregate, name);
plotdetails(2,1);

%% rise wave 3
temp = mean(NetTE_phi2_Itot)';
aggregate = [];
name = [];
aggregate = temp(rise(3,1):rise(3,2));
name = repmat(["\phi_2"],[numel(rise(3,1):rise(3,2)),1]);

temp = mean(NetTE_xi1_Itot)';
aggregate = [aggregate; temp(rise(3,1):rise(3,2))];
name = [name; repmat(["\xi_1"],[numel(rise(3,1):rise(3,2)),1]) ];

temp = mean(NetTE_sigma_Itot)';
aggregate = [aggregate; temp(rise(3,1):rise(3,2))];
name = [name; repmat(["\sigma"],[numel(rise(3,1):rise(3,2)),1]) ];
subplot(2,2,3)
vs = violinplot(aggregate, name);
plotdetails(3,1);

%% fall wave 1

figure(2)
temp = mean(NetTE_phi1_Itot)';
aggregate = [];
name = [];
aggregate = temp(fall(2,1):fall(1,2));
name = repmat(["\phi_1"],[numel(fall(1,1):fall(1,2)),1]);

temp = mean(NetTE_xi2_Itot)';
aggregate = [aggregate; temp(fall(1,1):fall(1,2))];
name = [name; repmat(["\xi_2"],[numel(fall(1,1):fall(1,2)),1]) ];

temp = mean(NetTE_alpha_Itot)';
aggregate = [aggregate; temp(fall(1,1):fall(1,2))];
name = [name; repmat(["\alpha"],[numel(fall(1,1):fall(1,2)),1]) ];
subplot(2,2,1)
vs = violinplot(aggregate, name);
plotdetails(1,2);

%% fall wave 2
temp = mean(NetTE_phi1_Itot)';
aggregate = [];
name = [];
aggregate = temp(fall(2,1):fall(2,2));
name = repmat(["\phi_1"],[numel(fall(2,1):fall(2,2)),1]);

temp = mean(NetTE_xi2_Itot)';
aggregate = [aggregate; temp(fall(2,1):fall(2,2))];
name = [name; repmat(["\xi_2"],[numel(fall(2,1):fall(2,2)),1]) ];

temp = mean(NetTE_alpha_Itot)';
aggregate = [aggregate; temp(fall(2,1):fall(2,2))];
name = [name; repmat(["\alpha"],[numel(fall(2,1):fall(2,2)),1]) ];
subplot(2,2,2)
vs = violinplot(aggregate, name);
plotdetails(2,2);

%% fall wave 3
temp = mean(NetTE_phi1_Itot)';
aggregate = [];
name = [];
aggregate = temp(fall(3,1):fall(3,2));
name = repmat(["\phi_1"],[numel(fall(3,1):fall(3,2)),1]);

temp = mean(NetTE_xi2_Itot)';
aggregate = [aggregate; temp(fall(3,1):fall(3,2))];
name = [name; repmat(["\xi_2"],[numel(fall(3,1):fall(3,2)),1]) ];

temp = mean(NetTE_alpha_Itot)';
aggregate = [aggregate; temp(fall(3,1):fall(3,2))];
name = [name; repmat(["\alpha"],[numel(fall(3,1):fall(3,2)),1]) ];
subplot(2,2,3)
vs = violinplot(aggregate, name);
plotdetails(3,2);


function plotdetails(n,m)
if (m == 1)
    title(sprintf('Rise wave %d',round(n)));
end
if (m==2)
    title(sprintf('Fall wave %d',round(n)));
end
%ylabel('Fuel Economy in MPG ');
xlim([-4, 8]);
grid on;
set(gca, 'color', 'none');
xtickangle(0);
%fprintf('Test %02.0f passed ok! \n ',n);
end

