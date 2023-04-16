clearvars

% data obtained as follows:
% > go to https://www.healthdata.org/node/7425, 
% > click on each reference scenario
% There are total 1062 samples in
% data_download_file_reference_YYYY_ILLINOIS.csv (newer set) but it doesnt contain data
% for mobility and mask use. So I created an augmented dataset with
% mobility and mask use taken from reference_YYYY.csv(older set) files.
 
data2020 = readtable('../data/data_download_file_reference_2020_ILLINOIS_augmented.csv');
data2021 = readtable('../data/data_download_file_reference_2021_ILLINOIS_augmented.csv');
data2022 = readtable('../data/data_download_file_reference_2022_ILLINOIS_augmented.csv');



allData = [data2020;data2021;data2022];
% idx = allData.location_id == 536;
ILdata = allData; % Dataset for illinois state

Infection = ILdata.inf_mean; % inf_mean
Vaccinated = ILdata.cumulative_all_fully_vaccinated; % cumulative_all_fully_vaccinated
Vaccinated(isnan(Vaccinated)) = 0;
Masked = ILdata.mask_use_mean;  % Taken from reference2020, 2021, 2022.csv mask_use_mean
Masked(isnan(Masked)) = 0;
Death = ILdata.daily_deaths; % daily_deaths
Mobility = ILdata.mobility_mean; % Taken from reference2020, 2021, 2022.csv, mobility_mean
Mobility(isnan(Mobility)) = 0;

dayStops = [1 332 697 1002];
popChicagoMetro = [9684738 9601605 9509934 9433330];
days = 1:1:1002;
popDays = interp1(dayStops,popChicagoMetro,days);


% plot & verify

figure(1); gcf;clf;
subplot(2,3,1);
h1 = plot(Infection);
hold on;
inf_data_current=csvread('../data/infectiousIllinois_ci.csv');
h2 = plot(inf_data_current(:,2));
legend([h1,h2],"New","Old")
title("Infection")

subplot(2,3,2);
h1 = plot(cumsum(Death));
hold on;
inf_data_current=csvread('../data/deathIllinois.csv');
h2 = plot(cumsum(inf_data_current));
legend([h1,h2],"New","Old")
title("Cumulative Deaths")


subplot(2,3,3);
h1 = plot(Vaccinated);
hold on;
inf_data_current=csvread('../data/vaccinatedIllinois.csv');
inf_data_current(isnan(inf_data_current)) = 0;
h2 = plot(inf_data_current);
legend([h1,h2],"New","Old")
title("Vaccination")


subplot(2,3,4);
h1 = plot(Mobility);
hold on;
inf_data_current=csvread('../data/mobilityIllinois.csv');
inf_data_current(isnan(inf_data_current)) = 0;
h2 = plot(inf_data_current);
legend([h1,h2],"New","Old")
title("Mobility")


subplot(2,3,5);
h1 = plot(Masked);
hold on;
inf_data_current=csvread('../data/maskIllinois.csv');
inf_data_current(isnan(inf_data_current)) = 0;
h2 = plot(inf_data_current);
legend([h1,h2],"New","Old")
title("Mask use")


% writematrix(Infection,'infectious.csv');
% writematrix(Vaccinated,'vaccinated.csv');
% writematrix(Masked,'mask.csv');
% writematrix(Death,'death.csv');
% writematrix(Mobility,'mobility.csv');
%writematrix(Infection,'infectiousIllinois.csv');
% writematrix(Vaccinated,'vaccinatedIllinois.csv');
% writematrix(Masked,'maskIllinois.csv');
% writematrix(Death,'deathIllinois.csv');
% writematrix(Mobility,'mobilityIllinois.csv');