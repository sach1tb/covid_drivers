clearvars

% data obtained as follows:
% > go to https://www.healthdata.org/node/7425, 
% > click on each reference scenario
 data2020 = readtable('../data/data_download_file_reference_2020_ILLINOIS.csv');
 data2021 = readtable('../data/data_download_file_reference_2021_ILLINOIS.csv');
 data2022 = readtable('../data/data_download_file_reference_2022_ILLINOIS.csv');

% data2020 = readtable('../data/reference2020.csv');
% data2021 = readtable('../data/reference2021.csv');
% data2022 = readtable('../data/reference2022.csv');
% 


allData = [data2020;data2021;data2022];
% idx = allData.location_id == 536;
ILdata = allData;
% instead write, Infection=ILdata.inf_mean
Infection = ILdata.inf_mean; % inf_mean
Vaccinated = table2array(ILdata(:,60)); % ??
Masked = table2array(ILdata(:,51)); % mask_use_mean, was 57 from reference
Death = table2array(ILdata(:,42)); % daily_deaths
Mobility = table2array(ILdata(:,49)); % ??
% days2020 = find(data2020.location_id == 536);
% days2021 = find(data2021.location_id == 536);
% days2022 = find(data2022.location_id == 536);
% plot(Vaccinated);

dayStops = [1 332 697 1002];
popChicagoMetro = [9684738 9601605 9509934 9433330];
days = 1:1:1002;
popDays = interp1(dayStops,popChicagoMetro,days);


% plot & verify

figure(1); gcf;clf;
subplot(2,3,1);
plot(Infection, 'r', 'linewidth', 2);
hold on;
inf_data_current=csvread('../data/infectiousIllinois_ci.csv');
plot(inf_data_current(:,2), 'w--');


% writematrix(Infection,'infectious.csv');
% writematrix(Vaccinated,'vaccinated.csv');
% writematrix(Masked,'mask.csv');
% writematrix(Death,'death.csv');
% writematrix(Mobility,'mobility.csv');
% writematrix(Infection,'infectiousIllinois.csv');
% writematrix(Vaccinated,'vaccinatedIllinois.csv');
% writematrix(Masked,'maskIllinois.csv');
% writematrix(Death,'deathIllinois.csv');
% writematrix(Mobility,'mobilityIllinois.csv');