clearvars

close all
clc
% data obtained from https://www.healthdata.org/node/7425 , click on
% reference scenario
 data2020 = readtable('reference2020.csv');
 data2021 = readtable('reference2021.csv');
 data2022 = readtable('reference2022.csv');
% data2021 = csvread('reference2021.csv');
% data2022 = csvread('reference2022.csv');
% 
allData = [data2020;data2021;data2022];
idx = allData.location_id == 536;
USdata = allData(idx,:);
Infection = table2array(USdata(:,5));
Vaccinated = table2array(USdata(:,60));
Masked = table2array(USdata(:,57));
Death = table2array(USdata(:,42));
Mobility = table2array(USdata(:,49));
days2020 = find(data2020.location_id == 536);
days2021 = find(data2021.location_id == 536);
days2022 = find(data2022.location_id == 536);
plot(Vaccinated);

dayStops = [1 332 697 1002];
popChicagoMetro = [9684738 9601605 9509934 9433330];
days = 1:1:1002;
popDays = interp1(dayStops,popChicagoMetro,days);
% writematrix(Infection,'infectious.csv');
% writematrix(Vaccinated,'vaccinated.csv');
% writematrix(Masked,'mask.csv');
% writematrix(Death,'death.csv');
%writematrix(Mobility,'mobility.csv');
% writematrix(Infection,'infectiousIllinois.csv');
% writematrix(Vaccinated,'vaccinatedIllinois.csv');
% writematrix(Masked,'maskIllinois.csv');
% writematrix(Death,'deathIllinois.csv');
% writematrix(Mobility,'mobilityIllinois.csv');