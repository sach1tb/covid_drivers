clearvars

close all
clc

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
plot(Vaccinated);
% writematrix(Infection,'infectious.csv');
% writematrix(Vaccinated,'vaccinated.csv');
% writematrix(Masked,'mask.csv');
% writematrix(Death,'death.csv');
%writematrix(Mobility,'mobility.csv');
writematrix(Infection,'infectiousIllinois.csv');
writematrix(Vaccinated,'vaccinatedIllinois.csv');
writematrix(Masked,'maskIllinois.csv');
writematrix(Death,'deathIllinois.csv');
writematrix(Mobility,'mobilityIllinois.csv');