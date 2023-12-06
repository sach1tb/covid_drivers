# covid_drivers
Source for paper: S. Butail, A. Bhattacharya, and M. Porfiri, "Estimating hidden relationships in dynamical systems: discovering drivers of infection rates of COVID-19"

## scripts
Run in the following order to reproduce the results
1. extract_data.m: takes data from different sources, filters for Illinois, and puts them in csv files in the data folder
2. ukfCovidRecon.m: runs the UKF with an updated ukfConstrained which keeps the estimates within bounds; also plots figures 3, 4, S1 and S2
3. mi_idot2param_sigma_shuffle.m: calculates mutual information between idot and parameters with shuffling of sigma points; also plots figures 5, 6, S3 and S4 (run first section and then go to plotting directly to avoid time consuming computations)
4. cmi_idot2param_sigma_shuffle.m: calculates conditional mutual information and plots supplementary figures S5-S8 and S9-S11 (run first section and then go to plotting directly to avoid time consuming computations)

## data
All data is in data folder. 

## other sources
boundedline folder is from https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m

