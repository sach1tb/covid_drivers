# covid_drivers
Source for paper: S. Butail, A. Bhattacharya, and M. Porfiri, "Estimating hidden relationships in dynamical systems: discovering drivers of infection rates of COVID-19"

Run in the following order to reproduce the results
1. extract_data.m: takes data from different sources, filters for Illinois, and puts them in csv files in the data folder
2. ukfCovidRecon.m: runs the UKF with an updated ukfConstrained which keeps the estimates within bounds
3. mi_idot2param_sigma_shuffle.m
4. cmi_idot2param_sigma_shuffle.m
5. 
