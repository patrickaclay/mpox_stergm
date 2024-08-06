This repository holds 7 R scripts to fit and run a stergm model to the 2022
mpox outbreak, as well as two files describing the the model methods, and one 
empty dataframe to use for model fitting. 

00_mpox_stergm_modules.R: This code contains functions to run the network model.
You will need to alter # of vaccines given per weel in param_msm.
but will also need to become familiar if you with to add/alter simulation modules. 

00_setup_functions.R: This code contains functions to build the sexual network. 
You should not need to alter it. 

01_MSM_network_setup_167l_1percent.R: This code builds the sexual network. 
You will want to go in and change the population size and rerun/check for stability.

02_param_set_create.R: This code build uses latin hypsercube sampling to create
prior parameter sets. May need to change prior distribution is having trouble fitting
model.

03_run_priors.R: This code runs the simulations for each prior parameter set. 
Will need to alter parameters such as % of vaccines that are given to MSM.
Additionally, 

04_fitting_priors.R: selects 10% best fitting prior paramter sets. 
Calls NC_data_through_jan25.csv, which is currently empty. 
Will need to enter daily cases, and then compute comulative daily cases, 
and a running 7 day mean. 
Will want to visualize fit sets vs empirical cases to make sure you are getting a good fit. 

05_run_mpox_stergm_scenarios: Runs a year long simulation for each fit parameter
set. Can alter whether vaccination happens, how much vaccination happens, etc.
(comments in code should give all options)

06_compile_posterior_runs.R: Will compile simulation output. Comments indicate
al the possible output paramters that you can choose to compile. 

main_stergm_methods.docx and supp_stergm_methods.docx are the methods sections
from a manuscript currently in review that will explain the model and methods. 

 

