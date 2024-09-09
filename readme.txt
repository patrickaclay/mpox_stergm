This repository holds 9 R scripts to fit and run a stergm model to the 2022
mpox outbreak, as well as two files describing the the model methods, and some
datasheets (which will need to be put into appropriate folders, examine read and save
functions in R scripts).

00_mpox_stergm_modules.R: This code contains functions to run the network model.

00_setup_functions.R: This code contains functions to build the sexual network. 

01_MSM_network_setup_167l_1percent.R: This code builds the sexual network. 

02_param_set_create.R: This code build uses latin hypsercube sampling to create
prior parameter sets. 

03_run_priors.R: This code runs the simulations for each prior parameter set. 

04_fitting_priors.R: selects 10% best fitting prior paramter sets. 
Calls NYC_data_through_jan25.csv. 

05_run_mpox_stergm_scenarios: Runs a year long simulation for each fit parameter
set. Can alter whether vaccination happens, how much vaccination happens, etc.
(comments in code should give all options)

06_compile_posterior_runs.R: Will compile simulation output. Comments indicate
al the possible output paramters that you can choose to compile. 

07_figure_code_for_git.R: Will create and save figures for manuscript. 

main_stergm_methods.docx and supp_stergm_methods.docx are the methods sections
from a manuscript currently in review that will explain the model and methods. 

If you run these files, will need to be in the same place as folders named "fits", "posterior_runs_real", and
"prior_runs_real" for code to save data into/read data from. 

 

