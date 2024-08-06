
#Here we run similations for one year, with different possible
#vaccination scenarios

#load in libraries
library(here)
library(parallel)
library(EpiModel)
library(EpiModelHIV)
library(doParallel)
library(foreach)
set.seed(12345)

#read in fit parameters, 100 paramter sets
posteriors <- read.csv("posterior_runs_real/parameter_posterior.csv")
posteriors <- posteriors[,3:9]


#set number of cores if running in parallel. 
#detectCores() will use all cores, I use 5 because using all 
#8 was causing memory problems
no_cores <- 5#detectCores()

# Setup cluster
registerDoParallel(makeCluster(no_cores))

#load the network
load(here("fits", "fit_167k_highestrisk_hiv_edp.rda"))

#input parameters
tx.seek.prob <- 0.8 #probability of seeking treatment
infectious.period <- 27 #days until natural recovery
vaccine.effectiveness.1 <- 0.75 #based on unpublished CDC model
vaccine.effectiveness.2 <- 0.89 #based on unpublished CDC model
vaccinate.time.til.effectiveness <- 14 #days after vaccine given til becomes effective
vaccine.multiplier <- 1 #if you want to multiply vaccines given by a set relative amount
vaccine.timing.adjust <- 0 #if you want to delay (positive number) or speed up (negative number) vaccines administered by set number of days
vaccine.delay <- 43 #time from start of simulation until time of first vaccines recorded (we do not include post-exposure doses given at beginning of outbreak)
init.infs <- 5 #number of initial infections
latent.period <- 4 #length of pre-symptomatic infectious period, based on Brosius et al 2023, J Med Virol
vaccination.proportion.param = 0.82 #proportion of vaccinations that were administered to MSM
vaccination.switch <- TRUE #set as false if want to remove vaccination
behavior.change.switch <- TRUE #set as false if want to remove behavioral adaptation


#For some reason need to run the model outside of the parallel process
#for a few time steps before the parallel track will accept things
#read in fit paramters
trans.prob <- posteriors$trans[1]
surge.time.param <- posteriors$eventtime[1]
surge.trans.param <- posteriors$eventtrans[1]
import.param <- posteriors$import[1]
behavior.change.amount.param = posteriors$behaveadapt[1]
source(here("00_mpox_stergm_modules.R"))
# load in all other parameters
params <- param_msm()
# sequence of modules to run model
controls <- control_msm(nsteps = 10, nsims=1, ncores=1)
# initial conditions
inits <- init_msm()
# run simulation
sim <- netsim(est, params, inits, controls)


##now run on parallel track
vaccination.switch <- TRUE
vaccine.multiplier <- 1

last.run <- 0
loops <- 20
runs.per.loop <- 5

source(here("00_mpox_stergm_modules.R"))


for(loops.i in 0:(loops-1)){
  
  mod_output_trad <- foreach(i=(last.run + (loops.i * runs.per.loop) + 1):(last.run + (loops.i * runs.per.loop) + runs.per.loop), .packages=c("here","EpiModelHIV","EpiModel")) %dopar% {
    
    trans.prob <- posteriors$trans[i]
    surge.time.param <- posteriors$eventtime[i]
    surge.trans.param <- posteriors$eventtrans[i]
    import.param <- posteriors$import[i]
    behavior.change.amount.param = posteriors$behaveadapt[i]
    
    # load in all other parameters
    params <- param_msm()
    # sequence of modules to run model
    controls <- control_msm(nsteps = 365, nsims=1, ncores=1)
    # initial conditions
    inits <- init_msm()
    # run simulation
    sim <- netsim(est, params, inits, controls)
    
  }
  
  saveRDS(mod_output_trad, file=paste(paste("posterior",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),sep="_"),"RData",sep="."))
  
}



