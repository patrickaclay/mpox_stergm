# Here, we run a our model forward for each prior parameter set

#load in libraries
library(here)
library(parallel)
library(EpiModel)
library(EpiModelHIV)
library(doParallel)
library(foreach)
set.seed(12345)



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
vaccine.effectiveness.1 <- 0.68 #based on unpublished CDC model
vaccine.effectiveness.2 <- 0.83 #based on unpublished CDC model
vaccinate.time.til.effectiveness <- 14 #days after vaccine given til becomes effective
vaccine.multiplier <- 1 #if you want to multiply vaccines given by a set relative amount
vaccine.timing.adjust <- 0 #if you want to delay (positive number) or speed up (negative number) vaccines administered by set number of days
vaccine.delay <- 43 #time from start of simulation until time of first vaccines recorded (we do not include post-exposure doses given at beginning of outbreak)
behave.delay  <- 8 #time from start of simulation until may 22 (behavioral adaptation begins)
init.infs <- 5 #number of initial infections
latent.period <- 4 #length of pre-symptomatic infectious period, based on Brosius et al 2023, J Med Virol
vaccination.proportion.param = 0.82 #proportion of vaccinations that were administered to MSM
vaccination.switch <- TRUE #set as false if want to remove vaccination
behavior.change.switch <- TRUE #set as false if want to remove behavioral adaptation
vacc.strategy <- 1 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri

#load in priors
params.set <- read.csv("params_set.csv")
params.set <- as.data.frame(params.set[,-1])

#we for some reason need to run things outside of the parallel track
#before the parallel track will accept things
trans.prob <- params.set$trans[1]
#pride.exposures <- round(params.set$event[1]/28)
surge.time.param <- params.set$eventtime[1]
surge.trans.param <- params.set$eventtrans[1]
import.param <- params.set$import[1]
behavior.change.amount.param = params.set$behaveadapt[1]


#source(here("02_modules_gradual_vacc.R"))
    source(here("00_mpox_stergm_modules.R"))
# load in all other parameters
params <- param_msm()
# sequence of modules to run model
controls <- control_msm(nsteps = 10, nsims=1, ncores=1)
# initial conditions
inits <- init_msm()
# run simulation
sim <- netsim(est, params, inits, controls)

#create looping paramters
last.run <- 0
loops <- 166
runs.per.loop <- 6



for(loops.i in 0:(loops-1)){
  
  mod_output_trad <- foreach(i=(last.run + (loops.i * runs.per.loop) + 1):(last.run + (loops.i * runs.per.loop) + runs.per.loop), .packages=c("here","EpiModelHIV","EpiModel")) %dopar% {
    
    trans.prob <- params.set$trans[i]
    surge.time.param <- params.set$eventtime[i]
    surge.trans.param <- params.set$eventtrans[i]
    import.param <- params.set$import[i]
    behavior.change.amount.param = params.set$behaveadapt[i]
        source(here("00_mpox_stergm_modules.R"))
    
    
    # load in all other parameters
    params <- param_msm()
    # sequence of modules to run model
    controls <- control_msm(nsteps = 250, nsims=1, ncores=1)
    # initial conditions
    inits <- init_msm()
    # run simulation
    sim <- netsim(est, params, inits, controls)
    
    
  }
  saveRDS(mod_output_trad, file=paste(paste("prior_runs_real/abc_nyc_output_low",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),sep="_"),"RData",sep="."))
  
  rm(mod_output_trad)

}


last.run <- 996
loops <- 1
runs.per.loop <- 4



for(loops.i in 0:(loops-1)){
  
  mod_output_trad <- foreach(i=(last.run + (loops.i * runs.per.loop) + 1):(last.run + (loops.i * runs.per.loop) + runs.per.loop), .packages=c("here","EpiModelHIV","EpiModel")) %dopar% {
    
    trans.prob <- params.set$trans[i]
    surge.time.param <- params.set$eventtime[i]
    surge.trans.param <- params.set$eventtrans[i]
    import.param <- params.set$import[i]
    behavior.change.amount.param = params.set$behaveadapt[i]
    source(here("00_mpox_stergm_modules.R"))
    
    
    # load in all other parameters
    params <- param_msm()
    # sequence of modules to run model
    controls <- control_msm(nsteps = 250, nsims=1, ncores=1)
    # initial conditions
    inits <- init_msm()
    # run simulation
    sim <- netsim(est, params, inits, controls)
    
    
  }
  saveRDS(mod_output_trad, file=paste(paste("prior_runs_real/abc_nyc_prior_output",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),sep="_"),"RData",sep="."))
  
  #rm(mod_output_trad)
  
}

   #stop the cluster
stopImplicitCluster() 
    

