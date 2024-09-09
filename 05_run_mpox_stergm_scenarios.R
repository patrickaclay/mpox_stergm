
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

posteriors <- read.csv("posterior_runs_real/parameter_posterior.csv")
posteriors <- posteriors[,3:9]


# Use the detectCores() function to find the number of cores in system
no_cores <- 5#detectCores()

# Setup cluster
registerDoParallel(makeCluster(no_cores))

#load(here("fits", "1k_fit_artnet.rda"))
load(here("fits", "fit_167k_highestrisk_hiv_edp.rda"))


tx.seek.prob <- 0.8 #probability of seeking treatment
infectious.period <- 27 #days until natural recovery
vaccine.effectiveness.1 <- 0.68#0.401
vaccine.effectiveness.2 <- 0.83#0.479
vaccinate.time.til.effectiveness <- 14
vaccine.multiplier <- 1
vaccine.timing.adjust <- 0
vaccine.delay <- 43 #time from start of simulation until June 26th (time of first vaccines recorded)
behave.delay <-8 #time from start of simulation untel May 22 (behavioral adaptation begins)
init.infs <- 5
latent.period <- 4
vaccination.proportion.param = 0.82
vaccination.switch <- TRUE
behavior.change.switch <- TRUE
vacc.strategy <- 1 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri



#we for some reason need to run things outside of the parallel track
#before the parallel track will accept things
trans.prob <- posteriors$trans[1]
surge.time.param <- posteriors$eventtime[1]
surge.trans.param <- posteriors$eventtrans[1]
import.param <- posteriors$import[1]
behavior.change.amount.param = posteriors$behaveadapt[1]
#source(here("02_modules_gradual_vacc.R"))
source(here("NYC_modules.R"))
# load in all other parameters
params <- param_msm()
# sequence of modules to run model
controls <- control_msm(nsteps = 10, nsims=1, ncores=1)
# initial conditions
inits <- init_msm()
# run simulation
sim <- netsim(est, params, inits, controls)


### Now run the full simulations for analysis

##1dpri Main

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 1 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_best_1dpri",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
}


##intermediate Main

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 2 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_best_int",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
  
}



##2 dose priority Main

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 3 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_best_2dpri",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
}





##No vaccination Main

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 0
  vacc.strategy <- 0 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_novacc",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
}




##1dpri 50% doses

for(iteration in 1:1){
  
vaccination.switch <- TRUE
vaccine.multiplier <- 0.5
vacc.strategy <- 1 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri

last.run <- 0
loops <- 20
runs.per.loop <- 5

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
  
  saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_limited_1dpri",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
  
}

}


##intermediate 50% doses

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 0.5
  vacc.strategy <- 2 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5


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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_limited_int",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  

}



##2 dose priority 50% doses

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 0.5
  vacc.strategy <- 3 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_limited_2dpri",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }

}




####Higher incremental vaccine

##1dpri high inc

vaccine.effectiveness.1 <- 0.64#0.401
vaccine.effectiveness.2 <- 0.85#0.479

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 1 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_best_1dpri_high_inc_VE",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
}


##intermediate high inc

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 2 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_best_int_high_inc_VE",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
  
}



##2 dose priority high inc

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 3 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_best_2dpri_high_inc_VE",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
}



####lower incremental vaccine

##1dpri low inc

vaccine.effectiveness.1 <- 0.72#0.401
vaccine.effectiveness.2 <- 0.81#0.479

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 1 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_best_1dpri_low_inc_VE",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
}


##intermediate low inc

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 2 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_best_int_low_inc_VE",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
  
}



##2 dose priority low inc

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 3 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_best_2dpri_low_inc_VE",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
}





######
#higher transmissability
##1dpri Main

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 1 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  for(loops.i in 0:(loops-1)){
    
    mod_output_trad <- foreach(i=(last.run + (loops.i * runs.per.loop) + 1):(last.run + (loops.i * runs.per.loop) + runs.per.loop), .packages=c("here","EpiModelHIV","EpiModel")) %dopar% {
      
      trans.prob <- max(1,posteriors$trans[i] * 1.2)
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_clade1_1dpri",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
}


##intermediate Main

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 2 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  
  for(loops.i in 0:(loops-1)){
    
    mod_output_trad <- foreach(i=(last.run + (loops.i * runs.per.loop) + 1):(last.run + (loops.i * runs.per.loop) + runs.per.loop), .packages=c("here","EpiModelHIV","EpiModel")) %dopar% {
      
      
      trans.prob <- max(1,posteriors$trans[i] * 1.2)
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_clade1_int",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
  
}



##2 dose priority Main

for(iteration in 1:1){
  
  vaccination.switch <- TRUE
  vaccine.multiplier <- 1
  vacc.strategy <- 3 #0 = no vacc, 1 = 1dpri, 2 = int, 3 = 2dpri
  
  last.run <- 0
  loops <- 20
  runs.per.loop <- 5
  
  
  for(loops.i in 0:(loops-1)){
    
    mod_output_trad <- foreach(i=(last.run + (loops.i * runs.per.loop) + 1):(last.run + (loops.i * runs.per.loop) + runs.per.loop), .packages=c("here","EpiModelHIV","EpiModel")) %dopar% {
      
      trans.prob <- max(1,posteriors$trans[i] * 1.2)
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
    
    saveRDS(mod_output_trad, file=paste(paste("posterior_runs_real/posterior_nyc_output_clade1_2dpri",(last.run + (loops.i * runs.per.loop) + 1),(last.run + (loops.i * runs.per.loop) + runs.per.loop),iteration,sep="_"),"RData",sep="."))
    
  }
  
}



