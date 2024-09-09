##This code contains modules for simulating mpox transmission
##through a population of MSM using a STERGM framework
##Written by Patrick Clay and Emily Pollock
##last updated 8/5/2024

# param_msm enters all model parameters
# control_msm runs all modules in a particular order, determines number and length of runs
# initialize_msm builds the initial network and assigns attributes to nodes
# simnet_msm controls the breaking and creating of sexual contact at each timestep
# edges_correct_msm determines behavioral adaptation
# infect_msm determines new infections at each timestep
# discord_edgelist and discord_edgelist_asym determine susceptible/infected pairs
# progress_msm determines state transition for infected individuals
# prevalence_msm - records summary statistics (e.g. new infections)
# vaccinate_msm - vaccinates individuals
# init_tergmLite - I overwrite part of the epimodel base functions
#because they were throwing errors



# Parameters ####
param_msm <- function(init.inf = 5,             # initial number of infected people 
                      behavior.change = behavior.change.switch,
                      vaccination = vaccination.switch,
                      behavior.change.amount = behavior.change.amount.param,
                      surge.time = surge.time.param,
                      surge.trans = surge.trans.param,
                      import = import.param,
                      behavior.delay = behave.delay,
                      
                      vaccine.multiply = vaccine.multiplier,
                      vaccination.proportion.msm = vaccination.proportion.param,
                      strategy = vacc.strategy,

                      act.rate.main = 1.54/7, 
                      act.rate.casual = 0.96/7, 
                      act.rate.instant = 1, 
                      e.to.a.rate = 1/(7.6 - latent.period),     # transition rate from latent (e) to asymptomatically infectious (a)
                      a.to.i.rate = 1/latent.period,       # transition rate from a to symptomatic. Currently have it set so that a does not transmit
                      testing.rate = 1/5,    # will take 8 days on average (median) for symptomatic person to seek testing
                      
                      inf.prob = trans.prob,       # probability of infection upon sexual contact
                      testing.prob = tx.seek.prob, # prob of inf. person seeking testing at all
                      i.to.r.rate = 1/infectious.period,  # natural clearance rate
                      
                      nodal.tx = TRUE, 
                      vaccinate.coverage.1.1 = 478*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.2 = 1573*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.3 = 2729*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.4 = 7298*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.5 = 11342*vaccine.multiplier*vaccination.proportion.param,
                      vaccinate.coverage.1.6 = 15542*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.7 = 20139*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.8 = 12273*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.9 = 6398*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.10 = 6118*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.11 = 3056*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.12 = 2490*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.13 = 1946*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.14 = 1748*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.15 = 1255*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.16 = 1090*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.17 = 845*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.18 = 706*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.19 = 444*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.20 = 393*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.21 = 359*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.22 = 159*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.23 = 222*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.24 = 148*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.25 = 150*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.26 = 158*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.27 = 115*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.28 = 118*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.29 = 111*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.30 = 98*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.31 = 95*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.32 = 101*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.33 = 87*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.34 = 74*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.35 = 48*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.36 = 53*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.37 = 66*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.1.38 = 66*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.1 = 6*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.2 = 9*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.3 = 24*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.4 = 43*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.5 = 65*vaccine.multiplier*vaccination.proportion.param, #sum of first five weeks, to account for 28 day gap
                      vaccinate.coverage.2.6 = 96*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.7 = 216*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.8 = 287*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.9 = 213*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.10 = 263*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.11 = 2340*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.12 = 7193*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.13 = 14588*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.14 = 9510*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.15 = 5390*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.16 = 2264*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.17 = 1717*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.18 = 1284*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.19 = 853*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.20 = 1040*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.21 = 705*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.22 = 227*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.23 = 342*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.24 = 269*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.25 = 276*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.26 = 214*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.27 = 134*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.28 = 187*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.29 = 182*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.30 = 146*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.31 = 89*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.32 = 117*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.33 = 100*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.34 = 88*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.35 = 78*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.36 = 101*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.37 = 86*vaccine.multiplier*vaccination.proportion.param, 
                      vaccinate.coverage.2.38 = 70*vaccine.multiplier*vaccination.proportion.param,
                      
                      vaccinate.timeline.1 = 0 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.2 = 7 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.3 = 14 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.4 = 21 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.5 = 28 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.6 = 35 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.7 = 42 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.8 = 49 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.9 = 56 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.10 = 63 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.11 = 70 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.12 = 77 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.13 = 84 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.14 = 91 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.15 = 98 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.16 = 105 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.17 = 112 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.18 = 119 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.19 = 126 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.20 = 133 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.21 = 140 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.22 = 147 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.23 = 154 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.24 = 161 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.25 = 168 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.26 = 175 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.27 = 182 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.28 = 189 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.29 = 196 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.30 = 203 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.31 = 210 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.32 = 217 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.33 = 224 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.34 = 231 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.35 = 238 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.36 = 245 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.37 = 252 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.38 = 259 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.timeline.over = 266 + vaccine.timing.adjust + vaccine.delay,  
                      vaccinate.interval = vaccinate.time.til.effectiveness,  # time after vaccine that dose becomes effective
                      vacc.effect.1 = vaccine.effectiveness.1,       # two weeks after first dose
                      vacc.effect.2 = vaccine.effectiveness.2,       # two weeks after second dose
                      rel.behave.change.1 = 0,
                      rel.behave.change.2 = 0.015,
                      rel.behave.change.3 = 0.113,
                      rel.behave.change.4 = 0.270,
                      rel.behave.change.5 = 0.379,
                      rel.behave.change.6 = 0.478,
                      rel.behave.change.7 = 0.496,
                      rel.behave.change.8 = 0.558,
                      rel.behave.change.9 = 0.838,
                      rel.behave.change.10 = 1,
                      rel.behave.change.11 = 0.964,
                      rel.behave.change.12 = 0.862,
                      rel.behave.change.13 = 0.791,
                      rel.behave.change.14 = 0.779,
                      rel.behave.change.15 = 0.739,
                      rel.behave.change.16 = 0.656,
                      rel.behave.change.17 = 0.589,
                      rel.behave.change.18 = 0.533,
                      rel.behave.change.19 = 0.415,
                      rel.behave.change.20 = 0.328,
                      rel.behave.change.21 = 0.355,
                      rel.behave.change.22 = 0.413,
                      rel.behave.change.23 = 0.431,
                      rel.behave.change.24 = 0.339,
                      rel.behave.change.25 = 0.347,
                      rel.behave.change.26 = 0.423,
                      rel.behave.change.27 = 0.483,
                      rel.behave.change.28 = 0.367,
                      rel.behave.change.29 = 0.211,
                      rel.behave.change.30 = 0.252,
                      rel.behave.change.31 = 0.301,
                      rel.behave.change.32 = 0.296,
                      rel.behave.change.33 = 0.185,


                      
                      ...) {
  
  ## Process parameters
  
  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))
  
  class(p) <- "param.net"
  return(p)
  
}

# Control Settings ####
control_msm <- function(simno = 1,
                        nsteps = 200,
                        start = 1,
                        nsims = 1,
                        ncores = 4,
                        cumulative.edgelist = FALSE,
                        truncate.el.cuml = 0,
                        initialize.FUN = initialize_msm,
                        progress.FUN = progress_msm,
                        vaccination.FUN = vaccinate_msm,
                        infection.FUN = infect_msm,
                        resim_nets.FUN = simnet_msm,
                        prev.FUN = prevalence_msm,
                        verbose.FUN = verbose.net,
                        module.order = NULL,
                        save.nwstats = FALSE,
                        save.other = c("el", "attr"),
                        tergmLite = TRUE,
                        tergmLite.track.duration = FALSE, # CPN2
                        set.control.ergm = control.simulate.formula(MCMC.burnin = 2e5),
                        set.control.stergm = control.simulate.network(),
                        verbose = TRUE,
                        skip.check = TRUE,
                        ...) {
  
  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)
  
  p$skip.check <- TRUE
  p$save.transmat <- FALSE
  
  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)
  p[["f.names"]] <- c(p[["bi.mods"]], p[["user.mods"]])
  p$save.other <- c("attr", "temp", "el")
  
  p$save.network <- FALSE
  if (is.null(p$verbose.int)) {
    p$verbose.int <- 1
  }
  
  p <- set.control.class("control.net", p) # CPN2
  return(p)
}


##################### Critical Network Modules #####################################

#### Initialize Module #####
# gets called once at the beginning of each simulation to construct networks
# and master dat object 

initialize_msm <- function(x, param, init, control, s) { #So this is what sets up the network that is then changed by simnet
  
  dat <- create_dat_object(param, init, control) #Ah, this is why everything depends on init
  
  #### Network Setup ####
  # Initial network setup
  # Simulate each network on the basis of their fits
  # Add to dat object 
  
  dat[["nw"]] <- list()
  nnets <- 3
  for (i in 1:nnets) {
    dat[["nw"]][[i]] <- simulate( 
      x[[i]][["fit"]],
      basis = x[[i]][["fit"]][["newnetwork"]],
      dynamic = FALSE
    )
  }
  nw <- dat[["nw"]]
  
  # Pull Network parameters
  dat[["nwparam"]] <- list()
  for (i in 1:nnets) {
    dat[["nwparam"]][i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }
  
  # Convert to tergmLite method
  dat <- init_tergmLite(dat)
  
  #### Nodal Attributes Setup ####
  dat[["attr"]] <- param[["netstats"]][["attr"]]
  
  num <- network.size(nw[[1]])
  dat <- append_core_attr(dat, 1, num)
  
  # Pull in attributes on network. 
  # We pull from the one-time network because this has both deg.main and deg.pers attributes 
  # (see estimation script)
  nwattr.all <- names(nw[[3]][["val"]][[3]])
  nwattr.use <- nwattr.all[!nwattr.all %in% c("na", "vertex.names")]
  for (i in seq_along(nwattr.use)) {
    dat$attr[[nwattr.use[i]]] <- get.vertex.attribute(nw[[3]], nwattr.use[i])
  }
  
  # Add other attributes 
  
  # First we need some initial conditions parameters
  init.inf <- get_param(dat, "init.inf")
  
  # Generate status vector based on nums init vaccinated and init infected (non-overlapping groups)
  # starting values s / i / v
  # initial infections among highest-activity groups unless init size is > than size of those groups 
  
  status <- rep("s", num)
  riskg <- get_attr(dat, "riskg")
  
  # initial infecteds, groups 4-6
  risk.high <- which(status == "s" & (riskg == "4" | riskg == "5" | riskg == "6"))
  
  if (length(risk.high) > init.inf) {
    infected <- sample(risk.high, init.inf, replace=FALSE)
  }  else {infected <- sample(which(status=="s"), init.inf, replace=FALSE)}
  
  status[infected] <- "a"
  
  # vaccinated time vector (NA until vaccinated)
  vaccTime1 <- rep(NA, num)
  vaccTime2 <- rep(NA, num)

  # infection time vector (NA until infected)
  infTime <- rep(NA, num)
  infTime[infected] <- 1
  
  # secondary infections from node vector (NA until an infection is transmitted from node to partner)
  secInfs <- rep(NA, num)
  secInfs[infected] <- 0
  
  # which generation is this infection
  infGen <- rep(NA, num)
  infGen[infected] <- 1
  
  # Pull sqrt age to get age 
  sqrt.age <- get_attr(dat, "sqrt.age")
  age <- sqrt.age^2
  
  # Set nodal attribute for tx seeking if nodal.tx=TRUE
  if (dat$param$nodal.tx){
    tx.prob <- get_param(dat, "testing.prob")
    tx.seek <- rep(0, num)
    seek <- sample(1:num, tx.prob*num, replace=FALSE)
    tx.seek[seek] <- 1
    dat <- set_attr(dat, "tx.seek", tx.seek)
  }
  
  # set attributes
  dat <- set_attr(dat, "secInfs", secInfs)     # record how many infections per infection
  dat <- set_attr(dat, "infTime", infTime)     # record time step of infection
  dat <- set_attr(dat, "infGen", infGen)       # generation of infection
  dat <- set_attr(dat, "vaccTime1", vaccTime1)   # record time step of vaccination
  dat <- set_attr(dat, "vaccTime2", vaccTime2)   # record time step of vaccination
  dat <- set_attr(dat, "status", status)       # initial status
  dat <- set_attr(dat, "age", sqrt.age^2)      # age attribute to accompany the sqrt.age attribute
  dat <- set_attr(dat, "recTime", rep(NA,num)) # infection recovery time 
  dat <- set_attr(dat, "testTime", rep(NA,num)) # infection testing time 
  dat <- set_attr(dat, "pre.riskg", riskg) #holds their pre-infection riskgroup as a reference. 

  dat$num.nw <- 3
  
  #### Other Setup ####
  dat[["stats"]] <- list()
  dat[["stats"]][["nwstats"]] <- list()
  dat[["temp"]] <- list()
  dat[["epi"]] <- list()
  
  dat <- set_epi(dat, "num", at = 1,  num)
  dat <- set_epi(dat, "cuml.infs", at = 1, init.inf)
  dat <- set_epi(dat, "cuml.cases", at = 1, 0)
  dat <- set_epi(dat, "prev", at = 1, init.inf)
  
  
  # Setup Partner List for all 3 networks 
  # (only gets updated if tracking turned on in control function)
#  for (n_network in seq_len(3)) {
#    dat <- update_cumulative_edgelist(dat, n_network)
#  }
  
  # Network statistics
  if (dat[["control"]][["save.nwstats"]]) {
    for (i in seq_along(x)) {
      nwL <- networkLite(dat[["el"]][[i]], dat[["attr"]])
      nwstats <- summary(
        dat[["control"]][["nwstats.formulas"]][[i]],
        basis = nwL,
        term.options = dat[["control"]][["mcmc.control"]][[i]][["term.options"]],
        dynamic = i < 3
      )
      
      dat[["stats"]][["nwstats"]][[i]] <- matrix(
        nwstats, 
        nrow = 1, ncol = length(nwstats),
        dimnames = list(NULL, names(nwstats))
      )
      
      dat[["stats"]][["nwstats"]][[i]] <- 
        as.data.frame(dat[["stats"]][["nwstats"]][[i]])
    }
  }
  
  class(dat) <- "dat"
  return(dat)
  
}



#### Simnet -- Simulate Networks & Update Coefs based on pop size changes ###############
simnet_msm <- function(dat, at) {
  
  #########################
  ## Nest the edges_correct function
  # (workaround for some functions not being accessible when running in parallel even if in global environ)
  
  edges_correct_msm <- function(dat, at) {
    
      behavior.change         <- get_param(dat, "behavior.change")
      behavior.delay          <- get_param(dat, "behavior.delay")

    if(behavior.change == TRUE){
    
      rel.behave.change.1         <- get_param(dat, "rel.behave.change.1")
      rel.behave.change.2         <- get_param(dat, "rel.behave.change.2")
      rel.behave.change.3         <- get_param(dat, "rel.behave.change.3")
      rel.behave.change.4         <- get_param(dat, "rel.behave.change.4")
      rel.behave.change.5         <- get_param(dat, "rel.behave.change.5")
      rel.behave.change.6         <- get_param(dat, "rel.behave.change.6")
      rel.behave.change.7         <- get_param(dat, "rel.behave.change.7")
      rel.behave.change.8         <- get_param(dat, "rel.behave.change.8")
      rel.behave.change.9         <- get_param(dat, "rel.behave.change.9")
      rel.behave.change.10         <- get_param(dat, "rel.behave.change.10")
      rel.behave.change.11         <- get_param(dat, "rel.behave.change.11")
      rel.behave.change.12         <- get_param(dat, "rel.behave.change.12")
      rel.behave.change.13         <- get_param(dat, "rel.behave.change.13")
      rel.behave.change.14         <- get_param(dat, "rel.behave.change.14")
      rel.behave.change.15         <- get_param(dat, "rel.behave.change.15")
      rel.behave.change.16         <- get_param(dat, "rel.behave.change.16")
      rel.behave.change.17         <- get_param(dat, "rel.behave.change.17")
      rel.behave.change.18         <- get_param(dat, "rel.behave.change.18")
      rel.behave.change.19         <- get_param(dat, "rel.behave.change.19")
      rel.behave.change.20         <- get_param(dat, "rel.behave.change.20")
      rel.behave.change.21         <- get_param(dat, "rel.behave.change.21")
      rel.behave.change.22         <- get_param(dat, "rel.behave.change.22")
      rel.behave.change.23         <- get_param(dat, "rel.behave.change.23")
      rel.behave.change.24         <- get_param(dat, "rel.behave.change.24")
      rel.behave.change.25         <- get_param(dat, "rel.behave.change.25")
      rel.behave.change.26         <- get_param(dat, "rel.behave.change.26")
      rel.behave.change.27         <- get_param(dat, "rel.behave.change.27")
      rel.behave.change.28         <- get_param(dat, "rel.behave.change.28")
      rel.behave.change.29         <- get_param(dat, "rel.behave.change.29")
      rel.behave.change.30         <- get_param(dat, "rel.behave.change.30")
      rel.behave.change.31         <- get_param(dat, "rel.behave.change.31")
      rel.behave.change.32         <- get_param(dat, "rel.behave.change.32")
      rel.behave.change.33         <- get_param(dat, "rel.behave.change.33")
      rel.behave.change <- c(rel.behave.change.1,rel.behave.change.2,rel.behave.change.3,
                             rel.behave.change.4,rel.behave.change.5,rel.behave.change.6,
                             rel.behave.change.7,rel.behave.change.8,rel.behave.change.9,
                             rel.behave.change.10,rel.behave.change.11,rel.behave.change.12,
                             rel.behave.change.13,rel.behave.change.14,rel.behave.change.15,
                             rel.behave.change.16,rel.behave.change.17,rel.behave.change.18,
                             rel.behave.change.19,rel.behave.change.20,rel.behave.change.21,
                             rel.behave.change.22,rel.behave.change.23,rel.behave.change.24,
                             rel.behave.change.25,rel.behave.change.26,rel.behave.change.27,
                             rel.behave.change.28,rel.behave.change.29,rel.behave.change.30,
                             rel.behave.change.31,rel.behave.change.32,rel.behave.change.33)
      
      behavior.change.amount <- get_param(dat, "behavior.change.amount")

      if(((at-behavior.delay)-1) < 2 | ((at-behavior.delay)-1) > 230){week.yesterday = 33}
      if(((at-behavior.delay)-1) >= 2 & ((at-behavior.delay)-1) <= 230){week.yesterday = floor((((at-behavior.delay)-1)-2)/7)+1}
      if((at-behavior.delay) < 2 | (at-behavior.delay) > 230){week.today = 33}
      if((at-behavior.delay) >= 2 & (at-behavior.delay) <= 230){week.today = floor(((at-behavior.delay)-2)/7)+1}
    old.num <- sum(dat$attr$active == 1, na.rm = TRUE) * (1 - behavior.change.amount * rel.behave.change[week.yesterday])                     #sets number of nodes at prior timestep
    new.num <- sum(dat$attr$active == 1, na.rm = TRUE) * (1 - behavior.change.amount * rel.behave.change[week.today]) #sets number of nodes at this timestep
    adjust <- log(new.num) - log(old.num)               #calculates log difference between those two
    
    
    for (i in length(dat$nwparam)) {
      
      coef.form1 <- get_nwparam(dat, network = i)$coef.form #get formation coefficient
      coef.form1[1] <- coef.form1[1] + adjust               #says to increase or decrease formation of edges
      dat$nwparam[[i]]$coef.form <- coef.form1              #re-records this number
    }
    }
    

    return(dat)
    
    
    
 }
  

  
  ##end nesting##
  #########################
  
  
  ## Grab parameters from dat object 
  cumulative.edgelist <- get_control(dat, "cumulative.edgelist") # are we tracking the cumulative edgelist (T/F)
  truncate.el.cuml <- get_control(dat, "truncate.el.cuml")       # how long in the past do we keep edgelist
  ## Grab parameters from dat object 
  set.control.stergm <- get_control(dat, "set.control.stergm")   # specific control settings for network simulation
  set.control.ergm <- get_control(dat, "set.control.ergm")       # specific control settings for network simulation
  save.nwstats <- get_control(dat, "save.nwstats")
  nwstats.formulas <- get_control(dat, "nwstats.formulas")
  
    ## Main network
    for (i in 1:length(dat$el)) {    #I believe this is where loops through overlapping networks
      nwparam <- EpiModel::get_nwparam(dat, network = i)   #get parameters of this network
      isTERGM <- ifelse(nwparam$coef.diss$duration > 1, TRUE, FALSE) #set isTERGM to true is relationships non-instantaneous
      
      nwL <- networkLite(dat[["el"]][[i]], dat[["attr"]])
      
      
#      if (get_control(dat, "tergmLite.track.duration") == TRUE) { #figure out how long relationships have lasted
#        nwL %n% "time" <- dat[["nw"]][[i]] %n% "time"
#        nwL %n% "lasttoggle" <- dat[["nw"]][[i]] %n% "lasttoggle"
#      }
      
      if(i==1){
        # update pers degree (nodal attribute based on network status)
        deg.pers <- get_attr(dat, "deg.pers")
        deg.pers <- get_degree(dat[["el"]][[2]])
        dat <- set_attr(dat, "deg.pers", deg.pers)
      }
      
      if(i==2){
        # update main degree (nodal attribute based on network status)
        deg.main <- get_attr(dat, "deg.main")
        deg.main <- get_degree(dat[["el"]][[1]])
        dat <- set_attr(dat, "deg.main", deg.main)
      }
      
      if (isTERGM == TRUE) { #The following happens only if tergm (relationships last)
        dat[["nw"]][[i]] <- simulate( #forms, dissolves contacts. generic function. simulate distribution corresponding to fitted model object
          nwL, #what are the attributes of each node
          formation = nwparam[["formation"]],
          dissolution = nwparam[["coef.diss"]][["dissolution"]],
          coef.form = nwparam[["coef.form"]],
          coef.diss = nwparam[["coef.diss"]][["coef.adj"]],
          constraints = nwparam[["constraints"]],
          time.start = at - 1,
          time.slices = 1,
          time.offset = 1,
          control = set.control.stergm,
          output = "final"
        )
      } else {  #The following happens if tergm is false, e.g. instantaneous relationships
        dat[["nw"]][[i]] <- simulate( #forms contacts
          basis = nwL,
          object = nwparam[["formation"]],
          coef = nwparam[["coef.form"]],
          constraints = nwparam[["constraints"]],
          control = set.control.ergm,
          dynamic = FALSE,
          nsim = 1,
          output = "network"
        )
      }
      
      dat[["el"]][[i]] <- as.edgelist(dat[["nw"]][[i]])
  
    
    if (get_control(dat, "save.nwstats") == TRUE) {
      term.options <- if (isTERGM == TRUE) {
        set.control.stergm$term.options
      } else {
        set.control.ergm$term.options
      }
      dat$stats$nwstats[[i]] <- rbind(dat$stats$nwstats[[i]],
                                      summary(dat$control$nwstats.formulas[[i]],
                                              basis = nwL,
                                              term.options = term.options,
                                              dynamic = isTERGM))
    }
    
  }
  
#  if (get_control(dat, "cumulative.edgelist") == TRUE) {
#    for (n_network in seq_len(3)) {
#      dat <- update_cumulative_edgelist(dat, n_network, truncate.el.cuml)
#    }
#  }
  
  return(dat)
}






################### DISEASE RELATED MODULES ############################
# Infection ####
##Infection
infect_msm <- function(dat, at) {
  

  # model-specific discordant edgelist function
  discord_edgelist_mpx <- function (dat, at, network){ 
    status <- get_attr(dat, "status")
    active <- get_attr(dat, "active")
    el <- get_edgelist(dat, network)
    del <- NULL
    if (nrow(el) > 0) {
      el <- el[sample(1:nrow(el)), , drop = FALSE]
      stat <- matrix(status[el], ncol = 2)
      isInf <- matrix(stat %in% c("i", "im", "a"), ncol = 2)
      isSus <- matrix(stat %in% c("s", "v1", "v2"), ncol = 2)
      SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
      ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
      pairs <- rbind(SIpairs, ISpairs[, 2:1])
      if (nrow(pairs) > 0) {
        sus <- pairs[, 1]
        inf <- pairs[, 2]
        del <- data.frame(at, sus, inf)
        keep <- rowSums(matrix(c(active[del$sus], active[del$inf]), 
                               ncol = 2)) == 2
        del <- del[keep, ]
        if (nrow(del) < 1) {
          del <- NULL
        }
      }
    }
    return(del)
  }
  
  
  #### Setup ####
  
  # Get attributes and parameters 
  active    <- get_attr(dat, "active")
  status    <- get_attr(dat, "status")
  secInfs   <- get_attr(dat, "secInfs")
  infTime   <- get_attr(dat, "infTime")
  infGen    <- get_attr(dat, "infGen")
  riskgroup <- get_attr(dat, "riskg")
  pre.riskgroup <- get_attr(dat, "pre.riskg") 
  vaccTime1  <- get_attr(dat, "vaccTime1")
  vaccTime2  <- get_attr(dat, "vaccTime2")
  testTime <- get_attr(dat, "testTime")
  
  inf.prob           <- get_param(dat, "inf.prob")
  act.rate.main      <- get_param(dat, "act.rate.main")
  act.rate.casual    <- get_param(dat, "act.rate.casual")
  act.rate.instant   <- get_param(dat, "act.rate.instant")
  vacc.effect.1        <- get_param(dat, "vacc.effect.1")
  vacc.effect.2        <- get_param(dat, "vacc.effect.2")
  vaccinate.interval <- get_param(dat, "vaccinate.interval")
  
  # Set up trackers 
  nInf <- 0
  idsInf <- NULL
  
  # By risk group 
  nInf.q1 <- 0
  nInf.q2 <- 0
  nInf.q3 <- 0
  nInf.q4 <- 0
  nInf.q5 <- 0
  nInf.q6 <- 0
  
  # New infections in each network
  nInfsMain <- 0
  nInfsPers <- 0
  nInfsInst   <- 0
  
  #### Transmissions ####
  # Loop through discordant edgelists in each network
  for(nw.transmit in 1:3){ 
    
    if(nw.transmit == 1){act.rate <- act.rate.main}
    if(nw.transmit == 2){
      behavior.change         <- get_param(dat, "behavior.change")
      behavior.delay         <- get_param(dat, "behavior.delay")
      
      if(behavior.change == FALSE){act.rate <- act.rate.casual}
      
      if(behavior.change == TRUE){
        rel.behave.change.1         <- get_param(dat, "rel.behave.change.1")
        rel.behave.change.2         <- get_param(dat, "rel.behave.change.2")
        rel.behave.change.3         <- get_param(dat, "rel.behave.change.3")
        rel.behave.change.4         <- get_param(dat, "rel.behave.change.4")
        rel.behave.change.5         <- get_param(dat, "rel.behave.change.5")
        rel.behave.change.6         <- get_param(dat, "rel.behave.change.6")
        rel.behave.change.7         <- get_param(dat, "rel.behave.change.7")
        rel.behave.change.8         <- get_param(dat, "rel.behave.change.8")
        rel.behave.change.9         <- get_param(dat, "rel.behave.change.9")
        rel.behave.change.10         <- get_param(dat, "rel.behave.change.10")
        rel.behave.change.11         <- get_param(dat, "rel.behave.change.11")
        rel.behave.change.12         <- get_param(dat, "rel.behave.change.12")
        rel.behave.change.13         <- get_param(dat, "rel.behave.change.13")
        rel.behave.change.14         <- get_param(dat, "rel.behave.change.14")
        rel.behave.change.15         <- get_param(dat, "rel.behave.change.15")
        rel.behave.change.16         <- get_param(dat, "rel.behave.change.16")
        rel.behave.change.17         <- get_param(dat, "rel.behave.change.17")
        rel.behave.change.18         <- get_param(dat, "rel.behave.change.18")
        rel.behave.change.19         <- get_param(dat, "rel.behave.change.19")
        rel.behave.change.20         <- get_param(dat, "rel.behave.change.20")
        rel.behave.change.21         <- get_param(dat, "rel.behave.change.21")
        rel.behave.change.22         <- get_param(dat, "rel.behave.change.22")
        rel.behave.change.23         <- get_param(dat, "rel.behave.change.23")
        rel.behave.change.24         <- get_param(dat, "rel.behave.change.24")
        rel.behave.change.25         <- get_param(dat, "rel.behave.change.25")
        rel.behave.change.26         <- get_param(dat, "rel.behave.change.26")
        rel.behave.change.27         <- get_param(dat, "rel.behave.change.27")
        rel.behave.change.28         <- get_param(dat, "rel.behave.change.28")
        rel.behave.change.29         <- get_param(dat, "rel.behave.change.29")
        rel.behave.change.30         <- get_param(dat, "rel.behave.change.30")
        rel.behave.change.31         <- get_param(dat, "rel.behave.change.31")
        rel.behave.change.32         <- get_param(dat, "rel.behave.change.32")
        rel.behave.change.33         <- get_param(dat, "rel.behave.change.33")
        rel.behave.change <- c(rel.behave.change.1,rel.behave.change.2,rel.behave.change.3,
                               rel.behave.change.4,rel.behave.change.5,rel.behave.change.6,
                               rel.behave.change.7,rel.behave.change.8,rel.behave.change.9,
                               rel.behave.change.10,rel.behave.change.11,rel.behave.change.12,
                               rel.behave.change.13,rel.behave.change.14,rel.behave.change.15,
                               rel.behave.change.16,rel.behave.change.17,rel.behave.change.18,
                               rel.behave.change.19,rel.behave.change.20,rel.behave.change.21,
                               rel.behave.change.22,rel.behave.change.23,rel.behave.change.24,
                               rel.behave.change.25,rel.behave.change.26,rel.behave.change.27,
                               rel.behave.change.28,rel.behave.change.29,rel.behave.change.30,
                               rel.behave.change.31,rel.behave.change.32,rel.behave.change.33)
        
        behavior.change.amount <- get_param(dat, "behavior.change.amount")
        
        if((at-behavior.delay) < 2 | (at-behavior.delay) > 230){week.today = 33}
        if((at-behavior.delay) >= 2 & (at-behavior.delay) <= 230){week.today = floor(((at-behavior.delay)-2)/7)+1}
        
        act.rate <- act.rate.casual * (1-(behavior.change.amount * rel.behave.change[week.today]))
      }}
    if(nw.transmit == 3){act.rate <- act.rate.instant}
    
    # get discordant edgelist
    del <- discord_edgelist_mpx(dat, at, network = nw.transmit)
    
    if (!(is.null(del))) {
      
      # add status column for sus 
      del$statusSus <- status[del$sus]
      # add column for vaccTime
      del$vaccTime1Sus <- vaccTime1[del$sus]
      del$vaccTime1Sus[which(is.na(del$vaccTime1Sus))] <- 0
      del$vaccTime2Sus <- vaccTime2[del$sus]
      del$vaccTime2Sus[which(is.na(del$vaccTime2Sus))] <- 0
      
      #add a column for status of I
      del$statusInf <- status[del$inf]
      #add a column for change in risk group
      del$riskchange <- pre.riskgroup[del$inf] - riskgroup[del$inf]
      
      # alter inf prob by vaccination status
      del$transProb <- inf.prob
      del$transProb[del$statusSus=="v1" & at > (del$vaccTime1Sus+vaccinate.interval)] <- inf.prob * (1-vacc.effect.1)
      del$transProb[del$statusSus=="v2" & at <= (del$vaccTime2Sus+vaccinate.interval)] <- inf.prob * (1-vacc.effect.1)
      del$transProb[del$statusSus=="v2" & at > (del$vaccTime2Sus+vaccinate.interval)] <- inf.prob * (1-vacc.effect.2)
      
      # act rate altered by status
      del$actRate <- act.rate
      
      if(nw.transmit == 1){
        del$actRate[del$statusInf == "i" | del$statusInf == "im"] <- act.rate/2
        del$actRate[del$statusInf == "i" & !is.na(del$testTime)] <- 0.0035/inf.prob}
      
      if(nw.transmit == 2){
        del$actRate[del$statusInf == "i" | del$statusInf == "im"] <- act.rate/2
        del$actRate[del$statusInf == "i" & !is.na(del$testTime)] <- 0}
      
      
      # final transmission probability 
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      
      # filter down to which pairs transmitted infection
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      
      idsNewInf <- unique(del$sus)
      idsOldInf <- unique(del$inf)
      
      if (length(idsNewInf) > 0) {
        
        
        status[idsNewInf]  <- "e"
        secInfs[idsNewInf] <- 0
        infTime[idsNewInf] <- at
        
        # in case any id infections more than 1 person during this time step
        if (length(idsOldInf) < length(idsNewInf)) {
          
          infGen[del$sus]  <- infGen[del$inf] + 1
          
          # count how many secondary infections per infecting id 
          tab <- table(del$inf)
          ids <- as.numeric(names(tab))
          infs <- as.numeric(tab)
          secInfs[ids] <- secInfs[ids] + infs
          
        } 
        
        # in case any susceptible ids get infected by more than one person at this time step
        # we pick the first time they show up in the edgelist
        if (length(idsNewInf) < length(idsOldInf)) {
          
          tab <- table(del$sus)
          ids <- as.numeric(names(which(tab>1)))
          
          for (i in 1:length(ids)){
            #if(length(ids)>1) {          browser()}
            d <- del[del$sus %in% ids[i],]
            if(i==1){ 
              first <- d[1,]
            } else {
              first <- rbind(first, d[1,])
            }
          }
          
          single_del <- del[-c(which(del$sus %in% ids)),]
          newdel <- rbind(single_del, first)
          
          idsNewInf <- newdel$sus
          idsOldInf <- newdel$inf
          
          secInfs[idsOldInf] <- secInfs[idsOldInf] + 1 
          infGen[idsNewInf]  <- infGen[idsOldInf] + 1 
          
        }
        
        
        if (length(idsNewInf) == length(idsOldInf)){ 
          
          secInfs[idsOldInf] <- secInfs[idsOldInf] + 1 
          infGen[idsNewInf]  <- infGen[idsOldInf] + 1 
        }
        
        if(nw.transmit == 1){nInfsMain <- length(idsNewInf)}
        if(nw.transmit == 2){nInfsPers <- length(idsNewInf)}
        if(nw.transmit == 3){nInfsInst <- length(idsNewInf)}
        
        idsInf <- c(idsInf, idsNewInf)
        
        nInf <- nInf + length(idsNewInf)
        infected_riskgrp <- riskgroup[idsNewInf]
        nInf.q1 <- nInf.q1 + length(infected_riskgrp[which(infected_riskgrp == 1)])
        nInf.q2 <- nInf.q2 + length(infected_riskgrp[which(infected_riskgrp == 2)])
        nInf.q3 <- nInf.q3 + length(infected_riskgrp[which(infected_riskgrp == 3)])
        nInf.q4 <- nInf.q4 + length(infected_riskgrp[which(infected_riskgrp == 4)])
        nInf.q5 <- nInf.q5 + length(infected_riskgrp[which(infected_riskgrp == 5)])
        nInf.q6 <- nInf.q6 + length(infected_riskgrp[which(infected_riskgrp == 6)])
        
      }
    } 
  }
  
  # Update Attrs and Trackers ####
  
  # nodal attributes 
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "infTime", infTime)
  dat <- set_attr(dat, "infGen", infGen)
  dat <- set_attr(dat, "secInfs", secInfs)
  
  # update epi trackers 
  dat <- set_epi(dat, "se.flow", at, nInf)
  
  prev.infs <- get_epi(dat, "cuml.infs", at=at-1)
  dat <- set_epi(dat, "cuml.infs", at, prev.infs + nInf)
  
  dat <- set_epi(dat, "se.flow.q1", at, nInf.q1)
  dat <- set_epi(dat, "se.flow.q2", at, nInf.q2)
  dat <- set_epi(dat, "se.flow.q3", at, nInf.q3)
  dat <- set_epi(dat, "se.flow.q4", at, nInf.q4)
  dat <- set_epi(dat, "se.flow.q5", at, nInf.q5)
  dat <- set_epi(dat, "se.flow.q6", at, nInf.q6)
  
  dat <- set_epi(dat, "se.flow.main", at, nInfsMain)
  dat <- set_epi(dat, "se.flow.pers", at, nInfsPers)
  dat <- set_epi(dat, "se.flow.ot", at, nInfsInst)
  
  ####superspreader event
  dat <- set_epi(dat, "superspreader.event", at, 0)

  if(at > 5){
  
    surge.time <- get_param(dat, "surge.time")
    surge.trans <- get_param(dat, "surge.trans")
    import <- get_param(dat, "import")
    inf_hr <- get_epi(dat,"inf_hr", at-1)
    sus_hr <- get_epi(dat,"sus_hr", at-1)
    N_hr <- get_epi(dat,"N_hr", at-1)
    
    
  if(at < surge.time){ #median 17 days from 5 "a" individuals to 5 cases, another 5 between five cases and pride.
    super.spreader <- sus_hr * (inf_hr/N_hr) * surge.trans + floor(import) + rbinom(1,1,(import-floor(import)))
    status <- get_attr(dat, "status")
    riskg <- get_attr(dat, "riskg")
  # initial infecteds
  risk.high <- which(status == "s" & (riskg == "4" | riskg == "5" | riskg == "6"))
  
  if (length(risk.high) < super.spreader) {super.spreader <- length(risk.high)}
  if (length(risk.high) >= super.spreader) {
    infected <- sample(risk.high, super.spreader, replace=FALSE)
  }  #else {infected <- sample(which(status=="s"), super.spreader, replace=FALSE)}
  
  status[infected] <- "e"
  
  old.infs <- get_epi(dat,"cuml.infs", at)
  dat <- set_epi(dat, "cuml.infs", at, old.infs + super.spreader)
  
  dat <- set_attr(dat, "status", status)
  dat <- set_epi(dat, "superspreader.event", at, 1)
  }}
  
  return(dat)
}


# Disease Progression ####
progress_msm <- function(dat, at) {
  
  # e  - latent class
  # a  - asymptomatic and infectious
  # i  - symptomatic and infectious
  # im - infectious but won't seek testing (not symptomatic enough, no access to care, etc) 
  # r  - recovered and immune 
  
  ## Get Params & Attributes 
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  recTime <- get_attr(dat, "recTime")
  testTime <- get_attr(dat, "testTime")
  riskg <- get_attr(dat, "riskg")
  pre.riskg <- get_attr(dat, "pre.riskg")
  tx.seek <- get_attr(dat, "tx.seek")
  
  e.to.a.rate <-     get_param(dat, "e.to.a.rate")
  a.to.i.rate <-     get_param(dat, "a.to.i.rate")
  i.to.r.rate <-     get_param(dat, "i.to.r.rate")
  testing.rate <-  get_param(dat, "testing.rate")
  if(at < 58){testing.rate <- 1/(1/testing.rate + (52 - at)*(10/52))}
  testing.prob <-  get_param(dat, "testing.prob")
    
    
    # natural recovery among infected 
    
    ## Natural Recovery among i & im 
    nRec       <- 0
    idsEligRec <- which(active == 1 & (status == "i" | status == "im"))
    nEligRec   <- length(idsEligRec)
    
    if (nEligRec > 0) {
      vecRec <- which(rbinom(nEligRec, 1, i.to.r.rate) == 1) # vector of who recovers without testing
      if (length(vecRec) > 0) {
        idsRec <- idsEligRec[vecRec]
        nRec   <- length(idsRec)
        status[idsRec] <- "r"
        recTime[idsRec] <- at
        riskg[idsRec] <- pre.riskg[idsRec]
        
      }
    }
    
    ## testing
    nTest      <- 0
    #idsEligTest <- which(active == 1 & status == "i") ## NEEDS TO ALSO FILTER BY TREATMENT TIME = NA
    #half the ascertainment rate before July first (aka at = 47)
    if(at < 47){
      idsEligTest <- which(active == 1 & status == "i" & is.na(testTime) & runif(length(active)) > 0.5)}
    if(at >= 47){
      idsEligTest <- which(active == 1 & status == "i" & is.na(testTime))}
    nEligTest   <- length(idsEligTest)
    
    
    if (nEligTest > 0) {
      vecTest <- which(rbinom(nEligTest, 1, testing.rate) == 1) # vector of who is tested
      if (length(vecTest) > 0) {
        idsTest <- idsEligTest[vecTest]
        nTest   <- length(idsTest)
        testTime[idsTest] <- at
        riskg[idsTest] <- 1
      }
    }
    
    ## Asympt. to Infectious but will not seek testing
    nInf_m       <- 0
    idsEligInf <- which(active == 1 & status == "a" & tx.seek == 0)
    nEligInf   <- length(idsEligInf)
    
    if (nEligInf > 0) {
      vecInf <- which(rbinom(nEligInf, 1, a.to.i.rate) == 1) 
      if (length(vecInf) > 0) {
        idsInf <- idsEligInf[vecInf]
        nInf_m  <- length(idsInf)
        status[idsInf] <- "im"
        riskg[idsInf] <- pmax(1,riskg[idsInf] - 1)
      }
    }
    
    ## Asympt to Infectious 
    nInf <- 0
    idsEligInf <- which(active == 1 & status == "a" & tx.seek == 1)
    nEligInf <- length(idsEligInf)
    
    if (nEligInf > 0) {
      vecInf <- which(rbinom(nEligInf, 1, a.to.i.rate) == 1) 
      if (length(vecInf) > 0) {
        idsInf <- idsEligInf[vecInf]
        nInf <- length(idsInf)
        status[idsInf] <- "i"
        riskg[idsInf] <- pmax(1,riskg[idsInf] - 1)
      }
    }
    
  
  
  ## Latent to Asymptomatically Infectious 
  nAsy       <- 0
  idsEligAsy <- which(active == 1 & status == "e")
  nEligAsy   <- length(idsEligAsy)
  
  if (nEligAsy > 0) {
    vecAsy <- which(rbinom(nEligAsy, 1, e.to.a.rate) == 1)
    if (length(vecAsy) > 0) {
      idsAsy <- idsEligAsy[vecAsy]
      nAsy   <- length(idsAsy)
      status[idsAsy] <- "a"
    }
  }
  
  
  # Update attributes
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "recTime", recTime)
  dat <- set_attr(dat, "testTime", testTime)
  dat <- set_attr(dat, "riskg", riskg)
  
  # Update epidemiology trackers 
  dat <- set_epi(dat, "ea.flow", at, nAsy)
  dat <- set_epi(dat, "ai.flow", at, nInf + nInf_m)
  dat <- set_epi(dat, "ir.flow", at, nRec)
  dat <- set_epi(dat, "test.flow", at, nTest)
  
  prev.cases <- get_epi(dat, "cuml.cases", at-1)
  new.cases <- prev.cases + nTest
  dat <- set_epi(dat, "cuml.cases", at, new.cases)
  
  return(dat)
}


# Track Prevalence & Other Metrics ####
prevalence_msm <- function(dat, at) {
  
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  riskg <- get_attr(dat, "riskg")
  tx.seek <- get_attr(dat, "tx.seek")
  
  
  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)
  
  dat <- set_epi(dat, "num",   at, sum(active == 1))
  dat <- set_epi(dat, "num.e", at, sum(active == 1 & status == "e", na.rm = TRUE))
  dat <- set_epi(dat, "num.a", at, sum(active == 1 & status == "a", na.rm = TRUE))
  dat <- set_epi(dat, "num.i", at, sum(active == 1 & (status == "i" | status == "im"), na.rm = TRUE))
  dat <- set_epi(dat, "num.r", at, sum(active == 1 & status == "r", na.rm = TRUE))
  dat <- set_epi(dat, "num.r.1.tx", at, sum(active == 1 & status == "r" & riskg == 1 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.2.tx", at, sum(active == 1 & status == "r" & riskg == 2 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.3.tx", at, sum(active == 1 & status == "r" & riskg == 3 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.4.tx", at, sum(active == 1 & status == "r" & riskg == 4 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.5.tx", at, sum(active == 1 & status == "r" & riskg == 5 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.6.tx", at, sum(active == 1 & status == "r" & riskg == 6 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.1.notx", at, sum(active == 1 & status == "r" & riskg == 1 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.2.notx", at, sum(active == 1 & status == "r" & riskg == 2 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.3.notx", at, sum(active == 1 & status == "r" & riskg == 3 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.4.notx", at, sum(active == 1 & status == "r" & riskg == 4 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.5.notx", at, sum(active == 1 & status == "r" & riskg == 5 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.r.6.notx", at, sum(active == 1 & status == "r" & riskg == 6 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1", at, sum(active == 1 & status == "v1", na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.1.tx", at, sum(active == 1 & status == "v1" & riskg == 1 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.2.tx", at, sum(active == 1 & status == "v1" & riskg == 2 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.3.tx", at, sum(active == 1 & status == "v1" & riskg == 3 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.4.tx", at, sum(active == 1 & status == "v1" & riskg == 4 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.5.tx", at, sum(active == 1 & status == "v1" & riskg == 5 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.6.tx", at, sum(active == 1 & status == "v1" & riskg == 6 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.1.notx", at, sum(active == 1 & status == "v1" & riskg == 1 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.2.notx", at, sum(active == 1 & status == "v1" & riskg == 2 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.3.notx", at, sum(active == 1 & status == "v1" & riskg == 3 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.4.notx", at, sum(active == 1 & status == "v1" & riskg == 4 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.5.notx", at, sum(active == 1 & status == "v1" & riskg == 5 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v1.6.notx", at, sum(active == 1 & status == "v1" & riskg == 6 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2", at, sum(active == 1 & status == "v2", na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.1.tx", at, sum(active == 1 & status == "v2" & riskg == 1 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.2.tx", at, sum(active == 1 & status == "v2" & riskg == 2 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.3.tx", at, sum(active == 1 & status == "v2" & riskg == 3 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.4.tx", at, sum(active == 1 & status == "v2" & riskg == 4 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.5.tx", at, sum(active == 1 & status == "v2" & riskg == 5 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.6.tx", at, sum(active == 1 & status == "v2" & riskg == 6 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.1.notx", at, sum(active == 1 & status == "v2" & riskg == 1 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.2.notx", at, sum(active == 1 & status == "v2" & riskg == 2 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.3.notx", at, sum(active == 1 & status == "v2" & riskg == 3 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.4.notx", at, sum(active == 1 & status == "v2" & riskg == 4 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.5.notx", at, sum(active == 1 & status == "v2" & riskg == 5 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.v2.6.notx", at, sum(active == 1 & status == "v2" & riskg == 6 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.s", at, sum(active == 1 & status == "s", na.rm = TRUE))
  dat <- set_epi(dat, "num.s.1.tx", at, sum(active == 1 & status == "s" & riskg == 1 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.2.tx", at, sum(active == 1 & status == "s" & riskg == 2 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.3.tx", at, sum(active == 1 & status == "s" & riskg == 3 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.4.tx", at, sum(active == 1 & status == "s" & riskg == 4 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.5.tx", at, sum(active == 1 & status == "s" & riskg == 5 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.6.tx", at, sum(active == 1 & status == "s" & riskg == 6 & tx.seek == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.1.notx", at, sum(active == 1 & status == "s" & riskg == 1 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.2.notx", at, sum(active == 1 & status == "s" & riskg == 2 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.3.notx", at, sum(active == 1 & status == "s" & riskg == 3 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.4.notx", at, sum(active == 1 & status == "s" & riskg == 4 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.5.notx", at, sum(active == 1 & status == "s" & riskg == 5 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "num.s.6.notx", at, sum(active == 1 & status == "s" & riskg == 6 & tx.seek == 0, na.rm = TRUE))
  dat <- set_epi(dat, "prev",  at, sum(active == 1 & status=="e") + 
                   sum(active == 1 & status=="a") + 
                   sum(active == 1 & status=="i") + 
                   sum(active == 1 & status=="im"))
  dat <- set_epi(dat, "inf_hr",  at, sum(active == 1 & 
                                                (status=="a" | status=="i" | status=="im") & 
                                                (riskg == 6 | riskg == 5 | riskg == 4)))
  dat <- set_epi(dat, "sus_hr",  at, sum(active == 1 & 
                                                (status=="s") & 
                                                (riskg == 6 | riskg == 5 | riskg == 4)))
  dat <- set_epi(dat, "N_hr",  at, sum(active == 1 & 
                                           (riskg == 6 | riskg == 5 | riskg == 4)))
  
  # growth rate & doubling time based on actual infections 
  
  nt <- get_epi(dat,"cuml.infs", at)
  n0 <- get_epi(dat, "prev", 1)
  
  r <- log(nt/n0)/at
  dtime <- log(2)/r
  
  dat <- set_epi(dat, "growth.rate", at, r)
  dat <- set_epi(dat, "doubling.time", at, dtime)
  
  # growth rate & doubling time based on Tested cases 
  # nRec is the number of cases that recover per day

  nt2 <- get_epi(dat, "cuml.cases", at)
  n02 <- 1
  
  if (nt2==0) {
    r2=0
    dtime2=NA
  } else {
    r2 <- log(nt2/n02)/at
    if (r2==0) {dtime2=NA} else {dtime2 <- log(2)/r2}
  }
  
  dat <- set_epi(dat, "growth.rate.diag", at, r2)
  dat <- set_epi(dat, "doubling.time.diag", at, dtime2)
  
  # effective reproduction number R(t)
  # calculated as the mean number of secondary cases among those infections that cleared at time t
  
  recTime <- get_attr(dat, "recTime")
  secInfs <- get_attr(dat, "secInfs")
  status <- get_attr(dat, "status")
  
  # pull those who recovered this time step 
  recs <- which(recTime==at)
  
  # how many did they on average infect
  if (length(recs)>0){
    meanSec <- mean(secInfs[recs], na.rm=T)
  } else { meanSec <- NA }
  
  dat <- set_epi(dat, "rt", at, meanSec)
  
  # 5-day moving average
  vec <- c("rt", "growth.rate", "growth.rate.diag", "doubling.time", "doubling.time.diag")
  
  if (at > 5){
    #browser()
    time <- (at-4):at
    
    for (i in 1:length(vec)){
      x <- get_epi(dat, vec[i], time)
      dat <- set_epi(dat, paste0(vec[i], ".avg"), at, mean(x, na.rm=T))
    }
    
  } #else(
  #for (i in length(vec)){
  #  #x <- get_epi(dat, vec[i], at)
  #  dat <- set_epi(dat, paste0(vec[i], ".avg"), at, x)
  #}
  
  # )
  
  return(dat)
}


# Vaccinatation ####
vaccinate_msm <- function(dat, at) {
  
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  vaccTime1 <- get_attr(dat, "vaccTime1")
  vaccTime2 <- get_attr(dat, "vaccTime2")
  riskg <- get_attr(dat, "riskg")
  tx.seek <- get_attr(dat, "tx.seek")
  
  
  vaccination     <- get_param(dat, "vaccination")
  vaccine.multiply    <- get_param(dat, "vaccine.multiply")
  vaccination.proportion.msm    <- get_param(dat, "vaccination.proportion.msm")
  strategy <- get_param(dat, "strategy")
  
  if(strategy == 1 & vaccination == TRUE){
  vaccinate.coverage.1.1 = 478*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.2 = 1573*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.3 = 2729*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.4 = 7298*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.5 = 11342*vaccine.multiply*vaccination.proportion.msm
  vaccinate.coverage.1.6 = 15542*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.7 = 20139*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.8 = 12273*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.9 = 6398*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.10 = 6118*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.11 = 3056*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.12 = 2490*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.13 = 1946*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.14 = 1748*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.15 = 1255*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.16 = 1090*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.17 = 845*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.18 = 706*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.19 = 444*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.20 = 393*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.21 = 359*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.22 = 159*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.23 = 222*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.24 = 148*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.25 = 150*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.26 = 158*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.27 = 115*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.28 = 118*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.29 = 111*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.30 = 98*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.31 = 95*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.32 = 101*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.33 = 87*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.34 = 74*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.35 = 48*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.36 = 53*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.37 = 66*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.1.38 = 66*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.1 = 6*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.2 = 9*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.3 = 24*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.4 = 43*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.5 = 65*vaccine.multiply*vaccination.proportion.msm #sum of first five weeks to account for 28 day gap
  vaccinate.coverage.2.6 = 96*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.7 = 216*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.8 = 287*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.9 = 213*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.10 = 263*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.11 = 2340*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.12 = 7193*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.13 = 14588*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.14 = 9510*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.15 = 5390*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.16 = 2264*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.17 = 1717*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.18 = 1284*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.19 = 853*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.20 = 1040*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.21 = 705*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.22 = 227*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.23 = 342*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.24 = 269*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.25 = 276*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.26 = 214*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.27 = 134*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.28 = 187*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.29 = 182*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.30 = 146*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.31 = 89*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.32 = 117*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.33 = 100*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.34 = 88*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.35 = 78*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.36 = 101*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.37 = 86*vaccine.multiply*vaccination.proportion.msm 
  vaccinate.coverage.2.38 = 70*vaccine.multiply*vaccination.proportion.msm
  }
  
  if(strategy == 2 & vaccination == TRUE){
    vaccinate.coverage.1.1 = 484*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.2 = 1582*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.3 = 2753*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.4 = 7341*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.5 = 10923*vaccine.multiply*vaccination.proportion.msm
    vaccinate.coverage.1.6 = 14056*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.7 = 17602*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.8 = 5219*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.9 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.10 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.11 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.12 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.13 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.14 = 8063*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.15 = 6645*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.16 = 3354*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.17 = 2562*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.18 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.19 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.20 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.21 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.22 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.23 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.24 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.25 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.26 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.27 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.28 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.29 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.30 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.31 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.32 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.33 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.34 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.35 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.36 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.37 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.38 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.1 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.2 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.3 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.4 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.5 = 484*vaccine.multiply*vaccination.proportion.msm #sum of first five weeks to account for 28 day gap
    vaccinate.coverage.2.6 = 1582*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.7 = 2753*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.8 = 7341*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.9 = 6611*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.10 = 6381*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.11 = 5396*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.12 = 9683*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.13 = 16534*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.14 = 3195*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.15 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.16 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.17 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.18 = 1990*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.19 = 1308*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.20 = 1433*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.21 = 1064*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.22 = 386*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.23 = 564*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.24 = 417*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.25 = 426*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.26 = 372*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.27 = 249*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.28 = 305*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.29 = 293*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.30 = 244*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.31 = 184*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.32 = 218*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.33 = 187*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.34 = 162*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.35 = 126*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.36 = 154*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.37 = 152*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.38 = 136*vaccine.multiply*vaccination.proportion.msm 
    
  }
  
  if(strategy == 3 & vaccination == TRUE){
    vaccinate.coverage.1.1 = 242*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.2 = 791*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.3 = 1377*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.4 = 3671*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.5 = 5704*vaccine.multiply*vaccination.proportion.msm
    vaccinate.coverage.1.6 = 7819*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.7 = 10178*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.8 = 6280*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.9 = 3306*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.10 = 3191*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.11 = 2698*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.12 = 4842*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.13 = 8267*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.14 = 5629*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.15 = 3323*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.16 = 1677*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.17 = 1281*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.18 = 995*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.19 = 654*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.20 = 717*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.21 = 532*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.22 = 193*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.23 = 282*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.24 = 209*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.25 = 213*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.26 = 186*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.27 = 125*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.28 = 153*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.29 = 147*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.30 = 122*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.31 = 92*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.32 = 109*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.33 = 94*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.34 = 81*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.35 = 63*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.36 = 77*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.37 = 76*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.1.38 = 68*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.1 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.2 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.3 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.4 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.5 = 0*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.6 = 242*vaccine.multiply*vaccination.proportion.msm #sum of first five weeks, to account for 28 day gap
    vaccinate.coverage.2.7 = 791*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.8 = 1377*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.9 = 3671*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.10 = 5704*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.11 = 7819*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.12 = 10178*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.13 = 6280*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.14 = 3306*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.15 = 3191*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.16 = 2698*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.17 = 4842*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.18 = 8267*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.19 = 5629*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.20 = 3323*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.21 = 1677*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.22 = 1281*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.23 = 995*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.24 = 654*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.25 = 717*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.26 = 532*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.27 = 193*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.28 = 282*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.29 = 209*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.30 = 213*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.31 = 186*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.32 = 125*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.33 = 153*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.34 = 147*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.35 = 122*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.36 = 92*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.37 = 109*vaccine.multiply*vaccination.proportion.msm 
    vaccinate.coverage.2.38 = 94*vaccine.multiply*vaccination.proportion.msm 

  }
    
  

  vaccinate.timeline.1    <- get_param(dat, "vaccinate.timeline.1")
  vaccinate.timeline.2    <- get_param(dat, "vaccinate.timeline.2")
  vaccinate.timeline.3    <- get_param(dat, "vaccinate.timeline.3")
  vaccinate.timeline.4    <- get_param(dat, "vaccinate.timeline.4")
  vaccinate.timeline.5    <- get_param(dat, "vaccinate.timeline.5")
  vaccinate.timeline.6    <- get_param(dat, "vaccinate.timeline.6")
  vaccinate.timeline.7    <- get_param(dat, "vaccinate.timeline.7")
  vaccinate.timeline.8    <- get_param(dat, "vaccinate.timeline.8")
  vaccinate.timeline.9    <- get_param(dat, "vaccinate.timeline.9")
  vaccinate.timeline.10    <- get_param(dat, "vaccinate.timeline.10")
  vaccinate.timeline.11    <- get_param(dat, "vaccinate.timeline.11")
  vaccinate.timeline.12    <- get_param(dat, "vaccinate.timeline.12")
  vaccinate.timeline.13    <- get_param(dat, "vaccinate.timeline.13")
  vaccinate.timeline.14    <- get_param(dat, "vaccinate.timeline.14")
  vaccinate.timeline.15    <- get_param(dat, "vaccinate.timeline.15")
  vaccinate.timeline.16    <- get_param(dat, "vaccinate.timeline.16")
  vaccinate.timeline.17    <- get_param(dat, "vaccinate.timeline.17")
  vaccinate.timeline.18    <- get_param(dat, "vaccinate.timeline.18")
  vaccinate.timeline.19    <- get_param(dat, "vaccinate.timeline.19")
  vaccinate.timeline.20    <- get_param(dat, "vaccinate.timeline.20")
  vaccinate.timeline.21    <- get_param(dat, "vaccinate.timeline.21")
  vaccinate.timeline.22    <- get_param(dat, "vaccinate.timeline.22")
  vaccinate.timeline.23    <- get_param(dat, "vaccinate.timeline.23")
  vaccinate.timeline.24    <- get_param(dat, "vaccinate.timeline.24")
  vaccinate.timeline.25    <- get_param(dat, "vaccinate.timeline.25")
  vaccinate.timeline.26    <- get_param(dat, "vaccinate.timeline.26")
  vaccinate.timeline.27    <- get_param(dat, "vaccinate.timeline.27")
  vaccinate.timeline.28    <- get_param(dat, "vaccinate.timeline.28")
  vaccinate.timeline.29    <- get_param(dat, "vaccinate.timeline.29")
  vaccinate.timeline.30    <- get_param(dat, "vaccinate.timeline.30")
  vaccinate.timeline.31    <- get_param(dat, "vaccinate.timeline.31")
  vaccinate.timeline.32    <- get_param(dat, "vaccinate.timeline.32")
  vaccinate.timeline.33    <- get_param(dat, "vaccinate.timeline.33")
  vaccinate.timeline.34    <- get_param(dat, "vaccinate.timeline.34")
  vaccinate.timeline.35    <- get_param(dat, "vaccinate.timeline.35")
  vaccinate.timeline.36    <- get_param(dat, "vaccinate.timeline.36")
  vaccinate.timeline.37    <- get_param(dat, "vaccinate.timeline.37")
  vaccinate.timeline.38    <- get_param(dat, "vaccinate.timeline.38")
  vaccinate.timeline.over <- get_param(dat, "vaccinate.timeline.over")
  
  duration <- vaccinate.timeline.2 - vaccinate.timeline.1  #time of each coverage period
  
  
  if(at < vaccinate.timeline.1 | (at >= vaccinate.timeline.over) | vaccination == FALSE){
    coverage.1 <- 0
    coverage.2 <- 0
    
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.1 & at < vaccinate.timeline.2){
    coverage.1 <- vaccinate.coverage.1.1
    coverage.2 <- vaccinate.coverage.2.1
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.2 & at < vaccinate.timeline.3){
    coverage.1 <- vaccinate.coverage.1.2
    coverage.2 <- vaccinate.coverage.2.2
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.3 & at < vaccinate.timeline.4){
    coverage.1 <- vaccinate.coverage.1.3
    coverage.2 <- vaccinate.coverage.2.3
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.4 & at < vaccinate.timeline.5){
    coverage.1 <- vaccinate.coverage.1.4
    coverage.2 <- vaccinate.coverage.2.4
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.5 & at < vaccinate.timeline.6){
    coverage.1 <- vaccinate.coverage.1.5
    coverage.2 <- vaccinate.coverage.2.5
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.6 & at < vaccinate.timeline.7){
    coverage.1 <- vaccinate.coverage.1.6
    coverage.2 <- vaccinate.coverage.2.6
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.7 & at < vaccinate.timeline.8){
    coverage.1 <- vaccinate.coverage.1.7
    coverage.2 <- vaccinate.coverage.2.7
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.8 & at < vaccinate.timeline.9){
    coverage.1 <- vaccinate.coverage.1.8
    coverage.2 <- vaccinate.coverage.2.8
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.9 & at < vaccinate.timeline.10){
    coverage.1 <- vaccinate.coverage.1.9
    coverage.2 <- vaccinate.coverage.2.9
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.10 & at < vaccinate.timeline.11){
    coverage.1 <- vaccinate.coverage.1.10
    coverage.2 <- vaccinate.coverage.2.10
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.11 & at < vaccinate.timeline.12){
    coverage.1 <- vaccinate.coverage.1.11
    coverage.2 <- vaccinate.coverage.2.11
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.12 & at < vaccinate.timeline.13){
    coverage.1 <- vaccinate.coverage.1.12
    coverage.2 <- vaccinate.coverage.2.12
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.13 & at < vaccinate.timeline.14){
    coverage.1 <- vaccinate.coverage.1.13
    coverage.2 <- vaccinate.coverage.2.13
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.14 & at < vaccinate.timeline.15){
    coverage.1 <- vaccinate.coverage.1.14
    coverage.2 <- vaccinate.coverage.2.14
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.15 & at < vaccinate.timeline.16){
    coverage.1 <- vaccinate.coverage.1.15
    coverage.2 <- vaccinate.coverage.2.15
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.16 & at < vaccinate.timeline.17){
    coverage.1 <- vaccinate.coverage.1.16
    coverage.2 <- vaccinate.coverage.2.16
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.17 & at < vaccinate.timeline.18){
    coverage.1 <- vaccinate.coverage.1.17
    coverage.2 <- vaccinate.coverage.2.17
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.18 & at < vaccinate.timeline.19){
    coverage.1 <- vaccinate.coverage.1.18
    coverage.2 <- vaccinate.coverage.2.18
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.19 & at < vaccinate.timeline.20){
    coverage.1 <- vaccinate.coverage.1.19
    coverage.2 <- vaccinate.coverage.2.19
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.20 & at < vaccinate.timeline.21){
    coverage.1 <- vaccinate.coverage.1.20
    coverage.2 <- vaccinate.coverage.2.20
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.21 & at < vaccinate.timeline.22){
    coverage.1 <- vaccinate.coverage.1.21
    coverage.2 <- vaccinate.coverage.2.21
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.22 & at < vaccinate.timeline.23){
    coverage.1 <- vaccinate.coverage.1.22
    coverage.2 <- vaccinate.coverage.2.22
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.23 & at < vaccinate.timeline.23){
    coverage.1 <- vaccinate.coverage.1.23
    coverage.2 <- vaccinate.coverage.2.23
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.23 & at < vaccinate.timeline.24){
    coverage.1 <- vaccinate.coverage.1.23
    coverage.2 <- vaccinate.coverage.2.23
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.24 & at < vaccinate.timeline.over){
    coverage.1 <- vaccinate.coverage.1.24
    coverage.2 <- vaccinate.coverage.2.24
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.25 & at < vaccinate.timeline.26){
    coverage.1 <- vaccinate.coverage.1.25
    coverage.2 <- vaccinate.coverage.2.25
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.26 & at < vaccinate.timeline.27){
    coverage.1 <- vaccinate.coverage.1.26
    coverage.2 <- vaccinate.coverage.2.26
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.27 & at < vaccinate.timeline.28){
    coverage.1 <- vaccinate.coverage.1.27
    coverage.2 <- vaccinate.coverage.2.27
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.28 & at < vaccinate.timeline.29){
    coverage.1 <- vaccinate.coverage.1.28
    coverage.2 <- vaccinate.coverage.2.28
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.29 & at < vaccinate.timeline.30){
    coverage.1 <- vaccinate.coverage.1.29
    coverage.2 <- vaccinate.coverage.2.29
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.30 & at < vaccinate.timeline.31){
    coverage.1 <- vaccinate.coverage.1.30
    coverage.2 <- vaccinate.coverage.2.30
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.31 & at < vaccinate.timeline.32){
    coverage.1 <- vaccinate.coverage.1.31
    coverage.2 <- vaccinate.coverage.2.31
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.32 & at < vaccinate.timeline.33){
    coverage.1 <- vaccinate.coverage.1.32
    coverage.2 <- vaccinate.coverage.2.32
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.33 & at < vaccinate.timeline.34){
    coverage.1 <- vaccinate.coverage.1.33
    coverage.2 <- vaccinate.coverage.2.33
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.34 & at < vaccinate.timeline.35){
    coverage.1 <- vaccinate.coverage.1.34
    coverage.2 <- vaccinate.coverage.2.34
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.35 & at < vaccinate.timeline.36){
    coverage.1 <- vaccinate.coverage.1.35
    coverage.2 <- vaccinate.coverage.2.35
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.36 & at < vaccinate.timeline.37){
    coverage.1 <- vaccinate.coverage.1.36
    coverage.2 <- vaccinate.coverage.2.36
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.37 & at < vaccinate.timeline.38){
    coverage.1 <- vaccinate.coverage.1.37
    coverage.2 <- vaccinate.coverage.2.37
  }
  
  if(vaccination == TRUE & at >= vaccinate.timeline.38 & at < vaccinate.timeline.over){
    coverage.1 <- vaccinate.coverage.1.38
    coverage.2 <- vaccinate.coverage.2.38
  }
  
  vaccinate.interval <- get_param(dat, "vaccinate.interval")

  
  ## vaccination first dose
  vacc.first.dose.now <- round(coverage.1/duration)
  nVacc.1 <- 0
  if(at < vaccinate.timeline.5){idsEligVacc <- which(active == 1 & 
      tx.seek == 1 & (riskg == 4 | riskg == 5 | riskg == 6) & (status == "s" | 
      status == "e" | status == "a") & is.na(vaccTime1))}
  if(at >= vaccinate.timeline.5 & at < vaccinate.timeline.7){idsEligVacc <- which(active == 1 & 
      tx.seek == 1 & (riskg == 3 | riskg == 4 | riskg == 5 | 
      riskg == 6) & (status == "s" | status == "e" | status == "a") & 
        is.na(vaccTime1))}
  if(at >= vaccinate.timeline.7){idsEligVacc <- which(active == 1 & 
      tx.seek == 1 & (riskg == 2 | riskg == 3 | riskg == 4 | riskg == 5 | 
      riskg == 6) & (status == "s" | status == "e" | status == "a") & 
        is.na(vaccTime1))}
  
  nEligVacc <- length(idsEligVacc)

  if (nEligVacc > vacc.first.dose.now) {
    vecVacc <- sample(idsEligVacc, vacc.first.dose.now, replace=FALSE) #vector of who is vaccinated
    nVacc.1 <- length(vecVacc)
    vaccTime1[vecVacc] <- at
    ## a and e could take vaccine slots, but will not receive vaccines
    vecVaccSuccess <- which(vaccTime1 == at & status == "s")
    status[vecVaccSuccess] <- "v1"
  }
  
  if (nEligVacc > 0 & nEligVacc < vacc.first.dose.now) {
    nVacc.1 <- nEligVacc
    vaccTime1[idsEligVacc] <- at
    ## a and e could take vaccine slots, but will not receive vaccines
    vecVaccSuccess <- which(vaccTime1 == at & status == "s")
    status[vecVaccSuccess] <- "v1"
  }
  
  ## vaccination second dose
  vacc.second.dose.now <- round(coverage.2/duration)
  nVacc.2 <- 0
  idsEligVacc <- which(active == 1 & (vaccTime1 > 0 & vaccTime1 < (at - 27)))
  nEligVacc <- length(idsEligVacc)
  
  if (nEligVacc > vacc.second.dose.now) {
    vecVacc <- sample(idsEligVacc, vacc.second.dose.now, replace=FALSE) #vector of who is vaccinated
    nVacc.2 <- length(vecVacc)
    vaccTime2[vecVacc] <- at
    ## a and e could take vaccine slots, but will not receive vaccines
    vecVaccSuccess <- which(vaccTime2 == at & status == "v1")
    status[vecVaccSuccess] <- "v2"
  }
  
  if (nEligVacc > 0 & nEligVacc < vacc.second.dose.now) {
    nVacc.2 <- nEligVacc
    vaccTime2[idsEligVacc] <- at
    ## a and e could take vaccine slots, but will not receive vaccines
    vecVaccSuccess <- which(vaccTime2 == at & status == "v1")
    status[vecVaccSuccess] <- "v2"
  }
  
  
  
  
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "vaccTime1", vaccTime1)
  dat <- set_attr(dat, "vaccTime2", vaccTime2)
  
  dat <- set_epi(dat, "vacc.1.flow", at, nVacc.1)
  dat <- set_epi(dat, "vacc.2.flow", at, nVacc.2)
  
  
  return(dat)
}


# create function to extract median and IQR at each time step for all outputs
extract_med_iqr <- function(sim){
  epi <- sim$epi
  timesteps <- sim$control$nsteps
  vars <- names(epi)
  x <- NULL
  d <- matrix(NA, nrow=timesteps, ncol=3)
  
  for (i in vars) {
    #browser()
    for (time in 1:timesteps){
      t <- summary(as.numeric(epi[[i]][time,]))
      vals <- cbind(t[3], t[2], t[5])
      d[time,] <- vals
    }
    #browser()
    colnames(d) <- c(paste0(i, ".med"), paste0(i, ".iqr1"), paste0(i, ".iqr3"))
    x <- cbind(x, d)
    d <- matrix(NA, nrow=timesteps, ncol=3)
  }
  
  x <- x[-1,]
  x <- as.data.frame(x)
  x$Day <- 1:nrow(x)
  
  return(x)
}


# create function to extract mean and IQR at each time step for all outputs
extract_mean_iqr <- function(sim){
  epi <- sim$epi
  timesteps <- sim$control$nsteps
  vars <- names(epi)
  x <- NULL
  d <- matrix(NA, nrow=timesteps, ncol=3)
  
  for (i in vars) {
    #browser()
    for (time in 1:timesteps){
      t <- summary(as.numeric(epi[[i]][time,]))
      vals <- cbind(t[4], t[2], t[5])
      d[time,] <- vals
    }
    #browser()
    colnames(d) <- c(paste0(i, ".mean"), paste0(i, ".iqr1"), paste0(i, ".iqr3"))
    x <- cbind(x, d)
    d <- matrix(NA, nrow=timesteps, ncol=3)
  }
  
  x <- x[-1,]
  x <- as.data.frame(x)
  x$Day <- 1:nrow(x)
  
  return(x)
}



#############
#############
#############

init_tergmLite <- function (dat) {
  num_nw <- ifelse(inherits(dat$nw, "network"), 1, length(dat$nw))
  dat$el <- list()
  dat$control$mcmc.control <- list()
  dat$control$nwstats.formulas <- list()
  for (i in 1:num_nw) {
    nwp <- dat$nwparam[[i]]
    is_tergm <- all(nwp$coef.diss$duration > 1)
    if (is(dat$nw[[i]], "networkDynamic")) {
      nw <- as.networkLite(network.collapse(dat$nw[[i]], 
                                            at = 1))
    }
    else {
      nw <- as.networkLite(dat$nw[[i]])
    }
    dat$el[[i]] <- as.edgelist(nw)
    attributes(dat$el[[i]])$vnames <- NULL
    nwstats_formula_name <- paste(c("nwstats.formula", if (num_nw > 
                                                           1) i), collapse = ".")
    nwstats_formula <- NVL(dat$control[[nwstats_formula_name]], 
                           trim_env(~.))
    if (identical(nwstats_formula, "formation")) 
      nwstats_formula <- nwp$formation
    dat$control$nwstats.formulas[[i]] <- nwstats_formula
    if (is_tergm) {
      mcmc_control_name <- paste(c("mcmc.control.tergm", 
                                   if (num_nw > 1) i), collapse = ".")
      dat$control$mcmc.control[[i]] <- check.control.class("simulate.formula.tergm", 
                                                           "init_tergmLite", NVL(dat$control[[mcmc_control_name]], 
                                                                                 control.simulate.formula.tergm()))
  #    if (dat$control$tergmLite.track.duration == TRUE) {
  #      nw %n% "time" <- 1
  #      nw %n% "lasttoggle" <- cbind(as.edgelist(nw), 
  #                                   1)
  #    }
    }
    else {
      mcmc_control_name <- paste(c("mcmc.control.ergm", 
                                   if (num_nw > 1) i), collapse = ".")
      dat$control$mcmc.control[[i]] <- check.control.class("simulate.formula", 
                                                           "init_tergmLite", NVL(dat$control[[mcmc_control_name]], 
                                                                                 control.simulate.formula()))
    }
    dat$nw[[i]] <- nw
  }
  return(dat)
}
