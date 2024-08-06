##set up network

library(EpiModel)

## Network characteristics ---------------------------------------------------

# Original model was built to have ability to have different behavior among White and Black MSM in Atlanta
# While there are parameter names that for W and B characteristics, they are the same, race not actually in the model 
# Differences in behavior among come from sexual activity groups not defined by race 

#You will want to change the population size 
#if make too small, might get instability
#in that case, will want to combine top two activity groups
pop.size <- 167000

#### params ####
num.B <- num.W <- pop.size/2

# mean/pers degree distributions matrices.
deg.mp.B <- deg.mp.W <-
  (matrix(c(0.297, 0.150, 0.132,
            0.307, 0.059, 0.055), byrow = TRUE, nrow = 2))



# Instant rates
mdeg.inst.B <- mdeg.inst.W <-
  (matrix(c(0.1092682/7, 0.1920175/7, 0.374595/7,
            0.06154875/7, 0.2049258/7, 0.3939091/7), byrow = TRUE, nrow = 2))

top5 <- TRUE 
# Instant rates (with additional partners based on highest risk group)
# We add to each category weighted on original dist of instant rates by other partnerships 
if (top5==TRUE) {
  mdeg.inst.B <- mdeg.inst.W <-
    (matrix(c(0.1092682/7, 0.1920175/7, 0.374595/7,
              0.06154875/7, 0.2049258/7, 0.3939091/7), byrow = TRUE, nrow = 2))
}
# Quintile distribution of overall AI rates
# adding risk group for top 5% based on Weiss et al. 2020
qnts.W <- qnts.B <- c(0.0000,
                      0.007/7,
                      0.038/7,
                      0.071/7,
                      0.221/7)

if (top5==TRUE) {
  qnts.W <- qnts.B <- c(0.005285245/7,
                        0.0582199/7,
                        0.1894323/7,
                        0.4456711/7,
                        1.103787/7,
                        2.598077/7)
}
# Proportion in same-race partnerships (main, casl, inst)
prop.hom.mpi.B <- prop.hom.mpi.W <- (c(0.9484, 0.9019, 0.9085) +
                                       c(0.9154, 0.8509, 0.8944))/2

# Mean sqrt age diffs (main, casl, inst)
sqrt.adiff.BB <- c(0.519, 0.875, 0.801)
sqrt.adiff.BW <- c(0.519, 0.875, 0.801)
sqrt.adiff.WW <- c(0.519, 0.875, 0.801)

# Mean durations
rates.main <- mean(c(1/407,
                     1/407,
                     1/407))
rates.pers <- mean(c(1/166,
                     1/166,
                     1/166))

durs.main <- 1/rates.main 
durs.pers <- 1/rates.pers 

# Age-sex-specific mortality rates
# Assuming leave pool at 40 then
# But the mortality rate exists but is 0 for ages 1-17 (just for indexing purposes in simulation)
# demography not currently implemented 

ages <- 18:39
asmr.B <- c(rep(0, 17),
            1-(1-c(rep(0.00159, 7),
                   rep(0.00225, 10),
                   rep(0.00348, 5)))^(1/(365/1)), 1)

asmr.W <- c(rep(0, 17),
            1-(1-c(rep(0.00103, 7),
                   rep(0.00133, 10),
                   rep(0.00214, 5)))^(1/(365/1)), 1)

# I, R, V role frequencies
role.B.prob <- role.W.prob <- (c(0.242, 0.321, 0.437))

time.unit = 1 # days

#mdeg.hiv.pers = c(0.5614588, 0.807601)
#mdeg.hiv.inst = c(0.02180057, 0.03882583)


## Create Target Statistics -------------------------------------------
source("00_setup_functions.R")

st <- calc_nwstats_msm(
  method = 1, #1 for 1 race models, 2 for two race models (so use 2)
  top5 = top5,
  time.unit = time.unit,
  num.B = num.B,
  num.W = num.W,
  deg.mp.B = deg.mp.B,
  deg.mp.W = deg.mp.W,
  mdeg.inst.B = mdeg.inst.B,
  mdeg.inst.W = mdeg.inst.W,
  qnts.B = qnts.B,
  qnts.W = qnts.W,
  prop.hom.mpi.B = prop.hom.mpi.B,
  prop.hom.mpi.W = prop.hom.mpi.W,
  balance = "mean",
  sqrt.adiff.BB = sqrt.adiff.BB,
  sqrt.adiff.WW = sqrt.adiff.WW,
  sqrt.adiff.BW = sqrt.adiff.BW,
  diss.main = ~offset(edges),
  diss.pers = ~offset(edges),
  durs.main = durs.main,
  durs.pers = durs.pers,
  ages = ages,
  asmr.B = asmr.B,
  asmr.W = asmr.W,
  role.B.prob = role.B.prob,
  role.W.prob = role.W.prob)#,
  #mdeg.hiv.pers = mdeg.hiv.pers,
  #mdeg.hiv.inst = mdeg.hiv.inst)


# 1. Main Model -----------------------------------------------------------

# Initialize network
nw.main <- base_nw_msm(st)

# Assign degree
nw.main <- assign_degree(nw.main, deg.type = "pers", nwstats = st)

# Formulas
formation.m <- ~edges +
  nodefactor(~deg.pers>0) +
  absdiff("sqrt.age") +
  offset(nodematch("role.class", diff = TRUE, levels = 1:2))

new.stats.m <- c(st$stats.m[1],
                 sum(st$stats.m[2:3]),
                 st$stats.m[4])
# Fit model
fit.m <- netest(nw.main,
                formation = formation.m,
                coef.form = c(-Inf, -Inf),
                target.stats = new.stats.m,
                coef.diss = st$coef.diss.m,
                constraints = ~bd(maxout = 1),
                set.control.ergm = control.ergm(MCMC.interval=2048, #2x the regular interval length
                                                MCMC.samplesize=2048))


saveRDS(fit.m, "fits/main.rds")

# 2. Casual Model ---------------------------------------------------------

# Initialize network
nw.pers <- nw.main

# Assign degree
nw.pers %v% "deg.main" <- get_degree(fit.m$newnetwork)

# Formulas
formation.p <- ~edges +
  nodefactor("deg.main") +
  absdiff("sqrt.age") +
  #nodefactor("hiv", levels=2) +
  degree(0) + 
  offset(nodematch("role.class", diff = TRUE, levels = 1:2))


deg.0.cas <- colSums(deg.mp.B)[1] * pop.size
new.stats.p <- c(st$stats.p[-3], deg.0.cas) # take out concurrent term, add deg0

# Fit model
fit.p <- netest(nw.pers,
                formation = formation.p,
                coef.form = c(-Inf, -Inf),
                target.stats = new.stats.p,
                coef.diss = st$coef.diss.p,
                constraints = ~bd(maxout = 3) + sparse,
                set.control.ergm = control.ergm(MCMC.interval=2048*2,
                                                MCMC.samplesize=2048*2
                  ))

saveRDS(fit.p, "fits/pers.rds")

# Fit inst model ----------------------------------------------------------

# Initialize network
nw.inst <- nw.pers

# Update pers degree
nw.inst %v% "deg.pers" <- get_degree(fit.p$newnetwork)


# adjust nodefactor for pers degree
tab <- table(nw.inst %v% "deg.main", nw.inst %v% "deg.pers") / 167000
weights.deg.pers.2.3 <- tab[,2:3] / rowSums(tab[,2:3])
inst.rates <- pop.size * deg.mp.W * mdeg.inst.W * time.unit
inst.rates.2 <- inst.rates[,3] * weights.deg.pers.2.3
inst.rates.new <- cbind(inst.rates[,1:2], inst.rates.2)

#new.stats.i <- c(new.stats.i[1],
#                 inst.rates.new[-1],
#                 new.stats.i[7:13])

new.stats.i <- c(st$stats.i[1],
                 inst.rates.new[-1],
                 st$stats.i[7:12])



# Formulas
formation.i <- ~edges +
  nodefactor(c("deg.main", "deg.pers")) +
  nodefactor("riskg", levels = -2) +
  absdiff("sqrt.age") +
  #nodefactor("hiv", levels=2) +
  offset(nodematch("role.class", diff = TRUE, levels = 1:2))

# Fit model
fit.i <- netest(nw.inst,
                formation = formation.i,
                target.stats = new.stats.i,
                coef.form = c(-Inf, -Inf),
                coef.diss = dissolution_coefs(~offset(edges), 1),
                constraints = ~sparse,
                # stochastic-approx: faster than MCMLE or MPLE for estimation on big networks
                # esp with lots of terms/constraints
                # (we still use MCMLE/MPLE for initializating the starting network)
                # can underestimate SEs, but we use here because it's SO.MUCH.FASTER.
                set.control.ergm = control.ergm(#main.method = "Stochastic-Approximation",
                                                MCMC.interval = 2048*2,
                                                MCMC.samplesize = 2048*2))

saveRDS(fit.i, "fits/inst.rds")

est <- list(fit.m, fit.p, fit.i) 



save(est, file = here("fits/fit_167k_highestrisk_hiv_edp.rda"))
#save(est, file = "~/OneDrive - CDC/mpox_nyc/fit_167k_highestrisk_hiv_edp.rda")

#load(here("fits", "fit_167k_highestrisk_hiv.rda"))
#fit.m <- est[[1]]
#fit.p <- est[[2]]
#fit.i <- est[[3]]

#check fit
m.dx <- netdx(fit.m, dynamic=TRUE, nsims=5, nsteps=500)
p.dx <- netdx(fit.p, dynamic=TRUE, nsims=5, nsteps=500)
i.dx <- netdx(fit.i, dynamic=FALSE, nsims=500)
i.dx2 <- netdx(fit.i, dynamic=FALSE, nsims=100, nwstats.formula = ~nodefactor("riskg", levels=1:6))
