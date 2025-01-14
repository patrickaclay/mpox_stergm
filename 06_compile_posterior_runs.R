#Here we compile results
#specifically we compile sheets for incident cases and cumulative cases by day
# (epi$cuml.cases, epi$test.flow)
# other options are epi$:

# num = number of agents
# num.x = number of individuals in state x (x = e, a, i ,r)
# num.r.x.tx = number of recovered treatment seeking individuals in activity group x (x = 1,2,3,4,5,6)
# num.r.x.notx = number of recovered non-treatment seeking individuals in activity group x (x = 1,2,3,4,5,6)
# num.v1.x.tx = number of single dosed treatment seeking individuals in activity group x (x = 1,2,3,4,5,6)
# num.v1.x.notx = number of single dosed non-treatment seeking individuals in activity group x (x = 1,2,3,4,5,6)
# num.v2.x.tx = number of double dosed treatment seeking individuals in activity group x (x = 1,2,3,4,5,6)
# num.v2.x.notx = number of double dosed non-treatment seeking individuals in activity group x (x = 1,2,3,4,5,6)
# num.s.x.tx = number of susceptible treatment seeking individuals in activity group x (x = 1,2,3,4,5,6)
# num.s.x.notx = number of susceptible non-treatment seeking individuals in activity group x (x = 1,2,3,4,5,6)
# prev = infection prevalence
# inf_hr = number of infectious high activity (groups 4-6) individuals
# sus_hr = number of susceptible high activity (groups 4-6) individuals
# N_hr = number high activity (groups 4-6) individuals
# growth.rate = growth rate of infections
# doubling.time = doubling time of infections
# growth.rate.diag = growth rate of diagnosed cases
# doubling.time.diag = doubling time of diagnosed cases
# cuml.infs = cumulative infections (including undiagnosed)
# se.flow = daily infections
# ea.flow = daily transitions from exposed to asymptomatic-infectious
# ai.flow = daily transitions from asymptomatic-infectious to symptomatic-infectious
# ir/flow = daily transitions from symptomatic-infectious to recovered
# se.flow.qx = daily infections by activity group (x = 1,2,3,4,5,6)
# se.flow.main = daily infections resulting from main partnerships
# se.flow.pers = daily infections resulting from casual relationships
# se.flow.ot = daily infections resultsing from one time relationships
# 

library(here)


#1dpri

last.run <- 0
loops <- 20
runs.per.loop <- 5

  incsheet <- c()
  cumsheet <- c()

for(loops.i in 0:(loops-1)){

current.data.set <- readRDS(paste(
  paste("posterior_runs_real/posterior_nyc_output_best_1dpri",
        (last.run + (loops.i * runs.per.loop) + 1),
        (last.run + (loops.i * runs.per.loop) + runs.per.loop),
        sep = '_'),"RDATA",sep="."))

for(i in 1:runs.per.loop){
  cumsheet<- cbind(cumsheet,current.data.set[[i]]$epi$cuml.cases[,1])
  incsheet<- cbind(incsheet,current.data.set[[i]]$epi$test.flow[,1])
}

}
  
  cumsheet <- as.data.frame(cumsheet)
  incsheet <- as.data.frame(incsheet)
  cumsheet$day <- seq(1,365)
  incsheet$day <- seq(1,365)
  

write.csv(cumsheet, "posterior_runs_real/cum_cases_best_1dpri.csv")
write.csv(incsheet, "posterior_runs_real/inc_cases_best_1dpri.csv")
