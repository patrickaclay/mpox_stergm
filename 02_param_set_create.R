#gerate fit parameter priors

library(lhs)
library(here)

h <- 1000 #choose number of points
lhs<-maximinLHS(h,5) #latin hypercube sampling

#fit min and max for uniform priors
#per exposure transmission probability
trans.min <- 0.1
trans.max <- 0.9
#transmission rate in extra-network super spreader events
eventtrans.min <- 0
eventtrans.max <- 0.5
#maximum reduction in sexual partnership rates due to behavioral adaptation
behaveadapt.min <- 0
behaveadapt.max <- 0.9
#how long does the period of importations and superspreader events last
eventtime.min <- 49
eventtime.max <- 111+16 #end of peak + infection interval
#importation rate per day
import.min <- 0
import.max <- 3



#Now we can generate a "parameter set" by rescaling our simulated 
#latin hypercube sample.
params.set <- cbind(
  trans = lhs[,1]*(trans.max-trans.min)+trans.min,
  eventtrans = lhs[,2]*(eventtrans.max-eventtrans.min)+eventtrans.min,
  behaveadapt = lhs[,3]*(behaveadapt.max-behaveadapt.min)+behaveadapt.min,
  eventtime = lhs[,4]*(eventtime.max-eventtime.min)+eventtime.min,
  import = lhs[,5]*(import.max-import.min)+import.min)

params.set <- as.data.frame(params.set)

write.csv(params.set, "params_set.csv")

