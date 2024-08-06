##Here, we fit prior stergm outputs to case data 

 
#First, lets load in libraries
library(here)
library(dplyr)
library(tidyverse)
library(janitor)
library(lubridate)
library(readxl)
library(here)
library(httr)
set.seed(12345)

ncdata <- read.csv("NC_data_through_jan25.csv")


#We fit via two procedures. In the first, we eliminate runs that differ
#in final cumulative cases by more than 50 percent 
#compared to empirical data at end of the for the fitting timeline
#so set those limits here

lower_cum_case <- 2300
upper_cum_case <- 4300


#OK, Now load in simulations
prior.index <- 0
runs.per <- 6
datasets <- 166

likelihood.results <- matrix(NA,1000,2)
likelihood.results <- as.data.frame(likelihood.results)
colnames(likelihood.results) <- c("set","likelihood")

row.index <- 1

all.cases <- c()

for(data.index in 1:datasets){
    
  current.data.set <- readRDS(paste(
                              paste("prior_runs_real/abc_nc_prior_output",
                                    prior.index + ((data.index-1) * runs.per) + 1,
                                    prior.index + ((data.index-1) * runs.per) + runs.per,
                                    sep = '_'),"RDATA",sep="."))
  for(run.index in 1:runs.per){
    
    param.index = prior.index + ((data.index-1) * runs.per) + run.index
    
    model.cases <- current.data.set[[run.index]]$epi$test.flow$sim1[1:length(ncdata$mean.cases)]
    model.cases[1] <- 0
    
    disp <- mean(model.cases)/(var(model.cases)/mean(model.cases) - 1)
    
    likelihood.results$set[row.index] <- param.index
    likelihood.results$likelihood[row.index] <- sum((gamma(disp + ncdata$mean.cases)/(factorial(ncdata$mean.cases) * gamma(disp))) * 
                                                    ((model.cases/(model.cases + disp))^ncdata$mean.cases) * 
                                                    (1 + model.cases/disp)^(-disp))
    if(sum(model.cases) < lower_cum_case | sum(model.cases) > upper_cum_case){
      likelihood.results$likelihood[row.index] <- 0
    }
  
  all.cases <- cbind(all.cases,model.cases)  
    
  row.index <- row.index + 1
  
 }
}


prior.index <- 996
runs.per <- 4
datasets <- 1


for(data.index in 1:datasets){
  
  current.data.set <- readRDS(paste(
    paste("prior_runs_real/abc_nc_prior_output",
          prior.index + ((data.index-1) * runs.per) + 1,
          prior.index + ((data.index-1) * runs.per) + runs.per,
          sep = '_'),"RDATA",sep="."))
  for(run.index in 1:runs.per){
    
    param.index = prior.index + ((data.index-1) * runs.per) + run.index
    
    model.cases <- current.data.set[[run.index]]$epi$test.flow$sim1[1:length(ncdata$mean.cases)]
    model.cases[1] <- 0
    
    disp <- mean(model.cases)/(var(model.cases)/mean(model.cases) - 1)
    
    likelihood.results$set[row.index] <- param.index
    likelihood.results$likelihood[row.index] <- sum((gamma(disp + ncdata$mean.cases)/(factorial(ncdata$mean.cases) * gamma(disp))) * 
                                                      ((model.cases/(model.cases + disp))^ncdata$mean.cases) * 
                                                      (1 + model.cases/disp)^(-disp))
    if(sum(model.cases) < lower_cum_case | sum(model.cases) > upper_cum_case){
      likelihood.results$likelihood[row.index] <- 0
    }
    
    all.cases <- cbind(all.cases,model.cases)  
    
    row.index <- row.index + 1
    
  }
}




write.csv(all.cases,paste(
  paste("incident_cases",
        1,
        prior.index + datasets * runs.per,
        sep = '_'),"csv",sep="."))

write.csv(likelihood.results,paste(
  paste("likelihoods",
        1,
        prior.index + datasets * runs.per,
        sep = '_'),"csv",sep="."))


###############

Likelihood_data <- read.csv("likelihoods_1_1000.csv")
Likelihood_data <- Likelihood_data[,-1]

parameter_posterior <- Likelihood_data[sample(seq_len(nrow(Likelihood_data)),100,prob=Likelihood_data$likelihood,replace = TRUE),]
parameter_posterior$sample <- 0
parameter_posterior$trans <- 0
parameter_posterior$eventtrans <- 0
parameter_posterior$behaveadapt <- 0
parameter_posterior$eventtime <- 0
parameter_posterior$import <- 0

parameter_priors <- read.csv("params_set.csv")

for(i in 1:length(parameter_posterior$trans)){
  parameter_posterior$sample[i] <- i
  parameter_posterior$trans[i] <- parameter_priors$trans[parameter_posterior$set[i]]
  parameter_posterior$eventtrans[i] <- parameter_priors$eventtrans[parameter_posterior$set[i]]
  parameter_posterior$behaveadapt[i] <- parameter_priors$behaveadapt[parameter_posterior$set[i]]
  parameter_posterior$eventtime[i] <- parameter_priors$eventtime[parameter_posterior$set[i]]
  parameter_posterior$import[i] <- parameter_priors$import[parameter_posterior$set[i]]
  
}
  

write.csv(parameter_posterior,"posterior_runs_real/parameter_posterior.csv")
