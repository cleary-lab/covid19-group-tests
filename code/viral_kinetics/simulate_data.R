library(tidyverse)
library(ggplot2)
## devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)
source("simulation_functions.R")
source("model_funcs.R")

## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 6900000

## Sample size
n <- 10000

## Epidemic growth rate
growth_rate <- 0.1

## Duration of epidemic so far in days
times <- 0:120

## Viral kinetics pars
## Change wd to local path
parTab <- read.csv("pars/partab.csv",stringsAsFactors=FALSE)
chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains", parTab, FALSE, 1, burnin=200000,multi=TRUE)
chain <- as.data.frame(chains[[2]])

pars <- list(viral_peak_mean=mean(chain$viral_peak_mean), viral_peak_sd=mean(chain$viral_peak_sd),
             wane_mean=mean(chain$wane_mean), wane_sd=mean(chain$wane_sd),
             tp_last_day=7)

## Simulate the epidemic process
epidemic_process <- simulate_epidemic_process(population_n,0.1,times)
epidemic_process$plot

## Simulate infection times
infection_times <- simulate_infection_times(n, epidemic_process$overall_prob_infection, 
                                            epidemic_process$incidence)
#incidence <- rep(1, 120)/120
#infection_times <- simulate_infection_times(n, 0.05, 
#                                            incidence)

## Simulate symptom onset times using default incubation period distribution
onset_times <- simulate_symptom_onsets(infection_times)

## Simulate viral loads for the sample population
viral_loads <- simulate_viral_loads(infection_times, onset_times, times, pars)

## Get viral loads of individuals that were infected
tmp <- viral_loads[which(infection_times > 0),]
par(mfrow=c(1,1))
## Look at their viral loads from time of infection onward
hist(apply(tmp, 1, function(x) mean(x[x > 0])))
## If incidence is growing, we have a skewed distribution because more and more individuals are 
## early in their viral loads

## Have a look to see the window of detection for each simulated individual
viral_loads_bin <- apply(viral_loads, 1, function(x) x >= 2)
#image(viral_loads_bin)

par(mfrow=c(3,2))
wow <- seq(60,120,by=10)
hist(viral_loads[viral_loads[,wow[1]] > 0, wow[1]],xlim=c(0,10),breaks=10)
hist(viral_loads[viral_loads[,wow[2]] > 0, wow[2]],xlim=c(0,10),breaks=10)
hist(viral_loads[viral_loads[,wow[3]] > 0, wow[3]],xlim=c(0,10),breaks=10)
hist(viral_loads[viral_loads[,wow[4]] > 0, wow[4]],xlim=c(0,10),breaks=10)
hist(viral_loads[viral_loads[,wow[5]] > 0, wow[5]],xlim=c(0,10),breaks=10)
hist(viral_loads[viral_loads[,wow[6]] > 0, wow[6]],xlim=c(0,10),breaks=10)
par(mfrow=c(1,1))

