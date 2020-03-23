library(tidyverse)
library(ggplot2)
## devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)
source("simulation_functions.R")
source("model_funcs.R")

## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 6900000

## Sample size
n <- 1000

## Epidemic growth rate
growth_rate <- 0.1

## Duration of epidemic so far in days
times <- 0:100

## Viral kinetics pars
## Change wd to local path
parTab <- read.csv("partab.csv",stringsAsFactors=FALSE)
chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains", parTab, FALSE, 1, burnin=200000,multi=TRUE)
chain <- as.data.frame(chains[[2]])

pars <- list(viral_peak_mean=mean(chain$viral_peak_mean), viral_peak_sd=mean(chain$viral_peak_sd),
             wane_mean=mean(chain$wane_mean), wane_sd=mean(chain$wane_sd),
             tp_last_day=7)

## Simulate the epidemic process
epidemic_process <- simulate_epidemic_process(population_n,0.1,times)
epidemic_process$plot

## Simulate infection times
infection_times <- simulate_infection_times(n, epidemic_process$overall_prob_infection, epidemic_process$incidence)

## Simulate symptom onset times using default incubation period distribution
onset_times <- simulate_symptom_onsets(infection_times)

## Simulate viral loads for the sample population
viral_loads <- simulate_viral_loads(infection_times, onset_times, times, pars)

## Have a look to see the window of detection for each simulated individual
viral_loads_bin <- apply(viral_loads, 1, function(x) x >= 2)
image(viral_loads_bin)
