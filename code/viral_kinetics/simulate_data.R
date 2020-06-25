library(tidyverse)
library(ggplot2)
library(deSolve)
library(ggpubr)
## devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)
source("simulation_functions.R")
source("model_funcs.R")

set.seed(0)
## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 6900000

## Sample size
n <- 10000

## Duration of epidemic in days
times <- 0:365

## Viral kinetics pars
## Change wd to local path
parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
#chains <- load_mcmc_chains("~/Documents/GitHub/covid19-group-tests/chains", parTab, FALSE, 1, burnin=200000,multi=TRUE)
chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains_multivariate_sputum/", parTab, FALSE, 10, burnin=1000000,multi=TRUE)
chain <- as.data.frame(chains[[2]])

pars <- list(viral_peak_mean=mean(chain$viral_peak_mean), viral_peak_sd=mean(chain$viral_peak_sd),
             wane_mean=mean(chain$wane_mean), wane_sd=mean(chain$wane_sd),
             tp_last_day=7)

## Simulate the epidemic process
#epidemic_process <- simulate_epidemic_process(population_n,0.1,times)
#epidemic_process$plot
seir_pars <- c("R0"=2.5,"gamma"=1/7,"sigma"=1/6.4,"I0"=100,"recovered0"=0)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)
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
simulated_data <- simulate_viral_loads_hinge(infection_times, times, chain, parTab)

simulated_data <- simulate_viral_loads(infection_times, onset_times, times, pars)
simulated_data <- simulate_viral_loads_by_inf2(infection_times[1:nrow(delays)], times,delays,0,3)


viral_loads <- simulated_data$viral_loads
viral_loads <- simulated_data
viral_loads_tmp <- viral_loads[1:100,]

viral_loads_melted <- reshape2::melt(viral_loads_tmp)
colnames(viral_loads_melted) <- c("i","t","viral_load")


viral_loads_melted$ct <- 40 - log2(10)*(viral_loads_melted$viral_load-LOD) + rnorm(nrow(viral_loads_melted),0, 5)
viral_loads_melted[viral_loads_melted$viral_load <= 0, "ct"] <- 40
viral_loads_melted[viral_loads_melted$ct > 40, "ct"] <- 40
viral_loads_melted %>% filter(ct < 40) %>% ggplot() + geom_histogram(aes(x=ct),binwidth=1)

p <- ggplot(viral_loads_melted) + 
  geom_tile(aes(x=t,y=i,fill=ct)) + 
  scale_fill_viridis_c(trans="reverse") +
  theme_bw() +
  ylab("Individual") + xlab("Time")
#pdf("viral_loads.pdf",height=5,width=6)
#p
#dev.off()

viral_loads <- simulated_data$viral_loads
viral_loads <- as_tibble(viral_loads)
colnames(viral_loads) <- times
#write_csv(viral_loads,"seir_viral_loads_no_tail.csv")

#simulated_data$kinetics_pars[which(wow > 100),]

## Get viral loads of individuals that were infected
tmp <- viral_loads[which(infection_times > 0),]
par(mfrow=c(1,1))
## Look at their viral loads from time of infection onward
hist(apply(tmp, 1, function(x) mean(x[x > 0])))
## If incidence is growing, we have a skewed distribution because more and more individuals are 
## early in their viral loads

