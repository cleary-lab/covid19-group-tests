library(tidyverse)
library(data.table)
library(coda)
library(MASS)
library(doParallel)
library(ggpubr)
library(rethinking)
library(extraDistr)
library(patchwork)
library(lazymcmc)
library(deSolve)
source("simulation_functions.R")
source("model_funcs.R")
source("model_funcs_multivariate_hinge.R")

set.seed(0)
## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 10000000

## Sample size
n <- 10000000

## Duration of epidemic in days
times <- 0:365

run_name <- "swab"
shift_twane <- 0

## Viral kinetics pars
## Change wd to local path
parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
chains <- load_mcmc_chains(paste0("~/Google Drive/nCoV/pool_samples/chains_",run_name), parTab, 
                           FALSE, 100, burnin=1000000,multi=TRUE)
chain <- as.data.frame(chains[[2]])
chain$wane_mean <- chain$wane_mean + shift_twane

## Simulate the epidemic process
seir_pars <- c("R0"=2.5,"gamma"=1/7,"sigma"=1/6.4,"I0"=100,"recovered0"=0)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)
epidemic_process$plot

## Simulate infection times
infection_times <- simulate_infection_times(n, epidemic_process$overall_prob_infection, 
                                            epidemic_process$incidence)

## Simulate viral loads for the sample population
simulated_data <- simulate_viral_loads_hinge(infection_times, times, chain, parTab)

viral_loads <- simulated_data$viral_loads
viral_loads_tmp <- viral_loads[1:1000,]
viral_loads_melted <- reshape2::melt(viral_loads_tmp)
colnames(viral_loads_melted) <- c("i","t","viral_load")

obs_vl <- simulated_data$obs
obs_vl_tmp <- obs_vl[1:1000,]
obs_vl_melted <- reshape2::melt(obs_vl_tmp)
colnames(obs_vl_melted) <- c("i","t","obs")

p_obs <- ggplot(obs_vl_melted) + 
  geom_tile(aes(x=t,y=i,fill=obs)) + 
  scale_fill_viridis_c() +
  theme_bw() +
  ylab("Individual") + xlab("Time")
pdf(paste0(run_name,"_observed_vls.pdf"),height=5,width=6)
p_obs
dev.off()
p_vl <- ggplot(viral_loads_melted) + 
  geom_tile(aes(x=t,y=i,fill=viral_load)) + 
  scale_fill_viridis_c() +
  theme_bw() +
  ylab("Individual") + xlab("Time")
pdf(paste0(run_name,"_viral_loads.pdf"),height=5,width=6)
p_vl
dev.off()

viral_loads <- simulated_data$viral_loads
viral_loads <- as_tibble(viral_loads)
colnames(viral_loads) <- times
write_csv(viral_loads,paste0("seir_viral_loads_",run_name,".csv"))

obs_dat <- simulated_data$obs
obs_dat <- as_tibble(obs_dat)
colnames(obs_dat) <- times
write_csv(obs_dat,paste0("seir_obs_viral_loads_",run_name,".csv"))

obs_vl <- simulated_data$obs
obs_vl_melted <- reshape2::melt(obs_vl)
colnames(obs_vl_melted) <- c("i","t","obs")
obs_vl_melted <- as_tibble(obs_vl_melted)
obs_vl_melted <- obs_vl_melted %>% mutate(ct = 40 - log2(10)*obs,
                                          ct = pmax(0, ct),
                                          t = t - 1)

plot_t <- c(80, 120, 160)
obs_vl_melted_tmp <- obs_vl_melted %>% 
  filter(t %in% plot_t) %>%
  filter(ct < 40)

t_label <- c("80"="Early","120"="Peak","160"="Late")
obs_vl_melted_tmp$t <- as.character(obs_vl_melted_tmp$t)
obs_vl_melted_tmp$t <- t_label[obs_vl_melted_tmp$t]
obs_vl_melted_tmp$t <- factor(obs_vl_melted_tmp$t, levels=c("Early","Peak","Late"))
pC <- obs_vl_melted_tmp %>% ggplot() +
  geom_histogram(aes(x=ct),binwidth=1,fill="grey70",col="black") +
  facet_wrap(~t,scales="free_y",ncol=1) +
  scale_x_continuous(breaks=seq(0,40,by=5),trans="reverse") +
  geom_hline(yintercept=0) +
  ylab("Count") +
  xlab("Ct") +
  theme_pubr() +
  theme(plot.tag=element_text(vjust=-3)) +
  labs(tag="C")

load(paste0("p_",run_name,".RData"))
load(paste0("p_onset_",run_name,".RData"))

pdf(paste0(run_name,"_sim.pdf"),height=6,width=8)
p_comb | pC
dev.off()

pdf(paste0(run_name,"_onset_sim.pdf"),height=6,width=8)
p_comb_onset | pC
dev.off()

