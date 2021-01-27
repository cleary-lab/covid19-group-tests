#################################################################
## SCRIPT 4: SIMULATE MANY SPUTUM VIRAL LOAD TRAJECTORIES WITH AND WITHOUT OBSERVATION ERROR
#################################################################
## 1. Reads in the MCMC fits to the Wolfel et al. SPUTUM data
## 2. Simulates a stochastic SEIR model with specified population size
## 3. Simulates viral load trajectories with observation error for each individual, and writes these to disk
## 4. Plots the distributions from the simulations
## NOTE: commented out code near the top is to allow multiple simulations to be run in parallel, if required.
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
library(Matrix)
library(odin)
setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")
source("functions/simulation_functions.R")
source("functions/model_funcs.R")
source("functions/odin_funcs.R")
source("functions/model_funcs_multivariate_hinge.R")

set.seed(0)

## Manage MCMC runs and parallel runs
nsims <- 1
simno <- 1
n_clusters <- 1
#cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
#registerDoParallel(cl)

## Entire population size for simulation - large enough to generate enough viral loads to overcome stochastic effects
population_n <- 12500000

## Sample size
n <- 12500000
## Duration of epidemic in days
times <- 0:365

run_name <- "sputum"
shift_twane <- 0 ## Can make the waning duration average shift_twane days longer
savewd <- "~/Google Drive/nCoV/sims_for_brian/"

## Viral kinetics pars
## Change wd to local path
parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
chains <- load_mcmc_chains(paste0("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/chains/chains_sputum"), parTab, 
                           FALSE, 100, burnin=1000000,multi=TRUE)
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)
chain$wane_mean <- chain$wane_mean + shift_twane

## Simulate the epidemic process
seir_pars <- c("R0"=2.5,"gamma"=1/7,"sigma"=1/6.4,"I0"=100,"recovered0"=0)


## UNCOMMENT THE COMMENTED LINES UP TO L84 TO RUN MULTIPLE SIMULATIONS IN PARALLEL
#res <- foreach(simno=1:nsims,.packages = c("tidyverse","rethinking","odin","data.table")) %dopar% {
# source("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/functions/odin_funcs.R")
par_file <- paste0(savewd,run_name,"/used_pars_",run_name,"_",simno,".csv")
vl_file <- paste0(savewd,run_name,"/seir_viral_loads_",run_name,"_",simno,".csv")
obs_file <- paste0(savewd,run_name,"/seir_obs_viral_loads_",run_name,"_",simno,".csv")

epidemic_process <- simulate_seir_process(population_n,seir_pars,times,ver="normal",beta_smooth=0.5,stochastic = TRUE)

## Simulate infection times
infection_times <- simulate_infection_times(n, epidemic_process$overall_prob_infection, 
                                            epidemic_process$incidence)
infection_times_dat <- tibble(i=seq_along(infection_times), inf_time=infection_times)
## Simulate viral loads for the sample population

## <0.2% of simulated viral loads are 11 or higher
simulated_data <- simulate_viral_loads_hinge(infection_times, times, chain, parTab,save_during=TRUE,
                                             save_block=50000,vl_file=vl_file,obs_file=obs_file,par_file=par_file,
                                             add_noise=TRUE,max_vl=11,simno=simno)
#list(simulated_data, infection_times_dat, epidemic_process)
#}

#simulated_data <- res[[1]][[1]]
#infection_times_dat <- res[[1]][[2]]
#epidemic_process <- res[[1]][[3]]
viral_loads <- simulated_data$viral_loads
viral_loads_tmp <- viral_loads[1:100,]
viral_loads_melted <- reshape2::melt(viral_loads_tmp)
colnames(viral_loads_melted) <- c("i","t","viral_load")
viral_loads_melted <- as_tibble(viral_loads_melted) %>% left_join(infection_times_dat)
viral_loads_melted <- viral_loads_melted %>% filter(inf_time >= 0) %>% arrange(-inf_time,i, t)
viral_loads_melted$i <- match(viral_loads_melted$i, unique(viral_loads_melted$i))


obs_dat <- simulated_data$obs
obs_vl_tmp <- obs_dat[1:100,]
obs_vl_melted <- reshape2::melt(obs_vl_tmp)
colnames(obs_vl_melted) <- c("i","t","obs")
obs_vl_melted <- as_tibble(obs_vl_melted) %>% left_join(infection_times_dat)
obs_vl_melted <- obs_vl_melted %>% filter(inf_time >= 0) %>% arrange(-inf_time,i, t)
obs_vl_melted$i <- match(obs_vl_melted$i, unique(obs_vl_melted$i))


p_obs_traj <- ggplot(obs_vl_melted %>% mutate(t = t-inf_time) %>% filter(t >0 & t < 100)%>%
                       filter(i %in% sample(unique(obs_vl_melted$i), pmin(length(unique(obs_vl_melted$i)),10)))) + 
  geom_line(aes(x=t,y=obs,group=i,col=as.factor(i)),size=0.25) + 
  theme_bw() +
  ylab("log10 viral load") + xlab("Time")
p_obs_traj

p_obs <- ggplot(obs_vl_melted) + 
  geom_tile(aes(x=t,y=i,fill=obs)) + 
  scale_fill_viridis_c() +
  theme_bw() +
  ylab("Individual") + xlab("Time")

pdf(paste0(savewd,run_name,"_observed_vls.pdf"),height=5,width=6)
p_obs
dev.off()
png(paste0(savewd,run_name,"_observed_vls.png"),height=5,width=6,res=300,units="in")
p_obs
dev.off()

p_vl <- ggplot(viral_loads_melted) + 
  geom_tile(aes(x=t,y=i,fill=viral_load)) + 
  scale_fill_viridis_c() +
  theme_bw() +
  ylab("Individual") + xlab("Time")

obs_vl_melted <- reshape2::melt(simulated_data$obs[seq(1, min(1000000, nrow(simulated_data$obs)),by=1),])
colnames(obs_vl_melted) <- c("i","t","obs")
obs_vl_melted <- as_tibble(obs_vl_melted)
obs_vl_melted <- obs_vl_melted %>% mutate(ct = 40 - log2(10)*obs,ct = pmax(0, ct),t = t - 1)

plot_t <- c(97, 127, 157)
obs_vl_melted_tmp <- obs_vl_melted %>% 
  filter(t %in% plot_t) %>%
  filter(ct < 40)

t_label <- c("97"="Early","127"="Peak","157"="Late")
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

## Source fit_data_multivariate_hinge_swab first
load(paste0("sims/p_",run_name,".RData"))
load(paste0("sims/p_onset_",run_name,".RData"))

pdf(paste0(savewd,run_name,"_sim.pdf"),height=6,width=8)
p_comb | pC
dev.off()

pdf(paste0(savewd,run_name,"_onset_sim.pdf"),height=6,width=8)
p_comb_onset | pC
dev.off()

p1 <- epidemic_process$incidence_plot + geom_vline(xintercept=plot_t,linetype="dashed")

ggsave( paste0(savewd,run_name,"_observed_trajectories.png"),plot=p_obs_traj,height=5,width=6,dpi=300,units="in")
ggsave(paste0(savewd,run_name,"_inc_plot.png"),plot=p1,height=4,width=7,dpi=300,units="in")


