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
setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")
source("functions/simulation_functions.R")
source("functions/model_funcs.R")
source("functions/model_funcs_multivariate_hinge.R")

set.seed(0)
## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 12500000

## Sample size
n <- 12500000
n <- 100000

## Duration of epidemic in days
times <- 0:365

run_name <- "swab"
shift_twane <- 0
par_file <- paste0("sims/used_pars_",run_name,".csv")
vl_file <- paste0("sims/seir_viral_loads_",run_name,".csv")
obs_file <- paste0("sims/seir_obs_viral_loads_",run_name,".csv")

## Viral kinetics pars
## Change wd to local path
parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
chains <- load_mcmc_chains(paste0("~/Google Drive/nCoV/pool_samples/chains_swab"), parTab, 
                           FALSE, 100, burnin=1000000,multi=TRUE)
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)
chain$wane_mean <- chain$wane_mean + shift_twane

## Simulate the epidemic process
seir_pars <- c("R0"=2.5,"gamma"=1/7,"sigma"=1/6.4,"I0"=100,"recovered0"=0)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)
epidemic_process$plot

## Simulate infection times
infection_times <- simulate_infection_times(n, epidemic_process$overall_prob_infection, 
                                            epidemic_process$incidence)
infection_times_dat <- tibble(i=seq_along(infection_times), inf_time=infection_times)
## Simulate viral loads for the sample population
simulated_data <- simulate_viral_loads_hinge(infection_times, times, chain, parTab,save_during=FALSE,
                                             save_block=10000,vl_file=vl_file,obs_file=obs_file,par_file=par_file,
                                             add_noise=TRUE)

obs_vl_melted <- reshape2::melt(simulated_data$obs[seq(1, min(1000000, nrow(simulated_data$obs)),by=1),])
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
  theme(plot.tag=element_text(vjust=-1), 
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.text=element_text(size=6)) +
  labs(tag="D")
pC


load("swab_sim_for_plot.RData")
load("swab_sim_means_for_plot.RData")
load("swab_obs_for_plot.RData")

load("sputum_sim_for_plot.RData")
load("sputum_sim_means_for_plot.RData")
load("sputum_obs_for_plot.RData")

prop_onset_true_swab <- prop_onset_true_swab %>% mutate(samp="Swab")
prop_onset_true_sputum <- prop_onset_true_sputum %>% mutate(samp = "Sputum")

prop_onset_true <- bind_rows(prop_onset_true_swab, prop_onset_true_sputum)

## Swab onset viral load draws
p_draws_onset_swab <- dat_all_onset_swab %>% 
  filter(samp %in% 1:50) %>%
  ggplot() + 
  geom_line(aes(x=t_shifted,y=y,group=samp),size=0.1) +
  geom_line(data=mean_line_onset_swab,aes(x=t_shifted,y=mean_line,col="Mean"),size=0.75) +
  geom_line(data=mean_line_onset_swab,aes(x=t_shifted,y=median_line,col="Median"),size=0.75) +
  coord_cartesian(ylim=c(0,12)) +
  scale_y_continuous(breaks=seq(0,12,by=1)) +
  scale_x_continuous(limits=c(-15,40), expand=c(0,0),breaks=seq(-15,40,by=5)) +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=0,linetype="dashed",col="red") +
  geom_hline(yintercept=2,linetype="dashed") +
  theme_pubr() +
  theme(plot.tag=element_text(vjust=-1), 
        legend.position=c(0.8,0.8),
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.text=element_text(size=6),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()
        )+
  scale_color_manual(values=c("blue","darkorange")) +
  guides(color=guide_legend(title=NULL)) +
  xlab("Days since symptom onset") +
  ylab("log10 RNA copies / swab") +
  labs(tag="A")

p_draws_onset_sputum <- dat_all_onset_sputum %>% 
  filter(samp %in% 1:50) %>%
  ggplot() + 
  geom_line(aes(x=t_shifted,y=y,group=samp),size=0.1) +
  geom_line(data=mean_line_onset_sputum,aes(x=t_shifted,y=mean_line,col="Mean"),size=0.75) +
  geom_line(data=mean_line_onset_sputum,aes(x=t_shifted,y=median_line,col="Median"),size=0.75) +
  coord_cartesian(ylim=c(0,12)) +
  scale_y_continuous(breaks=seq(0,12,by=1)) +
  scale_x_continuous(limits=c(-15,40), expand=c(0,0),breaks=seq(-15,40,by=5)) +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=0,linetype="dashed",col="red") +
  geom_hline(yintercept=2,linetype="dashed") +
  theme_pubr() +
  theme(plot.tag=element_text(vjust=-1), 
        legend.position=c(0.8,0.8),
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.text=element_text(size=6),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()
  )+
  scale_color_manual(values=c("blue","darkorange")) +
  guides(color=guide_legend(title=NULL)) +
  xlab("Days since symptom onset") +
  ylab("log10 RNA copies / ml") +
  labs(tag="B")


prop_detect_onset <-  prop_onset_true %>% ggplot() +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_line(aes(x=t_shifted,y=prop,col=samp),size=1) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1),expand=c(0,0)) +
  scale_x_continuous(limits=c(-15,40), expand=c(0,0),breaks=seq(-15,40,by=5)) +
  xlab("Days since symptom onset") +
  ylab("Proportion detectable") +
  scale_color_lancet() +
  guides(color=guide_legend(title=NULL)) +
  theme_pubr() +
  theme(plot.tag=element_text(vjust=-1),
        legend.position=c(0.8,0.8),
        axis.title=element_text(size=8),
        axis.text=element_text(size=6),
        legend.text=element_text(size=6)
  )+
  labs(tag="C")

p_main <- (p_draws_onset_swab / p_draws_onset_sputum / prop_detect_onset) | pC
p_main

ggsave("figs/FigS4.pdf", p_main,height=7,width=8,units="in")

