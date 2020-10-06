#devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)
library(tidyverse)
library(data.table)
library(coda)
library(MASS)
library(doParallel)
library(ggpubr)
library(deSolve)
library(patchwork)
library(ggthemes)
library(ggsci)
library(extraDistr)
library(rethinking)

setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")
source("functions/simulation_functions.R")
source("functions/model_funcs.R")
source("functions/model_funcs_multivariate_hinge.R")


n_clusters <- 5
cl <- makeCluster(n_clusters,setup_strategy = "sequential")
registerDoParallel(cl)

n_samps <- 100

parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)

## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 12500000

## Sample size
n <- 100000

## Duration of epidemic so far in days
times <- 0:365
## Viral kinetics pars
## Change wd to local path
chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains_sputum/", parTab, FALSE, 10, burnin=1000000,multi=TRUE)
chain <- as.data.frame(chains[[2]])

pars <- list(viral_peak_mean=mean(chain$viral_peak_mean), viral_peak_sd=mean(chain$viral_peak_sd),
             wane_mean=mean(chain$wane_mean), wane_sd=mean(chain$wane_sd),
             tp_last_day=7)

set.seed(123)
## Simulate the epidemic process
#epidemic_process <- simulate_epidemic_process(population_n,0.1,times)
#epidemic_process$plot
seir_pars <- c("R0"=2.5,"gamma"=1/7,"sigma"=1/6.4,"I0"=100,"recovered0"=0)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)

## Simulate infection times
infection_times <- simulate_infection_times(n, epidemic_process$overall_prob_infection, 
                                            epidemic_process$incidence)
infection_times_dat <- tibble(i=seq_along(infection_times), inf_time=infection_times)

## Simulate symptom onset times using default incubation period distribution
onset_times <- simulate_symptom_onsets(infection_times)

## Simulate viral loads for the sample population
simulated_data <- simulate_viral_loads_hinge(infection_times, times, chain, parTab,save_during=FALSE,
                                             save_block=100000,vl_file=vl_file,obs_file=obs_file,par_file=par_file,
                                             add_noise=TRUE)

viral_loads <- simulated_data$obs %>% as_tibble
colnames(viral_loads) <- times
viral_prev <- apply(viral_loads[,1:ncol(viral_loads)], 2, function(x) length(x[x > 2])/length(x))
pA_dat2 <- data.frame(prev=viral_prev,t=unique(epidemic_process$seir_outputs$time))
pA_dat2 <- pA_dat2 %>% mutate(growth=ifelse(t < pA_dat2 %>% filter(prev == max(prev)) %>% pull(t), 
                                            "Epidemic growth","Epidemic decline"))
viral_loads <- viral_loads %>% mutate(i = 1:n()) %>% pivot_longer(-i) 
viral_loads <- viral_loads %>% mutate(name = as.numeric(name))
viral_loads <- viral_loads %>% rename(t=name,log_vl=value)

## For each t
## Number of repeats

res <- foreach(b=c(30, 15, 10, 5, 3, 1),.packages = c("tidyverse")) %dopar% {
  all_res <- NULL
  all_tests <- NULL
  for(run in 1:n_samps) {
    ## Subset individuals to sample and find if actually positive
    tmp_sample <- viral_loads %>% 
      group_by(t) %>% 
      sample_n(30) %>%
      mutate(pool_id = 1:n()) %>%
      mutate(pool = ceiling(pool_id/b)) %>%
      mutate(vl = 10^log_vl - 1) %>%
      mutate(pool_vl = rpois(n(), vl/b)) %>%
      mutate(is_pos = log_vl >= 2,
             is_true_pos = log_vl > 0,
             b=b,
             run=run) 
    
    ## Generate pools and see if pool tests positive
    tmp_pool <- tmp_sample %>%
      group_by(t, pool) %>%
      summarize(pool_vl_tot = sum(pool_vl), .groups="keep") %>%
      mutate(pool_is_pos = pool_vl_tot> 100) %>%
      mutate(n_tests = 1 + as.numeric(pool_is_pos)*b*(b!=1))
    
    tmp_res <- tmp_sample %>% left_join(tmp_pool, by=c("t","pool")) %>%
      mutate(found = is_pos & pool_is_pos)
    
    total_tests <- tmp_pool %>% 
      group_by(t) %>% 
      summarize(n_tests_tot = sum(n_tests),.groups="keep") %>%
      mutate(run=run, b=b)
    
    all_tests[[run]] <- total_tests
    all_res[[run]] <- tmp_res
  }
  all_tests_comb <- do.call("bind_rows", all_tests)
  all_res_comb <- do.call("bind_rows", all_res)
  list(all_tests_comb, all_res_comb)
}
all_tests_res <- lapply(res, function(x) x[[1]])
all_tests_res <- do.call("bind_rows", all_tests_res)

p1 <- all_tests_res %>% group_by(t, b) %>% 
  summarize(mean_tests=mean(n_tests_tot)) %>% 
  left_join(pA_dat2) %>%
  filter(prev < 0.1) %>%
  ggplot() + 
  geom_point(aes(x=prev,y=mean_tests, col=as.factor(b)),size=0.01) + 
  geom_smooth(aes(x=prev,y=mean_tests, col=as.factor(b)),se=FALSE) + 
  geom_hline(yintercept=30) +
  facet_wrap(~growth, ncol=1) +
  theme_bw() + 
  theme(legend.position="bottom") +
  ylab("Number of tests used") +
  xlab("Prevalence")

all_res_res <- lapply(res, function(x) x[[2]])
all_res_res <- do.call("bind_rows", all_res_res)

p2 <- all_res_res %>% filter(is_pos == TRUE) %>%
  group_by(t, run,b) %>%
  summarize(sensitivity=sum(pool_is_pos)/sum(is_true_pos)) %>%
  group_by(t,b) %>%
  summarize(sensitivity=mean(sensitivity)) %>%
  left_join(pA_dat2) %>%
  filter(prev < 0.1) %>%
  ggplot() + 
  geom_point(aes(x=prev,y=sensitivity, col=as.factor(b)),size=0.05) + 
  geom_smooth(aes(x=prev,y=sensitivity, col=as.factor(b)),se=FALSE) + 
  facet_wrap(~growth, ncol=1) +
  scale_y_continuous(limits=c(0.5,1)) +
  theme_bw() +
  theme(legend.position="bottom") +
  ylab("Sensitivity")+
  xlab("Prevalence")

p3 <- all_res_res %>% filter(is_pos == TRUE) %>%
  group_by(t, b) %>%
  summarize(vl=median(log_vl),
            lower=quantile(log_vl, 0.025),
            upper=quantile(log_vl, 0.975)) %>%
  left_join(pA_dat2) %>%
  filter(prev < 0.1) %>%
  ggplot() + 
  geom_point(aes(x=prev,y=vl, col=as.factor(b)),size=0.05) + 
  geom_smooth(aes(x=prev,y=vl, col=as.factor(b)),se=FALSE) + 
  geom_smooth(aes(x=prev,y=lower, col=as.factor(b)),se=FALSE,linetype="dashed",size=0.25) + 
  geom_smooth(aes(x=prev,y=upper, col=as.factor(b)),se=FALSE,linetype="dashed",size=0.25) + 
  facet_wrap(~growth, ncol=1) +
  scale_y_continuous(limits=c(0,10),breaks=seq(0,10,by=1)) +
  theme_bw() +
  theme(legend.position="bottom") +
  ylab("Viral load (median and 95% quantiles)")+
  xlab("Prevalence")

main_p <- (p1 | p2  | p3)
main_p
ggsave("~/Documents/GitHub/covid19-group-tests/madrid/compare_pools_noise_sputum.png",width=8,height=5)
  