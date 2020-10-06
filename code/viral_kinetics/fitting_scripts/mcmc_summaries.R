library(tidyverse)
library(data.table)
library(coda)
library(ggpubr)
library(patchwork)
library(lazymcmc)
setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")

## Viral kinetics pars
## Change wd to local path
parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
chains <- load_mcmc_chains(paste0("~/Google Drive/nCoV/pool_samples/chains_sputum"), parTab, 
                           FALSE, 100, burnin=1000000,multi=TRUE)

mcmc_summary <- summary(chains[[1]])
get_summary_metric <- function(mcmc_sum,par_name){
  par_mean <- mcmc_sum$statistics[which(rownames(mcmc_sum$statistics) == par_name),"Mean"]
  par_lower <- mcmc_sum$quantiles[which(rownames(mcmc_sum$quantiles) == par_name),"2.5%"]
  par_median <- mcmc_sum$quantiles[which(rownames(mcmc_sum$quantiles) == par_name),"50%"]
  par_upper <- mcmc_sum$quantiles[which(rownames(mcmc_sum$quantiles) == par_name),"97.5%"]
  paste0(signif(par_median,3)," (", signif(par_lower,3), "-",signif(par_upper,3),")")
}
get_summary_metric(mcmc_summary,"viral_peak_mean")
get_summary_metric(mcmc_summary,"viral_peak_sd")

get_summary_metric(mcmc_summary,"wane_mean")
get_summary_metric(mcmc_summary,"wane_sd")

get_summary_metric(mcmc_summary,"rho_viral_wane")
get_summary_metric(mcmc_summary,"sd")

