library(lazymcmc)
library(tidyverse)
library(data.table)
source("model_funcs.R")
n_indiv <- 9

times <- seq(0,30,by=2)

lod <- 2
sd <- 2
tp_lower <- 0
tp_upper <- 7

viral_peak_mean <- 8
viral_peak_sd <- 2
viral_pars <- gamma_pars_from_mean_sd(viral_peak_mean, viral_peak_sd^2)
viral_par1 <- viral_pars[[1]]
viral_par2 <- viral_pars[[2]]

wane_mean <- 0.5
wane_sd <- 0.2
wane_pars <- gamma_pars_from_mean_sd(wane_mean, wane_sd^2)
wane_par1 <- wane_pars[[1]]
wane_par2 <- wane_pars[[2]]

all_tp <- NULL
all_viral_pars <- NULL
all_wanes <- NULL

viral_loads <- NULL
for(i in seq_len(n_indiv)){
  tp <- runif(1, tp_lower, tp_upper)
  viral_peak <- rgamma(1, viral_par1, scale=viral_par2)
  wane <- rgamma(1, wane_par1, scale=wane_par2)
  
  all_tp <- c(all_tp, tp)
  all_viral_pars <- c(all_viral_pars, viral_peak)
  all_wanes <- c(all_wanes, wane)
  
  y <- model_func(times, tp, viral_peak, wane)
  obs <- y + rnorm(length(y), 0, sd)
  obs[obs < lod] <- lod
  obs[times < tp] <- NA
  y[times < tp] <- NA
  
  viral_loads[[i]] <- data.frame(t=times,y=y,obs=obs, i=i)
}


viral_loads_dat <- do.call("rbind", viral_loads)
ggplot(viral_loads_dat) + 
  geom_line(aes(x=t,y=y),col="blue") + 
  geom_point(aes(x=t,y=obs),col="orange") + 
  geom_line(aes(x=t,y=obs),col="orange") + 
  geom_hline(yintercept=lod,linetype="dashed")+
  coord_cartesian(ylim=c(0,15))+
  facet_wrap(~i)

parTab <- read.csv("partab.csv",stringsAsFactors=FALSE)
f <- create_func_indivs(parTab, viral_loads_dat,ver="posterior")
f(parTab$values)

real_pars <- c("viral_peak"=viral_peak, "tp"=tp, "wane"=wane)

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=50000,"save_block"=1000)
startTab <- generate_start_tab(parTab)
x <- run_MCMC(startTab, viral_loads_dat, mcmcPars=mcmcPars1,filename="test1",CREATE_POSTERIOR_FUNC = create_func_indivs,PRIOR_FUNC=NULL,ver="posterior")
chain <- read.csv(x$file)
plot(coda::as.mcmc(chain[chain$sampno > 50000,]))