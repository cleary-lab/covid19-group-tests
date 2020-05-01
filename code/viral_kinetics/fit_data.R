#devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)
library(tidyverse)
library(data.table)
library(coda)
library(MASS)
library(doParallel)
library(ggpubr)

setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")

source("model_funcs.R")

n_clusters <- 3
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

## Choose working directory
#setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")

filenames <- paste0("drosten_data_", 1:3)

n_indiv <- 9

parTab <- read.csv("pars/partab.csv",stringsAsFactors=FALSE)
## Read in extracted Wolfel data
viral_loads_dat <- read.csv("data/drosten_data.csv")
ggplot(viral_loads_dat) + geom_line(aes(x=t,y=obs)) + facet_wrap(~i)
viral_loads_dat <- viral_loads_dat[,c("t","obs","i")]

## Run MCMC
mcmcPars1 <- c("iterations"=200000,"popt"=0.44,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=100000,"save_block"=1000)
mcmcPars2 <- c("iterations"=200000,"popt"=0.234,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=100000,"save_block"=1000)

parTab[parTab$names == "viral_peak_sd", "values"] <- 0.9
parTab[parTab$names == "viral_peak_sd", "fixed"] <- 1

prior_func <- function(pars){
  names(pars) <- parTab$names
  to_test <- pars["viral_peak_sd"]
  prior <- dnorm(to_test, 1, 0.25, TRUE)
  prior
}

f <- create_func_indivs(parTab, viral_loads_dat,ver="posterior")
f(parTab$values)
res <- foreach(i=seq_along(filenames),.packages = "lazymcmc") %dopar% {
  setwd("/Users/james/Google Drive/nCoV/pool_samples/chains")
  startTab <- generate_start_tab(parTab)
  ## Generate random starting conditions for the chain
  ## Note that a list of starting tables is created,
  ## one for each temperature chain

  x <- run_MCMC(startTab, viral_loads_dat, mcmcPars=mcmcPars1,filename=filenames[i],
                CREATE_POSTERIOR_FUNC = create_func_indivs,PRIOR_FUNC=prior_func,ver="posterior",
                OPT_TUNING=0.2)
  ## Check convergence
  chain <- read.csv(x$file)
  best_pars <- get_best_pars(chain)
  chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
  covMat <- cov(chain)
  mvrPars <- list(covMat,1,w=0.8)
  
  ## Start from best location of previous chain
  startTab$values <- best_pars
  ## Run second chain
  output <- run_MCMC(startTab, viral_loads_dat, mcmcPars=mcmcPars2,filename=filenames[i],
                     CREATE_POSTERIOR_FUNC = create_func_indivs,PRIOR_FUNC=prior_func,ver="posterior",
                     OPT_TUNING=0.2, mvrPars=mvrPars)
}

chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains", parTab, FALSE, 1, burnin=200000,multi=TRUE)
#plot(chains[[2]])
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)

## Simulate posterior draws over 100 days
## Note that t0 is going to be time of infection
times <- seq(0,100,by=0.1)
dat_fake <- data.frame(t=rep(times, n_indiv),obs=0,i=rep(seq_len(n_indiv),each=length(times)))

## 1000 posterior draws
nsamp <- 10000
samps <- sample(unique(chain$sampno), nsamp)
store_all <- matrix(nrow=nsamp,ncol=nrow(dat_fake))
store_all_obs <- matrix(nrow=nsamp,ncol=nrow(dat_fake))

f_model <- create_func_indivs(parTab, dat_fake, ver="model",remove_pre_tp = FALSE,draw_incubation_period = TRUE)

indivs <- unique(dat_fake$i)
ts <- unique(dat_fake$t)
## Solve model for each draw
for(i in seq_along(samps)){
  pars <- get_index_par(chain, samps[i])
  pred <- f_model(pars)
  pred[which(pred < 0)] <- 0
  ## Also simulate observations (ie. draw Normal noise)
  noises <- rnorm(length(pred),0, pars["sd"])
  obs <- pred + noises
  #obs <- pred
  #obs[which(obs > 0)] <- obs[which(obs>0)] + noises[which(obs > 0)]
  store_all[i,] <- pred
  store_all_obs[i,] <- obs
}

## Generate quantiles
quants <- t(apply(store_all,2, function(x) quantile(x,c(0.025,0.5,0.975),na.rm=TRUE)))
colnames(quants) <- c("lower","median","upper")
quants <- as_tibble(quants)
dat_fake <- as_tibble(dat_fake)
quants <- bind_cols(dat_fake, quants)

quants_obs <- t(apply(store_all_obs,2, function(x) quantile(x,c(0.025,0.5,0.975),na.rm=TRUE)))
colnames(quants_obs) <- c("lower","median","upper")
quants_obs <- as_tibble(quants_obs)
dat_fake <- as_tibble(dat_fake)
quants_obs <- bind_cols(dat_fake, quants_obs)

rect_dat <- data.frame(x1=-7,x2=14,y1=-1,y2=11)

use_indivs <- c(2,3,6,7)

ggplot() + 
  geom_rect(data=rect_dat,aes(ymin=y1,ymax=y2,xmin=x1,xmax=x2),fill="grey70",alpha=0.5) +
  geom_ribbon(data=quants_obs[quants_obs$t >= 14 & quants_obs$i %in% use_indivs,], aes(x=t,ymin=lower,ymax=upper),fill="#0072B2",alpha=0.4) +
  geom_ribbon(data=quants[quants$t >= 14 & quants$i %in% use_indivs,],aes(x=t,ymin=lower,ymax=upper),fill="#0072B2",alpha=0.6) +
  geom_line(data=quants[quants$t >= 14 & quants$i %in% use_indivs,], aes(x=t,y=median),col="#0072B2") +
  geom_ribbon(data=quants_obs[quants_obs$t <= 14 & quants_obs$i %in% use_indivs,], aes(x=t,ymin=lower,ymax=upper),fill="#009E73",alpha=0.4) +
  geom_ribbon(data=quants[quants$t <= 14 & quants$i %in% use_indivs,],aes(x=t,ymin=lower,ymax=upper),fill="#009E73",alpha=0.6) +
  geom_line(data=quants[quants$t <= 14 & quants$i %in% use_indivs,], aes(x=t,y=median),col="#009E73") +
  geom_point(data=viral_loads_dat[viral_loads_dat$i %in% use_indivs, ],aes(x=t+14,y=obs)) + 
  coord_cartesian(ylim=c(0,10),xlim=c(0,50)) +
  scale_y_continuous(breaks=seq(0,10,by=1)) +
  scale_x_continuous(breaks=seq(0,7*7,by=7),labels=seq(-14,7*5,by=7)) +
  geom_vline(xintercept=14,linetype="dashed") +
  geom_hline(yintercept=2,linetype="dashed") +
  theme_pubr() +
  theme(panel.grid.major = element_line(color="grey80",size=0.25)) +
  ylab("log10 RNA titre") +
  xlab("Days post symptom onset") +
  facet_wrap(~i, ncol=1)
  

chain_subset <- chain[,colnames(chain) %like% "mean" | colnames(chain) %like% "sd"]





## Simulate posterior draws over 100 days
## Note that t0 is going to be time of infection
times <- seq(0,100,by=0.1)
dat_fake <- data.frame(t=rep(times, n_indiv),obs=0,i=rep(seq_len(n_indiv),each=length(times)))

## 1000 posterior draws
nsamp <- 10000
samps <- sample(unique(chain$sampno), nsamp)
store_all <- matrix(nrow=nsamp,ncol=nrow(dat_fake))
store_all_obs <- matrix(nrow=nsamp,ncol=nrow(dat_fake))

f_model <- create_func_indivs(parTab, dat_fake, ver="model",remove_pre_tp = FALSE,draw_incubation_period = FALSE)

indivs <- unique(dat_fake$i)
ts <- unique(dat_fake$t)
## Solve model for each draw
for(i in seq_along(samps)){
  pars <- get_index_par(chain, samps[i])
  pred <- f_model(pars)
  ## Also simulate observations (ie. draw Normal noise)
  noises <- rnorm(length(pred),0, pars["sd"])
  obs <- pred + noises
  #pred[which(pred < 0)] <- 0
  #obs <- pred
  #obs[which(obs > 0)] <- obs[which(obs>0)] + noises[which(obs > 0)]
  store_all[i,] <- pred
  store_all_obs[i,] <- obs
}

## Generate quantiles
quants <- t(apply(store_all,2, function(x) quantile(x,c(0.025,0.5,0.975),na.rm=TRUE)))
colnames(quants) <- c("lower","median","upper")
quants <- as_tibble(quants)
dat_fake <- as_tibble(dat_fake)
quants <- bind_cols(dat_fake, quants)

quants_obs <- t(apply(store_all_obs,2, function(x) quantile(x,c(0.025,0.5,0.975),na.rm=TRUE)))
colnames(quants_obs) <- c("lower","median","upper")
quants_obs <- as_tibble(quants_obs)
dat_fake <- as_tibble(dat_fake)
quants_obs <- bind_cols(dat_fake, quants_obs)

mu_labels <- apply(chain[,colnames(chain) == "viral_peak"], 2, function(x){
  res <- signif(quantile(x,c(0.5,0.025,0.975)),3)
  final <- paste0("α = ",res[1], " (",res[2],"-",res[3],")\n")
  final
})
 

tp_labels <- apply(chain[,colnames(chain) == "tp"], 2, function(x){
  res <- signif(quantile(x,c(0.5,0.025,0.975)),3)
  final <- paste0("tp = ",res[1], " (",res[2],"-",res[3],")\n")
  final
})

wane_labels <- apply(chain[,colnames(chain) == "wane"], 2, function(x){
  res <- signif(quantile(x,c(0.5,0.025,0.975)),3)
  final <- paste0("ω = ",res[1], " (",res[2],"-",res[3],")")
  final
})

wow <- paste0(mu_labels,tp_labels, wane_labels)
label_dat <- data.frame(label=wow, i=1:length(wow))

figS1 <- ggplot() + 
  geom_ribbon(data=quants_obs, aes(x=t,ymin=lower,ymax=upper),fill="#0072B2",alpha=0.4) +
  geom_ribbon(data=quants,aes(x=t,ymin=lower,ymax=upper),fill="#0072B2",alpha=0.6) +
  geom_line(data=quants, aes(x=t,y=median),col="#0072B2") +
  geom_point(data=viral_loads_dat,aes(x=t,y=obs)) + 
  geom_text(data=label_dat,aes(x=21,y=8,label=label),parse=FALSE,size=3) +
  coord_cartesian(ylim=c(0,10),xlim=c(0,30)) +
  scale_y_continuous(breaks=seq(0,10,by=1)) +
  scale_x_continuous(breaks=seq(0,7*5,by=7)) +
  geom_hline(yintercept=2,linetype="dashed") +
  theme_pubr() +
  theme(panel.grid.major = element_line(color="grey80",size=0.25),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8)) +
  ylab("log10 RNA titre") +
  xlab("Days post symptom onset") +
  facet_wrap(~i, ncol=3)

cairo_pdf("FigS1.pdf",height=7,width=8)
figS1
dev.off()
