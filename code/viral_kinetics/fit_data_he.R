#devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)
library(tidyverse)
library(data.table)
library(coda)
library(MASS)
library(doParallel)
library(ggpubr)
library(rethinking)
library(extraDistr)

setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")

source("model_funcs.R")
source("model_funcs_onsets.R")

n_clusters <- 3
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

## Choose working directory
setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")

filenames <- paste0("he_data_", 1:3)

ct_dat <- read.csv("data/lau.csv",stringsAsFactors = FALSE)
ct_dat <- ct_dat %>% mutate(oro.N.ct = ifelse(oro.N.ct == "neg",40,oro.N.ct),
                            oro.N.ct = as.numeric(oro.N.ct))
ct_dat <- ct_dat %>% 
  mutate(viral_load = ((oro.N.ct-40)/-log2(10)))
colnames(ct_dat) <- c("i","t","obs","viral_load")
#ct_dat <- ct_dat %>% filter(i %in% 1:25)
n_indiv <- length(unique(ct_dat$i))

parTab <- read.csv("pars/partab_ct.csv",stringsAsFactors=FALSE)

## Run MCMC
mcmcPars1 <- c("iterations"=50000,"popt"=0.44,"opt_freq"=10000,
               "thin"=10,"adaptive_period"=50000,"save_block"=1000)
mcmcPars2 <- c("iterations"=50000,"popt"=0.234,"opt_freq"=10000,
               "thin"=10,"adaptive_period"=50000,"save_block"=1000)

parTab[parTab$names == "viral_peak_sd", "values"] <- 0.9
parTab[parTab$names == "viral_peak_sd", "fixed"] <- 0

parTab[parTab$names == "sd","fixed"] <- 1
parTab[parTab$names == "sd","values"] <- 1
parTab[parTab$names == "lod","values"] <- 40
parTab[parTab$names == "max_titre","values"] <- 0
parTab[parTab$names == "tp","values"] <- 0

## Enumerate out par tab rows
tmp_partab <- parTab[parTab$indiv == 1,]
for(i in 2:n_indiv){
  tmp_partab$indiv <- i
  parTab <- rbind(parTab, tmp_partab)
}


f <- create_func_indivs_new(parTab, ct_dat,ver="model")
ct_dat$predicted <- f(parTab$values)
ggplot(ct_dat) + geom_point(aes(x=t,y=obs)) + geom_line(aes(x=t,y=predicted),col="red") + facet_wrap(~i)


prior_func <- function(pars){
  names(pars) <- parTab$names
  to_test <- pars["viral_peak_sd"]
  prior <- dnorm(to_test, 1, 0.25, TRUE)
  prior
}
prior_func <- NULL
f <- create_func_indivs_new(parTab, ct_dat,ver="posterior")
f(parTab$values)



#res <- foreach(i=seq_along(filenames),.packages = "lazymcmc") %dopar% {
i <- 1
  setwd("/Users/james/Google Drive/nCoV/pool_samples/chains_new_model")
  startTab <- generate_start_tab(parTab)
  ## Generate random starting conditions for the chain
  ## Note that a list of starting tables is created,
  ## one for each temperature chain

  x <- run_MCMC(parTab, ct_dat, mcmcPars=mcmcPars1,filename=filenames[i],
                CREATE_POSTERIOR_FUNC = create_func_indivs_new,PRIOR_FUNC=prior_func,ver="posterior",
                OPT_TUNING=0.2)
  ## Check convergence
  chain <- read.csv(x$file)
  best_pars <- get_best_pars(chain)
  
  ct_dat_expanded <- expand.grid(i=unique(ct_dat$i), t=seq(-10,30,by=0.1))
  ct_dat_expanded <- ct_dat_expanded %>% arrange(i, t)
  f <- create_func_indivs_new(parTab, ct_dat_expanded,ver="model")
  ct_dat_expanded$predicted <- transform_to_ct(f(best_pars),40)
  ct_dat_expanded$predicted <- f(best_pars)
  ggplot(ct_dat_expanded) + 
    geom_point(data=ct_dat,aes(x=t,y=obs)) + 
    geom_line(aes(x=t,y=predicted),col="red") + 
    geom_vline(xintercept=0,linetype="dashed") +
    facet_wrap(~i)
  
  
  chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
  covMat <- cov(chain)
  mvrPars <- list(covMat,1,w=0.8)
  
  ## Start from best location of previous chain
  startTab$values <- best_pars
  ## Run second chain
  output <- run_MCMC(startTab, ct_dat, mcmcPars=mcmcPars2,filename=filenames[i],
                     CREATE_POSTERIOR_FUNC = create_func_indivs_new,PRIOR_FUNC=prior_func,ver="posterior",
                     OPT_TUNING=0.2, mvrPars=mvrPars)
#}

chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains", parTab, FALSE, 1, burnin=200000,multi=TRUE)
#plot(chains[[2]])
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)

## Simulate posterior draws over 100 days
## Note that t0 is going to be time of infection
times <- seq(-10,30,by=0.1)
dat_fake <- data.frame(t=rep(times, n_indiv),obs=0,i=rep(seq_len(n_indiv),each=length(times)))

## 1000 posterior draws
nsamp <- 1000
samps <- sample(unique(chain$sampno), nsamp)
store_all <- matrix(nrow=nsamp,ncol=nrow(dat_fake))
store_all_obs <- matrix(nrow=nsamp,ncol=nrow(dat_fake))

f_model <- create_func_indivs_new(parTab, dat_fake, ver="model",remove_pre_tp = FALSE,draw_incubation_period = FALSE)

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

use_indivs <- unique(ct_dat$i)

p1 <- ggplot() + 
  geom_ribbon(data=quants[quants$i %in% use_indivs,],aes(x=t,ymin=lower,ymax=upper),fill="#0072B2",alpha=0.25) +
  geom_line(data=quants[quants$i %in% use_indivs,], aes(x=t,y=median),col="#009E73") +
  geom_point(data=ct_dat[ct_dat$i %in% use_indivs, ],aes(x=t,y=obs)) + 
  scale_y_continuous(trans="reverse") +
  geom_vline(xintercept=0,linetype="dashed") +
  facet_wrap(~i) +
  ylab("Ct value") +
  xlab("Days since infection (including before)") +
  theme_bw()

pdf("new_fits.pdf",width=12,height=12)
p1
dev.off()


## 1000 posterior draws
nsamp <- 1000
samps <- sample(unique(chain$sampno), nsamp)
ts_subset <- seq(-14,30,by=0.1)

store_all <- matrix(nrow=nsamp,ncol=length(ts_subset))
pars_store <- matrix(nrow=nsamp,ncol=4)

## Solve model for each draw
for(i in seq_along(samps)){
  pars <- get_index_par(chain, samps[i])
  
  mus <- c(pars["viral_peak_mean"], pars["wane_mean"], pars["tp_mean"])
  rhos <- c(pars["rho_viral_wane"],pars["rho_viral_tp"],pars["rho_wane_tp"])
  sds <- c(pars["viral_peak_sd"], pars["wane_sd"], pars["tp_sd"])
  
  R <- diag(length(rhos))
  R[upper.tri(R)] <- R[lower.tri(R)]<- rhos
  
  tmp <- rmvnorm2(1, mus, sds, R)
  td <- runif(1, -6, 0)
  tmp <- c(tmp, td)
  
  pars_store[i,] <- tmp
  viral_peak <- tmp[1]
  wane <- tmp[2]
  tp <- tmp[3]
  y0 <- 40
  #y <- numeric(length(ts_subset))
  #y[ts_subset < tp] <- y0
  #y[ts_subset >= tp] <- viral_peak - wane*(ts_subset[ts_subset >= tp]-tp) + y0
  #

  
  y <- model_func_tonset_fit(ts_subset, tp, td, viral_peak, wane, y0)
  y[y > y0] <- y0
  store_all[i,] <- y
}

store_all_binary <- apply(store_all, 1, function(x) as.numeric(x < 40))
plot(ts_subset, rowSums(store_all_binary)/nsamp,type='l',ylim=c(0,1),ylab="Proportion detectable",xlab="Days since symptom onset")
abline(v=0,col="red")

## Generate quantiles
quants <- t(apply(store_all,2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)))
colnames(quants) <- c("lower","lower_mid","median","upper_mid","upper")
quants <- as_tibble(quants)
quants$t <- ts_subset

ggplot(quants) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid),alpha=0.5) + 
  geom_line(aes(x=t, y=median)) + 
  scale_y_continuous(trans="reverse")


