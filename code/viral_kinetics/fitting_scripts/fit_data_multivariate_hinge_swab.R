#################################################################
## SCRIPT 1: FIT VIRAL KINETICS MODEL TO SWAB DATA FROM WOLFEL ET AL. 2020
#################################################################
## 1. Reads in the Wolfel et al. SWAB data
## 2. Fits the single-hinge viral kinetics model and plots model fits
## 3. Simulates many viral kinetics curves from the posterior distribution for the kinetics parameters
## 4. Generates plots for the viral kinetics trajectories and proportion detectable WRT infection/symptom onset

#devtools::install_github("jameshay218/lazymcmc") ## Must install this package to re-fit chains
library(lazymcmc)
library(tidyverse)
library(data.table)
library(coda)
library(MASS)
library(doParallel)
library(ggpubr)
library(rethinking)
library(extraDistr)
library(patchwork)

## Choose working directory
setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")

## Load model functions
source("functions/model_funcs.R")
source("functions/model_funcs_multivariate_hinge.R")


run_name <- "swab"
run_name_chain <- "swab"
n_clusters <- 5
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

## Rerun MCMC chains or load existing?
rerun_chains <- FALSE
shift_twane <- 0
filenames <- paste0("drosten_data_", 1:n_clusters)
n_indiv <- 9

## Parameter control table
parTab <- read.csv("pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
## Read in extracted Wolfel data
viral_loads_dat <- read.csv("data/drosten_data.csv")
ggplot(viral_loads_dat) + geom_line(aes(x=t,y=obs)) + facet_wrap(~i)
viral_loads_dat <- viral_loads_dat[,c("t","obs","i")]
viral_loads_dat$t <- viral_loads_dat$t

## Run MCMC
mcmcPars1 <- c("iterations"=1000000,"popt"=0.44,"opt_freq"=1000,
               "thin"=100,"adaptive_period"=500000,"save_block"=1000)
mcmcPars2 <- c("iterations"=4000000,"popt"=0.234,"opt_freq"=1000,
               "thin"=100,"adaptive_period"=1000000,"save_block"=1000)

## Prior on incubation periods
prior_func <- function(pars){
  names(pars) <- parTab$names
  to_test <- pars[which(names(pars) == "incu")]
  prior <- dlnorm(to_test, 1.621, 0.418, TRUE)
  sum(prior)
}
#prior_func <- NULL
parTab[parTab$names == "tp","upper_bound"] <- 100
parTab[parTab$names == "tp","lower_bound"] <- -100
f <- create_func_indivs_multivariate_hinge(parTab, viral_loads_dat,ver="posterior",PRIOR_FUNC=prior_func)
f(parTab$values)
f_model <- create_func_indivs_multivariate_hinge(parTab, viral_loads_dat,ver="model",PRIOR_FUNC=prior_func)
f_model(parTab$values)

if(rerun_chains){
  res <- foreach(i=seq_along(filenames),.packages = c("lazymcmc","rethinking","extraDistr")) %dopar% {
    #setwd(paste0("/Users/james/Google Drive/nCoV/pool_samples/chains_",run_name))
    setwd(paste0("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/chains/chains_",run_name))
    startTab <- generate_start_tab(parTab)
    ## Generate random starting conditions for the chain
    ## Note that a list of starting tables is created,
    ## one for each temperature chain
  
    x <- run_MCMC(parTab, viral_loads_dat, mcmcPars=mcmcPars1,filename=filenames[i],
                  CREATE_POSTERIOR_FUNC = create_func_indivs_multivariate_hinge,PRIOR_FUNC=prior_func,ver="posterior",
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
                       CREATE_POSTERIOR_FUNC = create_func_indivs_multivariate_hinge,PRIOR_FUNC=prior_func,ver="posterior",
                       OPT_TUNING=0.2, mvrPars=mvrPars)
  }
}
## Read in pre-computed MCMC chains for diagnostics
chains <- load_mcmc_chains(paste0("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/chains/chains_",run_name_chain), 
                           parTab, TRUE, 10, burnin=1000000,multi=TRUE)
gelman.diag(chains[[1]])
effectiveSize(chains[[1]])

## Read in pre-computed MCMC chains for analyses
chains <- load_mcmc_chains(paste0("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/chains/chains_",run_name_chain), 
                           parTab, FALSE, 10, burnin=1000000,multi=TRUE)
#plot(chains[[2]])
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)

## Simulate posterior draws over 100 days
## Note that t0 is going to be time of infection
## 1000 posterior draws
nsamp <- 1000
samps <- sample(unique(chain$sampno), nsamp)
store_all <- NULL

times <- seq(0,100,by=0.1)
dat_fake <- data.frame(t=rep(times, n_indiv),obs=0,i=rep(seq_len(n_indiv),each=length(times)))

f_model <- create_func_indivs_multivariate_hinge(parTab, dat_fake,ver="model",PRIOR_FUNC=prior_func,for_plot=TRUE)
f_model(parTab$values)

indivs <- unique(dat_fake$i)
ts <- unique(dat_fake$t)

tmax <- 50
tmin <- -30
fill_numbers <- expand_grid(t=seq(-30,tmax),y=0,i=seq_len(n_indiv))

t_round <- 1
t_by <- 0.1

## Solve model for each draw
for(j in seq_along(samps)){
  if(j %% 100 == 0) print(j)
  pars <- get_index_par(chain, samps[j])
  
  ## Extract incubation periods
  incus <- round(pars[which(names(pars) == "incu")],t_round)
  incu_tab <- tibble(i=seq_along(incus),t_start=0,t_end=incus + tmax)
  
  ## Solve model from -incu_period to tmax
  dat_fake_tmp <- incu_tab %>%
    gather(var, t, -1) %>%
    dplyr::select(-var) %>%
    group_by(i) %>%
    tidyr::expand(t=full_seq(t,t_by)) %>%
    mutate(obs=0) %>% as.data.frame
  
  ## Create model function and solve
  f_model <- create_func_indivs_multivariate_hinge(parTab, dat_fake_tmp,ver="model",PRIOR_FUNC=prior_func,for_plot=TRUE)
  pred <- f_model(pars)
  ## Round to precision we're interested in
  pred <- pred %>% mutate(t = round(t, t_round))
  
  ## Find time period before viral loads are detectable
  ## VLs are detectable tshift days after infection, infection is -incu days from 0
  tshifts <- round(pars[which(names(pars) == "tshift")],t_round)
  tshifts <- tibble(i = seq_len(n_indiv),tshift=tshifts, incu=incus)
  
  pred <- pred %>% 
    group_by(i) %>%
    summarize(min_t=min(t),.groups="keep") %>% ## Get first observation for each individual
    bind_rows(expand_grid(i=seq_len(n_indiv),min_t=tmin)) %>% 
    group_by(i) %>%
    tidyr::expand(t=full_seq(min_t,t_by)) %>% ## Create time sequence from -30 to infection time
    ungroup() %>%
    mutate(y=0) %>% ## Pre-infection has VL of 0
    bind_rows(pred) %>%
    arrange(i, t) %>%
    mutate(samp = j) %>%
    left_join(tshifts,by="i") %>% 
    mutate(vl_start = - incu + tshift, ## When do VLs start being detectable?
           vl_noise = vl_start > t) %>% ## Boolean to indicate when vls are detectable
    mutate(noise = rnorm(n(), 0, pars["sd"])) %>% ## Draw noise
    as.data.table
  
  store_all[[j]] <- pred
}

store_all <- rbindlist(store_all)
store_all$t <- round(store_all$t, t_round)
store_all <- store_all %>% mutate(obs = y + noise,
                     obs=ifelse(vl_noise, y, obs), ## If before infection, just 0
                     obs=ifelse(obs < pars["lod"], pars["lod"], obs), ## If < LOD, then observed as LOD
                     ## If between LOD and limit quantification, would be observed as uniform
                     obs=ifelse(obs > pars["lod"] & obs < pars["limit_quantification"],  
                                runif(n(), pars["lod"],pars["limit_quantification"]),
                                obs))

quants <- store_all %>% group_by(i, t) %>%
  summarize(lower=quantile(y,0.025),
            mean=median(y),
            upper=quantile(y,0.975))

quants_obs <- store_all %>% group_by(i, t) %>%
  summarize(lower=quantile(obs,0.025),
            mean=median(obs),
            upper=quantile(obs,0.975))

rect_dat <- data.frame(x1=-30,x2=0,y1=-5,y2=13)
rect_dat2 <- data.frame(x1=-30,x2=50,y1=-5,y2=0)

## Look at model fits
p_swab <- ggplot(quants) +
  geom_rect(data=rect_dat,aes(ymin=y1,ymax=y2,xmin=x1,xmax=x2),fill="grey70",alpha=0.5) +
  geom_rect(data=rect_dat2,aes(ymin=y1,ymax=y2,xmin=x1,xmax=x2),fill="grey70",alpha=0.5) +
  geom_ribbon(data=quants_obs[quants_obs$t >= 0,], aes(x=t,ymin=lower,ymax=upper),fill="#0072B2",alpha=0.4) +
  geom_ribbon(data=quants[quants$t >= 0,], aes(x=t ,ymin=lower,ymax=upper),fill="#0072B2",alpha=0.6) +
  geom_line(data=quants[quants$t >= 0,],aes(x=t,y=mean),col="#0072B2") +
  
  geom_ribbon(data=quants_obs[quants_obs$t <= 0,], aes(x=t,ymin=lower,ymax=upper),fill="#009E73",alpha=0.4) +
  geom_ribbon(data=quants[quants$t <= 0,],aes(x=t,ymin=lower,ymax=upper),fill="#009E73",alpha=0.6) +
  geom_line(data=quants[quants$t <= 0,], aes(x=t,y=mean),col="#009E73") +
  
  #geom_line(data=store_all %>% filter(vl_noise == FALSE), aes(x=t,y=y,group=samp)) +
  
  coord_cartesian(xlim=c(-14,35),ylim=c(-2,12)) +
  geom_point(data=viral_loads_dat,aes(x=t,y=obs)) +
  scale_y_continuous(breaks=seq(0,12,by=2)) +
  scale_x_continuous(breaks=seq(0,35+14,by=7) - 14) +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=2,linetype="dashed",col="red") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_pubr() +
  theme(panel.grid.major = element_line(color="grey80",size=0.25),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8)) +
  ylab("log10 RNA copies / swab") +
  xlab("Days post symptom onset") +
  facet_wrap(~i, ncol=3)

png(paste0("plots/",run_name,"_fit.png"),width=8,height=6,res=300,units="in")
p_swab
dev.off()


###################################
## Look at posterior draws
###################################
## 100 posterior draws
nsamp <- min(10000, length(unique(chain$sampno)))
samps <- sample(unique(chain$sampno), nsamp,replace=TRUE)
dat_fake <- tibble(t=seq(0,100,by=0.1),i=1) %>% as.data.frame
dat_all <- NULL
parTab1 <- parTab %>% filter(indiv %in% c(0,1))

## Solve model for each draw
for(i in seq_along(samps)){
  if(i %% 1000 == 0) print(i)
  
  choose_indiv <- sample(seq_len(n_indiv),1)
  
  pars <- get_index_par(chain, samps[i])
  
  rhos <- c(pars["rho_viral_wane"])
  mus <- c(pars["viral_peak_mean"], pars["wane_mean"]+shift_twane)
  sds <- c(pars["viral_peak_sd"], pars["wane_sd"])
  
  R <- diag(2)
  R[upper.tri(R)] <- R[lower.tri(R)]<- rhos
  
  tmp <- rmvnorm2(1, mus, sds, R)
  
  pars_use <- pars[which(parTab$indiv %in% c(0,choose_indiv))]
  
  while(tmp[1,2] < 0){
    print("error - sampled waning rate with wrong sign")
    print(tmp)
    tmp <- rmvnorm2(1, mus, sds, R)
  }
  
  dat_fake$i <- choose_indiv
  pars_use["viral_peak"] <- tmp[1,1]
  pars_use["t_wane"] <- tmp[1,2]
  
  f_model <- create_func_indivs_multivariate_hinge(parTab[parTab$indiv %in% c(0,choose_indiv),], dat_fake,ver="model",PRIOR_FUNC=prior_func,for_plot=FALSE)
  pred <- f_model(pars_use)
  
  obs <- pred + rnorm(length(pred), 0, pars_use["sd"])
  
  obs_observable <- which(obs > pars_use["lod"])
  obs_unquantifiable <- which(obs < pars_use["limit_quantification"] & obs > pars_use["lod"])
  obs_below_lod <- which(obs < pars_use["lod"])
  
  obs[obs_below_lod] <- pars_use["lod"]
  obs[obs_unquantifiable] <- runif(length(obs[obs_unquantifiable]), pars_use["lod"], pars_use["limit_quantification"])
  
  dat_fake$t_shifted <- round(dat_fake$t - pars_use["incu"],1)
  #pred[pred < pars_use["lod"]] <- pars_use["lod"]
  dat_fake$y <- pred
  dat_fake$samp <- i
  dat_fake$obs <- obs
  
  dat_all[[i]] <- as_tibble(dat_fake)
}
dat_all <- rbindlist(dat_all)
dat_all <- as_tibble(dat_all) %>% dplyr::select(-i)

## Adding zeros for quantiles
extra_rows <- dat_all %>% group_by(samp) %>% 
  summarize(min_t=min(t_shifted),.groups="keep") %>% 
  bind_rows(expand_grid(samp=unique(dat_all$samp),min_t=-20)) %>%
  group_by(samp) %>%
  tidyr::expand(t_shifted=full_seq(min_t,0.1)) %>% ## Create time sequence from -30 to infection time
  ungroup() %>%
  mutate(y=0,obs=0,t_shifted=round(t_shifted,1)) ## Pre-infection has VL of 0

dat_all_onset <- dat_all %>% dplyr::select(samp, t_shifted, y, obs) %>%
  bind_rows(extra_rows) %>% arrange(samp, t_shifted)

mean_line_onset <- dat_all_onset %>% group_by(t_shifted) %>%
  summarize(median_line=median(y),
            mean_line=mean(y),
            .groups="keep")

prop_onset_true <- dat_all_onset %>% 
  #mutate(detectable = y > pars_use["lod"]) %>%
  mutate(detectable = y > 2) %>%
  group_by(t_shifted) %>%
  summarize(prop=sum(detectable)/n())

prop_onset_obs <- dat_all_onset %>% 
  #mutate(detectable = obs > pars_use["lod"]) %>%
  mutate(detectable = obs > 2) %>%
  group_by(t_shifted) %>%
  summarize(prop=sum(detectable)/n())

prop_onset_true_swab <- prop_onset_true
dat_all_onset_swab <- dat_all_onset
mean_line_onset_swab <- mean_line_onset

save(dat_all_onset_swab,file="sims/swab_sim_for_plot.RData")
save(mean_line_onset_swab,file="sims/swab_sim_means_for_plot.RData")
save(prop_onset_true_swab,file="sims/swab_obs_for_plot.RData")

quant_1 <- prop_onset_true %>% filter(t_shifted > 0) %>% mutate(diff1=abs(0.25-prop)) %>% filter(diff1==min(diff1))
quant_2 <- prop_onset_true %>% filter(t_shifted > 0) %>% mutate(diff1=abs(0.5-prop)) %>% filter(diff1==min(diff1))
quant_3 <- prop_onset_true %>% filter(t_shifted > 0) %>% mutate(diff1=abs(0.75-prop)) %>% filter(diff1==min(diff1))
quant_4 <- prop_onset_true %>% filter(t_shifted > 0) %>% mutate(diff1=abs(0.05-prop)) %>% filter(diff1==min(diff1))


prop_detect_onset <-  prop_onset_obs %>% ggplot() +
  geom_vline(xintercept=0,linetype="dashed") +
  
  ## 50th quantile
  geom_segment(data=quant_2, aes(x=-15,xend=t_shifted,y=0.5,yend=0.5),linetype="dotted",col="grey30") +
  geom_segment(data=quant_2, aes(x=t_shifted,xend=t_shifted,y=0,yend=0.5),linetype="dotted",col="grey30") +

  ## 95th quantile
  geom_segment(data=quant_4, aes(x=-15,xend=t_shifted,y=0.05,yend=0.05),linetype="dotted",col="grey30") +
  geom_segment(data=quant_4, aes(x=t_shifted,xend=t_shifted,y=0,yend=0.05),linetype="dotted",col="grey30") +
  
  geom_line(aes(x=t_shifted,y=prop),col="blue",size=1) +
  geom_line(data=prop_onset_true,aes(x=t_shifted,y=prop),col="red",size=1) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1),expand=c(0,0)) +
  scale_x_continuous(limits=c(-15,35), expand=c(0,0),breaks=seq(-15,35,by=5)) +
  xlab("Days since symptom onset") +
  ylab("Proportion detectable") +
  theme_pubr() +
  theme(panel.grid.major = element_line(color="grey70",size=0.1),
        plot.tag=element_text(vjust=-3)
        )+
  labs(tag="B")


p_draws_onset <- dat_all_onset %>% 
  filter(samp %in% 1:50) %>%
  ggplot() + 
  geom_line(aes(x=t_shifted,y=y,group=samp),size=0.1) +
  geom_line(data=mean_line_onset,aes(x=t_shifted,y=mean_line),col="blue",size=0.75) +
  geom_line(data=mean_line_onset,aes(x=t_shifted,y=median_line),col="darkorange",size=0.75) +
  coord_cartesian(ylim=c(0,12)) +
  scale_y_continuous(breaks=seq(0,12,by=1)) +
  scale_x_continuous(limits=c(-15,35), expand=c(0,0),breaks=seq(-15,35,by=5)) +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=0,linetype="dashed",col="red") +
  geom_hline(yintercept=2,linetype="dashed") +
  theme_pubr() +
  theme(plot.tag=element_text(vjust=-3))+
  xlab("Days since symptom onset") +
  ylab("log10 RNA copies / swab") +
  labs(tag="A")

p_comb_onset <- (p_draws_onset + prop_detect_onset) + plot_layout(ncol=1)
p_comb_onset

## By date of infection
mean_line <- dat_all %>% group_by(t) %>%
  summarize(median_line=median(y),
            mean_line=mean(y),
            .groups="keep")
prop_detect_dat <- dat_all %>% 
  mutate(detectable = y > pars_use["lod"]) %>%
  group_by(t) %>%
  summarize(prop=sum(detectable)/n())

prop_detect_dat <- dat_all %>%  
  mutate(detectable = obs > pars_use["lod"]) %>%
  group_by(t) %>%
  summarize(prop=sum(detectable)/n())

prop_detect <-  ggplot() +
  geom_vline(xintercept=5,linetype="dashed") +
  #geom_line(aes(x=t,y=prop),col="blue",size=1) +
  geom_line(data=prop_detect_dat,
            aes(x=t,y=prop),col="blue",size=1) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1)) +
  #scale_x_continuous(limits=c(0,50), expand=c(0,0),breaks=seq(0,50,by=5)) +
  xlab("Days since infection") +
  ylab("Proportion detectable") +
  theme_pubr() +
  theme(panel.grid.major = element_line(color="grey70",size=0.1),
        plot.tag=element_text(vjust=-3)
  )+
  labs(tag="B")

p_draws <- dat_all %>% 
  filter(samp %in% 1:50) %>%
  ggplot() + 
  geom_line(aes(x=t,y=y,group=samp),size=0.1) +
 # geom_line(data=mean_line,aes(x=t,y=mean_line),col="blue",size=0.75) +
  geom_line(data=mean_line,aes(x=t,y=median_line),col="blue",size=0.75) +
  coord_cartesian(ylim=c(0,12)) +
  scale_y_continuous(breaks=seq(0,12,by=1)) +
  scale_x_continuous(limits=c(0,35),expand=c(0,0),breaks=seq(0,35,by=5)) +
  geom_vline(xintercept=5,linetype="dashed") +
  geom_hline(yintercept=0,linetype="dashed",col="red") +
  geom_hline(yintercept=2,linetype="dashed") +
  theme_pubr() +
  theme(plot.tag=element_text(vjust=-3))+
  xlab("Days since infection") +
  ylab("log10 RNA copies / swab") +
  labs(tag="A")

library(patchwork)
p_comb <- (p_draws + prop_detect) + plot_layout(ncol=1)

save(p_comb,file=paste0("sims/p_",run_name,".RData"))
save(p_comb_onset,file=paste0("sims/p_onset_",run_name,".RData"))

png(paste0("plots/",run_name,"_draws.png"),width=8,height=8,res=300,units="in")
p_comb
dev.off()
png(paste0("plots/",run_name,"_onset_draws.png"),width=8,height=8,res=300,units="in")
p_comb_onset
dev.off()

quant_1
quant_2
quant_3
quant_4

png(paste0("plots/swab_sims.png"),width=8,height=8,res=300,units="in")
p_draws/prop_detect_onset
dev.off()


pdf(paste0("plots/swab_sims.pdf"),width=8,height=8)
p_draws/prop_detect_onset
dev.off()
