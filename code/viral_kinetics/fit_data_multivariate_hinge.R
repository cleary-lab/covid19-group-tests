#devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)
#devtools::load_all("~/Documents/GitHub/lazymcmc_jh/")
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
source("model_funcs_multivariate_hinge.R")

n_clusters <- 5
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

## Choose working directory
#setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")

filenames <- paste0("drosten_data_", 1:n_clusters)

n_indiv <- 9

parTab <- read.csv("pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
## Read in extracted Wolfel data
viral_loads_dat <- read.csv("data/drosten_sputum.csv")
ggplot(viral_loads_dat) + geom_line(aes(x=t,y=obs)) + facet_wrap(~i)
viral_loads_dat <- viral_loads_dat[,c("t","obs","i")]
viral_loads_dat$t <- viral_loads_dat$t


## Run MCMC
mcmcPars1 <- c("iterations"=1000000,"popt"=0.44,"opt_freq"=1000,
               "thin"=100,"adaptive_period"=500000,"save_block"=1000)
mcmcPars2 <- c("iterations"=400000,"popt"=0.234,"opt_freq"=1000,
               "thin"=100,"adaptive_period"=1000000,"save_block"=1000)


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

res <- foreach(i=seq_along(filenames),.packages = c("lazymcmc","rethinking","extraDistr")) %dopar% {
  setwd("/Users/james/Google Drive/nCoV/pool_samples/chains_multivariate_sputum")
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

chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains_multivariate_sputum/", parTab, FALSE, 10, burnin=1000000,multi=TRUE)
plot(chains[[2]])
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)

## Simulate posterior draws over 100 days
## Note that t0 is going to be time of infection
times <- seq(0,100,by=0.1)
dat_fake <- data.frame(t=rep(times, n_indiv),obs=0,i=rep(seq_len(n_indiv),each=length(times)))

## 1000 posterior draws
nsamp <- 1000
samps <- sample(unique(chain$sampno), nsamp)
store_all <- matrix(nrow=nsamp,ncol=nrow(dat_fake))
store_all_obs <- matrix(nrow=nsamp,ncol=nrow(dat_fake))
store_all <- NULL
store_all_obs <- NULL

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
for(i in seq_along(samps)){
  if(i %% 10 == 0) print(i)
  pars <- get_index_par(chain, samps[i])
  
  ## Extract incubation periods
  incus <- round(pars[which(names(pars) == "incu")],t_round)
  incu_tab <- tibble(i=seq_along(incus),t_start=0,t_end=incus + tmax)
  
  ## Solve model from -incu_period to tmax
  dat_fake_tmp <- incu_tab %>%
    gather(var, t, -1) %>%
    dplyr::select(-var) %>%
    group_by(i) %>%
    expand(t=full_seq(t,t_by)) %>%
    mutate(obs=0) %>%
    as.data.frame
  
  ## Create model function and solve
  f_model <- create_func_indivs_multivariate_hinge(parTab, dat_fake_tmp,ver="model",PRIOR_FUNC=prior_func,for_plot=TRUE)
  pred <- f_model(pars)
  
  ## Round to precision we're interested in
  pred$t <- round(pred$t, t_round)
  
  ## Find time period before viral loads are detectable
  ## VLs are detectable tshift days after infection, infection is -incu days from 0
  tshifts <- round(pars[which(names(pars) == "tshift")],t_round)
  tshifts <- tibble(i = seq_len(n_indiv),tshift=tshifts, incu=incus)
  
  pred <- pred %>% 
    group_by(i) %>%
    summarize(min_t=min(t),.groups="keep") %>% ## Get first observation for each individual
    bind_rows(expand_grid(i=seq_len(n_indiv),min_t=tmin)) %>% 
    group_by(i) %>%
    expand(t=full_seq(min_t,t_by)) %>% ## Create time sequence from -30 to infection time
    ungroup() %>%
    mutate(y=0) %>% ## Pre-infection has VL of 0
    bind_rows(pred) %>%
    arrange(i, t) %>%
    mutate(samp = i) %>%
    left_join(tshifts,by="i") %>% 
    mutate(vl_start = - incu + tshift, ## When do VLs start being detectable?
           vl_noise = vl_start >= t) %>% ## Boolean to indicate when vls are detectable
    mutate(noise = rnorm(n(), 0, pars["sd"])) ## Draw noise
  
  store_all <- bind_rows(store_all, pred)
}

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
            mean=mean(y),
            upper=quantile(y,0.975))

quants_obs <- store_all %>% group_by(i, t) %>%
  summarize(lower=quantile(obs,0.025),
            mean=mean(obs),
            upper=quantile(obs,0.975))

rect_dat <- data.frame(x1=-30,x2=0,y1=-5,y2=13)

p_sputum <- ggplot(quants) +
  geom_rect(data=rect_dat,aes(ymin=y1,ymax=y2,xmin=x1,xmax=x2),fill="grey70",alpha=0.5) +
  geom_ribbon(data=quants_obs[quants_obs$t >= 0,], aes(x=t,ymin=lower,ymax=upper),fill="#0072B2",alpha=0.4) +
  geom_ribbon(data=quants[quants$t >= 0,], aes(x=t ,ymin=lower,ymax=upper),fill="#0072B2",alpha=0.6) +
  geom_line(data=quants[quants$t >= 0,],aes(x=t,y=mean),col="#0072B2") +
  
  geom_ribbon(data=quants_obs[quants_obs$t <= 0,], aes(x=t,ymin=lower,ymax=upper),fill="#009E73",alpha=0.4) +
  geom_ribbon(data=quants[quants$t <= 0,],aes(x=t,ymin=lower,ymax=upper),fill="#009E73",alpha=0.6) +
  geom_line(data=quants[quants$t <= 0,], aes(x=t,y=mean),col="#009E73") +
  
  coord_cartesian(xlim=c(-20,50),ylim=c(-2,12)) +
  geom_point(data=viral_loads_dat,aes(x=t,y=obs)) +
  scale_y_continuous(breaks=seq(0,12,by=1)) +
  scale_x_continuous(breaks=seq(-20,50,by=7)) +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=2,linetype="dashed") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_pubr() +
  theme(panel.grid.major = element_line(color="grey80",size=0.25)) +
  ylab("log10 RNA copies / ml (sputum)") +
  xlab("Days post symptom onset") +
  facet_wrap(~i, ncol=3)

png("sputum_fit.png",width=8,height=7,res=300,units="in")
p_sputum
dev.off()


## 1000 posterior draws
nsamp <- 100
samps <- sample(unique(chain$sampno), nsamp)
dat_fake <- tibble(t=seq(0,50,by=0.1),i=1) %>% as.data.frame
dat_all <- NULL
parTab1 <- parTab %>% filter(indiv %in% c(0,1))
## Solve model for each draw
for(i in seq_along(samps)){
  if(i %% 10 == 0) print(i)
  pars <- get_index_par(chain, samps[i])
  
  rhos <- c(pars["rho_viral_wane"])
  
  mus <- c(pars["viral_peak_mean"], pars["wane_mean"])
  sds <- c(pars["viral_peak_sd"], pars["wane_sd"])
  
  R <- diag(2)
  R[upper.tri(R)] <- R[lower.tri(R)]<- rhos
  
  tmp <- rmvnorm2(1, mus, sds, R)
  
  pars_use <- pars[which(parTab$indiv %in% c(0,1))]
  
  pars_use["viral_peak"] <- tmp[1,1]
  pars_use["t_wane"] <- tmp[1,2]
  
  f_model <- create_func_indivs_multivariate_hinge(parTab[parTab$indiv %in% c(0,1),], dat_fake,ver="model",PRIOR_FUNC=prior_func,for_plot=FALSE)
  pred <- f_model(pars_use)
  dat_fake$y <- pred
  dat_fake$samp <- i
  dat_all <- bind_rows(dat_all, dat_fake)
}

dat_all %>% 
  ggplot() + geom_line(aes(x=t,y=y,group=samp)) +
  coord_cartesian(ylim=c(0,12)) +
  scale_y_continuous(breaks=seq(0,12,by=1)) +
  geom_vline(xintercept=5,linetype="dashed") +
  geom_hline(yintercept=0) +
  geom_hline(yintercept=2)


chain_subset <- chain[,colnames(chain) %like% "mean" | colnames(chain) %like% "sd"]


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
