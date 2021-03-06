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
source("functions/odin_funcs.R")
source("functions/model_funcs_multivariate_hinge.R")

cols <- ggsci::pal_lancet()(10)
n_indiv <- 9

set.seed(123)

parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
## Read in extracted Wolfel data

viral_loads_dat <- read.csv("~/Google Drive/nCoV/pool_samples/drosten_data.csv")
ggplot(viral_loads_dat) + geom_line(aes(x=t,y=obs)) + facet_wrap(~i)
viral_loads_dat <- viral_loads_dat[,c("t","obs","i")]

## Run MCMC
chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains_swab/", parTab, FALSE, 10, burnin=1000000,multi=TRUE)
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)

## Simulate posterior draws over 100 days
## Note that t0 is going to be time of infection
times <- seq(0,100,by=0.1)
dat_fake <- data.frame(t=rep(times, n_indiv),obs=0,i=rep(seq_len(n_indiv),each=length(times)))

## 10000 posterior draws
nsamp <- 10000
samps <- sample(unique(chain$sampno), nsamp)
#store_all <- matrix(nrow=nsamp,ncol=nrow(dat_fake))
#store_all_obs <- matrix(nrow=nsamp,ncol=nrow(dat_fake))

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
rerun_pc <- TRUE
if(rerun_pc){

  store_all <- NULL
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
}

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

use_indivs <- c(2,3,7)

## Revision: plot bar for points when below limit of quantification
viral_loads_dat_bars <- viral_loads_dat[viral_loads_dat$i %in% use_indivs, ]
viral_loads_dat_bars <- viral_loads_dat_bars %>% filter(obs < 2)


pC <- ggplot() + 
  geom_rect(data=rect_dat,aes(ymin=y1,ymax=y2,xmin=x1,xmax=x2),fill="grey70",alpha=0.5) +
  geom_rect(data=rect_dat2,aes(ymin=y1,ymax=y2,xmin=x1,xmax=x2),fill="grey70",alpha=0.5) +
  geom_ribbon(data=quants_obs[quants_obs$t >= 0 & quants_obs$i %in% use_indivs,], aes(x=t,ymin=lower,ymax=upper),fill="#ED0000FF",alpha=0.4) +
  geom_ribbon(data=quants[quants$t >= 0 & quants_obs$i %in% use_indivs,], aes(x=t ,ymin=lower,ymax=upper),fill="#ED0000FF",alpha=0.6) +
  geom_line(data=quants[quants$t >= 0 & quants_obs$i %in% use_indivs,],aes(x=t,y=mean),col="#ED0000FF") +
  
  geom_ribbon(data=quants_obs[quants_obs$t <= 0 & quants_obs$i %in% use_indivs,], aes(x=t,ymin=lower,ymax=upper),fill="#0099B4FF",alpha=0.4) +
  geom_ribbon(data=quants[quants$t <= 0 & quants_obs$i %in% use_indivs,],aes(x=t,ymin=lower,ymax=upper),fill="#0099B4FF",alpha=0.6) +
  geom_line(data=quants[quants$t <= 0 & quants_obs$i %in% use_indivs,], aes(x=t,y=mean),col="#0099B4FF") +
  
  geom_point(data=viral_loads_dat[viral_loads_dat$i %in% use_indivs, ] %>% filter(obs > 2),aes(x=t,y=obs,shape="Positive, quantified"),size=0.5) + 
  geom_segment(data=viral_loads_dat_bars,aes(x=t,xend=t,y=0,yend=2,col="Positive, unquantified")) +
  scale_color_manual(values=c("Positive, unquantified"="black")) +
  scale_shape_manual(values=c("Positive, quantified"=19)) +
  
  coord_cartesian(xlim=c(-14,35),ylim=c(0,12)) +
  scale_y_continuous(breaks=seq(0,12,by=2), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,35+14,by=7) - 14) +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=2,linetype="dashed") +
  #geom_hline(yintercept=0,linetype="dashed") +
  theme_pubr() +
  theme(panel.grid.major = element_line(color="grey80",size=0.25),
        axis.text=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.text = element_text(family="sans",size=8),
        legend.spacing.y = unit(0.01, "cm"),
        legend.margin = margin(0,0,0,0,unit="pt"),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.box="vertical",
        legend.title=element_blank(),
        strip.text=element_blank(),
        panel.spacing=unit(0.5,"lines"),
        plot.margin = margin(0,0,-2,0,"mm")) +
  ylab("LOG10 RNA copies per swab") +
  xlab("Days post symptom onset") +
  facet_wrap(~i, ncol=1)+
  labs(tag="B")

## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 12500000

## Sample size
n <- 1000

## Duration of epidemic so far in days
times <- 0:365

## Viral kinetics pars
## Change wd to local path
chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains_swab/", parTab, FALSE, 10, burnin=1000000,multi=TRUE)
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
                                             add_noise=FALSE)

viral_loads <- simulated_data$viral_loads %>% as_tibble
test_times <- seq(1,350,by=1)
wow <- sapply(test_times, function(x){
  tmp <- x - infection_times
  tmp <- tmp[tmp > 0]
  tmp[!is.na(tmp)]
})
res <- lapply(wow, function(x) quantile(x,c(0.025,0.5,0.975)))
res <- do.call("rbind",res)
colnames(res) <- c("lower","median","upper")
res <- as.data.frame(res)
res$time <- test_times
pAb <- ggplot(res) + geom_ribbon(aes(ymin=lower,ymax=upper,x=time),alpha=0.5) + 
  geom_line(aes(x=time,y=median)) + theme_bw() +
  ylab("Time since infection") +
  xlab("Date")

viral_loads <- viral_loads %>% mutate(indiv=1:n())

n_samps <- 500
viral_loads_tmp <- viral_loads[sample(1:nrow(viral_loads), n_samps),]

viral_loads_melted <- reshape2::melt(viral_loads_tmp,id.vars="indiv")
colnames(viral_loads_melted) <- c("i","t","Viral load")
viral_loads_melted$t <- as.numeric(viral_loads_melted$t)
viral_loads_melted <- as_tibble(viral_loads_melted) %>% left_join(infection_times_dat)
viral_loads_melted <- viral_loads_melted %>% arrange(i, t)
#viral_loads_melted <- viral_loads_melted %>% filter(inf_time >= 0) %>% arrange(-inf_time,i, t)
#viral_loads_melted$i <- match(viral_loads_melted$i, unique(viral_loads_melted$i))

viral_loads_melted1 <- viral_loads_melted
viral_loads_melted1[viral_loads_melted$`Viral load` <= 0, "Viral load"] <- NA
viral_loads_melted1[viral_loads_melted$`Viral load` > 11, "Viral load"] <- 11

pB <- ggplot(viral_loads_melted1) + 
  geom_tile(aes(x=t,y=as.factor(i),fill=`Viral load`)) + 
  geom_rect(aes(xmin=101,xmax=181,ymin=224,ymax=251),fill=NA,col="grey10") +
  #scale_fill_viridis_c(limits=c(0,10),breaks=seq(0,10,by=2)) +
  scale_fill_gradientn(colours = c("white","#00468bff","#00468bff","#ffdd55ff","#ff751a","#ed0000ff","#ed0000ff"),
                       breaks=seq(0,11,by=2),limits=c(0,11.5),na.value="white")+
  scale_x_continuous(limits=c(50,250), expand=c(0,0), breaks=seq(0,365,by=50),labels=seq(0,365,by=50)) +
  theme_pubr() +
  theme(axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.position=c(0.9,0.8),
        legend.text = element_text(family="sans",size=6),
        legend.title=element_text(family="sans",size=8),
        plot.margin = margin(0,0,0,-10,"mm"),
        plot.background=element_blank()) +
  ylab("") + xlab("Time (days)")+
  labs(tag="E")

indivs <- unique(viral_loads_melted1$i)
i_subsamp <- indivs[225:250]
pB1 <- ggplot(viral_loads_melted1[viral_loads_melted1$i %in% i_subsamp,]) + 
  geom_raster(aes(x=t,y=as.factor(i),fill=`Viral load`)) + 
  #scale_fill_viridis_c(limits=c(0,10),breaks=seq(0,10,by=2)) +
  #scale_fill_gradient2(low="white",mid="#4dbbd5ff",high="#E64B35FF",midpoint=5,limits=c(0,10)) +
  #scale_fill_viridis_c(limits=c(0,10),breaks=seq(0,10,by=2)) +
  scale_fill_gradientn(colours = c("white","#00468bff","#00468bff","#ffdd55ff","#ff751a","#ed0000ff","#ed0000ff"),
                       breaks=seq(0,11,by=2),limits=c(0,11.5),na.value="white")+
  scale_x_continuous(limits=c(100,180), expand=c(0,0), breaks=seq(0,365,by=25),labels=seq(0,365,by=25)) +
  theme_pubr() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.text = element_text(family="sans",size=6),
        legend.title=element_text(family="sans",size=8),
        legend.position="none",
        plot.margin = margin(0,-10,0,0,"mm"),
        plot.background = element_blank()) +
  ylab("Individual") + xlab("Time (days)")+
  labs(tag="C")

viral_prev <- apply(viral_loads[,1:(ncol(viral_loads)-1)], 2, function(x) length(x[x > 2])/length(x))

pA_dat1 <- data.frame(y=epidemic_process$incidence,x=unique(epidemic_process$seir_outputs$time),Variable="Incidence (new infections)")
pA_dat2 <- data.frame(y=viral_prev,x=unique(epidemic_process$seir_outputs$time),Variable="Prevalence (viral load > 2)")
pA_dat <- bind_rows(pA_dat1, pA_dat2)

pA <- ggplot(pA_dat) + geom_line(aes(x=x,y=y,col=Variable)) + 
  theme_pubr() +
  scale_x_continuous(expand=c(0,0), limits=c(50,250),breaks=seq(0,300,by=50)) +
  scale_color_manual(values=c("#ed0000ff","#00466bff")) +
  scale_y_continuous(limits=c(0,0.3),breaks=seq(0,0.3,by=0.05)) +
  theme(legend.position=c(0.8,0.7),
        panel.grid.major = element_line(color="grey80",size=0.25),
        axis.text=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.text = element_text(family="sans",size=6),
        legend.title=element_text(family="sans",size=6),
        plot.margin = margin(0,5,0,0,"mm")) +
  ylab("Per capita infected") +
  xlab("")+
  labs(tag="D")

viral_loads <- simulated_data$viral_loads %>% as_tibble
viral_loads <- viral_loads %>% mutate(indiv=1:n())
viral_loads_melted2 <- reshape2::melt(viral_loads,id.vars="indiv")
colnames(viral_loads_melted2) <- c("i","t","Viral load")
viral_loads_melted2$t <- as.numeric(viral_loads_melted2$t)
viral_loads_melted2 <- as_tibble(viral_loads_melted2)
viral_loads_for_hist <- viral_loads_melted2 %>% filter(t == 140 & `Viral load` > 0) %>% sample_n(10000)
p_hist <- viral_loads_for_hist %>% ggplot() + geom_histogram(aes(x=`Viral load`),fill="grey40",col="black",binwidth=0.25) +
  theme_pubr() + 
  scale_y_continuous(limits=c(0,400),expand = c(0,0)) + 
  scale_x_continuous(breaks=seq(0,12,by=2),limits=c(-0.5,12)) +
  xlab("log10 RNA copies per swab") +
  ylab("Count") +
  theme(axis.text=element_text(family="sans",size=8),
      axis.title=element_text(family="sans",size=8),
      legend.text = element_text(family="sans",size=6),
      legend.title=element_text(family="sans",size=8),
      legend.position="none",
      plot.margin = margin(0,0,0,0,"mm"),
      plot.background = element_blank()) +
  labs(tag="C")
p_hist
## Go along the viral loads matrix and plot the distribution
tmp_viral_loads <- viral_loads[,seq(51,251,by=25)]
colnames(tmp_viral_loads) <- seq(50,250,by=25)

melted_viral_loads_tmp <- reshape2::melt(tmp_viral_loads)
colnames(melted_viral_loads_tmp) <- c("time","value")
melted_viral_loads_tmp$time <- as.numeric(as.character(melted_viral_loads_tmp$time))

subset_melted <- melted_viral_loads_tmp %>% filter(value > 0) %>% group_by(time)# %>% sample_n(100)

tmp_p <- ggplot(melted_viral_loads_tmp[melted_viral_loads_tmp$value > 0,])+ 
  geom_violin(aes(y=value,x=time,group=time),draw_quantiles = c(0.025,0.5,0.975),trim=TRUE,fill="grey90") +
  geom_jitter(data=subset_melted,aes(x=time,y=value),width=3,height=0,size=0.1) +
  #scale_x_continuous(limits=c(50,250), expand=c(0,0)) +
  theme_pubr()+
  xlab("Time (days)") +
  ylab("LOG10 viral load")+ 
  theme(legend.position=c(0.9,0.8),
        panel.grid.major = element_line(color="grey80",size=0.25),
        axis.text=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.text = element_text(family="sans",size=6),
        legend.title=element_text(family="sans",size=8),
        plot.margin = margin(0,-5,0,0,"mm"),
        panel.background = element_blank(),
        plot.background = element_blank())+
  labs(tag="E")

right_p <- pA/pB + plot_layout(heights=c(1,2.5))
left_p_bot <- ((pC+theme(legend.position=c(0.8,0.8))) + plot_spacer() + plot_layout(ncol=1,heights=c(20,1))) | pB1
left_p_bot
left_p <- (plot_spacer() / left_p_bot) + plot_layout(heights=c(1,1.5))

#png("Fig2.png",height=6,width=8,res=300,units="in")
#left_p | right_p + plot_layout(widths=c(1.3,1))
#dev.off()


#pdf("Fig2_base_revision2.pdf",height=6,width=8)
#left_p | right_p + plot_layout(widths=c(1.3,1))
#dev.off()

chains1 <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains_swab/", parTab, TRUE, 10, burnin=1000000,multi=TRUE,chainNo=FALSE)

gelman.diag(chains1[[1]])
gelman_diagnostics(chains1[[1]])

chains1 <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains_swab/", parTab, TRUE, 100, burnin=1000000,multi=TRUE,chainNo=TRUE)

chains_tmp <- lapply(chains1[[1]], function(x){
  tmp <- as.data.frame(x)
  tmp$sampno <- 1:nrow(tmp)
  tmp
})


chains_comb <- do.call("rbind",chains_tmp)
chains_comb <- reshape2::melt(chains_comb, id.vars=c("sampno","chain"))
chains_comb$variable <- as.character(chains_comb$variable)
use_pars <- c("sd","viral_peak_mean","viral_peak_sd",
              "wane_mean","wane_sd","tshift","desired_mode","incu")
par_key <- c("sd"="sigma", 
             "viral_peak_mean"="bar(alpha)","viral_peak_sd"="sigma^alpha",
             "wane_mean"="bar(omega)","wane_sd"="sigma^omega",
             "desired_mode"="t[p]","tshift"="t[g]","incu"="t[inc]")

chains_comb <- chains_comb %>% filter(variable %in% use_pars)
chains_comb$variable <- par_key[chains_comb$variable]

p1 <- ggplot(chains_comb) + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~variable,scales="free",labeller=label_parsed,nrow=2) +
  scale_color_lancet() + 
  xlab("Iteration") +
  ylab("Value") +
  theme_pubr()+
  theme(legend.position="none",
        text=element_text(family="sans",size=8)) +
  labs(tag="A")

p2 <- ggplot(chains_comb) + 
  geom_density(aes(x=value,fill=as.factor(chain)),alpha=0.5) + 
  facet_wrap(~variable,scales="free",labeller=label_parsed,ncol=1) +
  scale_color_lancet() +
  xlab("Value") +
  ylab("Density") +
  theme_pubr() +
  theme(legend.position="none",
        text=element_text(family="sans",size=8),axis.text=element_text(size=6)) +
  labs(tag="B")

#pdf("FigS2.pdf",height=7,width=5)
# p1 | p2
#dev.off()

## Same think but for sputum
chains2 <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains_sputum/", parTab, TRUE, 100, burnin=1000000,multi=TRUE,chainNo=TRUE)

chains_tmp2 <- lapply(chains2[[1]], function(x){
  tmp <- as.data.frame(x)
  tmp$sampno <- 1:nrow(tmp)
  tmp
})


chains_comb2 <- do.call("rbind",chains_tmp2)
chains_comb2 <- reshape2::melt(chains_comb2, id.vars=c("sampno","chain"))
chains_comb2$variable <- as.character(chains_comb2$variable)
use_pars <- c("sd","viral_peak_mean","viral_peak_sd",
              "wane_mean","wane_sd","tshift","desired_mode","incu")
par_key <- c("sd"="sigma", 
             "viral_peak_mean"="bar(alpha)","viral_peak_sd"="sigma^alpha",
             "wane_mean"="bar(omega)","wane_sd"="sigma^omega",
             "desired_mode"="t[p]","tshift"="t[g]","incu"="t[inc]")

chains_comb2 <- chains_comb2 %>% filter(variable %in% use_pars)
chains_comb2$variable <- par_key[chains_comb2$variable]

p3 <- ggplot(chains_comb2) + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~variable,scales="free",labeller=label_parsed,nrow=2) +
  scale_color_lancet() + 
  xlab("Iteration") +
  ylab("Value") +
  theme_pubr()+
  theme(legend.position="none",
        text=element_text(family="sans",size=8),
        axis.text=element_text(size=6)) +
  labs(tag="B")

pdf("FigS2.pdf",height=5,width=7)
p1 / p3
dev.off()

chain1 <- chains$chain

format_estimates <- function(x){
  tmp <- signif(quantile(x, c(0.5,0.025,0.975)),3)
  res <- paste0(tmp[1], " (", tmp[2],"-",tmp[3],")")
  res
}
format_estimates(chain2[,"viral_peak_mean"])

## Figure S1
mu_labels <- apply(chain[,colnames(chain) %like% "viral_peak"], 2, function(x){
  res <- signif(quantile(x,c(0.5,0.025,0.975)),3)
  final <- paste0("α = ",res[1], " (",res[2],"-",res[3],")\n")
  final
})
mu_labels <- mu_labels[3:11]

wane_labels <- apply(chain[,colnames(chain) %like% "wane"], 2, function(x){
  res <- signif(quantile(x,c(0.5,0.025,0.975)),3)
  final <- paste0("tw = ",res[1], " (",res[2],"-",res[3],")")
  final
})
wane_labels <- wane_labels[4:12]

wow <- paste0(mu_labels, wane_labels)
label_dat <- data.frame(label=wow, i=1:length(wow))


use_indivs <- 1:9


## Revision: plot bar for points when below limit of quantification
viral_loads_dat_bars_supp <- viral_loads_dat %>% filter(obs < 2)

pS1 <- ggplot() + 
  geom_rect(data=rect_dat,aes(ymin=y1,ymax=y2,xmin=x1,xmax=x2),fill="grey70",alpha=0.5) +
  geom_rect(data=rect_dat2,aes(ymin=y1,ymax=y2,xmin=x1,xmax=x2),fill="grey70",alpha=0.5) +
  geom_ribbon(data=quants_obs[quants_obs$t >= 0 & quants_obs$i %in% use_indivs,], aes(x=t,ymin=lower,ymax=upper),fill="#ED0000FF",alpha=0.4) +
  geom_ribbon(data=quants[quants$t >= 0 & quants_obs$i %in% use_indivs,], aes(x=t ,ymin=lower,ymax=upper),fill="#ED0000FF",alpha=0.6) +
  geom_line(data=quants[quants$t >= 0 & quants_obs$i %in% use_indivs,],aes(x=t,y=mean),col="#ED0000FF") +
  
  geom_ribbon(data=quants_obs[quants_obs$t <= 0 & quants_obs$i %in% use_indivs,], aes(x=t,ymin=lower,ymax=upper),fill="#0099B4FF",alpha=0.4) +
  geom_ribbon(data=quants[quants$t <= 0 & quants_obs$i %in% use_indivs,],aes(x=t,ymin=lower,ymax=upper),fill="#0099B4FF",alpha=0.6) +
  geom_line(data=quants[quants$t <= 0 & quants_obs$i %in% use_indivs,], aes(x=t,y=mean),col="#0099B4FF") +
  
  geom_point(data=viral_loads_dat[viral_loads_dat$i %in% use_indivs, ] %>% filter(obs > 2),aes(x=t,y=obs,shape="Positive, quantified"),size=0.5) + 
  geom_segment(data=viral_loads_dat_bars_supp,aes(x=t,xend=t,y=0,yend=2,col="Positive, unquantified")) +
  scale_color_manual(values=c("Positive, unquantified"="black")) +
  scale_shape_manual(values=c("Positive, quantified"=19)) +
  
  geom_text(data=label_dat,aes(x=21,y=10,label=label),parse=FALSE,size=3) +
  coord_cartesian(xlim=c(-14,35),ylim=c(0,12)) +
  scale_y_continuous(breaks=seq(0,12,by=1), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,35+14,by=7) - 14) +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=2,linetype="dashed") +
  #geom_hline(yintercept=0,linetype="dashed") +
  
  theme_pubr() +
  theme(panel.grid.major = element_line(color="grey80",size=0.25),
        axis.text=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.text = element_text(family="sans",size=8),
        legend.title=element_blank(),
        legend.position="bottom",
        strip.background = element_blank(),
        strip.text=element_blank(),
        plot.margin = margin(0,0,0,0,"mm")) +
  ylab("LOG10 RNA copies per swab") +
  xlab("Days post symptom onset") +
  facet_wrap(~i, ncol=3)+
  labs(tag="")

cairo_pdf("FigS_fits.pdf",height=5,width=8)
pS1
dev.off()
