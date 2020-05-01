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

setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/")
source("simulation_functions.R")
source("model_funcs.R")

cols <- ggsci::pal_npg()(10)
n_indiv <- 9

set.seed(1)

parTab <- read.csv("~/Google Drive/nCoV/pool_samples/partab.csv",stringsAsFactors=FALSE)
## Read in extracted Wolfel data

viral_loads_dat <- read.csv("~/Google Drive/nCoV/pool_samples/drosten_data.csv")
ggplot(viral_loads_dat) + geom_line(aes(x=t,y=obs)) + facet_wrap(~i)
viral_loads_dat <- viral_loads_dat[,c("t","obs","i")]

## Run MCMC
mcmcPars1 <- c("iterations"=200000,"popt"=0.44,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=100000,"save_block"=1000)
mcmcPars2 <- c("iterations"=200000,"popt"=0.234,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=100000,"save_block"=1000)

parTab[parTab$names == "viral_peak_sd", "values"] <- 0.9
parTab[parTab$names == "viral_peak_sd", "fixed"] <- 1

chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains", parTab, FALSE, 1, burnin=200000,multi=TRUE)
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)

## Simulate posterior draws over 100 days
## Note that t0 is going to be time of infection
times <- seq(0,100,by=0.1)
dat_fake <- data.frame(t=rep(times, n_indiv),obs=0,i=rep(seq_len(n_indiv),each=length(times)))

## 10000 posterior draws
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
  #pred[which(pred < 0)] <- 0
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

rect_dat <- data.frame(x1=-14,x2=21,y1=-1,y2=11)

use_indivs <- c(2,3,6,7)

pC <- ggplot() + 
  geom_rect(data=rect_dat,aes(ymin=y1,ymax=y2,xmin=x1,xmax=x2),fill="grey70",alpha=0.5) +
  geom_ribbon(data=quants_obs[quants_obs$t >= 21 & quants_obs$i %in% use_indivs,], aes(x=t,ymin=lower,ymax=upper),fill="#DC0000FF",alpha=0.4) +
  geom_ribbon(data=quants[quants$t >= 21 & quants$i %in% use_indivs,],aes(x=t,ymin=lower,ymax=upper),fill="#DC0000FF",alpha=0.6) +
  geom_line(data=quants[quants$t >= 21 & quants$i %in% use_indivs,], aes(x=t,y=median),col="#DC0000FF") +
  #geom_ribbon(data=quants_obs[quants_obs$t <= 21 & quants_obs$i %in% use_indivs,], aes(x=t,ymin=lower,ymax=upper),fill="#00a087ff",alpha=0.4) +
  geom_ribbon(data=quants[quants$t <= 21 & quants$i %in% use_indivs,],aes(x=t,ymin=lower,ymax=upper),fill="#ffaa00",alpha=0.6) +
  geom_line(data=quants[quants$t <= 21 & quants$i %in% use_indivs,], aes(x=t,y=median),col="#de7e00") +
  geom_point(data=viral_loads_dat[viral_loads_dat$i %in% use_indivs, ],aes(x=t+21,y=obs),size=0.5) + 
  coord_cartesian(ylim=c(0,10),xlim=c(0,50)) +
  scale_y_continuous(breaks=seq(0,10,by=2)) +
  scale_x_continuous(breaks=seq(0,7*9,by=7),labels=seq(-21,7*6,by=7)) +
  geom_vline(xintercept=21,linetype="dashed") +
  geom_hline(yintercept=2,linetype="dashed") +
  theme_pubr() +
  theme(panel.grid.major = element_line(color="grey80",size=0.25),
        axis.text=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.text = element_text(family="sans",size=8),
        strip.background = element_blank(),
        strip.text=element_blank(),
        plot.margin = margin(0,0,0,0,"mm")) +
  ylab("LOG10 RNA copies/ml, g, swab") +
  xlab("Days post symptom onset") +
  facet_wrap(~i, ncol=1)+
  labs(tag="B")



pC
## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 6900000

## Sample size
n <- 500000

## Epidemic growth rate
growth_rate <- 0.1

## Duration of epidemic so far in days
times <- 0:365

## Viral kinetics pars
## Change wd to local path
parTab <- read.csv("pars/partab.csv",stringsAsFactors=FALSE)
chains <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains", parTab, FALSE, 1, burnin=200000,multi=TRUE)
chain <- as.data.frame(chains[[2]])

pars <- list(viral_peak_mean=mean(chain$viral_peak_mean), viral_peak_sd=mean(chain$viral_peak_sd),
             wane_mean=mean(chain$wane_mean), wane_sd=mean(chain$wane_sd),
             tp_last_day=7)

## Simulate the epidemic process
#epidemic_process <- simulate_epidemic_process(population_n,0.1,times)
#epidemic_process$plot
seir_pars <- c("R0"=2.2,"gamma"=1/7,"sigma"=1/6.4,"I0"=100)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)
pA <- epidemic_process$incidence_plot

## Simulate infection times
infection_times <- simulate_infection_times(n, epidemic_process$overall_prob_infection, 
                                            epidemic_process$incidence)
#incidence <- rep(1, 120)/120
#infection_times <- simulate_infection_times(n, 0.05, 
#                                            incidence)

## Simulate symptom onset times using default incubation period distribution
onset_times <- simulate_symptom_onsets(infection_times)

## Simulate viral loads for the sample population
#simulated_data <- simulate_viral_loads(infection_times, onset_times, times, pars)

#viral_loads <- simulated_data$viral_loads
viral_loads <- read_csv("seir_viral_loads.csv")

inf_times <- apply(viral_loads, 1, function(x) which(x > 0)[1])

test_times <- seq(1,350,by=1)
wow <- sapply(test_times, function(x){
  tmp <- x - inf_times
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

pA/pAb

viral_loads <- viral_loads %>% mutate(indiv=1:n())


n_samps <- 500
viral_loads_tmp <- viral_loads[sample(1:nrow(viral_loads), n_samps),]

viral_loads_melted <- reshape2::melt(viral_loads_tmp,id.vars="indiv")
colnames(viral_loads_melted) <- c("i","t","Viral load")
viral_loads_melted$t <- as.numeric(viral_loads_melted$t)
viral_loads_melted1 <- viral_loads_melted
viral_loads_melted1[viral_loads_melted$`Viral load` <= 0, "Viral load"] <- NA
viral_loads_melted1 <- viral_loads_melted1 %>% arrange(i, t)

pB <- ggplot(viral_loads_melted1) + 
  geom_tile(aes(x=t,y=as.factor(i),fill=`Viral load`)) + 
  geom_rect(aes(xmin=100,xmax=200,ymin=224,ymax=251),fill=NA,col="grey10") +
  #scale_fill_viridis_c(limits=c(0,10),breaks=seq(0,10,by=2)) +
  scale_fill_gradientn(colours = c("white","#3c5488ff","#3c5488ff","#ffdd55ff","#ff751a","#dc000fff","#dc000fff"),
                       breaks=seq(0,10,by=2),limits=c(0,10.5),na.value="white")+
  scale_x_continuous(limits=c(0,365), expand=c(0,0), breaks=seq(0,365,by=50),labels=seq(0,365,by=50)) +
  theme_pubr() +
  theme(axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.text = element_text(family="sans",size=6),
        legend.title=element_text(family="sans",size=8),
        legend.position=c(0.9,0.8),
        plot.margin = margin(0,0,0,-5,"mm"),
        plot.background=element_blank()) +
  ylab("") + xlab("Time (days)")+
  labs(tag="E")
pB
indivs <- unique(viral_loads_melted1$i)
i_subsamp <- indivs[225:250]
pB1 <- ggplot(viral_loads_melted1[viral_loads_melted1$i %in% i_subsamp,]) + 
  geom_raster(aes(x=t,y=as.factor(i),fill=`Viral load`)) + 
  #scale_fill_viridis_c(limits=c(0,10),breaks=seq(0,10,by=2)) +
  #scale_fill_gradient2(low="white",mid="#4dbbd5ff",high="#E64B35FF",midpoint=5,limits=c(0,10)) +
  #scale_fill_viridis_c(limits=c(0,10),breaks=seq(0,10,by=2)) +
  scale_fill_gradientn(colours = c("white","#3c5488ff","#3c5488ff","#ffdd55ff","#ff751a","#dc000fff","#dc000fff"),
                       breaks=seq(0,10,by=2),limits=c(0,10.5),na.value="white")+
  scale_x_continuous(limits=c(100,200), expand=c(0,0), breaks=seq(0,365,by=25),labels=seq(0,365,by=25)) +
  theme_pubr() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.text = element_text(family="sans",size=6),
        legend.title=element_text(family="sans",size=8),
        legend.position="none",
        plot.margin = margin(0,-5,0,0,"mm"),
        plot.background = element_blank()) +
  ylab("Individual") + xlab("Time (days)")+
  labs(tag="C")
pB1

viral_prev <- apply(viral_loads[,1:(ncol(viral_loads)-1)], 2, function(x) length(x[x > 0])/length(x))

pA_dat1 <- data.frame(y=epidemic_process$incidence,x=unique(epidemic_process$seir_outputs$time),Variable="Incidence (new infections)")
pA_dat2 <- data.frame(y=viral_prev,x=unique(epidemic_process$seir_outputs$time),Variable="Prevalence (viral load > 1)")
pA_dat <- bind_rows(pA_dat1, pA_dat2)

pA <- ggplot(pA_dat) + geom_line(aes(x=x,y=y,col=Variable)) + 
  theme_pubr() +
  scale_x_continuous(expand=c(0,0), limits=c(0,365),breaks=seq(0,365,by=50)) +
  scale_color_manual(values=c("#dc000fff","#3c5488ff")) +
  scale_y_continuous(limits=c(0,0.31),breaks=seq(0,0.3,by=0.05)) +
  theme(legend.position=c(0.8,0.7),
        panel.grid.major = element_line(color="grey80",size=0.25),
        axis.text=element_text(family="sans",size=8),
        axis.title=element_text(family="sans",size=8),
        legend.text = element_text(family="sans",size=6),
        legend.title=element_text(family="sans",size=6),
        plot.margin = margin(-5,0,-5,-5,"mm")) +
  ylab("Per capita infected") +
  xlab("")+
  labs(tag="D")
pA

## Go along the viral loads matrix and plot the distribution
tmp_viral_loads <- viral_loads[,seq(51,251,by=25)]
colnames(tmp_viral_loads) <- seq(50,250,by=25)

melted_viral_loads_tmp <- reshape2::melt(tmp_viral_loads)
colnames(melted_viral_loads_tmp) <- c("time","value")
melted_viral_loads_tmp$time <- as.numeric(as.character(melted_viral_loads_tmp$time))

subset_melted <- melted_viral_loads_tmp %>% filter(value > 0) %>% group_by(time) %>% sample_n(100)

tmp_p <- ggplot(melted_viral_loads_tmp[melted_viral_loads_tmp$value > 0,])+ 
  geom_violin(aes(y=value,x=time,group=time),draw_quantiles = c(0.025,0.5,0.975),trim=TRUE,fill="grey90") +
  geom_jitter(data=subset_melted,aes(x=time,y=value),width=3,height=0,size=0.1) +
  scale_x_continuous(breaks=seq(0,365,by=50),limits=c(0,365),expand=c(0,0)) +
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
pA/tmp_p

right_p <- pA/pB + plot_layout(heights=c(1,2.5))
left_p_bot <- pC | pB1
left_p <- (plot_spacer() / left_p_bot) + plot_layout(heights=c(1,1.5))

png("Fig2.png",height=6,width=8,res=300,units="in")
left_p | right_p + plot_layout(widths=c(1.3,1))
dev.off()


pdf("Fig2_base.pdf",height=6,width=8)
left_p | right_p + plot_layout(widths=c(1.3,1))
dev.off()

parTab[parTab$names == "viral_peak_sd","fixed"] <- 1
chains1 <- load_mcmc_chains("~/Google Drive/nCoV/pool_samples/chains", parTab, TRUE, 10, 
                            burnin=200000,multi=TRUE,chainNo = TRUE)
gelman.diag(chains1[[1]])

gelman_diagnostics(chains1[[1]])

chains_tmp <- lapply(chains1[[1]], function(x){
  tmp <- as.data.frame(x)
  tmp$sampno <- 1:nrow(tmp)
  tmp
})


chains_comb <- do.call("rbind",chains_tmp)
chains_comb <- reshape2::melt(chains_comb, id.vars=c("sampno","chain"))

use_pars <- c("sd","viral_peak_mean","wane_mean","wane_sd","tp_mean","tp_sd")
par_key <- c("sd"="sigma", "viral_peak_mean"="bar(alpha)","wane_mean"="bar(omega)","wane_sd"="sigma^omega","tp_mean"="bar(tp)","tp_sd"="sigma^{tp}")

chains_comb <- chains_comb %>% filter(variable %in% use_pars)
chains_comb$variable <- par_key[chains_comb$variable]


p1 <- ggplot(chains_comb) + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~variable,scales="free_y",labeller=label_parsed,ncol=1) +
  scale_color_viridis_d() + 
  xlab("Iteration") +
  ylab("Value") +
  theme_pubr()+
  theme(legend.position="none",
        text=element_text(family="sans",size=8)) +
  labs(tag="A")

p2 <- ggplot(chains_comb) + 
  geom_density(aes(x=value,fill=as.factor(chain)),alpha=0.5) + 
  facet_wrap(~variable,scales="free",labeller=label_parsed,ncol=1) +
  scale_fill_viridis_d() +
  xlab("Value") +
  ylab("Density") +
  theme_pubr() +
  theme(legend.position="none",
        text=element_text(family="sans",size=8)) +
  labs(tag="B")

pdf("FigS2.pdf",height=7,width=5)
p1 | p2
dev.off()

format_estimates <- function(x){
  tmp <- signif(quantile(x, c(0.5,0.025,0.975)),3)
  res <- paste0(tmp[1], " (", tmp[2],"-",tmp[3],")")
  res
}
format_estimates(wow$viral_peak_mean)

