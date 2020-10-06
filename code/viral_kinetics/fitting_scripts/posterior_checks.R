library(coda)
library(bayesplot)
library(tidyverse)
library(ggplot2)
library(lazymcmc)
library(ggpubr)
library(patchwork)
library(rethinking)
library(data.table)
library(ggsci)

parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)

###################################
## SWAB
###################################
chains_swab <- load_mcmc_chains(paste0("~/Google Drive/nCoV/pool_samples/chains_swab/"), 
                           parTab, TRUE, 100, burnin=1000000,multi=TRUE,chainNo = FALSE)
gelman.diag(chains_swab[[1]])
effectiveSize(chains_swab[[1]])
chains_swab <- load_mcmc_chains(paste0("~/Google Drive/nCoV/pool_samples/chains_swab/"), 
                                parTab, TRUE, 100, burnin=1000000,multi=TRUE,chainNo = TRUE)
chains_swab_list <- chains_swab[[1]]
#mcmc_trace(chains_swab_list)
#chains_swab_list <- lapply(chains_swab_list, function(x) x[,c("chain","sd", "viral_peak_mean", "viral_peak_sd", "wane_mean", "wane_sd", 
#                                          "rho_viral_wane", "viral_peak", "t_wane", "desired_mode", "tshift","incu")])

chains_tmp <- lapply(chains_swab_list, function(x){
  tmp <- as.data.frame(x)
  tmp$sampno <- 1:nrow(tmp)
  tmp
})

chains_comb <- do.call("rbind",chains_tmp)

viral_peaks_samp <- numeric(nrow(chains_comb))
wane_samp <- numeric(nrow(chains_comb))
incu_samp <- numeric(nrow(chains_comb))
tshift_samp <- numeric(nrow(chains_comb))
desired_mode_samp <- numeric(nrow(chains_comb))

incu_cols <- which(colnames(chains_comb) %like% "incu")
tshift_cols <- which(colnames(chains_comb) %like% "tshift")
desired_mode_cols <- which(colnames(chains_comb) %like% "desired_mode")

for(i in 1:nrow(chains_comb)){
  choose_indiv <- sample(seq_len(9),1)
  viral_peak_mean <- chains_comb$viral_peak_mean[i]
  wane_mean <- chains_comb$wane_mean[i]
  viral_peak_sd <- chains_comb$viral_peak_sd[i]
  wane_sd <- chains_comb$wane_sd[i]
  rho_viral_wane <- chains_comb$rho_viral_wane[i]
  tmp <- rmvnorm2(1, c(viral_peak_mean,wane_mean),c(viral_peak_sd,wane_sd),
                  matrix(c(1,rho_viral_wane,rho_viral_wane,1),ncol=2))
  viral_peaks_samp[i] <- tmp[1,1]
  wane_samp[i] <- tmp[1,2]
  
  while(tmp[1,2] < 0){
    tmp <- rmvnorm2(1, c(viral_peak_mean,wane_mean),c(viral_peak_sd,wane_sd),
                    matrix(c(1,rho_viral_wane,rho_viral_wane,1),ncol=2))
    viral_peaks_samp[i] <- tmp[1,1]
    wane_samp[i] <- tmp[1,2]
  }
  
  incu_samp[i] <- chains_comb[i,incu_cols[choose_indiv]]
  tshift_samp[i] <- chains_comb[i,tshift_cols[choose_indiv]]
  desired_mode_samp[i] <- chains_comb[i,desired_mode_cols[choose_indiv]]
}
chains_comb <- cbind(chains_comb, viral_peaks_samp, wane_samp, incu_samp, tshift_samp, desired_mode_samp)

chains_comb <- reshape2::melt(chains_comb, id.vars=c("sampno","chain"))
chains_comb$variable <- as.character(chains_comb$variable)

use_pars <- c("sd","viral_peak_mean","viral_peak_sd","wane_mean","wane_sd","rho_viral_wane",
              "viral_peaks_samp", "wane_samp","desired_mode_samp","tshift_samp","incu_samp")

par_key <- c("sd"="sigma", "viral_peak_mean"="bar(alpha)","viral_peak_sd"="sigma^alpha","wane_mean"="bar(t[w])","wane_sd"="sigma^{t[w]}",
             "rho_viral_wane"="rho","incu_samp"="t[inc]",
             "viral_peaks_samp"="alpha", "wane_samp"="t[w]","desired_mode_samp"="t[p]","tshift_samp"="t[g]")

par_key <- c("sd"="sigma", "viral_peak_mean"="bar(alpha)","viral_peak_sd"="sigma^alpha","wane_mean"="bar(t[w])","wane_sd"="sigma^{t[w]}",
             "rho_viral_wane"="rho","incu_samp"="t[inc]",
             "viral_peaks_samp"="alpha", "wane_samp"="t[w]","desired_mode_samp"="t[p]","tshift_samp"="t[g]")

par_levels <- c("bar(alpha)","sigma^alpha","alpha",
                "bar(t[w])", "sigma^{t[w]}", "t[w]",
                "t[inc]","t[p]","t[g]",
                "sigma", "rho")

chains_comb <- chains_comb %>% filter(variable %in% use_pars)
chains_comb$variable <- par_key[chains_comb$variable]
chains_comb$Data <- "Swab"
chains_comb$variable <- factor(chains_comb$variable, levels=par_levels)

p1 <- ggplot(chains_comb) + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~variable,scales="free_y",labeller=label_parsed,ncol=3) +
  scale_color_lancet() + 
  xlab("Iteration") +
  ylab("Value") +
  theme_pubr()+
  theme(legend.position="none",
        text=element_text(family="sans",size=8),
        axis.text=element_text(size=6),
        strip.background=element_blank(),
        strip.text=element_text(size=8,face="bold"),
        plot.tag = element_text(face="bold"))+
  labs(tag="A")

p2 <- ggplot(chains_comb) + 
  geom_density(aes(x=value,fill=as.factor(chain)),alpha=0.5) + 
  facet_wrap(~variable,scales="free",labeller=label_parsed,ncol=3) +
  scale_fill_lancet() +
  xlab("Value") +
  ylab("Density") +
  theme_pubr() +
  theme(legend.position="none",
        text=element_text(family="sans",size=8)) +
  labs(tag="B")


###################################
## SPUTUM
###################################
chains_sputum <- load_mcmc_chains(paste0("~/Google Drive/nCoV/pool_samples/chains_sputum/"), 
                                parTab, TRUE, 100, burnin=1000000,multi=TRUE,chainNo = FALSE)
gelman.diag(chains_sputum[[1]])
effectiveSize(chains_sputum[[1]])
chains_sputum <- load_mcmc_chains(paste0("~/Google Drive/nCoV/pool_samples/chains_sputum/"), 
                                parTab, TRUE, 100, burnin=1000000,multi=TRUE,chainNo = TRUE)
chains_sputum_list <- chains_sputum[[1]]
#mcmc_trace(chains_sputum_list)
#chains_sputum_list <- lapply(chains_sputum_list, function(x) x[,c("chain","sd", "viral_peak_mean", "viral_peak_sd", "wane_mean", "wane_sd", 
#                                                                  "rho_viral_wane", "viral_peak", "t_wane", "desired_mode", "tshift","incu")])

chains_tmp <- lapply(chains_sputum_list, function(x){
  tmp <- as.data.frame(x)
  tmp$sampno <- 1:nrow(tmp)
  tmp
})

chains_comb_sputum <- do.call("rbind",chains_tmp)

viral_peaks_samp <- numeric(nrow(chains_comb_sputum))
wane_samp <- numeric(nrow(chains_comb_sputum))
incu_samp <- numeric(nrow(chains_comb_sputum))
tshift_samp <- numeric(nrow(chains_comb_sputum))
desired_mode_samp <- numeric(nrow(chains_comb_sputum))

incu_cols <- which(colnames(chains_comb_sputum) %like% "incu")
tshift_cols <- which(colnames(chains_comb_sputum) %like% "tshift")
desired_mode_cols <- which(colnames(chains_comb_sputum) %like% "desired_mode")

for(i in 1:nrow(chains_comb_sputum)){
  choose_indiv <- sample(seq_len(9),1)
  viral_peak_mean <- chains_comb_sputum$viral_peak_mean[i]
  wane_mean <- chains_comb_sputum$wane_mean[i]
  viral_peak_sd <- chains_comb_sputum$viral_peak_sd[i]
  wane_sd <- chains_comb_sputum$wane_sd[i]
  rho_viral_wane <- chains_comb_sputum$rho_viral_wane[i]
  tmp <- rmvnorm2(1, c(viral_peak_mean,wane_mean),c(viral_peak_sd,wane_sd),
                  matrix(c(1,rho_viral_wane,rho_viral_wane,1),ncol=2))
  viral_peaks_samp[i] <- tmp[1,1]
  wane_samp[i] <- tmp[1,2]
  
  while(tmp[1,2] < 0){
    tmp <- rmvnorm2(1, c(viral_peak_mean,wane_mean),c(viral_peak_sd,wane_sd),
                    matrix(c(1,rho_viral_wane,rho_viral_wane,1),ncol=2))
    viral_peaks_samp[i] <- tmp[1,1]
    wane_samp[i] <- tmp[1,2]
  }
  
  incu_samp[i] <- chains_comb_sputum[i,incu_cols[choose_indiv]]
  tshift_samp[i] <- chains_comb_sputum[i,tshift_cols[choose_indiv]]
  desired_mode_samp[i] <- chains_comb_sputum[i,desired_mode_cols[choose_indiv]]
}
chains_comb_sputum <- cbind(chains_comb_sputum, viral_peaks_samp, wane_samp, incu_samp, tshift_samp, desired_mode_samp)

chains_comb_sputum <- reshape2::melt(chains_comb_sputum, id.vars=c("sampno","chain"))
chains_comb_sputum$variable <- as.character(chains_comb_sputum$variable)


use_pars <- c("sd","viral_peak_mean","viral_peak_sd","wane_mean","wane_sd","rho_viral_wane",
              "viral_peaks_samp", "wane_samp","desired_mode_samp","tshift_samp","incu_samp")

par_key <- c("sd"="sigma", "viral_peak_mean"="bar(alpha)","viral_peak_sd"="sigma^alpha","wane_mean"="bar(t[w])","wane_sd"="sigma^{t[w]}",
             "rho_viral_wane"="rho","incu_samp"="t[inc]",
             "viral_peaks_samp"="alpha", "wane_samp"="t[w]","desired_mode_samp"="t[p]","tshift_samp"="t[g]")

par_key <- c("sd"="sigma", "viral_peak_mean"="bar(alpha)","viral_peak_sd"="sigma^alpha","wane_mean"="bar(t[w])","wane_sd"="sigma^{t[w]}",
             "rho_viral_wane"="rho","incu_samp"="t[inc]",
             "viral_peaks_samp"="alpha", "wane_samp"="t[w]","desired_mode_samp"="t[p]","tshift_samp"="t[g]")

par_levels <- c("bar(alpha)","sigma^alpha","alpha",
                "bar(t[w])", "sigma^{t[w]}", "t[w]",
                "t[inc]","t[p]","t[g]",
                "sigma", "rho")


chains_comb_sputum <- chains_comb_sputum %>% filter(variable %in% use_pars)
chains_comb_sputum$variable <- par_key[chains_comb_sputum$variable]
chains_comb_sputum$Data <- "Sputum"
chains_comb_sputum$variable <- factor(chains_comb_sputum$variable, levels=par_levels)

p3 <- ggplot(chains_comb_sputum) + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~variable,scales="free_y",labeller=label_parsed,ncol=3) +
  scale_color_lancet() + 
  xlab("Iteration") +
  ylab("Value") +
  theme_pubr()+
  theme(legend.position="none",
        text=element_text(family="sans",size=8),
        axis.text=element_text(size=6),
        strip.background=element_blank(),
        strip.text=element_text(size=8,face="bold"),
        plot.tag = element_text(face="bold"))+
  labs(tag="A")

p4 <- ggplot(chains_comb_sputum) + 
  geom_density(aes(x=value,fill=as.factor(chain)),alpha=0.5) + 
  facet_wrap(~variable,scales="free",labeller=label_parsed,ncol=3) +
  scale_fill_lancet() +
  xlab("Value") +
  ylab("Density") +
  theme_pubr() +
  theme(legend.position="none",
        text=element_text(family="sans",size=8)) +
  labs(tag="B")




###################################
## COMBINED
###################################
chains_comb_all <- bind_rows(chains_comb, chains_comb_sputum)

p5 <- ggplot(chains_comb_all) + 
  geom_density(aes(x=value,fill=Data),alpha=0.5) + 
  facet_wrap(~variable,scales="free",labeller=label_parsed,ncol=3) +
  scale_fill_lancet() +
  xlab("Value") +
  ylab("Density") +
  theme_pubr() +
  theme(text=element_text(family="sans",size=8),
        axis.text=element_text(size=6),
        strip.background=element_blank(),
        strip.text=element_text(size=8,face="bold"),
        plot.tag = element_text(face="bold"),
        legend.position=c(0.8,0.1))

pdf("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/figs/convergence_swab.pdf",height=8,width=8)
p1 / p2
dev.off()
pdf("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/figs/convergence_sputum.pdf",height=8,width=8)
p3 / p4
dev.off()

pdf("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/figs/convergence_both.pdf",height=8,width=8)
p1 / (p3 + labs(tag="B"))
dev.off()

pdf("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/figs/densities_all.pdf",height=8,width=8)
p5
dev.off()
