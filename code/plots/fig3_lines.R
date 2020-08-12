library(tidyverse)
library(patchwork)
library(ggsci)

setwd("~/Documents/GitHub/covid19-group-tests/code/plots/")

dat <- read_csv("summary.all_results.csv")
old_colnames <- colnames(dat)
colnames(dat) <- c("N","b","q","t","p_true","p_sample","p_mle")

dat$MSE <- (dat$p_mle - dat$p_true)^2
dat$rBias <- 100*abs(dat$p_mle - dat$p_true)/dat$p_true
dat$true_diff <- dat$p_true - dat$p_mle
dat$true_diff_pct <- abs(dat$true_diff)/dat$p_true
dat$label <- paste0(dat$N,"_",dat$b,"_",dat$q)
dat$prev_group <- cut(dat$p_true*100, c(0,0.01,0.05,0.1,0.5,1,5,10,100))
dat$sample_prev_group <- cut(dat$p_sample*100, c(0,0.01,0.05,0.1,0.5,1,5,10,100))
dat <- dat %>% filter(t < 137 & t > 20)

labels <- c("<0.01%", "(0.01%, 0.05%]", "(0.05%, 0.1%]", "(0.1%, 0.5%]", "(0.5%, 1%]", 
            "(1%, 5%]", "(5%, 10%]", ">10%")

p1 <- dat %>% group_by(N) %>% filter(b == min(b)) %>% 
  group_by(N,prev_group) %>% summarize(mean_rbias = mean(rBias),
                                              sd_rbias = sd(rBias)) %>%
  mutate(upper = mean_rbias + sd_rbias)  %>%
  ggplot() + 
  geom_line(aes(x=as.numeric(prev_group),y=mean_rbias,col=as.factor(N)),size=0.5) +
  geom_point(aes(x=as.numeric(prev_group),y=mean_rbias,col=as.factor(N)),size=1) +
  #geom_line(aes(x=as.numeric(prev_group),y=upper,col=as.factor(b))) +
  scale_color_lancet(guide=guide_legend(title="Number of samples",ncol = 2)) +
  scale_x_continuous(breaks=seq_along(levels(dat$prev_group)),labels=labels) +
  #scale_color_continuous(low="#00468BFF",high="#ED0000FF",guide=guide_colorbar(title="Number of samples")) +
  theme_light() +
  theme(#axis.text.x=element_text(angle=45,hjust=1),
    axis.text.y=element_text(size=6),
    axis.title.y=element_text(size=8),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.title=element_text(size=8),
    legend.text=element_text(size=6),
    legend.position=c(0.8,0.6),
        axis.title.x=element_blank()) +
  #scale_y_log10() +
  ylab("Average relative bias (%)") +
  labs(tag = "B")

p2 <- dat %>% group_by(N) %>% filter(b == min(b)) %>% 
  group_by(N,b,prev_group) %>% summarize(mean_rbias = sqrt(mean(MSE)),
                                           sd_rbias = sd(MSE)) %>%
  ggplot() + 
  geom_line(aes(x=as.numeric(prev_group),y=mean_rbias,col=as.factor(N)),size=0.5) +
  geom_point(aes(x=as.numeric(prev_group),y=mean_rbias,col=as.factor(N)),size=1) +
  scale_y_log10(labels=function(x) sprintf("%g", x)) +
  scale_x_continuous(breaks=seq_along(levels(dat$prev_group)),labels=labels) +
  #scale_color_continuous(low="#00468BFF",high="#ED0000FF",guide=guide_colorbar(title="Number of samples")) +
  scale_color_lancet(guide=guide_legend(title="Number of samples",byrow = TRUE)) +
  theme_light() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        legend.position="none") +
  xlab("Prevalence in population (%)") +
  ylab("Average root mean square error") +
  labs(tag = "C")

p3 <- dat %>% filter(label %in% c("576_6_96","2304_48_48","18432_48_384")) %>%
  select(N, b, q, t, p_mle, label, p_true) %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label=factor(label,levels=c("N=576, b=6, n=96","N=2304, b=48, n=48","N=18432, b=48, n=384"))) %>%
  #pivot_longer(c(p_true,p_sample),names_to="objective") %>%
  group_by(label,p_true) %>%
  summarize(mean_est=mean(p_mle),
            sd_est=sd(p_mle)) %>%
  mutate(upper=mean_est+sd_est,
         lower=mean_est-sd_est) %>%
  mutate(lower = pmax(2e-5,lower)) %>%
  ggplot() +
  geom_abline(slope=1,intercept=0,linetype="dashed",col="#ED0000FF") +
  geom_ribbon(aes(x=p_true,ymin=lower,ymax=upper),alpha=0.25) +
  geom_line(aes(x=p_true,y=mean_est)) +
  ylab("Estimated prevalence") +
  xlab("Prevalence in population (per capita)") +
  scale_y_log10(labels=function(x) sprintf("%g", x)) +
  scale_x_log10(labels=function(x) sprintf("%g", x)) +
  facet_wrap(~label,nrow=1) +
  theme_light() +
theme(axis.text.x=element_text(size=6),
      axis.text.y=element_text(size=6),
      axis.title.y=element_text(size=8),
      axis.title.x=element_text(size=8),
      legend.title=element_text(size=8),
      legend.text=element_text(size=6),
      strip.text = element_text(size=8,color="black"),
      strip.background = element_blank(),
      legend.position="bottom")
main_p <- (p3/p1/p2) + plot_layout(heights=c(1.5,2,2))
main_p
ggsave("fig3_test.png",main_p,height=8,width=7,units="in")
ggsave("fig3_simple.png",p3,height=2.5,width=7,units="in")



pS3 <- dat %>% filter(label %in% c("576_6_96","2304_48_48","18432_48_384")) %>%
  select(N, b, q, t, p_mle, label, p_sample,sample_prev_group) %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label=factor(label,levels=c("N=576, b=6, n=96","N=2304, b=48, n=48","N=18432, b=48, n=384"))) %>%
  mutate(est_diff = (p_sample - p_mle)/p_sample) %>%
  #pivot_longer(c(p_true,p_sample),names_to="objective") %>%
  group_by(label,N,t) %>%
  summarize(mean_est=mean(est_diff),
            sd_est=sd(est_diff)) %>%
  mutate(upper=mean_est+sd_est,
         lower=mean_est-sd_est) %>%
  mutate(lower = pmax(2e-5,lower)) %>%
  ggplot() +
  geom_abline(slope=1,intercept=0,linetype="dashed",col="#ED0000FF") +
  #geom_ribbon(aes(x=p_sample,ymin=lower,ymax=upper),alpha=0.25) +
  geom_line(aes(x=t,y=mean_est)) +
  ylab("Estimated prevalence") +
  xlab("Prevalence in population (per capita)") +
  #scale_y_log10(labels=function(x) sprintf("%g", x)) +
  #scale_x_log10(labels=function(x) sprintf("%g", x)) +
  geom_vline(aes(xintercept=1/N),linetype="dashed") +
  facet_wrap(~label,nrow=1) +
  theme_light() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        strip.text = element_text(size=8,color="black"),
        strip.background = element_blank(),
        legend.position="bottom")
pS3

p3_alt <- dat %>% 
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  group_by(N) %>%
  filter(b == min(b)) %>%
  mutate(label = factor(label, levels=c("N=288, b=6, n=48", 
                                           "N=576, b=6, n=96", 
                                           "N=1152, b=6, n=192", 
                                           "N=2304, b=6, n=384", 
                                           "N=4608, b=12, n=384", 
                                           "N=9216, b=24, n=384", 
                                           "N=18432, b=48, n=384"))) %>%
  select(N, b, q, t, p_mle, label, p_true) %>% 
  group_by(N, b,label,p_true) %>%
  summarize(mean_est=mean(p_mle),
            sd_est=sd(p_mle)) %>%
  mutate(upper=mean_est+sd_est,
         lower=mean_est-sd_est) %>%
  mutate(lower = pmax(2e-5,lower)) %>%
  ggplot() +
  geom_abline(slope=1,intercept=0,linetype="dashed",col="#ED0000FF") +
  geom_ribbon(aes(x=p_true,ymin=lower,ymax=upper),alpha=0.25) +
  geom_line(aes(x=p_true,y=mean_est)) +
  ylab("Estimated prevalence") +
  xlab("Prevalence in population (per capita)") +
  scale_y_log10(labels=function(x) sprintf("%g", x)) +
  scale_x_log10(labels=function(x) sprintf("%g", x)) +
  facet_wrap(~label) +
  theme_light() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        strip.text = element_text(size=8,color="black"),
        strip.background = element_blank(),
        legend.position="bottom") 

ggsave("fig3_alt.png",p3_alt,height=6,width=8,units="in")


prev_dat <- dat %>% select(t,p_true) %>% distinct()

dat$prev_group <- cut(dat$p_true*100,c(0,0.01,0.05,0.1,0.5,1,5,10,max(dat$p_true)*100))

p2 <- dat %>%  
    group_by(N) %>% filter(b == min(b)) %>%
    ungroup() %>%
    group_by(prev_group, N, b, q) %>%
    summarize(x=mean(true_diff),
              lower=quantile(true_diff,0.025),
              upper=quantile(true_diff,0.975)) %>%
    ggplot() +
    geom_ribbon(aes(x=as.numeric(prev_group),ymin=lower,ymax=upper),alpha=0.25)+
    geom_line(aes(x=as.numeric(prev_group),y=x)) +
    scale_x_continuous(labels=levels(dat$prev_group), breaks=seq_along(levels(dat$prev_group))) +
  scale_y_continuous(limits=c(-0.25,0.25)) +
  ylab("True prevalence - estimated prevalence") +
  xlab("True prevalence") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(~N)

p3 <- dat %>%  
  mutate(true_diff_pct = p_mle/p_true) %>%
  group_by(N) %>% filter(b == min(b)) %>%
  ungroup() %>%
  group_by(prev_group, N, b, q) %>%
  summarize(x=mean(true_diff_pct),
            lower=quantile(true_diff_pct,0.025),
            upper=quantile(true_diff_pct,0.975)) %>%
  ggplot() +
  geom_hline(yintercept=1,linetype="dashed",col="#ED0000FF") +
  geom_ribbon(aes(x=as.numeric(prev_group),ymin=lower,ymax=upper),alpha=0.25)+
  geom_line(aes(x=as.numeric(prev_group),y=x)) +
  scale_x_continuous(labels=levels(dat$prev_group), breaks=seq_along(levels(dat$prev_group))) +
  ylab("Estimate/true prevalence") +
  xlab("True prevalence") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(~N,scales="free_y")




dat$sample_diff <- dat$p_sample - dat$p_mle
dat$sample_diff_pct <- abs(dat$true_diff)/dat$p_sample
dat$samp_prev_group <- cut(dat$p_sample*100,c(0,0.01,0.05,0.1,0.5,1,5,10,max(dat$p_sample)*100))

p4 <- dat %>%  
  group_by(N) %>% filter(b == min(b)) %>%
  ungroup() %>%
  group_by(samp_prev_group, N, b, q) %>%
  summarize(x=mean(sample_diff),
            lower=quantile(sample_diff,0.025),
            upper=quantile(sample_diff,0.975)) %>%
  ggplot() +
  geom_ribbon(aes(x=as.numeric(samp_prev_group),ymin=lower,ymax=upper),alpha=0.25)+
  geom_line(aes(x=as.numeric(samp_prev_group),y=x)) +
  scale_x_continuous(labels=levels(dat$samp_prev_group), breaks=seq_along(levels(dat$samp_prev_group))) +
  scale_y_continuous(limits=c(-0.25,0.25)) +
  ylab("Prevalence in sample - estimated prevalence") +
  xlab("Prevalence in sample") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(~N)

p5 <- dat %>%  
  group_by(N) %>% filter(b == min(b)) %>%
  ungroup() %>%
  group_by(samp_prev_group, N, b, q) %>%
  summarize(x=mean(sample_diff_pct),
            lower=quantile(sample_diff_pct,0.025),
            upper=quantile(sample_diff_pct,0.975)) %>%
  ggplot() +
  geom_ribbon(aes(x=as.numeric(samp_prev_group),ymin=lower,ymax=upper),alpha=0.25)+
  geom_line(aes(x=as.numeric(samp_prev_group),y=x)) +
  scale_x_continuous(labels=levels(dat$samp_prev_group), breaks=seq_along(levels(dat$samp_prev_group))) +
  ylab("Relative difference between estimate and prevalence in sample") +
  xlab("Prevalence in sample") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  scale_y_continuous(limits=c(0,5)) +
  facet_wrap(~N,scales="free_y")

pSY <- dat %>% mutate(sample_diff = (p_sample - p_true)/p_true) %>%
  group_by(N) %>% filter(b == min(b)) %>%
  ungroup() %>%
  group_by(label, N, b, q, p_true) %>%
  summarize(mean_sample_err = mean(sample_diff),
            sd_sample_err = sd(sample_diff)) %>%
  mutate(upper=mean_sample_err + sd_sample_err,
         lower=mean_sample_err - sd_sample_err) %>%
  ggplot() +
  geom_ribbon(aes(x=p_true,ymin=lower,ymax=upper),alpha=0.25) +
  geom_line(aes(x=p_true,y=mean_sample_err)) +
  geom_vline(aes(xintercept=1/N),linetype="dashed") +  
  scale_x_log10() +
  xlab("Prevalence in population") +
  ylab("Relative difference between prevalence in sample and population") +
  theme_light() +
  scale_y_continuous(limits=c(-10,10)) +
  facet_wrap(~N)+
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        strip.text = element_text(size=8,color="black"),
        strip.background = element_blank(),
        legend.position="bottom") 


pSX <- dat %>%
  #group_by(N) %>% filter(b == min(b)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x=p_sample,y=p_mle),size=0.01,alpha=0.5) +
  geom_abline(slope=1,intercept=0,linetype="dashed",col="#ED0000FF") +
  scale_y_log10(limits=c(5e-6,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), labels=function(x) sprintf("%g", x)) +
  scale_x_log10(limits=c(5e-6,1),breaks=c(0.00001,0.0001,0.001,0.01,0.1,1),labels=function(x) sprintf("%g", x)) +
  geom_hline(aes(yintercept=1/N), linetype="longdash",col="grey50",size=0.5) +
  theme_light() +
  ylab("Estimated prevalence") +
  xlab("Prevalence in sample") +
  facet_wrap(~N)+
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        strip.text = element_text(size=8,color="black"),
        strip.background = element_blank(),
        legend.position="bottom") 


ggsave("figSx.png",pSX,width=8,height=7,units="in")
ggsave("figSY.png",pSY,width=8,height=7,units="in")


pSZ <- dat %>%
  #group_by(N) %>% filter(b == min(b)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x=p_true,y=p_mle),size=0.01,alpha=0.1,col="red") +
  geom_abline(slope=1,intercept=0,linetype="dashed",col="#ED0000FF") +
  scale_y_log10(limits=c(5e-6,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), labels=function(x) sprintf("%g", x)) +
  scale_x_log10(limits=c(5e-6,1),breaks=c(0.00001,0.0001,0.001,0.01,0.1,1),labels=function(x) sprintf("%g", x)) +
  geom_hline(aes(yintercept=1/N), linetype="longdash",col="blue",size=0.5) +
  geom_vline(aes(xintercept=1/N), linetype="longdash",col="blue",size=0.5) +
  theme_light() +
  ylab("Estimated prevalence") +
  xlab("Prevalence in population") +
  facet_wrap(~N)+
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        strip.text = element_text(size=8,color="black"),
        strip.background = element_blank(),
        legend.position="bottom") 
pSZ
ggsave("figSZ.png",pSZ,width=8,height=7,units="in")
