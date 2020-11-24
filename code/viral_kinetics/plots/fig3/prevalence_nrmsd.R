library(tidyverse)
library(patchwork)
library(ggpubr)

setwd("~/Google Drive/nCoV/pool_samples/plots/")
dat <- read_csv("~/Google Drive/nCoV/pool_samples/prevalence_data/AllResults.csv")
colnames(dat) <- c("b","m","t","trial","true_freq","sampled_freq","est_freq")
dat$N <- dat$b*dat$m

## Get SEIR results
prevalence <- unique(dat[,c("t","true_freq")])
## Epidemic switch point
peak_t <- prevalence %>% filter(true_freq == max(true_freq)) %>% pull(t)
## Group by growth or decline phase and buckets of prevalence
prevalence <- prevalence %>% 
  mutate(period=ifelse(t <= peak_t, "Growth","Decline")) %>%
  mutate(prev_group=cut(true_freq, c(0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,1),include.lowest=TRUE))
prevalence$period <- factor(prevalence$period, levels=c("Growth","Decline"))
## Get residuals between estimated prevalence and true prevalence
## Calculate ratios of estimates, log ratio and squared residuals
dat <- dat %>% mutate(residual_ratio = est_freq/true_freq,
                      log_residual_ratio = log(residual_ratio),
                      residual = est_freq - true_freq,
                      residual_sqrd = residual^2)

## Join with prevalence data
dat <- dat %>% left_join(prevalence)

## Distribution of residuals
p_residuals_density_2304_6 <- dat %>% filter(N == 2304 & b == 6) %>% ggplot() + 
  geom_density(aes(x=residual,y=..scaled..),fill="grey") + 
  geom_vline(xintercept=0,linetype="dashed") +
  ggtitle("N=2304, b = 6") +
  scale_x_continuous(limits=c(-0.25,0.25)) +
  xlab("Estimated prevalence - true prevalence") +
  ylab("Scaled density") +
  theme_pubr() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.position="right")+
  facet_grid(prev_group~period)
pdf("residual_densities_2304_6.pdf",width=5,height=7)
p_residuals_density_2304_6
dev.off()


p_residuals_density_2304_12 <- dat %>% filter(N == 2304 & b == 12) %>% ggplot() + 
  geom_density(aes(x=residual,y=..scaled..),fill="grey") + 
  geom_vline(xintercept=0,linetype="dashed") +
  ggtitle("N=2304, b = 12") +
  scale_x_continuous(limits=c(-0.25,0.25)) +
  xlab("Estimated prevalence - true prevalence") +
  ylab("Scaled density") +
  theme_pubr() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.position="right")+
  facet_grid(prev_group~period)
pdf("residual_densities_2304_12.pdf",width=5,height=7)
p_residuals_density_2304_12
dev.off()


p_residuals_density_2304_24 <- dat %>% filter(N == 2304 & b == 24) %>% ggplot() + 
  geom_density(aes(x=residual,y=..scaled..),fill="grey") + 
  geom_vline(xintercept=0,linetype="dashed") +
  ggtitle("N=2304, b = 24") +
  scale_x_continuous(limits=c(-0.25,0.25)) +
  xlab("Estimated prevalence - true prevalence") +
  ylab("Scaled density") +
  theme_pubr() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.position="right")+
  facet_grid(prev_group~period)
pdf("residual_densities_2304_24.pdf",width=5,height=7)
p_residuals_density_2304_24
dev.off()


p_residuals_density_2304_48 <- dat %>% filter(N == 2304 & b == 48) %>% ggplot() + 
  geom_density(aes(x=residual,y=..scaled..),fill="grey") + 
  geom_vline(xintercept=0,linetype="dashed") +
  ggtitle("N=2304, b = 48") +
  scale_x_continuous(limits=c(-0.25,0.25)) +
  xlab("Estimated prevalence - true prevalence") +
  ylab("Scaled density") +
  theme_pubr() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.position="right")+
  facet_grid(prev_group~period)
pdf("residual_densities_2304_48.pdf",width=5,height=7)
p_residuals_density_2304_48
dev.off()


## Distribution of ratios to true value
p_ratio_2304_6 <- dat %>% filter(N == 2304 & b == 6) %>% ggplot() + 
  geom_density(aes(x=log_residual_ratio,y=..scaled..),fill="grey") + 
  geom_vline(xintercept=0,linetype="dashed") +
  ggtitle("N=2304, b = 6") +
  #scale_x_continuous(limits=c(-0.25,0.25)) +
  xlab("log(estimated prevalence / true prevalence)") +
  ylab("Scaled density") +
  theme_pubr() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.position="right")+
  facet_grid(prev_group~period)
pdf("p_ratio_2304_6.pdf",width=5,height=7)
p_ratio_2304_6
dev.off()


## Distribution of ratios to true value
p_ratio_2304_12 <- dat %>% filter(N == 2304 & b == 12) %>% ggplot() + 
  geom_density(aes(x=log_residual_ratio,y=..scaled..),fill="grey") + 
  geom_vline(xintercept=0,linetype="dashed") +
  ggtitle("N=2304, b = 12") +
  #scale_x_continuous(limits=c(-0.25,0.25)) +
  xlab("log(estimated prevalence / true prevalence)") +
  ylab("Scaled density") +
  theme_pubr() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.position="right")+
  facet_grid(prev_group~period)
pdf("p_ratio_2304_12.pdf",width=5,height=7)
p_ratio_2304_12
dev.off()

## Distribution of ratios to true value
p_ratio_2304_24 <- dat %>% filter(N == 2304 & b == 24) %>% ggplot() + 
  geom_density(aes(x=log_residual_ratio,y=..scaled..),fill="grey") + 
  geom_vline(xintercept=0,linetype="dashed") +
  ggtitle("N=2304, b = 24") +
  #scale_x_continuous(limits=c(-0.25,0.25)) +
  xlab("log(estimated prevalence / true prevalence)") +
  ylab("Scaled density") +
  theme_pubr() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.position="right")+
  facet_grid(prev_group~period)
pdf("p_ratio_2304_24.pdf",width=5,height=7)
p_ratio_2304_24
dev.off()


p_ratio_2304_48 <- dat %>% filter(N == 2304 & b == 48) %>% ggplot() + 
  geom_density(aes(x=log_residual_ratio,y=..scaled..),fill="grey") + 
  geom_vline(xintercept=0,linetype="dashed") +
  ggtitle("N=2304, b = 48") +
  #scale_x_continuous(limits=c(-0.25,0.25)) +
  xlab("log(estimated prevalence / true prevalence)") +
  ylab("Scaled density") +
  theme_pubr() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.position="right")+
  facet_grid(prev_group~period)
pdf("p_ratio_2304_48.pdf",width=5,height=7)
p_ratio_2304_48
dev.off()

## Summarise to get RMSE/RMSE and normalised RMSE
summary <- dat %>% 
  group_by(b,m,prev_group, period, N) %>% 
  summarise(RMSE=sqrt(mean(residual_sqrd)),
            NRMSE = RMSE/mean(true_freq))

## Get factor levels and order by total samples then number of pools 
## (more samples better, more pools better)
summary_levels <- summary %>% 
  ungroup() %>%
  arrange(-N, -b) %>% 
  select(b, m, N) %>% 
  unique() %>%
  mutate(label=paste0(N, "=",b,"x",m))
summary <- summary %>% left_join(summary_levels)
summary$label <- factor(summary$label, levels=rev(summary_levels$label))

## Overall plot
p_main <- ggplot(summary) + 
  geom_tile(aes(x=prev_group,y=label,fill=NRMSE)) + 
  ylab("Number of samples = number of pools x size of pool") +
  xlab("True prevalence") +
  #scale_fill_gradient2(low="blue",mid="orange",high="red",midpoint=2,limits=c(0,2.5)) +
  scale_fill_viridis_c(limits=c(0,2.5)) + 
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.position="right")+
  facet_wrap(~period,ncol=1)
  #facet_grid(N~period)
pdf("overall_plot.pdf",height=5,width=7)
p_main
dev.off()


## Overall plot with RMSE
p_main_rmse <- ggplot(summary) + 
  geom_tile(aes(x=prev_group,y=label,fill=RMSE)) + 
  ylab("Number of samples = number of pools x size of pool") +
  xlab("True prevalence") +
  #scale_fill_gradient2(low="blue",mid="orange",high="red",midpoint=2,limits=c(0,2.5)) +
  scale_fill_viridis_c() + 
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.position="right")+
  facet_wrap(~period,ncol=1)
#facet_grid(N~period)
p_main_rmse
pdf("rmse_plot.pdf",height=5,width=7)
p_main_rmse
dev.off()

p_by_N <- ggplot(summary) + 
  geom_tile(aes(x=prev_group,y=as.factor(b),fill=NRMSE)) + 
  ylab("Number of pools") +
  xlab("True prevalence") +
  scale_fill_viridis_c(limits=c(0,2.5)) + 
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position="right")+
  facet_grid(period~N)

pdf("p_by_n.pdf",height=5,width=7)
p_by_N
dev.off()

p_288 <- ggplot(summary[summary$N == 288,]) + 
  geom_tile(aes(x=prev_group,y=as.factor(b),fill=NRMSE)) + 
  ylab("Number of pools") +
  xlab("True prevalence") +
  scale_fill_viridis_c() + 
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.position="right")+
  ggtitle("N=288") +
  facet_wrap(~period,ncol=1)

p_576<- ggplot(summary[summary$N == 576,]) + 
  geom_tile(aes(x=prev_group,y=as.factor(b),fill=NRMSE)) + 
  ylab("Number of pools") +
  xlab("True prevalence") +
  scale_fill_viridis_c() + 
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.position="right")+
  ggtitle("N=576") +
  facet_wrap(~period,ncol=1)


p_1152 <- ggplot(summary[summary$N == 1152,]) + 
  geom_tile(aes(x=prev_group,y=as.factor(b),fill=NRMSE)) + 
  ylab("Number of pools") +
  xlab("True prevalence") +
  scale_fill_viridis_c() + 
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.position="right")+
  ggtitle("N=1152") +
  facet_wrap(~period,ncol=1)


p_2304 <- ggplot(summary[summary$N == 2304,]) + 
  geom_tile(aes(x=prev_group,y=as.factor(b),fill=NRMSE)) + 
  ylab("Number of pools") +
  xlab("True prevalence") +
  scale_fill_viridis_c() + 
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.position="right")+
  ggtitle("N=2304") +
  facet_wrap(~period,ncol=1)

p_4608 <- ggplot(summary[summary$N == 4608,]) + 
  geom_tile(aes(x=prev_group,y=as.factor(b),fill=NRMSE)) + 
  ylab("Number of pools") +
  xlab("True prevalence") +
  scale_fill_viridis_c() + 
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.position="right")+
  ggtitle("N=4608") +
  facet_wrap(~period,ncol=1)

pdf("p_288.pdf",height=5,width=7)
p_288
dev.off()
pdf("p_576.pdf",height=5,width=7)
p_576
dev.off()
pdf("p_1152.pdf",height=5,width=7)
p_1152
dev.off()
pdf("p_2304.pdf",height=5,width=7)
p_2304
dev.off()
pdf("p_4608.pdf",height=5,width=7)
p_4608
dev.off()