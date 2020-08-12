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

dat1 <- dat %>% filter(label %in% c("576_6_96","2304_24_96","18432_48_384")) %>%
  select(N, b, q, t, p_mle, label, p_true) %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label=factor(label,levels=c("N=576, b=6, n=96","N=2304, b=24, n=96","N=18432, b=48, n=384")))
  
fig3 <- dat1 %>%
  ggplot() +
  geom_hline(aes(yintercept=1/N), linetype="dashed",col="#925E9FFF",size=0.5) +
  geom_vline(aes(xintercept=1/N), linetype="dashed",col="#925E9FFF",size=0.5) +
  geom_point(data=dat1,aes(x=p_true,y=p_mle),size=0.05,alpha=0.1,col="black",shape=19) +
  geom_abline(slope=1,intercept=0,linetype="dashed",col="#ED0000FF",size=0.5) +
  ylab("Estimated prevalence") +
  xlab("Prevalence in population (per capita)") +
  scale_y_log10(limits=c(5e-6,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), 
                labels=c("<0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  scale_x_log10(limits=c(5e-6,1),breaks=c(0.00001,0.0001,0.001,0.01,0.1,1),
                labels=c("0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  facet_wrap(~label,ncol=1) +
  theme_light() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        strip.text = element_text(size=8,color="white"),
        strip.background = element_rect(fill="grey70",color="black"),
        legend.position="bottom",
        plot.tag = element_text(size=10,face="bold")) +
  labs(tag="A")

ggsave("fig3_scatter.png",fig3,height=2.5,width=7,units="in")


dat2 <- dat
dat2$prev_group <- cut(dat2$p_true*100, c(0,0.01,0.1,1,10,100))
dat2 <- dat2 %>% filter(label %in% c("576_6_96","2304_24_96","18432_48_384")) %>%
  select(N, b, q, t, p_mle, label, p_true,p_sample,prev_group) %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label=factor(label,levels=c("N=576, b=6, n=96","N=2304, b=24, n=96","N=18432, b=48, n=384"))) %>%
  mutate(i=1:n()) %>%
  pivot_longer(c(p_mle,p_true,p_sample))

label_key <- c("p_true"="Prevalence in population","p_sample"="Prevalence in sample","p_mle"="Estimated prevalence")
dat2$name <- label_key[dat2$name]
dat2$name <- factor(dat2$name, levels=label_key)

legend_label <- tibble(label=factor("N=576, b=6, n=96",levels=c("N=576, b=6, n=96","N=2304, b=24, n=96","N=18432, b=48, n=384")),
                       prev_group="(0,0.01]",
                       text2=c("Prevalence in\npopulation","Prevalence in\nsample","Estimated\nprevalence"),
                       text=factor(c("Prevalence in population","Prevalence in sample","Estimated prevalence"),
                                   levels=c("Prevalence in population","Prevalence in sample","Estimated prevalence")))


fig3b <- dat2 %>% 
  ggplot() + 
  geom_point(aes(x=name,y=value),size=0.05) + 
  geom_line(aes(x=name,y=value,group=i),size=0.005,col="grey50") + 
  geom_text(data=legend_label,aes(x=text,label=text2),y=-1,col="#AD002AFF",
            size=1.5,angle=90) +
  facet_grid(label~prev_group,switch="x") +
  theme_light() +
  ylab("Prevalence per capita (log scale)") +
  xlab("Prevalence range") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.y = element_text(size=8,color="white"),
        strip.background.y = element_rect(fill="grey70",color="black"),
        strip.text.x = element_text(size=6,color="black"),
        strip.background.x = element_blank(),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        plot.tag = element_text(size=10,face="bold")) +
  scale_y_log10(limits=c(5e-6,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), 
                labels=c("<0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  labs(tag="B")
fig3b

fig3_overall <- (fig3 | wrap_elements(full = fig3b, clip = FALSE)) +plot_layout(widths=c(1,3))
ggsave("fig3_overall.png",fig3_overall,height=5,width=8,units="in")
ggsave("fig3_overall.pdf",fig3_overall,height=5,width=8,units="in")

dat %>%
  group_by(N) %>% filter(b == min(b)) %>% select(N, b, q) %>% distinct() %>% 
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>% pull(label) -> labels
dput(labels)
labels <- labels[c(1,6,2,7,5,3,4)]

figS3 <- dat %>%
  group_by(N) %>% filter(b == min(b)) %>%
  ungroup() %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label = factor(label, levels=labels)) %>% 
  select(N,b,q,p_mle,p_sample,p_true,label) %>%
  rename(`Prevalence in population` = p_true,
         `Prevalence in sample` = p_sample) %>%
  pivot_longer(c(`Prevalence in population`, `Prevalence in sample`)) %>%
  ggplot() +
  geom_hline(aes(yintercept=1/N), linetype="dashed",col="#925E9FFF",size=0.5) +
  geom_vline(aes(xintercept=1/N), linetype="dashed",col="#925E9FFF",size=0.5) +
  geom_point(aes(x=value,y=p_mle),size=0.05,alpha=0.1,col="black",shape=19) +
  geom_abline(slope=1,intercept=0,linetype="dashed",col="#ED0000FF",size=0.5) +
  ylab("Estimated prevalence") +
  xlab("Target prevalence") +
  scale_y_log10(limits=c(5e-6,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), 
                labels=c("<0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  scale_x_log10(limits=c(5e-6,1),breaks=c(0.00001,0.0001,0.001,0.01,0.1,1),
                labels=c("0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  facet_grid(label~name) +
  theme_light() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        strip.text = element_text(size=6,color="white"),
        strip.background = element_rect(fill="grey70",color="black"),
        legend.position="bottom")
figS3

ggsave("figS3_scatter.png",figS3,height=8,width=6,units="in")

labels1 <- paste0("N=",unique(dat$N))

figS4 <- dat %>% mutate(sample_diff = (p_sample - p_true)/p_true) %>%
  group_by(N) %>% filter(b == min(b)) %>%
  ungroup() %>%
  group_by(label, N, b, q, p_true) %>%
  summarize(mean_sample_err = mean(sample_diff),
            sd_sample_err = sd(sample_diff)) %>%
  mutate(upper=mean_sample_err + sd_sample_err,
         lower=mean_sample_err - sd_sample_err) %>%
  mutate(label=factor(paste0("N=",N), levels=labels1)) %>%
  ggplot() +
  geom_ribbon(aes(x=p_true,ymin=lower,ymax=upper),alpha=0.25) +
  geom_line(aes(x=p_true,y=mean_sample_err)) +
  geom_vline(aes(xintercept=1/N),linetype="dashed") +  
  scale_x_log10() +
  xlab("Prevalence in population") +
  ylab("Relative difference between prevalence in sample and population") +
  theme_light() +
  scale_y_continuous(limits=c(-10,10)) +
  facet_wrap(~label)+
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        strip.text = element_text(size=6,color="white"),
        strip.background = element_rect(fill="grey70",color="black"),
        legend.position="bottom") 

ggsave("figS4.png",figS4,height=6,width=8,units="in")



ggsave("fig3b.png",fig3b,height=6,width=8,units="in")
