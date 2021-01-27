library(tidyverse)
library(patchwork)
library(ggsci)
library(lemon)
AAAS_palette <- c("blue1"="#3B4992FF","red1"="#EE0000FF","green1"="#008B45FF",
                  "purple1"="#631879FF","teal1"="#008280FF","red2"="#BB0021FF",
                  "purple2"="#5F559BFF","purple3"="#A20056FF",
                  "grey1"="#808180FF","black"="#1B1919FF")

setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/plots/")

dat <- read_csv("fig3/summary.all_results.csv")
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
  geom_hline(aes(yintercept=1/N), linetype="dashed",col="grey50",size=0.5) +
  geom_vline(aes(xintercept=1/N), linetype="dashed",col="grey50",size=0.5) +
  geom_point(data=dat1,aes(x=p_true,y=p_mle),size=0.05,alpha=0.1,col="black",shape=19) +
  geom_abline(slope=1,intercept=0,linetype="dashed",col="#EE0000FF",size=0.5) +
  ylab("Estimated prevalence") +
  xlab("Prevalence in population (per capita)") +
  scale_y_log10(limits=c(5e-6,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), 
                labels=c("<0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  scale_x_log10(limits=c(5e-6,1),breaks=c(0.00001,0.0001,0.001,0.01,0.1,1),
                labels=c("0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  facet_wrap(~label,ncol=1) +
  theme_classic() +
  theme(axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5),
        axis.title.y=element_text(size=6),
        axis.title.x=element_text(size=6),
        legend.title=element_text(size=5),
        legend.text=element_text(size=5),
        panel.grid.major=element_line(size=0.1,color="grey50"),
        strip.text = element_text(size=5,color="black",face="bold"),
        #strip.background = element_rect(fill="grey70",color="black"),
        strip.background = element_rect(fill="white",color="white"),
        legend.position="bottom",
        plot.tag = element_text(size=10)) +
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
        strip.text.y = element_text(size=7,color="white"),
        #strip.background.y = element_rect(fill="grey70",color="black"),
        strip.background.y = element_rect(fill="grey70",color="black"),
        strip.text.x = element_text(size=6,color="black"),
        strip.background.x = element_blank(),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=7),
        axis.title.x=element_text(size=7),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        plot.tag = element_text(size=10)) +
  scale_y_log10(limits=c(5e-6,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), 
                labels=c("<0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  labs(tag="B")
fig3b


upper_lim <- 1.2
lower_lim <- 0.8

dat4 <- dat
dat_0.001 <- dat4 %>% filter(label %in% c("576_6_96","2304_24_96","18432_48_384")) %>%
  select(N, b, q, t, p_mle, label, p_true,p_sample,prev_group) %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label=factor(label,levels=c("N=576, b=6, n=96","N=2304, b=24, n=96","N=18432, b=48, n=384"))) %>%
  mutate(i=1:n(),
         group=0.00001) %>%
  ## Get +/- 10% of 0.001
  filter(p_true < 0.00001*upper_lim & p_true > 0.00001*lower_lim)

dat_0.01 <- dat4 %>% filter(label %in% c("576_6_96","2304_24_96","18432_48_384")) %>%
  select(N, b, q, t, p_mle, label, p_true,p_sample,prev_group) %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label=factor(label,levels=c("N=576, b=6, n=96","N=2304, b=24, n=96","N=18432, b=48, n=384"))) %>%
  mutate(i=1:n(),
         group=0.0001) %>%
  ## Get +/- 10% of 0.01
  filter(p_true < 0.0001*upper_lim & p_true > 0.0001*lower_lim)
dat_0.1 <- dat4 %>% filter(label %in% c("576_6_96","2304_24_96","18432_48_384")) %>%
  select(N, b, q, t, p_mle, label, p_true,p_sample,prev_group) %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label=factor(label,levels=c("N=576, b=6, n=96","N=2304, b=24, n=96","N=18432, b=48, n=384"))) %>%
  mutate(i=1:n(),group=0.001) %>%
  ## Get +/- 10% of 0.1
  filter(p_true < 0.001*upper_lim & p_true > 0.001*lower_lim)
dat_1 <- dat4 %>% filter(label %in% c("576_6_96","2304_24_96","18432_48_384")) %>%
  select(N, b, q, t, p_mle, label, p_true,p_sample,prev_group) %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label=factor(label,levels=c("N=576, b=6, n=96","N=2304, b=24, n=96","N=18432, b=48, n=384"))) %>%
  mutate(i=1:n(),
         group=0.01) %>%
  ## Get +/- 10% of 1
  filter(p_true < 0.01*upper_lim & p_true > 0.01*lower_lim)
dat_10 <- dat4 %>% filter(label %in% c("576_6_96","2304_24_96","18432_48_384")) %>%
  select(N, b, q, t, p_mle, label, p_true,p_sample,prev_group) %>%
  mutate(label = paste0("N=",N,", b=",b,", n=",q)) %>%
  mutate(label=factor(label,levels=c("N=576, b=6, n=96","N=2304, b=24, n=96","N=18432, b=48, n=384"))) %>%
  mutate(i=1:n(),
         group=0.1) %>%
  ## Get +/- 10% of 1
  filter(p_true < 0.1*upper_lim & p_true > 0.1*lower_lim)

dat5 <- bind_rows(dat_0.001, dat_0.01,dat_0.1,dat_1,dat_10) %>%
  pivot_longer(c(p_mle, p_true, p_sample))


label_key <- c("p_true"="True prevalence","p_sample"="Prevalence in sample","p_mle"="Estimated prevalence")
dat5$name <- label_key[dat5$name]
dat5$name <- factor(dat5$name, levels=label_key)

legend_label <- tibble(label=factor("N=576, b=6, n=96",levels=c("N=576, b=6, n=96","N=2304, b=24, n=96","N=18432, b=48, n=384")),
                       group="0.0001",
                       text2=c("True prevalence","Prevalence in\nsample","Estimated\nprevalence"),
                       text=factor(c("True prevalence","Prevalence in sample","Estimated prevalence"),
                                   levels=c("True prevalence","Prevalence in sample","Estimated prevalence")))

dat5$group <- sprintf("%g", dat5$group)

dat_medians <- dat5 %>% group_by(group, name, label) %>% 
  summarize(median_val=median(value)) %>%
  mutate(i=1:n())

inaccurate_levels <- dat5 %>% filter(as.numeric(group) < 1/N) %>% select(label, group) %>% distinct()

dodge_x <- 0.2
dat5 <- dat5 %>% mutate(x_dodged=as.numeric(name) + runif(n(), -dodge_x,dodge_x))

fig3b_alt <- ggplot() +
  geom_rect(data=inaccurate_levels, aes(fill="p < 1/N"), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.25) +
  geom_line(data=dat5,
            aes(x=x_dodged,y=value,group=i),size=0.001,col="grey30") + 
  geom_point(data=dat5,shape=20,alpha=0.5,
             aes(col=name,x=x_dodged,y=value),
             size=0.0001) + 
  geom_line(data=dat_medians,aes(x=name,y=median_val,group=i,linetype="Median estimate"),col="orange",size=0.5) +
  geom_hline(data=dat5,aes(yintercept=as.numeric(group)),col="#EE0000FF",size=0.25,linetype="dashed") +
  geom_rect(data=dat5 %>% filter(name=="Estimated prevalence"),aes(xmin=as.numeric(name)-0.5,xmax=as.numeric(name)+0.5,ymin=5e-6,ymax=1),
            fill="white",alpha=0,color="#008B45FF",size=0.25) +
  #geom_text(data=legend_label,aes(x=text,label=text2),y=-1,col="#AD002AFF",
  #          size=1.5,angle=90) +
  facet_wrap(label~paste0("True prevalence = ", group),ncol=4,dir="h")+#,switch="x") +
  #scale_color_lancet(guide=guide_legend(title=NULL,#"Testing strategy",
  #                                      override.aes = list(size=2,alpha=1),order=1)) +
  #scale_color_lancet(guide="none") +
  #scale_color_manual(values=c("True prevalence"="#EE0000FF","Prevalence in sample"="#3B4992FF","Estimated prevalence"="#008B45FF"),guide="none")+
  scale_color_manual(values=c("True prevalence"="#EE0000FF","Prevalence in sample"="black","Estimated prevalence"="black"),guide="none")+
  scale_linetype_discrete() +
  scale_fill_manual(values=c("grey70")) +
  guides(linetype=guide_legend(title=NULL,order=2,fill="white"),fill=guide_legend(title=NULL,order=99,fill="white")) +
  scale_x_discrete(labels=c("True prevalence","Prevalence in\nsample","Estimated\nprevalence"))+
  theme_classic() +
  ylab("Prevalence per capita (log scale)") +
  xlab("Prevalence in population (per capita)") +
  theme(#axis.text.x=element_blank(),
      axis.text.x=element_text(angle=45,hjust=1,size=6),
      panel.grid.major = element_line(colour="grey50",size=0.1),
        #axis.ticks.x=element_blank(),
      strip.text= element_text(
        margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt"),
        size=5,color="black",face="bold"),
        #strip.text.y = element_text(size=5,color="black",face="bold"),
        #strip.background.x = element_rect(size=0.1),
        strip.background = element_rect(size=0.1,color="white"),
        #strip.text.x = element_text(size=5,color="black",face='bold'),
        axis.text.y=element_text(size=5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6),
        #axis.title.x=element_text(size=8),
        legend.title=element_text(size=6),
        legend.text=element_text(size=5),
      legend.key=element_rect(fill="white",color="white"),
        legend.box.background = element_rect(fill="white",color="white"),
        legend.margin = margin(-5,-5,-5,-5),
        legend.background = element_rect(fill="white",color="white"),
        #legend.position="top",
        legend.position=c(0.87,0.1),
        plot.background = element_blank(),
        plot.tag = element_text(size=10,face="bold")) +
  scale_y_log10(limits=c(5e-6,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), 
                labels=c("<0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  labs(tag="B")
fig3b_alt

fig3_overall <- (fig3 | wrap_elements(full = fig3b_alt, clip = FALSE)) +plot_layout(widths=c(1,3))
ggsave("Fig3_revised.png",fig3_overall,height=5,width=8,units="in")
ggsave("Fig3_revised.pdf",fig3_overall,height=5,width=8,units="in")

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
  geom_hline(aes(yintercept=1/N), linetype="dashed",col="grey50",size=0.5) +
  geom_vline(aes(xintercept=1/N), linetype="dashed",col="grey50",size=0.5) +
  geom_point(aes(x=value,y=p_mle),size=0.05,alpha=0.1,col="black",shape=19) +
  geom_abline(slope=1,intercept=0,linetype="dashed",col="#ED0000FF",size=0.5) +
  ylab("Estimated prevalence") +
  xlab("Target prevalence") +
  scale_y_log10(limits=c(5e-6,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), 
                labels=c("<0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  scale_x_log10(limits=c(5e-6,1),breaks=c(0.00001,0.0001,0.001,0.01,0.1,1),
                labels=c("0.00001",sapply(c(0.0001,0.001,0.01,0.1,1), function(x) sprintf("%g",x)))) +
  facet_grid(label~name) +
  theme_classic() +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_text(size=6),
        legend.text=element_text(size=6),
        panel.grid.major=element_line(size=0.1,color="grey50"),
        strip.text = element_text(size=5,color="black",face="bold"),
        #strip.background = element_rect(fill="grey70",color="black"),
        strip.background = element_rect(fill="white",color="white"),
        legend.position="bottom",
        plot.tag = element_text(size=10)) 
figS3

ggsave("FigS5.png",figS3,height=8,width=6,units="in")
ggsave("FigS5.pdf",figS3,height=8,width=6,units="in")

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

