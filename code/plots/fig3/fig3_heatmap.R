library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggsci)
library(ggthemes)

setwd("~/Documents/GitHub/covid19-group-tests/code/plots/")

dat <- read_csv("summary.growth.nrmse.csv")

head(dat)

dat <- dat %>%
  mutate(`Total samples` = as.factor(`Total samples`),
         Pools = as.factor(Pools))

options(digits=5)

p1 <- ggplot(dat) +
  geom_line(aes(x=`Prevalence in population (lower end of range)`,
                y=NRMSE_pop,
                linetype=Pools)) +
  geom_hline(yintercept=1, linetype="solid", col="#AD002AFF",size=0.5) +
  facet_wrap(~`Total samples`,ncol=4) +
  ylab("Normalized root mean square error") +
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1),
                labels = function(x) sprintf("%g", x)) +
  scale_y_continuous(breaks=seq(0,10,by=2)) +
  theme_bw() +
  theme(
    ## Axis text and titles
    axis.text.x = element_text(size=8,family="sans"),
    axis.text.y=element_text(size=8,family="sans"),
    axis.title.x=element_text(size=8,family="sans",vjust=-1),
    axis.title.y=element_text(size=8,family="sans"),
    
    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    
    ## Title
    plot.title = element_text(family="sans",size=8,face="bold",hjust=0.5),
    
    ## Legends
    legend.position=c(0.9,0.2),
    legend.title=element_text(size=8,family="sans"),
    legend.text=element_text(size=8,family="sans"),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="sans"),
    strip.background=element_rect(fill="#f0f0f0"))


ggsave("fig3.png",p1, height=5,width=8,units="in")
ggsave("fig3.pdf",p1, height=5,width=8,units="in")


dat <- read_csv("summary.growth.nrmse.csv")
dat <- dat %>%  mutate(`Total samples` = as.factor(`Total samples`))
p2 <- ggplot(dat %>% group_by(`Total samples`)  %>% filter(Pools==min(Pools))) +
  geom_line(aes(x=`Prevalence in population (lower end of range)`,
                y=NRMSE_pop,
                col=`Total samples`)) +
  geom_hline(yintercept=1, linetype="dashed", col="grey40",size=0.5) +
  ylab("Normalized root mean square error") +
  scale_x_log10(        
    breaks=c(0.0001,0.001,0.01,0.1),
                labels = function(x) sprintf("%g", x)) +
  scale_y_continuous(breaks=seq(0,10,by=2)) +
  scale_color_lancet() +
  theme_bw() +
  theme(
    ## Axis text and titles
    axis.text.x = element_text(size=8,family="sans"),
    axis.text.y=element_text(size=8,family="sans"),
    axis.title.x=element_text(size=8,family="sans",vjust=-1),
    axis.title.y=element_text(size=8,family="sans"),
    
    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    
    ## Title
    plot.title = element_text(family="sans",size=8,face="bold",hjust=0.5),
    
    ## Legends
    legend.position=c(0.75,0.75),
    legend.title=element_text(size=8,family="sans"),
    legend.text=element_text(size=8,family="sans"),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="sans"),
    strip.background=element_rect(fill="#f0f0f0"))
p2

ggsave("fig3_alt.png",p2, height=4,width=6,units="in")
ggsave("fig3_alt.pdf",p2, height=4,width=6,units="in")
