library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggsci)
library(ggthemes)
library(ggpubr)

setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/plots/bimodal_results/")


## Incidence/prevalence plot
inc_dat <- read_csv("incidence.csv")
prev_dat <- read_csv("prev_vl2.csv")

## Prevalence range
prev_dat %>% filter(t >= 150) %>% summarize(range(prev*100))


prev_dat$Variable <- "Prevalence (viral load > 100)"
inc_dat$Variable <- "Incidence (new infections)"
colnames(prev_dat) <- c("t","value","Variable")
colnames(inc_dat) <- c("t","value","Variable")
plot_dat <- bind_rows(prev_dat,inc_dat)

p_dat <- ggplot(plot_dat) + geom_line(aes(x=t,y=value,col=Variable)) +
  theme_pubr() +
  #scale_x_continuous(expand=c(0,0), limits=c(0,200),breaks=seq(0,200,by=50)) +
  scale_color_manual(values=c("#ed0000ff","#00466bff")) +
  scale_y_continuous(limits=c(0,0.02),breaks=seq(0,0.02,by=0.005)) +
  theme(legend.position=c(0.5,0.7),
        panel.grid.major = element_line(color="grey80",size=0.25),
        axis.text=element_text(family="sans",size=10),
        axis.title=element_text(family="sans",size=10),
        legend.text = element_text(family="sans",size=8),
        legend.title=element_text(family="sans",size=8),
        plot.margin = margin(0,5,0,0,"mm")) +
  ylab("Per capita infected") +
  xlab("Days since epidemic start")

dat <- read_csv(file = "summary.resource_t0-150_t1-200.csv")
dat$`Average samples tested per day` <- floor(dat$`Design samples`*dat$`Design runs per day`)
dat <- dat %>% 
  mutate(version = ifelse(`Design sample split`== 1, 
                          ifelse(`Design pools` == 1, "Simple", "Individual"),
                          "Combinatorial"))
dat1 <- dat

dat$`Test budget` <- as.factor(dat$`Test budget`)
dat$`Sample budget` <- as.factor(dat$`Sample budget`)
dat$label <- paste0("N:",dat$`Design samples`,", ",
                    "b:",dat$`Design pools`,", ",
                    "q:",dat$`Design sample split`,"\n",
                    #"Runs: ", round(dat$`Design runs per day`, 1))
                    "Tested: ", floor(dat$`Design runs per day` * dat$`Design samples`))

dat <- dat %>% mutate(label1 = ifelse(version =="Combinatorial", label, 
                         ifelse(version == "Simple", paste0("Tested: ", 
                                                            floor(dat$`Design runs per day` * dat$`Design samples`)),NA)))

#dat <- dat %>% mutate(`Design effectiveness` = ifelse(version == "Combinatorial", `Design effectiveness`, NA))

indiv_test_boundary <- dat %>% 
  filter((version == "Individual" & `Test budget` == `Design samples`)) %>%
  mutate(i=n())



simple_pool_boundary <- dat %>% 
  filter(version == "Simple") %>%
  filter(`Design samples` == min(`Design samples`))

dat <- dat %>% mutate(version = ifelse(version == "Simple", paste0("Simple pools of ", `Design samples`), version))

simple_pool_boundary2 <- dat %>% 
  filter(version == "Simple pools of 8" &
           as.numeric(`Sample budget`) > 3)


dat <- dat %>% group_by(`Test budget`) %>% mutate(label1=ifelse(label1==lag(label1, 1) & version == "Combinatorial", NA,label1))
dat$version <- factor(dat$version, levels=c("Individual","Simple pools of 4","Simple pools of 8","Combinatorial"))

p1 <- ggplot(dat) + 
  geom_tile(aes(y=`Test budget`, x=`Sample budget`,
                fill=`Average samples tested per day`),col="grey90") +
  geom_point(aes(y=`Test budget`, x=`Sample budget`,
                 shape=`version`),size=7,col="grey50",alpha=0.1) +
  geom_text(aes(y=`Test budget`, x=`Sample budget`,
                label=label1,col=`Average samples tested per day` > 3500),size=1.25) +
  #geom_tile(data=dat %>% filter(version == "Simple" & `Design samples` == 8),
  #          aes(y=`Test budget`, x=`Sample budget`), fill="#925E9FFF") +
  #geom_tile(data=dat %>% filter(version == "Simple" & `Design samples` == 4),
  #          aes(y=`Test budget`, x=`Sample budget`), fill="#925E9F99") +
  #geom_text(x=4,y=7,label="Individual testing optimal \n (1 test per sample)",fontface="plain",
  #          size=4, family="sans") +
  geom_step(data=indiv_test_boundary,
             mapping=aes(x=`Sample budget`,y=`Test budget`,group=1),
             position = position_nudge(y = 0.5,x=0.5),size=0.5, col="black") +
     geom_segment(x=1.5,xend=1.5,y=0,yend=1.5, size=0.5, col="black") +
     #geom_step(data=simple_pool_boundary,
    #           mapping=aes(x=`Sample budget`,y=`Test budget`,group=1),
    #           position = position_nudge(y = 0.5,x=0.5),size=0.5, col="black") +
  geom_segment(x=3.5,xend=3.5,y=0,yend=1.5, size=0.5, col="black") +
  geom_segment(x=4.5,xend=4.5,y=1.5,yend=2.5, size=0.5, col="black") +
  geom_segment(x=3.5,xend=4.5,y=1.5,yend=1.5, size=0.5, col="black") +
  geom_step(data=simple_pool_boundary2,
            mapping=aes(x=`Sample budget`,y=`Test budget`,group=1),
            position = position_nudge(y = 0.5,x=0.5),size=0.5, col="black") +
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0)) +
  scale_color_manual(values=c("black","white"),guide="none") +
  scale_fill_continuous(low="white", high="#00468BFF",
                       limits=c(0,6144),
                         guide=guide_colorbar(direction="horizontal",
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              title.vjust=5,
                                              ticks=FALSE,
                                              title="Average number of samples tested",
                                              barheight=0.5,barwidth=8,frame.color="grey40",
                                              order=1),
                        na.value="white") +
                        #limits=c(dat %>% filter(version == "Combinatorial") %>% 
                        #           summarize(x=min(`Design effectiveness`)) %>% 
                        #           pull(x),
                        #         max(dat$`Design effectiveness`,na.rm=TRUE))) +
  scale_shape_manual(values=c(15,16,17,18),
                     guide=guide_legend(title=NULL,#"Testing strategy",
                                        title.position="bottom",
                                        title.hjust=0.5,
                                        order=0,
                                        #nrow=2,
                                        #byrow=TRUE,
                                        override.aes = list(size=3,alpha=1))) +
  guides() +
  ggtitle("Design maximizing positive sample identification for bimodal prevalence range 0.373%-1.82%") +
  ylab("Daily test capacity") +
  xlab("Daily sample collection capacity") +
  theme_classic() +
  theme(panel.grid=element_blank(),
        legend.position="bottom",
        plot.title=element_text(hjust=0.5,size=8),
        axis.ticks = element_line(color="grey40"),
        axis.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.spacing.x = unit(0.1, 'pt'))
p1



summary_dat <- dat1 %>% 
  filter(version != "Individual") %>%
  mutate(version2 = ifelse(version == "Simple", paste0("Simple pools of ", `Design samples`),
                           "Combinatorial")) %>%
  #mutate(version3 = ifelse(version2 == "Combinatorial", 
  #                         ifelse(`Design sample split` == 2, "Combinatorial q=2", "Combinatorial q=3 or 4"),version2)) %>%
  group_by(version2) %>%
  summarize(mean_effectiveness=mean(`Design effectiveness`),
            min_effectiveness=min(`Design effectiveness`),
            max_effectiveness=max(`Design effectiveness`)) %>%
  mutate(min_effectiveness = ifelse(version2 == "Combinatorial", min_effectiveness, NA),
         max_effectiveness = ifelse(version2 == "Combinatorial", max_effectiveness, NA))
summary_dat$version2 <- factor(summary_dat$version2, 
                               #              levels=rev(c("Simple pools of 4","Simple pools of 8","Combinatorial q=2",
                               #                          "Combinatorial q=3 or 4")))
                               levels = rev(c("Simple pools of 4","Simple pools of 8","Combinatorial")))

summary_dat$`Stage of outbreak` <- "Prevalence 0.373%-1.82%"


p3 <- ggplot(summary_dat, aes(y=version2)) + 
  geom_bar(aes(x=mean_effectiveness,fill=`Stage of outbreak`),stat="identity",col="black",size=0.25,
           position=position_dodge()) +
  geom_errorbar(aes(xmin=min_effectiveness,xmax=max_effectiveness,group=`Stage of outbreak`),
                width=0.15,position=position_dodge(.9)) +
  coord_cartesian(xlim=c(1,10.5)) +
  scale_x_continuous(breaks=c(1,seq(2,10,by=2)),expand=c(0,0)) +
  scale_fill_lancet(guide = guide_legend(reverse = TRUE))+
  ylab("") +
  #xlab("Effectiveness relative to individual testing") +
  xlab("Fold increase in positives identified \nrelative to individual testing") +
  theme_classic() +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=8),
        plot.title=element_text(hjust=0.5,size=8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        legend.position = c(0.6,0.8),
        plot.margin = margin(0,5,0,-5))


ggsave("heatmap.pdf",p1,height=5,width=6,units="in")
ggsave("heatmap.png",p1,height=5,width=6,units="in",dpi=300)


ggsave("both.pdf",(p1 | p3) + plot_layout(widths=c(3,1)),height=5,width=8,units="in")
ggsave("both.png",(p1 | p3) + plot_layout(widths=c(3,1)),height=5,width=8,units="in",dpi=300)


ggsave("barchart.pdf",p3,height=5,width=3,units="in")
ggsave("barchart.png",p3,height=5,width=3,units="in",dpi=300)


ggsave("incprev.pdf",p_dat,height=3,width=9,units="in")
ggsave("incprev.png",p_dat,height=3,width=9,units="in",dpi=300)

