library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggsci)

setwd("~/Documents/GitHub/covid19-group-tests/code/plots/fig4/")

dat <- read_csv(file = "summary.resource_t0-190_t1-250.tmp.csv")
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
  filter(`Design samples` == max(`Design samples`))

dat <- dat %>% mutate(version = ifelse(version == "Simple", paste0("Simple pools of ", `Design samples`), version))
dat <- dat %>% group_by(`Test budget`) %>% mutate(label1=ifelse(label1==lag(label1, 1) & version == "Combinatorial", NA,label1))
dat$version <- factor(dat$version, levels=c("Individual","Simple pools of 2", "Simple pools of 4",
                                            "Simple pools of 8","Simple pools of 16",
                                            "Combinatorial"))

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
     geom_step(data=simple_pool_boundary,
               mapping=aes(x=`Sample budget`,y=`Test budget`,group=1),
               position = position_nudge(y = 0.5,x=0.5),size=0.5, col="black") +
     geom_segment(x=4.5,xend=4.5,y=0,yend=1.5, size=0.5, col="black") +
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0)) +
  scale_color_manual(values=c("black","white"),guide="none") +
  scale_fill_continuous(low="white", high="#AD002AFF",
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
  scale_shape_manual(values=c(15,4,16,3,18),
                     guide=guide_legend(title=NULL,#"Testing strategy",
                                        title.position="bottom",
                                        title.hjust=0.5,
                                        nrow=2,
                                        order=0,
                                        #nrow=2,
                                        #byrow=TRUE,
                                        override.aes = list(size=3,alpha=1))) +
  guides() +
  ggtitle("Design maximizing positive sample identification for prevalence range 0.03%-2.21%\n during epidemic decline") +
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
                    levels = rev(c("Simple pools of 2","Simple pools of 4","Simple pools of 8",
                                   "Simple pools of 16","Combinatorial")))

summary_dat$`Stage of outbreak` <- "Prevalence 0.03%-2.21%"

summary_dat_all <- bind_rows(summary_dat)

summary_dat_all$`Stage of outbreak` <- factor(summary_dat_all$`Stage of outbreak`,
                                                 levels=rev(c("Prevalence 1.03%-9.90%","Prevalence 0.03%-2.21%")))


p3 <- ggplot(summary_dat_all, aes(y=version2)) + 
  geom_bar(aes(x=mean_effectiveness),fill="#AD002AFF",stat="identity",col="black",size=0.25,
           position=position_dodge()) +
  geom_errorbar(aes(xmin=min_effectiveness,xmax=max_effectiveness,group=`Stage of outbreak`),
                width=0.15,position=position_dodge(.9)) +
  coord_cartesian(xlim=c(1,20.5)) +
  scale_x_continuous(breaks=c(1,seq(5,20.5,by=5)),expand=c(0,0)) +
  #scale_fill_lancet(guide = guide_legend(reverse = TRUE))+
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
p_main <- (p1 | p3) +  plot_layout(widths=c(3,1))
p_main
ggsave("fig4_decline.png",p_main, height=5,width=8,units="in")
ggsave("fig4_decline.pdf",p_main, height=5,width=8,units="in")
