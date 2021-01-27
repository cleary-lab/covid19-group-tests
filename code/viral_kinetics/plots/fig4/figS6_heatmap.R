library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggsci)

setwd("~/Documents/GitHub/covid19-group-tests/code/plots/fig4/")

dat <- read_csv(file = "summary.resource_t0-80_t1-108.csv")
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
  filter(version == "Simple pools of 4" &
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
  geom_segment(x=4.5,xend=4.5,y=1.5,yend=3.5, size=0.5, col="black") +
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
  ggtitle("Design maximizing positive sample identification for prevalence range 1.03%-9.90%") +
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


ggsave("figS6.png",p1, height=5,width=6,units="in")
ggsave("figS6.pdf",p1, height=5,width=6,units="in")

