library(tidyverse)

setwd("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/data/")

dat1 <- read.csv("patient1.csv", header=FALSE) %>% mutate(i = 1)
dat2 <- read.csv("patient2.csv", header=FALSE) %>% mutate(i = 2)
dat3 <- read.csv("patient3.csv", header=FALSE) %>% mutate(i = 3)
dat4 <- read.csv("patient4.csv", header=FALSE) %>% mutate(i = 4)
dat5 <- read.csv("patient5.csv", header=FALSE) %>% mutate(i = 5)
dat6 <- read.csv("patient6.csv", header=FALSE) %>% mutate(i = 6)
dat7 <- read.csv("patient7.csv", header=FALSE) %>% mutate(i = 7)
dat8 <- read.csv("patient8.csv", header=FALSE) %>% mutate(i = 8)
dat9 <- read.csv("patient9.csv", header=FALSE) %>% mutate(i = 9)


dat <- bind_rows(dat1, dat2, dat3,
                 dat4, dat5, dat6, 
                 dat7, dat8, dat9)

colnames(dat) <- c("t","obs","i")

dat$t <- round(dat$t, 0)
dat$obs <- round(dat$obs, 2)
dat$obs[dat$obs < 0.5] <- 0

ggplot(dat) +
  geom_line(aes(x=t,y=obs))+
  facet_wrap(~i)

write_csv(dat, "drosten_sputum.csv")
