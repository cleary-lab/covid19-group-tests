setwd("~/Google Drive/share_sims_vl/swab/")
library(tidyverse)
dat <- read_csv("seir_viral_loads_swab.csv")

prev <- dat %>% group_by(t) %>% mutate(is_pos=vl > 0) %>% summarize(y=sum(is_pos)/12500000)

prev <- prev %>% 
  mutate(phase1 = "peak") %>%
  mutate(phase1 = ifelse(y < 0.1 & t < 150,"growth",phase1),
         phase1 = ifelse(y < 0.1 & t > 150,"decline",phase1))

dat <- dat %>% 
  mutate(phase="peak") %>%
  mutate(phase = ifelse(t <= 108, "growth", phase),
         phase = ifelse(t >= 168, "decline", phase))

prev1 <- dat %>% 
  group_by(t) %>% 
  mutate(is_pos=vl > 2) %>% 
  summarize(n_detectable=sum(is_pos),n=n()) %>%
  mutate(sens=n_detectable/n) %>%
  left_join(prev)

prev2 <- dat %>% 
  left_join(prev) %>%
  group_by(phase1) %>% 
  mutate(is_pos=vl > 2) %>% 
  summarize(n_detectable=sum(is_pos),n=n()) %>%
  mutate(sens=n_detectable/n)

prev3 <- dat %>%
  mutate(is_pos=vl > 2) %>% 
  summarize(n_detectable=sum(is_pos),n=n()) %>%
  mutate(sens=n_detectable/n)


prev1 %>% ggplot() + geom_line(aes(x=t,y=sens))
