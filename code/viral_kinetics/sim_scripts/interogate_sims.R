library(tidyverse)

## Swab viral loads
swab_vl <- read_csv("~/Google Drive/nCoV/sims_for_brian/swab/seir_viral_loads_swab_1.csv")
head(swab_vl)
range(swab_vl$vl)
prev_detectable <- swab_vl %>% group_by(t) %>% mutate(is_pos = vl > 2) %>% summarize(prev=sum(is_pos)/12500000)
prev_detectable %>% ggplot() + geom_line(aes(x=t,y=prev))
swab_vl[1:500,] %>% group_by(i) %>% mutate(t = t - min(t)) %>% ggplot() + geom_line(aes(x=t,y=vl,group=i))
swab_vl %>% filter(t %in% c(90,120,150)) %>% ggplot() + geom_histogram(aes(x=vl)) + facet_wrap(~t)

## Swab viral loads with observation error
swab_obs <- read_csv("~/Google Drive/nCoV/sims_for_brian/swab/seir_obs_viral_loads_swab_1.csv")
head(swab_obs)
range(swab_obs$vl)
prev_detectable_obs <- swab_obs %>% group_by(t) %>% mutate(is_pos = vl > 2) %>% summarize(prev=sum(is_pos)/12500000)
prev_detectable_obs %>% ggplot() + geom_line(aes(x=t,y=prev))
swab_obs[1:500,] %>% group_by(i) %>% mutate(t = t - min(t)) %>% ggplot() + geom_line(aes(x=t,y=vl,group=i))
swab_obs %>% filter(t %in% c(90,120,150)) %>% ggplot() + geom_histogram(aes(x=vl)) + facet_wrap(~t)

## Sputum viral loads
sputum_vl <- read_csv("~/Google Drive/nCoV/sims_for_brian/sputum/seir_obs_viral_loads_sputum_1.csv")
head(sputum_vl)
range(sputum_vl$vl)
prev_detectable_sputum <- sputum_vl %>% group_by(t) %>% mutate(is_pos = vl > 2) %>% summarize(prev=sum(is_pos)/12500000)
prev_detectable_sputum %>% ggplot() + geom_line(aes(x=t,y=prev))
sputum_vl[1:500,] %>% group_by(i) %>% mutate(t = t - min(t)) %>% ggplot() + geom_line(aes(x=t,y=vl,group=i))
sputum_vl %>% filter(t %in% c(90,120,150)) %>% ggplot() + geom_histogram(aes(x=vl)) + facet_wrap(~t)

## Sputum viral loads with observation error
sputum_obs <- read_csv("~/Google Drive/nCoV/sims_for_brian/sputum/seir_obs_viral_loads_sputum_1.csv")
head(sputum_obs)
range(sputum_obs$vl)
prev_detectable_obs_sputum <- sputum_obs %>% group_by(t) %>% mutate(is_pos = vl > 2) %>% summarize(prev=sum(is_pos)/12500000)
prev_detectable_obs_sputum %>% ggplot() + geom_line(aes(x=t,y=prev))
sputum_obs[1:500,] %>% group_by(i) %>% mutate(t = t - min(t)) %>% ggplot() + geom_line(aes(x=t,y=vl,group=i))
sputum_obs %>% filter(t %in% c(90,120,150)) %>% ggplot() + geom_histogram(aes(x=vl)) + facet_wrap(~t)


## Look at one of the switch model ones observations
switch_vl <- read_csv("~/Google Drive/nCoV/sims_for_brian/swab_switch_SEIR/seir_viral_loads_swab_switch_SEIR_2.csv")
switch_vl <- switch_vl %>% arrange(i, t)
used_pars <- read_csv("~/Google Drive/nCoV/sims_for_brian/swab_switch_SEIR/used_pars_swab_switch_SEIR_1.csv")
hist(switch_vl$vl)
head(switch_vl)
range(switch_vl$vl)
prev_detectable_switch <- switch_vl %>% group_by(t) %>% mutate(is_pos = vl > 2) %>% summarize(prev=sum(is_pos)/12500000)
prev_detectable_switch %>% ggplot() + geom_line(aes(x=t,y=prev))
switch_vl[1:200,] %>% group_by(i) %>% mutate(t = t - min(t)) %>% ggplot() + geom_line(aes(x=t,y=vl,group=i,col=as.factor(i))) + facet_wrap(~i)
used_pars1 <- used_pars[unique(switch_vl[1:200,]$i),]

switch_vl %>% filter(t %in% c(90,120,150)) %>% ggplot() + geom_histogram(aes(x=vl)) + facet_wrap(~t)


## Look at one of the switch model ones
switch_obs <- read_csv("~/Google Drive/nCoV/sims_for_brian/swab_switch_SEIR/seir_obs_viral_loads_swab_switch_SEIR_1.csv")
head(switch_obs)
range(switch_obs$vl)
prev_detectable_switch_obs <- switch_obs %>% group_by(t) %>% mutate(is_pos = vl > 2) %>% summarize(prev=sum(is_pos)/12500000)
prev_detectable_switch_obs %>% ggplot() + geom_line(aes(x=t,y=prev))
switch_obs[1:200,] %>% group_by(i) %>% mutate(t = t - min(t) - 5) %>% ggplot() + geom_line(aes(x=t,y=vl,group=i))
switch_obs %>% filter(t %in% c(90,120,150)) %>% ggplot() + geom_histogram(aes(x=vl)) + facet_wrap(~t)
