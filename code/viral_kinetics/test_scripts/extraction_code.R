library(data.table)
library(tidyverse)

dat <- fread("~/Downloads/used_pars_sputum.csv")
dat$i <- 1:nrow(dat)
nrow(dat)
nrow(dat[complete.cases(dat),])
dat <- dat %>% mutate(isna=is.na(infection_time))
dat[89200:89300,] %>% View
dat[1000001,]

###### SPUTUM
dat_sputum <- fread("~/Google Drive/share_sims_vl/sputum/used_pars_sputum.csv")
vl <- fread("~/Google Drive/share_sims_vl/sputum/seir_viral_loads_sputum.csv")
vl1 <- vl %>% filter(i < 100001) %>% arrange(i, t)
## How many rows
nrow(dat_sputum)

## How many rows with full non-NA values
nrow(dat_sputum[complete.cases(dat_sputum),])

## View 1 + 1 millionth row
dat_sputum[1000001,]

###### SWAB
dat_swab <- fread("~/Google Drive/share_sims_vl/swab/used_pars_swab.csv")
## How many rows
nrow(dat_swab)

## How many rows with full non-NA values
nrow(dat_swab[complete.cases(dat_swab),])

## View 1 + 1 millionth row
dat_swab[1000001,]

###### SWAB SHIFTED
dat_swab_shifted <- fread("~/Google Drive/share_sims_vl/swab_shifted/used_pars_swab_shifted.csv")
## How many rows
nrow(dat_swab_shifted)

## How many rows with full non-NA values
nrow(dat_swab_shifted[complete.cases(dat_swab_shifted),])

## View 1 + 1 millionth row
dat_swab_shifted[1000001,]
