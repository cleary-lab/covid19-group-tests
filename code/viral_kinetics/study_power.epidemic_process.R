library(tidyverse)
library(ggplot2)
library(deSolve)
library(ggpubr)
source("simulation_functions.R")
source("model_funcs.R")

## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 6900000

## Duration of epidemic in days
times <- 0:365

initial_incidence <- 0.005

# generate 2 incidence curves
# curve 1
r0 <- 1.5
recovered <- 0.1
seir_pars <- c("R0"=r0,"gamma"=1/7,"sigma"=1/6.4,"I0"= population_n*initial_incidence, "recovered0"=population_n*recovered)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)
pdf('../../examples/study_power.incidence.high.pdf')
epidemic_process$plot
dev.off()
incidence = viral_loads <- as_tibble(epidemic_process$incidence)
write_csv(incidence,"../../examples/study_power.incidence.high.csv")
# 
# curve 2
r0 <- 1.1
recovered <- 0.2
seir_pars <- c("R0"=r0,"gamma"=1/7,"sigma"=1/6.4,"I0"= population_n*initial_incidence, "recovered0"=population_n*recovered)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)
pdf('../../examples/study_power.incidence.low.pdf')
epidemic_process$plot
dev.off()
incidence = viral_loads <- as_tibble(epidemic_process$incidence)
write_csv(incidence,"../../examples/study_power.incidence.low.csv")