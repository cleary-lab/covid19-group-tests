library(ggplot2)
library(tidyverse)
library(patchwork)
library(extraDistr)
library(rethinking)

N <- 100000

## Pars to vary
day_undetectable_mean <- 15
day_undetectable_sd <- 3
viral_peak_sd <- 0.5
y0 <- 0
obs_sd <- 1.35
mu_wane_corr <- -0.5

lower <- c(0, 0, 0, -10, 0, -1)
upper <- c(50,10,10,0,10,1)

infer <- c(1,2,4,5)
lower <- lower[infer]
upper <- upper[infer]

pars <- c(day_undetectable_mean, day_undetectable_sd, viral_peak_sd, y0, obs_sd, mu_wane_corr)
pars <- pars[infer]

## Fixed pars
## Incubation period
lognorm_mu <- 1.621
lognorm_sd <- 0.418

tp_mean <- -2
viral_peak_mean <- 9

tp_sd <- 1

wane_min <- 0.01

delay_max <- 3
delay_min <- 0.5

tp_tinc_corr <- -0.8

LOD <- 3


ages <- c(-5,5,10,15,20,25,30,35,40)-5
desired_probs <- c(0, 1,0.9,0.72,0.49,0.29,0.14,0.04,0)
aim_prop <- tibble(t_alt=ages,prop_target=desired_probs)

ts <- seq(0,50,by=0.1)

f <- function(pars, N, return="cost") {
  day_undetectable_mean <- pars[1]
  day_undetectable_sd <- pars[2]
  #viral_peak_sd <- pars[3]
  y0 <- pars[3]
  obs_sd <- pars[4]
  #mu_wane_corr <- pars[6]
  
  mus <- c(tp_mean, 0,viral_peak_mean,day_undetectable_mean)
  sds <- c(tp_sd, 1,viral_peak_sd,day_undetectable_sd)
  
  R <- diag(length(mus))
  R[1,2] <- R[2,1] <- tp_tinc_corr
  R[3,4] <- R[4,3] <- mu_wane_corr
  colnames(R) <- c("tp","tinc","viral_peak","day_undetectable")
  rownames(R) <- c("tp","tinc","viral_peak","day_undetectable")
  
  delays <- rmvnorm2(N, mus, sds, R)
  delays <- as.data.frame(delays)
  colnames(delays) <- c("tp","tinc","viral_peak","day_undetectable")
  delays$tinc <- exp(lognorm_mu + lognorm_sd*delays$tinc)
  delays$window <- delays$tinc + delays$tp
  delays$window_min <- pmax(delays$window, delay_min)
  delays$tp1 <- delays$window_min - delays$tinc
  delays$window_max <- pmin(delays$window_min, delay_max)
  delays$tlatent <- runif(nrow(delays), delay_min, delays$window_max)
  delays <- as_tibble(delays)
  
  trajs <- matrix(nrow=nrow(delays),ncol=length(ts))
  delays$i <- 1:nrow(delays)
  
  ## Generate these trajectories
  for(i in 1:nrow(trajs)){
    tp <- delays$tp[i]
    mu <- delays$viral_peak[i]
    wane <- delays$day_undetectable[i]
    td <- delays$tlatent[i]
    tsympt <- delays$tinc[i]
    #trajs[i,] <- model_func_tinf2(ts, tsympt, tp, td, mu, wane, y0)
    trajs[i,] <- model_func_tinf2(ts, tsympt, tp, td, mu, wane, y0, LOD) + rnorm(length(ts),0,obs_sd)
  }
  
  ## Get version where trajectories are aligned to time wrt symptom onset
  trajs1 <- trajs
  trajs1 <- reshape2::melt(trajs)
  colnames(trajs1) <- c("i","t","vl")
  trajs1$t <- ts[trajs1$t]
  trajs1 <- left_join(trajs1, delays,by="i") %>% as_tibble
  trajs1$t_alt <- round(trajs1$t - trajs1$tinc,0)

  prop_detectable <- trajs1 %>% mutate(detectable=ifelse(vl >= LOD, 1, 0)) %>%
    group_by(t_alt) %>%
    summarise(prop_detectable=sum(detectable)/n())
  
  prop_detectable <- prop_detectable %>% right_join(aim_prop,by="t_alt")
  
  cost <- (prop_detectable$prop_detectable - prop_detectable$prop_target)^2
  print(sum(cost))
  if(return=="cost"){
    return(sum(cost))
  } else {
    return(trajs1)
  }
}

start_pars <- runif(length(pars), lower, upper)
start_pars <- c(15, 4, -5, 2)


#fit <- optim(pars,f, N=1000,method="L-BFGS-B",
#             lower=lower,upper=upper)
fit <- optim(start_pars,f, N=1000)
pars <- fit$par

trajs1 <- f(pars, N=N,return="traj")
trajs_summary <- trajs1 %>% group_by(t_alt) %>% summarise(lower=quantile(vl, 0.025),
                                               lower_mid=quantile(vl, 0.25),
                                         median=quantile(vl,0.5),
                                         upper_mid=quantile(vl,0.75),
                                         upper=quantile(vl,0.975))
p1 <- trajs_summary %>%  ggplot() + 
  geom_ribbon(aes(x=t_alt,ymin=lower,ymax=upper),fill="blue",alpha=0.25) +
  geom_ribbon(aes(x=t_alt,ymin=lower_mid,ymax=upper_mid),fill="blue",alpha=0.5) +
  geom_line(aes(x=t_alt,y=median)) +
  scale_x_continuous(limits=c(-20,30))+
  xlab("Time since symptom onset") +
  ylab("log10 viral load") +
  coord_cartesian(ylim=c(0,12))+
  scale_y_continuous(breaks=seq(0,12,by=1)) +
  geom_hline(yintercept=LOD,linetype="dashed")


quants <- trajs1 %>% group_by(t) %>%
  summarise(lower=quantile(vl, 0.025),
            lower_mid=quantile(vl, 0.25),
            median=quantile(vl,0.5),
            upper_mid=quantile(vl,0.75),
            upper=quantile(vl,0.975))
p2 <- ggplot(quants) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),fill="blue",alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid),fill="blue",alpha=0.5) +
  geom_line(aes(x=t,y=median)) +
  xlab("Time since infection") +
  ylab("log10 viral load") +
  coord_cartesian(ylim=c(0,12)) +
  scale_y_continuous(breaks=seq(0,12,by=1)) +
  geom_hline(yintercept=LOD,linetype="dashed")
        


p3 <- trajs1 %>% mutate(detectable=ifelse(vl >= LOD, 1, 0)) %>%
  group_by(t_alt) %>%
  summarise(prop_detectable=sum(detectable)/n()) %>%
  ggplot() +
  geom_line(aes(x=t_alt,y=prop_detectable)) +
  geom_point(data=aim_prop,aes(x=t_alt,y=prop_target),col="red") +
  xlab("Time since symptom onset") +
  ylab("Proportion detectable (>= 3)") +
scale_x_continuous(limits=c(-20,40)) +
  scale_y_continuous(breaks=seq(0,1,by=0.1))


(p1 | p2) / p3

delays_melted <- reshape2::melt(delays, id.vars="i")
#delays_melted %>% ggplot() +
#  geom_histogram(aes(x=value)) + facet_wrap(~variable, scales="free")

