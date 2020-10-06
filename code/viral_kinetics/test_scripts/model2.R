n_indiv <- 1000



tinc <- rlnorm(n_indiv,1.621,0.418)
tp <- runif(n_indiv, -2, 1)
tp_mod <- tp - pmin(0, tinc + tp - 0.5)
tg <- runif(n_indiv, 0.2, 3)
tg_mod <- tg + pmin(0, tinc + tp_mod - tg - 0.2)

peak_time <- tinc + tp_mod
start_time <- tinc + tp_mod - tg

viral_peaks <- rnorm(n_indiv, 8, 1)
twane <- rnorm(n_indiv,15,2)


t <- seq(0,30,by=0.1)
vls <- matrix(nrow=n_indiv,ncol=length(t))

pars <- cbind(tinc,tp_mod, tg_mod, viral_peaks, twane)

mus <- c(0, -2, -0.1,8, 15)
sigmas <- c(1,0.5,0.5,1,2)
R <- diag(length(sigmas))
colnames(R) <- rownames(R) <- c("tinc","tp","tg","vl","tw")

#R["tinc","tp"] <- R["tp","tinc"] <- -0.8
R["tinc","tg"] <- R["tg","tinc"] <- 0.8
#R["tp","tg"] <- R["tg","tp"] <- 0.8


pars <- rmvnorm2(n_indiv, mus, sigmas, R)
pars[,1] <- exp(1.621 + 0.418*pars[,1])
colnames(pars) <- c("tinc","tp","tg","vl","tw")
pars <- as.data.frame(pars)
pars$tg <- exp(pars$tg)

pairs(pars)


par(mfrow=c(2,3))
hist(pars$tinc)
hist(pars$tp)
hist(pars$tg)
hist(pars$vl)
hist(pars$tw)
par(mfrow=c(1,1))


for(i in 1:nrow(pars)){
  tinc1 <- pars[i,1]
  tp1 <- pars[i,2]
  tg1 <- pars[i,3]
  vl1 <- pars[i,4]
  tw1 <- pars[i,5]
  
  tp1 <- tp1 - min(0, tinc1 + tp1 - 0.5)
  tg1 <- tg1 + min(0, tinc1 + tp1 - tg1 - 0.2)
  
  wane_rate <- vl1/(tinc1 + tp1 + tw1)
  growth_rate <- vl1/tg1
  
  y <- numeric(length(t))
  
  latent_period <- (tinc1 + tp1 - tg1)
  
  y[t < latent_period] <- 0
  
  y[t >= latent_period & t < (tinc1 + tp1)] <- growth_rate * (t[t >= latent_period & t < (tinc1 + tp1)] - latent_period)
  
  y[t >= (tinc1 + tp1)] <- vl1 - wane_rate*(t[t >= (tinc1 + tp1)] - (tinc1 + tp1))
  
  vls[i,] <- y
}
vls[vls < 0] <- 0
plot(vls[1,],type='l',ylim=c(0,15))
apply(vls, 1, lines)


pars_mod <- pars
pars_mod$tp_mod <- pars_mod$tp - pmin(0, pars_mod$tinc + pars_mod$tp - 0.5)
pars_mod$tg_mod <- pars_mod$tg + pmin(0, pars_mod$tinc + pars_mod$tp_mod - pars_mod$tg - 0.2)

pairs(pars_mod)

pars_mod_melt <- reshape2::melt(pars_mod)

ggplot(pars_mod_melt) + geom_histogram(aes(x=value)) + facet_wrap(~variable,scales="free")
