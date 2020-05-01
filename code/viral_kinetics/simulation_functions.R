simulate_epidemic_process <- function(n_indivs, growth_rate, times){
  incidence <- exp(growth_rate*times)/n_indivs
  incidence <- incidence - incidence[1]
  incidence[length(incidence)] <- 0
  p_infected <- sum(incidence)
  print(paste0("Cumulative incidence to date: ", signif(p_infected,3)))
  p1 <- ggplot(tibble(time=times[1:(length(times)-1)],incidence=incidence[1:(length(times)-1)])) +
    ylab("Daily probability of infection") +
    xlab("Time") +
    theme_bw() +
    geom_line(aes(x=time,y=incidence))
  return(list(plot=p1, incidence=incidence, overall_prob_infection=p_infected))
}

seir_ode <- function(t,Y,par){
  S<-Y[1]
  E<-Y[2]
  I<-Y[3]
  R<-Y[4]
  inc<-Y[5]
  N <- sum(Y[1:4])
  
  beta<-par[1]
  sigma<-par[2]
  gamma<-par[3]
  
  dYdt<-vector(length=4)
  dYdt[1]= -beta*I*S/N 
  dYdt[2]= beta*I*S/N - sigma*E
  dYdt[3]= sigma*E - gamma*I
  dYdt[4]= gamma*I
  dYdt[5] = beta*I*S/N
  
  return(list(dYdt))
}

simulate_seir_process <- function(n_indivs, pars, times){
  # Set parameter values
  R0 <- pars["R0"]
  sigma<-pars["sigma"];
  gamma<-pars["gamma"];
  beta <- R0*gamma
  N <- n_indivs
  I0 <- pars["I0"]
  R0 <- pars["R0"]
  
  init <- c(N-I0-R0,0,I0,0,R0)
  t<-times
  par<-c(beta,sigma,gamma)
  # Solve system using lsoda
  sol<-lsoda(init,t,seir_ode,par)
  # Plot solution
  sol <- as.data.frame(sol)
  colnames(sol) <- c("time","S","E","I","R","cumulative_incidence")
  incidence <- diff(c(0,sol$cumulative_incidence/N))
  prevalence <- (sol$E + sol$I)/N
  sol <- reshape2::melt(sol, id.vars="time")
  sol$value <- sol$value/N
  
  p <- ggplot(sol) + 
    geom_line(aes(x=time,y=value,col=variable)) + 
    ylab("Per capita") + 
    xlab("Date") +
    theme_bw()
  p_inc <- ggplot(data.frame(x=t,y=incidence,y1=prevalence)) + geom_line(aes(x=x,y=y),col="red") +
    geom_line(aes(x=x,y=y1),col="blue") +
    ylab("Per capita incidence (red) and prevalence (blue)") + 
    xlab("Date") +
    theme_bw()
  
  return(list(plot=p, incidence_plot=p_inc, incidence=incidence, 
              seir_outputs=sol,prevalence=prevalence,
              overall_prob_infection=sum(incidence)))  
}

simulate_infection_times <- function(n, p_infected, incidence){
  scaled_incidence <- incidence/sum(incidence)
  are_infected <- numeric(n)
  infection_times <- numeric(n)
  for(i in 1:n){
    infection <- rbinom(1,1, p_infected)
    are_infected[i] <- infection
    if(infection == 1){
      t_inf <- sample(1:length(incidence), 1, prob=scaled_incidence)
      infection_times[i] <- t_inf
    } else {
      infection_times[i] <- -1
    }
  }
  return(infection_times)
}

simulate_symptom_onsets <- function(infection_times, incubation_period_par1=1.62, incubation_period_par2=0.418){
  onset_times <- numeric(length(infection_times))
  
  for(i in seq_along(onset_times)){
    if(infection_times[i] > 0){
      tmp_onset_time <- rlnorm(1, incubation_period_par1, incubation_period_par2)
      onset_times[i] <- infection_times[i] + tmp_onset_time
    } else {
      onset_times[i] <- -1
    }
  }
  return(onset_times)
}


simulate_viral_loads <- function(infection_times, onset_times, times, kinetics_pars){
  viral_peak_mean <- kinetics_pars$viral_peak_mean
  viral_peak_sd <- kinetics_pars$viral_peak_sd
  
  wane_mean <- kinetics_pars$wane_mean
  wane_sd <- kinetics_pars$wane_sd
  
  tp_last_day <- kinetics_pars$tp_last_day
  
  viral_loads <- matrix(0, nrow=length(infection_times), ncol=length(times))
  
  n <- length(infection_times)
  
  viral_peaks <- pmax(rnorm(n, viral_peak_mean, viral_peak_sd),0.01)
  tps <- runif(n,0,tp_last_day)
  wanes <- pmax(rnorm(n, wane_mean, wane_sd),0.001)
  
  for(i in seq_along(infection_times)){
    if(i %% 1000 == 0) print(i)
    if(infection_times[i] > 0){
      viral_peak <- viral_peaks[i]
      tp <- tps[i]
      wane <- wanes[i]
      
      incubation_period <- onset_times[i] - infection_times[i]
      
      viral_load <- model_func_tinf(times, infection_times[i], tp+incubation_period, viral_peak, wane)
      #if(infection_times[i] + 60 < length(viral_load)){
      #  viral_load[(infection_times[i] + 60):length(viral_load)] <- 0
      #}
      viral_load[viral_load < 0] <- 0
      viral_loads[i,] <- viral_load
      #viral_loads[i,infection_times[i]:ncol(viral_loads)] <- viral_load[seq_along(infection_times[i]:ncol(viral_loads))]
    }
  }
  return(list(viral_loads=viral_loads,kinetics_pars=data.frame(peaks=viral_peaks,tp=tps,wane=wanes)))
}

