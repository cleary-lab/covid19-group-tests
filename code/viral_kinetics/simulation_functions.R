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

simulate_symptom_onsets <- function(infection_times, incubation_period_par1=1.57, incubation_period_par2=0.65){
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
  
  for(i in seq_along(infection_times)){
    if(i %% 1000 == 0) print(i)
    if(infection_times[i] > 0){
      viral_peak <- max(rnorm(1, viral_peak_mean, viral_peak_sd),0.01)
      tp <- runif(1,0,tp_last_day)
      wane <- max(rnorm(1, wane_mean, wane_sd),0.001)
      incubation_period <- onset_times[i] - infection_times[i]
      
      viral_load <- model_func(times, tp+incubation_period, viral_peak, wane)
      viral_load[viral_load < 0] <- 0
      
      viral_loads[i,infection_times[i]:ncol(viral_loads)] <- viral_load[seq_along(infection_times[i]:ncol(viral_loads))]
    }
  }
  return(viral_loads)
}

