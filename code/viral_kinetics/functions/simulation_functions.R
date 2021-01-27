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

seir_ode_switch <- function(t,Y,par){
  S<-Y[1]
  E<-Y[2]
  I<-Y[3]
  R<-Y[4]
  inc<-Y[5]
  N <- sum(Y[1:4])
  
  t_switch1 <- par[1]
  t_switch2 <- par[2]
  
  if(t < t_switch1){
    beta <- par[3]
  } else if(t > t_switch1 & t < t_switch2) {
    beta <- par[4]
  } else {
    beta <- par[5]
  }
  
  sigma<-par[6]
  gamma<-par[7]
  
  dYdt<-vector(length=4)
  dYdt[1]= -beta*I*S/N 
  dYdt[2]= beta*I*S/N - sigma*E
  dYdt[3]= sigma*E - gamma*I
  dYdt[4]= gamma*I
  dYdt[5] = beta*I*S/N
  
  return(list(dYdt))
}


simulate_seir_process <- function(n_indivs, pars, times, ver="normal",beta_smooth=0.8,stochastic=TRUE){
  ####################################################
  ## Stochastic model
  ####################################################
  if(stochastic){
    if(ver == "switch"){
      gamma1 <- pars["gamma"]
      sigma1 <- pars["sigma"]
      beta1 <- pars["R0_1"]*gamma1
      beta2 <- pars["R0_2"]*gamma1
      beta3 <- pars["R0_3"]*gamma1
      I0 <- pars["I0"]
      N <- n_indivs
      ## Odin stochastic SEIR model generator
      betas <- rep(beta3, length(times))
      betas[which(times < pars["t_switch2"])] <- beta2
      betas[which(times < pars["t_switch1"])] <- beta1
      betas <- smooth.spline(betas,spar=beta_smooth)$y
      seir <- seir_generator_interpolate(betat=times,betay=betas,sigma=sigma1,gamma=gamma1,S_ini=n_indivs-I0,I_ini=I0)
      
    } else {
      gamma1 <- pars["gamma"]
      sigma1 <- pars["sigma"]
      beta1 <- pars["R0"]*gamma1
      I0 <- pars["I0"]
      N <- n_indivs
      ## Odin stochastic SEIR model generator
      seir <- seir_generator(beta=beta1,sigma=sigma1,gamma=gamma1,S_ini=n_indivs-I0,I_ini=I0)
    }
    ## Solve model
    res <- seir$run(times)
    ## Make sure we get a simulation with an outbreak - keep trying until it takes off
    while(max(res[,"I"]) <= I0) res <- seir$run(times)
    
    ## Get raw incidence and overall probability of infection
    incidence <- res[,"inc"]/n_indivs
    #overall_prob <- max(res[,"R"])/n_indivs
    res <- as.data.frame(res)
    colnames(res) <- c("time","S","E","I","R","inc")
  } else {
    # Set parameter values
    sigma<-pars["sigma"];
    gamma<-pars["gamma"];
    N <- n_indivs
    I0 <- pars["I0"]
    recovered0 <- pars["recovered0"]
    
    init <- c(N-I0-recovered0,0,I0,0,recovered0)
    t<-times
    
    if(ver == "switch"){
      par <- c(pars["t_switch1"],pars["t_switch2"],
               pars["R0_1"]*(pars["gamma"]),pars["R0_2"]*(pars["gamma"]),pars["R0_3"]*(pars["gamma"]),
               pars["sigma"],pars["gamma"]
      )
      sol<-lsoda(init,t,seir_ode_switch,par)
    } else {
      R0 <- pars["R0"]
      beta <- R0*gamma
      par<-c(beta,sigma,gamma)
      # Solve system using lsoda
      sol<-lsoda(init,t,seir_ode,par)
    }
    # Plot solution
    res <- as.data.frame(sol)
    colnames(res) <- c("time","S","E","I","R","cumulative_incidence")
    incidence <- diff(c(0,sol$cumulative_incidence/N))
  }
  
  prevalence <- (res$E + res$I)/N
  res <- reshape2::melt(res, id.vars="time")
  res$value <- res$value/N
  
  p <- ggplot(res) + 
    geom_line(aes(x=time,y=value,col=variable)) + 
    ylab("Per capita") + 
    xlab("Date") +
    theme_bw()
  p_inc <- ggplot(data.frame(x=times,y=incidence,y1=prevalence)) + geom_line(aes(x=x,y=y),col="red") +
    geom_line(aes(x=x,y=y1),col="blue") +
    ylab("Per capita incidence (red) and prevalence (blue)") + 
    xlab("Date") +
    theme_bw()
  
  return(list(plot=p, incidence_plot=p_inc, incidence=incidence, 
              seir_outputs=res,prevalence=prevalence,
              overall_prob_infection=sum(incidence)))  
}

simulate_infection_times <- function(n, p_infected, incidence){
  scaled_incidence <- incidence/sum(incidence)
  are_infected <- rbinom(n,1,p_infected)
  t_inf <- sample(1:length(incidence), n, prob=scaled_incidence,replace=TRUE)
  are_infected[are_infected < 1] <- -1
  infection_times <- t_inf * are_infected
  infection_times[infection_times < 0] <- -1
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


simulate_viral_loads_hinge <- function(infection_times, times, chain_input, parTab,
                                       save_during=FALSE, save_block=10000,
                                       vl_file=NULL, obs_file=NULL, par_file=NULL,
                                       add_noise=TRUE,max_vl=NULL,simno=NA){
  chain <- as.matrix(chain_input)
  
  ## Re-sample rows
  n <- length(infection_times)
  samps <- sample(1:nrow(chain), n,replace=TRUE)
  
  ## Storage for simulated data
  dat_fake <- tibble(t=times,i=1) %>% as.data.frame
  
  ## If saving sims in blocks, reduce storage object size
  n_save <- n
  if(save_during){
    n_save <- save_block
  }
  
  ## Storage of simulated trajectories and used parameters
  viral_loads <- matrix(0, nrow=n_save, ncol=length(times))
  obs_dat <- matrix(0, nrow=n_save, ncol=length(times))
  used_pars <- matrix(NA,nrow=n_save,ncol=7)
  
  ## For first save, write file and save column names
  append <- FALSE
  save_colnames <- TRUE
  
  index <- 1
  block_no <- 1
  
  ## Make model functions for each possible sampled individual
  model_funcs <- NULL
  for(i in 1:9){
    dat_fake$i <- i
    model_funcs[[i]] <- create_func_indivs_multivariate_hinge(parTab[parTab$indiv %in% c(0,i),], 
                                                     dat_fake,ver="model",PRIOR_FUNC=prior_func,for_plot=FALSE)
  }
  
  ## Some pre-computation before the big loop
  R <- diag(2)
  
  ## For each infection
  for(i in seq_along(infection_times)){
    ## If individual was ever infected
    if(infection_times[i] > 0){
      
      ## Use parameters randomly from one of the 9 individuals
      choose_indiv <- sample(seq_len(9),1)
      dat_fake$i <- choose_indiv
      
      ## Get posterior draw
      pars <- chain[samps[i],2:(ncol(chain)-1)]
      pars_use <- pars[which(parTab$indiv %in% c(0,choose_indiv))]
      
      ## Construct mvrnorm sample
      rhos <- c(pars["rho_viral_wane"])
      mus <- c(pars["viral_peak_mean"], pars["wane_mean"])
      sds <- c(pars["viral_peak_sd"], pars["wane_sd"])
      
      R[upper.tri(R)] <- R[lower.tri(R)]<- rhos
      tmp <- rmvnorm2(1, mus, sds, R)
      ## Keep resampling until non-negative waning duration
      ## and viral peak
      while(tmp[1,2] < 0 | tmp[1,1] < 0 |
            (pars_use["incu"] - (pars_use["desired_mode"] + pars_use["tshift"]) + tmp[1,2]) < 0
            ){
        tmp <- rmvnorm2(1, mus, sds, R)
      }
      
      pars_use["viral_peak"] <- tmp[1,1]
      pars_use["t_wane"] <- tmp[1,2]
      
      ## STORE PARAMETERS USED
      used_pars[index,1] <- pars_use["viral_peak"]
      used_pars[index,2] <- pars_use["t_wane"]
      used_pars[index,3] <- pars_use["desired_mode"]
      used_pars[index,4] <- pars_use["tshift"]
      used_pars[index,5] <- pars_use["incu"]
      used_pars[index,6] <- pars_use["sd"]
      used_pars[index,7] <- infection_times[i]
      
      ## Shift these parameters to account for infection time
      pars_use["tshift"] <- pars_use["tshift"] + infection_times[i]
      pars_use["incu"] <- pars_use["incu"] + infection_times[i]
      
      ## Create model function and solve
      pred <- model_funcs[[choose_indiv]](pars_use)
      
      if(add_noise){
        ## Add observation error
        obs <- pred + rnorm(length(pred), 0, pars_use["sd"])
        obs_observable <- which(obs > pars_use["lod"])
        obs_unquantifiable <- which(obs < pars_use["limit_quantification"] & obs > pars_use["lod"])
        obs_below_lod <- which(obs < pars_use["lod"])
        
        obs[obs_below_lod] <- pars_use["lod"]
        obs[obs_unquantifiable] <- runif(length(obs[obs_unquantifiable]), pars_use["lod"], pars_use["limit_quantification"])
        obs[which(times < pars_use["tshift"])] <- 0
        obs_dat[index,] <- obs
      } else {
        pred[pred < pars_use["lod"]] <- pars_use["lod"]
        obs_dat[index,] <- pred
      }
      pred[pred < pars_use["lod"]] <- pars_use["lod"]
      viral_loads[index,] <- pred
    }
    index <- index + 1
    
    if(i %% save_block == 0){
      print(i)
      if(save_during){
        ## Truncate viral loads and observations to be at most max_vl
        if(!is.null(max_vl)){
          viral_loads <- pmin(viral_loads, max_vl)
          obs_dat <- pmin(obs_dat, max_vl)
        }
        
        index <- 1
        
        ## Save viral loads made so far
        vl_to_save <- Matrix::Matrix(viral_loads, sparse = TRUE)
        vl_to_save <- as.data.frame(Matrix::summary(vl_to_save))
        vl_to_save$j <- vl_to_save$j - 1
        #vl_to_save$block <- block_no
        vl_to_save$i <- (block_no-1)*save_block + vl_to_save$i
        if(!is.na(simno)){
          vl_to_save$simno <- simno
          colnames(vl_to_save) <- c("i","t","vl","simno")
        } else {
          colnames(vl_to_save) <- c("i","t","vl")
        }
        
        data.table::fwrite(vl_to_save, file = vl_file,sep = ",",col.names = save_colnames, 
                           row.names = FALSE, append = append)
        
        ## Save observations made so far
        if(add_noise){
          obs_to_save <- Matrix::Matrix(obs_dat, sparse = TRUE)
          obs_to_save <- as.data.frame(Matrix::summary(obs_to_save))
          obs_to_save$j <- obs_to_save$j - 1
          #obs_to_save$block <- block_no
          obs_to_save$i <- (block_no-1)*save_block + obs_to_save$i
          if(!is.na(simno)){
            obs_to_save$simno <- simno
            colnames(obs_to_save) <- c("i","t","vl","simno")
          } else {
            colnames(obs_to_save) <- c("i","t","vl")
          }
          data.table::fwrite(obs_to_save, file = obs_file,sep = ",",col.names = save_colnames, 
                             row.names = FALSE, append = append)
        }
        
        ## Save parameters so far
        colnames(used_pars) <- c("viral_peak","t_wane","mode","tshift","incu","sd","infection_time")
        data.table::fwrite(as.data.table(used_pars), file = par_file,sep = ",",col.names = save_colnames, 
                           row.names = FALSE, append = append)
        
        append <- TRUE
        save_colnames <- FALSE
        
        block_no <- block_no + 1
        
        ## If not done, reset storage objects. Otherwise, keep and return
        if(i < length(infection_times)){
          viral_loads <- matrix(0, nrow=n_save, ncol=length(times))
          obs_dat <- matrix(0, nrow=n_save, ncol=length(times))
          used_pars <- matrix(NA,nrow=n_save,ncol=7)
        }
      }
    }
  }
  if(!is.null(max_vl)){
    viral_loads <- pmin(viral_loads, max_vl)
    obs_dat <- pmin(obs_dat, max_vl)
  }
  return(list(viral_loads=viral_loads,obs=obs_dat, pars=used_pars))
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


simulate_viral_loads_by_inf <- function(n_indivs, times, kinetics_pars, incubation_period_par1=1.62, incubation_period_par2=0.418){
  viral_peak_mean <- kinetics_pars$viral_peak_mean
  viral_peak_sd <- kinetics_pars$viral_peak_sd
  
  wane_mean <- kinetics_pars$wane_mean
  wane_sd <- kinetics_pars$wane_sd
  
  tp_last_day <- kinetics_pars$tp_last_day
  
  viral_loads <- matrix(0, nrow=n_indivs, ncol=length(times))
  
  viral_peaks <- pmax(rnorm(n_indivs, viral_peak_mean, viral_peak_sd),0.01)
  tps <- runif(n_indivs,0,tp_last_day)
  wanes <- pmax(rnorm(n_indivs, wane_mean, wane_sd),0.001)
  
  for(i in seq_len(n_indivs)){
    if(i %% 1000 == 0) print(i)
      viral_peak <- viral_peaks[i]
      tp <- tps[i]
      wane <- wanes[i]
      
      incubation_period <- rlnorm(1, incubation_period_par1, incubation_period_par2)
      
      viral_load <- model_func_tinf(times, 0, tp+incubation_period, viral_peak, wane)
      viral_load[viral_load < 0] <- 0
      viral_loads[i,] <- viral_load
  }
  return(list(viral_loads=viral_loads,kinetics_pars=data.frame(peaks=viral_peaks,tp=tps,wane=wanes)))
}


simulate_viral_loads_by_inf2 <- function(infection_times, times, kinetics_pars,y0,LOD){
  n_indivs <- length(infection_times)
  viral_loads <- matrix(0, nrow=n_indivs, ncol=length(times))
  
  for(i in seq_along(infection_times)){
    if(i %% 1000 == 0) print(i)
    if(infection_times[i] > 0){
      tp <- kinetics_pars$tp[i]
      mu <- kinetics_pars$viral_peak[i]
      wane <- kinetics_pars$day_undetectable[i]
      td <- kinetics_pars$tlatent[i]
      tsympt <- kinetics_pars$tinc[i]
      tinf <- infection_times[i]
      
      viral_load <- model_func_tinf2_sim(times, tinf, tsympt, tp, td, mu, wane, y0, LOD)
      viral_load[viral_load < y0] <- y0
      viral_loads[i,] <- viral_load
    }
  }
  return(viral_loads=viral_loads)
}


