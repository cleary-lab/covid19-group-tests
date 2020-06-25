viral_load_func_hinge <- function(pars, obs_t, convert_ct=TRUE){
  tshift <- pars["tshift"]
  desired_mode <- pars["desired_mode"]
  
  height <- pars["viral_peak"]
  obs_sd <- pars["obs_sd"]
  true_0 <- pars["true_0"]
  yintercept <- pars["intercept"]
  t_onset <- pars["incu"]
  tw <- pars["t_wane"]
  
  y <- numeric(length(obs_t))
  wane_rate <- height / (t_onset - desired_mode + tw)
  
  y[obs_t <= tshift] <- true_0
  y[obs_t > tshift & obs_t <= desired_mode] <- height * (obs_t[obs_t > tshift & obs_t <= desired_mode]-tshift) / (desired_mode-tshift)
  y[obs_t > desired_mode] <- height - wane_rate * (obs_t[obs_t > desired_mode] - desired_mode)
  
  ct <- y
  if(convert_ct){
    ct <- yintercept - log2(10) * (ct-pars["LOD"])
  }
  ct
}



## Create the posterior function for model fitting
create_func_indivs_multivariate <- function(parTab, dat, PRIOR_FUNC=NULL,ver="model", 
                               remove_pre_tp=FALSE, draw_incubation_period=FALSE,
                               incu_par1=1.57,incu_par2=0.65){
  
  ## We'll be creating a vector of predictions that is the same length
  ## as the vector of observations
  obs <- dat$obs
  ts <- dat$t
  par_names <- parTab$names
  
  n_indivs <- length(unique(dat$indiv))
  
  ## Keep track of which individual we're solving the model for
  indiv_ids <- parTab$indiv
  indiv_ids_unique <- unique(indiv_ids)
  indiv_ids_unique <- indiv_ids_unique[indiv_ids_unique != 0]
  
  f <- function(pars){
    names(pars) <- par_names
    
    ## Overall parameters
    viral_peaks <- pars[which(par_names == "viral_peak")]
    wane_pars <- pars[which(par_names == "wane")]
    tp_pars <- pars[which(par_names == "tp")]
    
    predicted <- NULL
    ## Solve model for each individual
    for(i in indiv_ids_unique){
      ts_subset <- dat[dat$i == i, "t"]
      
      tinf <- 0
      incu_period <- 0
      if(draw_incubation_period){
        incu_period <- rlnorm(1, incu_par1, incu_par2)
        tinf <- 21 - incu_period
        tinf <- max(0, tinf)
      }
      
      subset_pars <- pars[which(indiv_ids  == i)]
      viral_peak <- subset_pars["viral_peak"]
      tp <- subset_pars["tp"] + incu_period
      wane <- subset_pars["wane"]
      
      y <- numeric(length(ts_subset))
      growth <- viral_peak/tp
      y[ts_subset < tinf] <- 0
      y[ts_subset <= tp + tinf & ts_subset > tinf] <- growth*(ts_subset[ts_subset <= tp + tinf & ts_subset > tinf]-tinf)
      y[ts_subset > tp + tinf] <- viral_peak + wane*tp - wane*(ts_subset[ts_subset > tp + tinf] - tinf)
    
      ## Do we want to include the pre-peak predictions in the likelihood?
      if(remove_pre_tp){
        y[which(ts_subset < tp)] <- NA
      }
      predicted <- c(predicted, y)
      
    }
    
    ## Can use same function to return model predictions or likelihoods
    if(ver == "model") {
      return(predicted)
    } else {
      ## Likelihoods for observations
      lik <- likelihood(obs, predicted, pars["sd"],pars["max_titre"], pars["lod"])
      
      
      all_pars <- matrix(nrow=length(viral_peaks),ncol=3)
      all_pars[,1] <- viral_peaks
      all_pars[,2] <- wane_pars
      all_pars[,3] <- tp_pars
      
      mus <- c(pars["viral_peak_mean"], pars["wane_mean"], pars["tp_mean"])
      rhos <- c(pars["rho_viral_wane"],pars["rho_viral_tp"],pars["rho_wane_tp"])
      sds <- c(pars["viral_peak_sd"], pars["wane_sd"], pars["tp_sd"])
      
      R <- diag(length(rhos))
      R[upper.tri(R)] <- R[lower.tri(R)]<- rhos
      
      ## Correlation prior
      corr_prior <- dlkjcorr(R, eta=2, TRUE)
      
      ## Sd priors
      sd_prior_viral <- dhcauchy(pars["viral_peak_sd"], 1, TRUE)
      sd_prior_wane <- dhcauchy(pars["wane_sd"], 1, TRUE)
      sd_prior_tp <- dhcauchy(pars["tp_sd"], 1, TRUE)
      
      ## Kinetics parameters
      kinetics_prior <- sum(dmvnorm2(all_pars, mus, sds, R, log=TRUE))
      
      random_effects <- kinetics_prior + sd_prior_viral + sd_prior_wane + sd_prior_tp + corr_prior
      
      #lik <- sum(lik) + random_effects
      lik <- random_effects
      if(!is.null(PRIOR_FUNC)){
        lik <- lik + PRIOR_FUNC(pars)
      }
      lik
      
    }
  }
}