likelihood_vl <- function(obs, predicted, sd, min_vl=0, limit_quantification=2){
  
  normal_vl <- which(obs > limit_quantification)
  unquantified_vl <- which(obs <= limit_quantification & obs > min_vl)
  low_vl <- which(obs <= min_vl)
  
  #liks <- dnorm(obs, predicted, sd, TRUE)
  liks <- numeric(length(obs))
  liks[normal_vl] <- dnorm(obs[normal_vl],predicted[normal_vl],sd,TRUE)
  liks[unquantified_vl] <- log(pnorm(limit_quantification, predicted[unquantified_vl], sd, lower.tail=TRUE) - pnorm(min_vl, predicted[unquantified_vl], sd, lower.tail=TRUE))
  liks[low_vl] <- pnorm(min_vl, predicted[low_vl], sd, lower.tail=TRUE, TRUE)
  
  
  #liks[normal_vl] <- dgumbel(obs[normal_vl],predicted[normal_vl],sigma=sd,TRUE)
  #liks[low_vl] <- pgumbel(min_vl,predicted[low_vl],sigma=sd, lower.tail=TRUE, TRUE)
  
  liks
}


## Create the posterior function for model fitting
create_func_indivs_multivariate_hinge<- function(parTab, dat, PRIOR_FUNC=NULL,ver="model", 
                                                 for_plot=FALSE){
  
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
    wane_pars <- pars[which(par_names == "t_wane")]
    t_shifts <- pars[which(par_names == "tshift")]
    desired_modes <- pars[which(par_names == "desired_mode")]
    true_0 <- pars[which(par_names == "true_0")]
    t_onsets <- pars[which(par_names == "incu")]
    
    predicted <- NULL
    pred_dat <- NULL
    ## Solve model for each individual
    for(i in indiv_ids_unique){
      ## Extract this individual's parameters
      subset_pars <- pars[which(indiv_ids  == i)]
      viral_peak <- subset_pars["viral_peak"]
      tw <- subset_pars["t_wane"]
      tshift <- subset_pars["tshift"]
      desired_mode <- subset_pars["desired_mode"] + tshift
      t_onset <- subset_pars["incu"]
      
      ## Observed data is shifted by incubation period
      if(ver == "posterior"){
        obs_t <- dat[dat$i == i, "t"] + t_onset
      } else {
        obs_t <- dat[dat$i == i, "t"]
      }
      
      y <- rep(true_0, length(obs_t))
      wane_rate <- viral_peak / (t_onset - desired_mode + tw)
      
      ## Latent period
      y[obs_t <= tshift] <- true_0
      ## Growth period
      y[obs_t > tshift & obs_t <= desired_mode] <- viral_peak * (obs_t[obs_t > tshift & obs_t <= desired_mode]-tshift) / (desired_mode-tshift)
      ## Wane period
      y[obs_t > desired_mode] <- viral_peak - wane_rate * (obs_t[obs_t > desired_mode] - desired_mode)
      
      if(for_plot){
        use_t <- obs_t - t_onset
        pred_dat <- bind_rows(pred_dat, tibble(t=use_t,y=y,i=i))
      } else {
        predicted <- c(predicted, y)
      }
      
    }
    
    ## Can use same function to return model predictions or likelihoods
    if(ver == "model") {
      if(for_plot){
        return(pred_dat)
      } else {
        return(predicted)
      }
    } else {
      ## Likelihoods for observations
      #lik <- likelihood(obs, predicted, pars["sd"],pars["max_titre"], pars["lod"])
      lik <- likelihood_vl(obs, predicted, pars["sd"],pars["lod"], pars["limit_quantification"])
      all_pars <- matrix(nrow=length(viral_peaks),ncol=2)
      all_pars[,1] <- viral_peaks
      all_pars[,2] <- wane_pars
      #all_pars[,3] <- desired_modes
      
      mus <- c(pars["viral_peak_mean"], pars["wane_mean"])#, pars["tp_mean"])
      rhos <- c(pars["rho_viral_wane"])#,pars["rho_viral_tp"]),pars["rho_wane_tp"])
      sds <- c(pars["viral_peak_sd"], pars["wane_sd"])#, pars["tp_sd"])
      
      #R <- diag(length(rhos))
      R <- diag(2)
      R[upper.tri(R)] <- R[lower.tri(R)]<- rhos
      
      ## Correlation prior
      corr_prior <- dlkjcorr(R, eta=2, TRUE)
      
      ## Sd priors
      sd_prior_viral <- dhcauchy(pars["viral_peak_sd"], 1, TRUE)
      sd_prior_wane <- dhcauchy(pars["wane_sd"], 1, TRUE)
      #sd_prior_tp <- dhcauchy(pars["tp_sd"], 1, TRUE)
      
      ## Kinetics parameters
      kinetics_prior <- sum(dmvnorm2(all_pars, mus, sds, R, log=TRUE))
      #random_effects <- 0
      random_effects <- kinetics_prior + sd_prior_viral + sd_prior_wane + corr_prior # + sd_prior_tp
      
      lik <- sum(lik) + random_effects
      #lik <- random_effects
      if(!is.null(PRIOR_FUNC)){
        lik <- lik + PRIOR_FUNC(pars)
      }
      lik
      
    }
  }
}