## Viral kinetics function
model_func <- function(ts, tp, viral_peak, wane){
  y <- numeric(length(ts))
  growth <- viral_peak/tp
  y[ts <= tp] <- growth*ts[ts <= tp]
  y[ts > tp] <- viral_peak + wane*tp - wane*(ts[ts > tp])
  y
}
## Viral kinetics function
model_func_tinf <- function(ts, tinf, tp, viral_peak, wane){
  y <- numeric(length(ts))
  growth <- viral_peak/tp
  y[ts <= tinf] <- 0
  y[ts <= tp + tinf & ts > tinf] <- growth*(ts[ts <= tp + tinf & ts > tinf]-tinf)
  y[ts > tp + tinf] <- viral_peak + wane*tp - wane*(ts[ts > tp + tinf] - tinf)
  y
}
## Calculate likelihood of observations given model predicted viral loads
## Assuming normally distributed measurement error, with censoring
## at upper and lower limtis of detection
likelihood <- function(obs, predicted, sd, max_titre, lod){
  high_titres <- which(obs >= max_titre)
  normal_titres <- which(obs < max_titre & obs > lod)
  low_titres <- which(obs <= lod)
  
  liks <- numeric(length(obs))
  
  liks[high_titres] <- pnorm(max_titre, predicted[high_titres], sd, lower.tail=FALSE, TRUE)
  liks[normal_titres] <- dnorm(obs[normal_titres], predicted[normal_titres], sd, TRUE)
  liks[low_titres] <- pnorm(lod, predicted[low_titres], sd, lower.tail=TRUE, TRUE)
  
  liks
}

## Get gamma shape and scale from specified mean and variance
gamma_pars_from_mean_sd <- function(gamma_mean, gamma_var){
  scale <- gamma_var/gamma_mean
  shape <- gamma_mean/scale
  return(list("shape"=shape,"scale"=scale))
}


## Create the posterior function for model fitting
create_func_indivs <- function(parTab, dat, PRIOR_FUNC=NULL,ver="model", 
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

      #growth <- viral_peak/tp
      #pre_peak <- growth*ts_subset[ts_subset <= tp]
      #post_peak <- viral_peak + wane*tp - wane*(ts_subset[ts_subset > tp])
      
      #y <- c(pre_peak, post_peak)
      
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
      
      #viral_peak_gamma_pars <- gamma_pars_from_mean_sd(pars["viral_peak_mean"],pars["viral_peak_sd"]^2)
      #viral_peaks_random_effects <- sum(dgamma(viral_peaks, viral_peak_gamma_pars[[1]], scale=viral_peak_gamma_pars[[2]],log=TRUE))
      #wane_gamma_pars <- gamma_pars_from_mean_sd(pars["wane_mean"],pars["wane_sd"]^2)
      #wane_random_effects <- sum(dgamma(wane_pars, wane_gamma_pars[[1]], scale=wane_gamma_pars[[2]],log=TRUE))
      #tp_gamma_pars <- gamma_pars_from_mean_sd(pars["tp_mean"],pars["tp_sd"]^2)
      #tp_random_effects <- sum(dgamma(tp_pars, tp_gamma_pars[[1]], scale=tp_gamma_pars[[2]],log=TRUE))
      
      ## Normally distributed parameters
      viral_peaks_random_effects <- sum(dnorm(viral_peaks, pars["viral_peak_mean"], pars["viral_peak_sd"],log=TRUE))
      wane_random_effects <- sum(dnorm(wane_pars, pars["wane_mean"], pars["wane_sd"],log=TRUE))
      tp_random_effects <- sum(dnorm(tp_pars, pars["tp_mean"], pars["tp_sd"], log=TRUE))
      
      random_effects <- viral_peaks_random_effects + wane_random_effects + tp_random_effects
      lik <- sum(lik) + random_effects
     
      if(!is.null(PRIOR_FUNC)){
        lik <- lik + PRIOR_FUNC(pars)
      }
      lik
      
    }
  }
  
}

## Generate random starting parameter table
generate_start_tab <- function (par_tab) 
{
  for (i in 1:nrow(par_tab)) {
    if (par_tab[i, "fixed"] == 0) {
      par_tab[i, "values"] <- runif(1, par_tab[i, "lower_bound"], 
                                    par_tab[i, "upper_bound"])
    }
  }
  return(par_tab)
}

get_index_par <- function (chain, index) 
{
  par_names <- colnames(chain)[2:(ncol(chain) - 1)]
  par <- as.numeric(chain[chain$sampno == index, 2:(ncol(chain) - 
                                                      1)])
  names(par) <- par_names
  return(par)
}
