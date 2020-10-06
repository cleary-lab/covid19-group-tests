## Viral kinetics function
model_func <- function(ts, tp, viral_peak, wane){
  y <- numeric(length(ts))
  growth <- viral_peak/tp
  y[ts <= tp] <- growth*ts[ts <= tp]
  y[ts > tp] <- viral_peak + wane*tp - wane*(ts[ts > tp])
  y
}
## Viral kinetics function
model_func_tinf <- function(ts, tinf,tp, viral_peak, wane){
  y <- numeric(length(ts))
  growth <- viral_peak/tp
  y[ts <= (tinf)] <- 0
  y[ts > (tinf) & ts <= (tinf + tp)] <- growth*(ts[ts > (tinf) & ts <= (tinf +  tp)]-(tinf))
  y[ts > (tp + tinf)] <- viral_peak + wane*tp - wane*(ts[ts > (tp + tinf)] - (tinf))
  y
}

## Viral kinetics function
model_func_tinf2 <- function(ts, tsympt,tp,td, viral_peak, day_undetectable,y0,LOD=3){
  y <- numeric(length(ts))
  growth <- viral_peak/abs(td)
  wane_rate <- (viral_peak-LOD)/(day_undetectable-tp)
  
  viral_latent_period <- ts <= (tsympt + tp - td)
  viral_growth_period <- (ts > (tsympt + tp - td)) & (ts <= tsympt + tp)
  viral_wane_period <- ts >  tsympt + tp
  
  y[viral_latent_period] <- y0
  y[viral_growth_period] <- growth*(ts[viral_growth_period] - (tsympt + tp - td))
  y[viral_wane_period] <- viral_peak - wane_rate*(ts[viral_wane_period] - (tsympt + tp))
  y
}

## Viral kinetics function
model_func_tinf2_sim <- function(ts, tinf, tsympt,tp,td, viral_peak, day_undetectable,y0,LOD=3){
  y <- numeric(length(ts))
  growth <- viral_peak/abs(td)
  wane_rate <- (viral_peak-LOD)/(day_undetectable-tp)
  
  viral_latent_period <- ts <= (tinf + tsympt + tp - td)
  viral_growth_period <- (ts > (tinf + tsympt + tp - td)) & (ts <= tinf + tsympt + tp)
  viral_wane_period <- ts >  tinf + tsympt + tp
  
  y[viral_latent_period] <- y0
  y[viral_growth_period] <- growth*(ts[viral_growth_period] - (tinf + tsympt + tp - td))
  y[viral_wane_period] <- viral_peak - wane_rate*(ts[viral_wane_period] - (tinf + tsympt + tp))
  y
}



## Viral kinetics function - parameterised by how many days before symptom onset viral loads start increasing
## and then subsequent days to peak viral load
model_func_tonset <- function(ts, td, tp, viral_peak, wane, y0=0){
  y <- numeric(length(ts))
  growth <- viral_peak/tp 
  y[ts <= td] <- y0
  y[ts > td & ts <= (td + tp)] <- growth*(ts[ts > td & ts <= (td + tp)]-td) + y0
  y[ts > (tp + td)] <- viral_peak + wane*(tp + td) - wane*ts[ts > (tp + td)] + y0
  y
}

## Viral kinetics function - paramterised by how many days prior to symptom onset viral loads peak
## assume that we cannot observe the uptick phase
model_func_tonset_fit <- function(ts, tp, td=-4, viral_peak, wane, y0=0){
  y <- numeric(length(ts))
  growth <- viral_peak/abs(td)
  y[ts <= (td+tp)] <- y0
  y[ts > (td+tp) & ts < tp] <- growth*(ts[ts > (td+tp) & ts < tp]-(td+tp)) + y0
  y[ts >= tp] <- viral_peak - wane*(ts[ts >= tp]-tp) + y0
  y
}

## Viral kinetics function - paramterised by how many days prior to symptom onset viral loads peak
## assume that we cannot observe the uptick phase
model_func_tonset_fit2 <- function(ts, tp, viral_peak, wane, y0=0){
  y <- numeric(length(ts))
  y[ts < tp] <- y0
  y[ts >= tp] <- viral_peak - wane*(ts[ts >= tp]+tp) + y0
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

transform_to_ct <- function(vl, intercept){
  ct <- intercept - log2(10)*vl
  ct
}

## Calculate likelihood of observations given model predicted viral loads
## Assuming normally distributed measurement error, with censoring
## at upper and lower limtis of detection
likelihood_discrete <- function(obs, predicted, sd, min_ct=0, lod=40){
  #predicted <- transform_to_ct(predicted, lod)
  
  high_titres <- which(obs <= min_ct)
  normal_titres <- which(obs > min_ct & obs < lod)
  low_titres <- which(obs >= lod)
  
  liks <- numeric(length(obs))
  
  liks[high_titres] <- pnorm(min_ct, predicted[high_titres], sd, lower.tail=TRUE, TRUE)
  liks[normal_titres] <- log(pnorm(obs[normal_titres], predicted[normal_titres], sd, TRUE, FALSE) - pnorm(obs[normal_titres]-1, predicted[normal_titres], sd, TRUE, FALSE))
  liks[low_titres] <- pnorm(lod, predicted[low_titres], sd, lower.tail=FALSE, TRUE)
  
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
  par <- chain[chain$sampno == index, 2:(ncol(chain) - 1)]
  if(nrow(par) > 1){
    par <- par[sample(1:nrow(par),1),]
  }
  par <- as.numeric(par)
  names(par) <- par_names
  return(par)
}



