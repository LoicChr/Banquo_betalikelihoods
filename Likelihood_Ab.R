likelihoodCov <- function(obs, pred, phi, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL, psi = NULL){
  # This likelihood function actually calculates the likelihood of the mean covers according to a dirichlet.
  # It supports various likelihood parametrization option. 
  prob <- matrix(NA, nrow = nrow(pred), ncol= ncol(pred))
  x <- obs > 0
  if (llparam %in% c("beta_conv", "none")){
    psi = psi
  }else{
    psi = plogis(pred, scale = exp(a1), location =b1)
  }
  
  prob_pres <- dbinom(obs > 0, 1, psi, log = T)
  if (llparam %in% c("phi_param", "none")){
    mean_cover <- pred[x]
  }else{
    mean_cover <- plogis(pred, scale = exp(a2), location = b2)
  }
  alpha = mean_cover * phi #Conversion in alpha/beta shape parameters
  beta = (1-mean_cover) * phi
  j <- match(obs[x], cover_class$cover_class)
  
  prob_cat <- log(pbeta(cover_class$maxCov[j], alpha, beta) - pbeta(cover_class$minCov[j], alpha, beta))
  prob_ab <- prob_pres
  prob_ab[x] <- prob_ab[x]+prob_cat
  
  sum(prob_ab)
}


likelihoodAb <-  function(pars){
  names(pars) <- list_params
  # This likelihood function is used to choose the option to calculate species abundance across sites
  # Either just traitspace or banquo with different traits.

  if (traits_biotic == "none"){
      pred <- P_S_E_tr
    }
  if (traits_biotic %in% c("Height", "SLA")){
      P_S_E_interactions <- banquo(P_S_E_tr, as.data.frame(trait), avg.out = T, intercept = pars.tr["intercept"],
                                   mu = pars["mu"], sigma = pars["sigma"], intra = 1)
      pred <- as.matrix(P_S_E_interactions)
  }
  
  if (traits_biotic == "2trait"){
    P_S_E_interactions <- banquo(P_S_E_tr, tr= cbind(trait1, trait2), avg.out = T,
                                 intercept = pars["intercept"],
                                 mu = c(pars["mu1"], pars["mu2"]),
                                 sigma = c(pars["sigma1"],pars["sigma2"]),
                                 rho = 0,
                                 intra =1)
    pred <- as.matrix(P_S_E_interactions)
  }
  pred <- pred +  1e-4
  #Computation of the likelihood
  if ( llparam == "phi_param"){
    loglik <- likelihoodCov(comm.red, pred, phi = pars["phi"], a1 = pars["a1"], b1 = pars["b1"])
  }
  if ( llparam == "beta_conv"){
    loglik <- likelihoodCov(comm.red, pred, phi = pars["phi"], a2 = pars["a2"], b2 = pars["b2"], psi = pars["psi"])
  }
  if ( llparam == "both"){
    loglik <- likelihoodCov(comm.red, pred, phi = pars["phi"], a1 = pars["a1"], b1 = pars["b1"], a2 = pars["a2"], b2 = pars["b2"])
  }
  if ( llparam == "none"){
    loglik <- likelihoodCov(comm.red, pred, phi = pars["phi"], psi = pars["psi"])
  }
    return(loglik)
}