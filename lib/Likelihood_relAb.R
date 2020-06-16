likelihoodRelCov <- function(obs, pred, phi, a = NULL, b = NULL){
  # This likelihood function actually calculates the likelihood of the mean covers according to a dirichlet.
  if (llparam %in% c("beta_conv")){
    pred <- plogis(pred, scale = exp(a), location = b)
  }
  prob <- sapply(1:nrow(obs), function(j) ddirichlet(obs[j,],pred[j,]*phi,log = T))
  loglik <- sum(prob)
  return(loglik)
}
likelihoodRelAb <-  function(pars){
  # This likelihood function is used to choose the option to calculate species abundance across sites
  # Either just traitspace or banquo with different traits.
    names(pars) <- list_params
    if (traits_biotic == "none"){
      P_S_E_interactions <- as.matrix(sweep(P_S_E_tr, 1, rowSums(P_S_E_tr), "/"))
    }
    # Transformation of the parameters
    if (traits_biotic %in% c("Height", "SLA")){ 
      P_S_E_interactions <- banquo(P_S_E_tr, as.data.frame(trait), avg.out = T, intercept = pars["intercept"], 
                                   mu = pars["mu"], sigma = pars["sigma"], intra = 1, std = T)
    }
    if (traits_biotic %in% c("2tr")){ 
      P_S_E_interactions <- banquo(P_S_E_tr, tr= cbind(trait1, trait2), avg.out = T, 
                                   intercept = pars["intercept"], 
                                   mu = c(pars["mu1"], pars["mu2"]), 
                                   sigma = c(pars["sigma1"],pars["sigma2"]),
                                   rho = 0,
                                   intra =1, std = F)
    }
    # Add offset
    P_S_E_interactions <- as.matrix(sweep(P_S_E_interactions, 1, rowSums(P_S_E_interactions), "/"))
    P_S_E_interactions <- P_S_E_interactions  + 1e-4
    pred <- as.matrix(sweep(P_S_E_interactions, 1, rowSums(P_S_E_interactions), "/"))

    #Computation of the likelihood
    if (llparam == "beta_conv"){
      loglik <- likelihoodRelCov(comm.red, pred, phi = pars["phi"], a = pars["a"], b = pars["b"])
    }else{
      loglik <- likelihoodRelCov(comm.red, pred, phi = pars["phi"])
    }

    return(loglik)
}
