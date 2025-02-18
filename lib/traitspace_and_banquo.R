###############################################################################
#                                                                             #
#   Banquo and Traitspace functions                                           #
#   From the article: Integrating traits, competition and abiotic filtering   #
#                 improves biodiversity predictions                           #
#   Authors: Loïc Chalmandrier*, Daniel B. Stouffer, Daniel C. Laughlin        #
#   *contact                                                                  #
###############################################################################


traitspace <- function(trait.model, env, PT.Sk, N = 100, avg.out = TRUE){
  ##### Arguments.
  # trait.model: linear model predicting the trait values as a function of environmental variables.
  # env: data.frame containing the environmental values of the sites under study.
  # PT.Sk: a list containing the distribution of species traits. Here Mclust objects.
  # N: the number of random draws used to the MonteCarlo integration.
  # avg.out: should the random draws be averaged out for the final output?
  
  require(mclust)
  require(mvtnorm)
  
  if (is.null(row.names(env))){
    row.names(env) <- paste("Site", 1:nrow(env), "-")
  }else{
    row.names(env) <- gsub("_", "-", row.names(env))
  }
  
  #Extraction of the trait-environment prediction
  pred_multi <- predict(trait.model, newdata=env, interval="prediction", level=0.95)
  row.names(pred_multi) <- row.names(env)
  cov_mat <- cov(as.matrix(trait.model$residuals))
  
  # Construction of the simulated data
  ## STEP 2a: Drawing samples from P(T/E)
  env.to.fit <- data.frame(env[gl(nrow(env), N),])
  dimnames(env.to.fit) <- list(paste(row.names(env)[gl(nrow(env), N)], rep(1:N, nrow(env)), sep = "_"), colnames(env))
  tr_mean_pred <- pred_multi[gl(nrow(env), N), ]
  row.names(tr_mean_pred) <- paste(row.names(pred_multi)[gl(nrow(pred_multi), N)], rep(1:N, nrow(pred_multi)), sep = "_")
  
  if (all(dim(cov_mat) == c(1,1))) tr_mean_pred <- tr_mean_pred[, 1, drop = F]
  
  tr_sample <- as.data.frame(t(apply(tr_mean_pred, 1, function(x) rmvnorm(1, x, cov_mat))))
  if (all(dim(cov_mat) == c(1,1))) tr_sample <- as.data.frame(t(tr_sample))
  
  ## computing(P(T/E))
  P_T_E <-do.call(rbind, lapply(1:nrow(tr_sample), function(j) dmvnorm(tr_sample[j,],tr_mean_pred[j,], cov_mat)))
  row.names(P_T_E) <- row.names(tr_sample)
  
  ## Step 2b: Computing the likelihood P(T/Sk) using Mclust done earlier
  P_T_S <- lapply(PT.Sk, function(pdf){
    mclust::dens(pdf$modelName,tr_sample,parameters=pdf$parameters)
  })
  
  P_T_S <- do.call(cbind, P_T_S)
  colnames(P_T_S) <- names(PT.Sk)
  
  ## Step 2c: Computing posterior P(Sk/T,E) using Bayes theorem
  P_S_T_E <- exp(sweep(log(P_T_S), 1, log(rowSums(P_T_S)), "-"))
  
  ## Step 2d: Posterior P(Sk/T) by integrating out T's (with log)
  P_S_E_all <- exp(sweep(log(P_S_T_E), 1, log(P_T_E), "+"))
  
  if (!avg.out){
    colnames(P_S_E_all) <- names(PT.Sk)
    row.names(P_S_E_all) <- paste(row.names(env)[gl(nrow(env), N)],rep(1:N, nrow(env)) , sep ="_")
    return(P_S_E_all )
  }else{
    sites <- sapply(strsplit(row.names(P_S_E_all), "_"), function(x) x[1])
    ### Step 2d.1: Monte Carlo integration across trait samples
    P_S_E_all.sp <- split(as.data.frame(P_S_E_all), sites)
    P_S_E_traitspace_unnorm <- do.call(rbind, lapply(P_S_E_all.sp, function(x) apply(x, 2, mean)))
    #P_S_E_traitspace <- sweep(P_S_E_traitspace_unnorm, 1, rowSums(P_S_E_traitspace_unnorm) ,"/" )
    return(P_S_E_traitspace_unnorm)
  }
}

banquo <- function(P_S_E, tr = NULL, intercept=1, mu=-1.5, sigma =0, intra = NULL, rho = 0, std = F, alphas = NULL, avg.out = TRUE, out_alpha = FALSE){
  require(corpcor)
  ### Banquo computes species abundances based on the traitspace object, the trait that controls competition 
  #### and the parameters
  ##### Arguments.
  # P_S_E: the output of the traitspace function
  # tr : the trait(s) used to compute the interaction matrix
  # intercept, mu, sigma, intra, rho. The parameters to generate the interaction matrix
  # avg.out:  should the random draws be averaged out for the final output?
  # out_alpha : should the interaction matrix be saved in the output?
  if (is.null(alphas)){
    if (length(mu) ==1 & length(sigma) == 1 & ncol(tr) == 1 & ncol(P_S_E) == nrow(tr)){
      alphas <- interaction.matrix(scale(tr[,1]), intercept = intercept, mu = mu, sigma = sigma, intra = intra, std = std)
      
    } else  if (length(mu) ==2 & length(sigma) == 2 & ncol(tr) == 2  & ncol(P_S_E) == nrow(tr)){
      alphas <- interaction.matrix_bi(trait1 = scale(tr[,1]), trait2 = scale(tr[,2]), 
                                      intercept = intercept, mu1 = mu[1], sigma1 = sigma[1],
                                      mu2 = mu[2], sigma2 = sigma[2], rho = rho,
                                      intra = intra, std = std)
    }else{
      stop("tr, mu, sigma, P_S_E do not match")
    }
  }
  #Compute the inverse of the full alpha matrix once since it's the same throughout
  alphas.inv <- corpcor::pseudoinverse(alphas)
  
  P_S_E_all_interactions <- P_S_E
  
  integration.biotic <- function(Ks){
    Ns <- alphas.inv %*% as.matrix(Ks)
    # at least one species has been competitively excluded
    excluded <- c()
    if(any(Ns < 0)){
      A <- alphas
      Ks_red <- Ks
      excluded <- c(excluded, which.min(Ns))
      Ks_red <- Ks[-excluded]
      A <- A[-excluded, -excluded]
      Ns <- Ks*0
      Ns[-excluded] <- corpcor::pseudoinverse(A) %*% as.matrix(Ks_red)
      
      # unfortunately we need to keep going
      while(any(Ns < 0)){
        excluded <- c(excluded, which.min(Ns))
        A <- alphas
        Ks_red <- Ks
        Ks_red <- Ks[-excluded]
        A <- as.matrix(A[-excluded, -excluded])
        Ns <- Ks*0
        Ns[-excluded] <- corpcor::pseudoinverse(A) %*% as.matrix(Ks_red)
      }
    }
    return(Ns)
  }
   
  # Determine abundances for each trait sample
  for(k in 1:nrow(P_S_E)){
    Ns <- integration.biotic(Ks = P_S_E[k,])
    P_S_E_all_interactions[k,] <- Ns
   }
  dimnames(P_S_E_all_interactions) <- dimnames(P_S_E)
  if (!avg.out){
    return(P_S_E_all_interactions)
  }else{
    ### Step 2d.1: Monte Carlo integration across trait samples
    fac <- sapply(strsplit(row.names(P_S_E_all_interactions), '_'), function(x) x[1])
    P_S_E_all.sp <- split(as.data.frame(P_S_E_all_interactions), as.factor(fac))
    P_S_E_interactions_unnorm <- do.call(rbind, lapply(P_S_E_all.sp, function(x) apply(x, 2, mean)))
   # P_S_E_interactions_unnorm  <- sweep(P_S_E_interactions_unnorm , 1, rowSums(P_S_E_interactions_unnorm ) ,"/" )
    
    #MC normalisation
    if (any(P_S_E_all_interactions < 0)){
      P_S_E_interactions_unnorm[P_S_E_all_interactions < 0] <- 0
      }
     if (out_alpha){
      out <- list(P_S_E_interactions_unnorm, alphas)
    }else{
      out <- P_S_E_interactions_unnorm
    }
    return(out)
  }
}


