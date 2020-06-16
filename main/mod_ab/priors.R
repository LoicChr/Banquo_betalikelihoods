density <- function(pars){
  names(pars)<- list_params
  if (traits_biotic == "none"){
    dbio <- 0
  }
  if (traits_biotic %in% c("Height", "SLA")){
    d1 <- dlnorm(pars["intercept"], meanlog = log(5), sdlog = 0.8, log = T)
    d2 <- dnorm(pars["mu"], 0, 1.1, log = T) 
    d3 <- dlnorm(pars["sigma"], meanlog = log(0.9), sdlog = 0.9, log = T)
    dbio <- d1+d2+d3
  }
  if (traits_biotic == "2tr"){
    d1 <- dlnorm(pars["intercept"], meanlog = log(2.8), sdlog = 0.8, log = T)
    d2 <- dnorm(pars["mu1"], 0, 1.1, log = T) 
    d3 <- dlnorm(pars["sigma1"], meanlog = log(0.9), sdlog = 0.9, log = T)
    d4 <- dnorm(pars["mu1"], 0, 1.1, log = T) 
    d5 <- dlnorm(pars["sigma1"], meanlog = log(0.9), sdlog = 0.9, log = T)
    dbio <- d1+d2+d3+d4+d5
  }
  # Loglik params
  if (llparam == "phi_param"){
    da <- dunif(pars["a1"],-6, 0.5, log = T)
    db <- dunif(pars["b1"],0,1, log = T)
    dphi <- dhcauchy(pars["phi"], 4, log = TRUE)
    dll <- da+db+dphi
  } else if (llparam == "beta_conv") {
    dpsi <- dbeta(pars["psi"],1,1, log = TRUE)
    da2 <- dunif(pars["a2"],-6, 2, log = T)
    db2 <- dunif(pars["b2"],0,5, log = T)
    dphi <- dhcauchy(pars["phi"], 4, log = TRUE)
    dll <- dpsi+da2+db2+dphi
  } else if (llparam == "both"){
    da1 <- dunif(pars["a1"],-6, 2, log = T)
    db1 <- dunif(pars["b1"],0,2.5, log = T)
    da2 <- dunif(pars["a2"],-6, 2, log = T)
    db2 <- dunif(pars["b2"],0,5, log = T)
    dphi <- dhcauchy(pars["phi"], 4, log = TRUE)
    dll <- da1+db1+da2+db2+dphi
  }
  else{
    dpsi <- dbeta(pars["psi"],1,1, log = TRUE)
    dphi <- dhcauchy(pars["phi"], 4, log = TRUE)
    dll <- dpsi + dphi
  }
  dll + dbio
}

sampler <- function(n=1){
  if (traits_biotic == "none"){
    dbio <- NULL
  }
  if (traits_biotic %in% c("Height", "SLA")){
    d1 <- rlnorm(n, meanlog = log(5), sdlog = 0.8)
    d2 <- rnorm(n, 0, 1.1) 
    d3 <- rlnorm(n, meanlog = log(0.9), sdlog = 0.9)
    dbio <- cbind(d1, d2, d3)
  }
  if (traits_biotic == "2tr"){
    d1 <- rlnorm(n, meanlog = log(2.8), sdlog = 0.8)
    d2 <- rnorm(n, 0, 1.1) 
    d3 <- rlnorm(n, meanlog = log(0.9), sdlog = 0.9)
    d4 <- rnorm(n, 0, 1.1) 
    d5 <- rlnorm(n, meanlog = log(0.9), sdlog = 0.9)
    dbio <- cbind(d1, d2, d3, d4, d5)
  }
  if (llparam == "phi_param"){
    da <- runif(n,-6, 0.5)
    db <- runif(n,0,1)
    dphi <- rhcauchy(n, 4)
    
    dll <- cbind(da, db, dphi)
  } else if (llparam == "beta_conv") {
    dpsi <- rbeta(n,1,1)
    da2 <- runif(n,-6, 2)
    db2 <- runif(n,0,5)
    dphi <- rhcauchy(n, 4)
    dll <- cbind(dpsi,da2,db2,dphi)
  } else if (llparam == "both"){
    da1 <- runif(n,-6, 2)
    db1 <- runif(n,0,5)
    da2 <- runif(n,-6, 2)
    db2 <- runif(n,0,5)
    dphi <- rhcauchy(n, 4)
    dll <- cbind(da1,db1,da2,db2,dphi)
  }else{
    dpsi <- rbeta(n,1,1)
    dphi <- rhcauchy(n, 2)
    
    dll <- cbind(dpsi, dphi)
  }
  cbind(dbio, dll)
}

if (traits_biotic %in% c("Height", "SLA")){
  list_params_bio <- c("intercept", "mu", "sigma")
  bounds_bio <- data.frame(lower = c(0.4,-2.7,0.2), upper = c(15,2.7,4), row.names = list_params_bio)
} else if (traits_biotic == "2tr"){
  list_params_bio <- c("intercept", "mu1", "sigma1", "mu2", "sigma2")
  bounds_bio <- data.frame(lower = c(0,-2.7,0.2,-2.7,0.2), upper = c(20,2.7,4,2.7,4), row.names = list_params_bio)
}else{
  list_params_bio <- NULL
  bounds_bio <- NULL
}
if (llparam == "phi_param"){
  list_params_ll <- c("a1","b1", "phi")
  bounds_ll <- data.frame(lower = c(-6,0,0), upper = c(0.5,2.5,25), row.names = list_params_ll)
} else if (llparam == "both"){
  list_params_ll <- c("a1","b1","a2","b2", "phi")
  bounds_ll <- data.frame(lower = c(-6,0,-6,0,0), upper = c(5,2,5,2, 25), row.names = list_params_ll)
} else if (llparam == "beta_conv"){
  list_params_ll <- c("psi", "a2","b2", "phi")
  bounds_ll <- data.frame(lower = c(0,-6,0,0), upper = c(1,2,5,25), row.names = list_params_ll)
} else{
  list_params_ll <- c("psi","phi")
  bounds_ll <- data.frame(lower = c(0,0), upper = c(1,100), row.names = list_params_ll)
}
list_params <- c(list_params_bio, list_params_ll)
bounds <- rbind(bounds_bio, bounds_ll)
rm(list = c("bounds_bio", "bounds_ll", "list_params_ll", "list_params_bio"))