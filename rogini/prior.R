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
  # Likelihood parameters
  if (llparam == "beta_conv"){
    da <- dunif(pars["a"],-6, 2, log = T)
    db <- dunif(pars["b"],-2,5, log = T)
    dphi <- dhcauchy(pars["phi"], 4, log = TRUE)
    dll <- da+db+dphi
  }else if (llparam == "psi"){
    dpsi <- dunif(pars["psi"],0,1)
    dphi <- dhcauchy(pars["phi"], 4, log = TRUE)
    dll <- dpsi+dphi
  }else{
    dphi <- dhcauchy(pars["phi"], 4, log = TRUE)
    dll <- dphi
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
  if (llparam == "beta_conv"){
    da <- runif(n,-6, 2)
    db <- runif(n,-2,5)
    dphi <- rhcauchy(n, 4)
    dll <- cbind(da,db,dphi)
  }else if (llparam == "psi"){
    dpsi <- runif(n,0,1)
    dphi <- rhcauchy(n, 4)
    dll <- cbind(dpsi, dphi)
  }else{
    dphi <- rhcauchy(n, 4)
    dll <- cbind(dphi)
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
if (llparam == "beta_conv"){
  list_params_ll <- c("a","b", "phi")
  bounds_ll <- data.frame(lower = c(-6,-2,0), upper = c(0.5,1,25), row.names = list_params_ll)
}else if(llparam == "psi"){
  list_params_ll <- c("psi", "phi")
  bounds_ll <- data.frame(lower = c(0,0), upper = c(1,25), row.names = list_params_ll)
  
}else{
  list_params_ll <- c("phi")
  bounds_ll <- data.frame(lower = c(0), upper = c(100), row.names = list_params_ll)
}
list_params <- c(list_params_bio, list_params_ll)
bounds <- rbind(bounds_bio, bounds_ll)
rm(list = c("bounds_bio", "bounds_ll", "list_params_ll", "list_params_bio"))
