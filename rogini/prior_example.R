list_params <- c("mu", "sd")

density <- function(pars){
  names(pars)<- list_params
  d1 <- dnorm(pars["mu"], -2, 0.1, log = T)
  d2 <- dlnorm(pars["sd"], 2, 0.1, log = T)
  return(d1+d2)
}

sampler <- function(n){
  d1 <- rnorm(n, -2, 0.1)
  d2 <- rlnorm(n, 2, 0.1)
  cbind(d1, d2)
}

bounds <- data.frame(lower = c(-2.4,4), upper = c(-1.6,12), row.names = list_params)
