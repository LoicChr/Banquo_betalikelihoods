bounds
N = 3000
dat <- data.frame(intercept = runif(N,0,20), mu = runif(N,-3,3), sigma = runif(N, 0,6))

cortoId <- apply(dat, 1, function(pars){
  alphas = interaction.matrix(scale(trait), intercept = pars[1], mu = pars[2], sigma = pars[3], intra = 1, std = T)
  cor(as.numeric(alphas), as.numeric(diag(nrow(alphas))))
})
plot(dat[,1], cortoId)

prior <- apply(dat, 1, function(x){
  d1 <- dlnorm(x[1], meanlog = log(5), sdlog = 0.8, log = T)
  d2 <- dnorm(x[2], 0, 1.417593, log = T) 
  d3 <- dlnorm(x[3], meanlog = log(1.3), sdlog = 0.5, log = T)
  dbio <- d1+d2+d3
})
