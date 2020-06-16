rm(list=ls())

#Librairies
library(plyr)
library(ade4)
library(extraDistr)
library(BayesianTools)

#Load functions
source("lib/interactionMatrix.R")
source("lib/traitspace_and_banquo.R")
source("lib/Likelihood_Ab.R")

# Data
# Beginning preparation ##############
source("main/mod_ab/data_prep.R")

traits_biotic <-  c("none", "Height", "SLA", "2tr")[1] 
llparam <- c("none","phi_param" ,"beta_conv","both")[1]
source("main/mod_ab/priors.R")

if (traits_biotic == 'Height'){
  trait <- t.avg$maxht
}
if (traits_biotic == 'SLA'){
  trait <- t.avg$sla
}
if (traits_biotic == "2tr"){
  trait1 <- t.avg$Axis1
  trait2 <- t.avg$Axis2
}

prior <- createPrior(density = density, sampler = sampler, lower = bounds[,1], upper = bounds[,2])
bayesianSetup <- createBayesianSetup(likelihoodAb, prior, names = list_params)
# # settings for the sampler
settings <- list(iterations = 15000*3, nrChains = 1)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)
save(list = ls(), file = paste0(result_file, "/obs_chain", id, ".Rdata"))

## llparam : none ; -1901 vs. -1826.98
## llparam : phi_param ; -1864.532 vs. -1829.378
## llparam : beta_conv ; -1766.469  vs. -1765.748
## llparam : both ; -1728 vs. -1767.189

