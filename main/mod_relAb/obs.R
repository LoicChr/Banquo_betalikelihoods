rm(list=ls())

#Librairies
library(plyr)
library(ade4)
library(extraDistr)
library(BayesianTools)

#Load functions
source("lib/interactionMatrix.R")
source("lib/traitspace_and_banquo.R")
source("lib/Likelihood_relAb.R")

# Data
# Beginning preparation ##############
source("main/mod_relAb/data_prep.R")

llparam <-"beta_conv"
traits_biotic <-  "none"

source("main/mod_relAb/priors.R")

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
bayesianSetup <- createBayesianSetup(likelihoodRelAb, prior, names = list_params)
# # settings for the sampler
settings <- list(iterations = 15000*3, nrChains = 1)
# P_S_E_tr <- matrix(1/15, ncol = 15, nrow = 93)

out2 <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

save(list = ls(), file = paste0(result_file, "/obs_chain", id, ".Rdata"))
