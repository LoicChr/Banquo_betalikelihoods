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
# P_S_E_tr <- matrix(1/15, ncol = 15, nrow = 93) # To try H0

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)


