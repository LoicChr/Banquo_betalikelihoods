#############
# This piece of script is designed to generate a prior distribution for a MCMC sampler. 
# To that end, it generates a sampler function, a density function and a boundaries dataframe for each of its parameters
# The basic script to do that is in prior_simple.R, it is the barebone version of prior.R without modeling options

# Modeling "options" specified by the variables "traits_biotic" and "llparam"
# According to those variables, the density, sampler and bounds are changed. 

# Ideally, I'd like to simplify the script in prior.R so that if I need to change the hyperparameters of the prior distributions, 
# I don't have to change it manually in the sampler part, in the density part and in the bounds part.



### Charge one of the options
traits_biotic = c("none", "Height", "SLA", "2tr")[1]
llparam = c("beta_conv", "psi", "none")[1]

source("prior.R")
# The source creates the function density, sampler and create the bounds dataframe

# Check the functioning of the functions
pars = sampler(1)
density(pars)
