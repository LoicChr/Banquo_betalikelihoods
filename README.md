# Banquo_betalikelihoods

The script allows to play with different options. 

The main scripts are : (1) main/mod_ab/obs.R for the abundance/cover model, and (2) main/mod_relAb/obs.R for the abundance/cover model.
In their respective folder, there is a data preparation script data_prep.R and a script to set the priors and parameter bounds (priors.R). Both are called by their respective obs.R script.

The likelihood functions are in Likelihood_Ab.R (abundance model) and Likelihood_relAb.R (relative abundance model). These two scripts contains first the core likelihood function itself and a function wrapped around it that can call Banquo to modify species covers.

Each model comes with a serie of options to calculate the mean cover: 
(1) "none" - Traitspace model (no modif of the traitspace object stored in data).
(2) "Height" - Banquo model using height
(3) "SLA" - Banquo model using SLA
(4) "2tr" - Banquo model using height and SLA

A serie of option to calculate the likelihood
Abundance/cover model : at its core, it's a zero inflated beta likelihood with a phi "precision parameter" for the beta distribution. The option are the following
(1) llparam = "none": the probability of presence is equal for every species in every sites. And the mean cover used in the beta are the species mean cover predicted by Traitspace/Banquo
(2) llparam = "phi_param": very badly named, it actually refers to the Bernouilli parameter psi for the zeros, not the phi parameter ('Sorry'). In that case, the psi parameter is an expit function to parametrize of species mean cover predicted by Traitspace/Banquo
(3) llparam = "beta_conv": Mean cover abundances in the beta are an expit function to parametrize of species mean cover predicted by Traitspace/Banquo
(4) llparam = "both": Both the mean cover and the probability of presence/absence are expit functions to parametrize of species mean cover predicted by Traitspace/Banquo

Relative abundance model : at its core, it's a Dirichlet with a phi "precision parameter". Zeros and Ones are removed by adding an offset in both the observed community matrix and the predicted community matrix.
The option are the following
(1) llparam = "none": the probability of presence is equal for every species in every sites. And the mean cover used in the beta are the species mean cover predicted by Traitspace/Banquo
(2) llparam = "beta_conv": Mean cover abundances in the beta are an expit function to parametrize of species mean cover predicted by Traitspace/Banquo
