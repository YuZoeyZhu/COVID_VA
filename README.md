# COVID_VA

This repository contains the reproducible materials for the paper _Hierarchical Latent Class Models for Mortality Surveillance Using Partially Verified Verbal Autopsies_.


## Overview of the structure and code scripts 

+ model/: R functions for fitting the Bayesian latent variable models with structured priors 'RW', 'Indep' and 'Fixed', and the Baseline model without structured prior (BL). Rcpp dependencies are included under the Rcpp folder.

  * `fit_model_with_structured_prior.R` and `fit_model_BL.R`: fit models without aggregations.
  * `fit_model_aggToTime.R`: fit models with aggregation to the time level.
  * `fit_model_PopLevel.R`: fit model with aggregation to the population level.


+ simulation/: Simulation study with two types of verification mechanisms. In case (i), the verification process depends on observed symptoms and covariates $X$ and $D$. In case (ii), the verification process also depends on the unobserved cause of death $Y$.

  * `sim_with_sex_time_age.R`: includes the function of obtaining the simulated data for different verification mechanisms. 
  * `analysis_no_aggregate_sim_data.R`: for fitting the models without aggregations.
  * `analysis_aggregate_sim_data.R`: for fitting the models with aggregations to the time level and the population level.

  
  
+ visualization/: Includes R codes for figures and analysis results.
  * `CRPS_calculation.R`: for calculating the CRPS.
  * `figures_script.R`: including the figure visualization scripts of aggregation results comparisons, prevalence posterior mean and CI analysis, and simulated true prevalence visualization.


+ Reproducible_examples/: R Markdown for a reproducible simulation study example.

