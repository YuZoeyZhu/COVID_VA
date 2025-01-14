library(Rcpp)

# Load models and data generation function
sourceCpp("../model/Rcpp/fit_model_BL.cpp")
sourceCpp("../model/Rcpp/fit_model_with_structured_prior.cpp")

source("../model/fit_model_BL.R")
source("../model/fit_model_with_structured_prior.R")

source("../simulation/sim_with_sex_time_age.R")

# parameters set up
n_Sex = 2
n_Age = 8
n_Time = 10

n_subsize = 100 # sample size of each sub-population
q = 10
K = 10
DEPEND_ON_Y = FALSE # or TRUE, corresponding to case i and ii of whether sampling with dependency on Y

# generate pi_true as prevalence
p_S = c(0.2, 0.3)

x <- 1:n_Time
p_T <- (-6 + 10*x - x^2) / n_Time / 2

x <- 1:n_Age
p_A = (-6 + 10*x - x^2) / n_Age / 2

Nsta = n_Sex * n_Time * n_Age
pi_true = rep(NA, Nsta)
S_ind = c(rep(0, n_Time * n_Age), rep(1, n_Time*n_Age))
T_ind =  rep(rep(1:n_Time, each = n_Age), n_Sex)
A_ind =  rep(rep(1:n_Age,  n_Time), n_Sex)
for(i in 1:Nsta){
  pi_true[i] = round(expit(-1 + p_S[S_ind[i]+1] + p_T[T_ind[i]] + p_A[A_ind[i]]), 2) 
}

# sampling iterations
E = 2000
BURN_IN = 1000

for (i in 1:50){
  print(paste("it's iter", i))
  
  # generate simulated dataset
  sim_data = sim_with_sex_time_age(pi_true, n_subsize = n_subsize, q = q, n_Sex = n_Sex, n_Time = n_Time, n_Age = n_Age, K = K, DEPEND_ON_Y = DEPEND_ON_Y)
  
  
  # fit models
  fitted_model_BL = fit_model_BL(sim_data, E = E, BURN_IN = BURN_IN)
  fitted_model_Fixed = fit_model_with_structured_prior(sim_data, K = K, E = E, BURN_IN = BURN_IN, STRUCTURED_PRIOR = "Fixed")
  fitted_model_Indep = fit_model_with_structured_prior(sim_data, K = K, E = E, BURN_IN = BURN_IN, STRUCTURED_PRIOR = "Indep")
  fitted_model_RW = fit_model_with_structured_prior(sim_data, K = K, E = E, BURN_IN = BURN_IN, STRUCTURED_PRIOR = "RW")
  
  
  # save model results
  models_fit = list()
  models_fit$BL = fitted_model_BL
  models_fit$Fixed = fitted_model_Fixed
  models_fit$Indep = fitted_model_Indep
  models_fit$RW = fitted_model_RW
  models_fit$sim_data = sim_data
  
  save(models_fit, file = paste0("../simulation_fitted_model_itr_", i, "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))
  
  metrics <- list(
    lambda_bias_model_BL = models_fit$BL$bias,
    lambda_bias_model_Fixed = models_fit$Fixed$bias,
    lambda_bias_model_Indep = models_fit$Indep$bias,
    lambda_bias_model_RW = models_fit$RW$bias
  )
  
  save(metrics, file = paste0("../metric/simulation_metrics_itr_", i, "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))
  
}


