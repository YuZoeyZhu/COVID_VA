######### aggregation simulation analysis #########
library(Rcpp)

# Load models 
sourceCpp("../model/Rcpp/fit_model_BL_aggToTime.cpp")
sourceCpp("../model/Rcpp/fit_model_Indep_aggToTime.cpp")
sourceCpp("../model/Rcpp/fit_model_Fixed_aggToTime.cpp")
sourceCpp("../model/Rcpp/fit_model_RW_aggToTime.cpp")
sourceCpp("../model/Rcpp/fit_model_PopLevel.cpp")

source("../model/fit_model_PopLevel.R")
source("../model/fit_model_aggToTime.R")


# parameters set up
n_Age = 8
n_Time = 10
n_Sex = 2

DEPEND_ON_Y = FALSE # or TRUE, corresponding to case i and ii of whether sampling with dependency on Y
n_subsize = 100

E = 2000
BURN_IN = 1000


# save results
models_fit_noAge = list()
models_fit_noAge_bias = list()
agg_time_RW_noAge = agg_time_Indep_noAge = agg_time_Fixed_noAge = agg_time_BL_noAge = matrix(NA, nrow = 50, ncol = n_Data)
agg_overall_PopLevel = rep(NA, 50)

for(i in 1:50){
  # Load simulation data
  name1 <- paste0("../simulation_fitted_model_itr_", i, "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda")
  
  if(file.exists(name1)){
    print(i)
    load(name1)
    
    # Fit models with aggregation to time level
    RW_agg_time = fit_model_aggToTime(models_fit$sim_data, K = K, E = E, BURN_IN = BURN_IN, STRUCTURED_PRIOR = "RW")
    Indep_agg_time = fit_model_aggToTime(models_fit$sim_data, K = K, E = E, BURN_IN = BURN_IN, STRUCTURED_PRIOR = "Indep")
    Fixed_agg_time = fit_model_aggToTime(models_fit$sim_data, K = K, E = E, BURN_IN = BURN_IN, STRUCTURED_PRIOR = "Fixed")
    BL_agg_time = fit_model_aggToTime(models_fit$sim_data, K = K, E = E, BURN_IN = BURN_IN, STRUCTURED_PRIOR = "BL")
    
    # Fit model with aggregation to the population level
    PopLevel_agg_overall = fit_model_PopLevel(models_fit$sim_data, K = K, E = E, BURN_IN = BURN_IN)
    
    models_fit_noAge$RW = RW_agg_time
    models_fit_noAge$Indep = Indep_agg_time
    models_fit_noAge$Fixed = Fixed_agg_time
    models_fit_noAge$BL = BL_agg_time
    
    model_fit_PopLevel = PopLevel_agg_overall
    
    save(models_fit_noAge, file = paste0("../simulation_fitted_model_noAge_itr_", i, "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))
    save(model_fit_PopLevel, file = paste0("../simulation_fitted_model_PopLevel_itr_", i, "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))
  }
}
