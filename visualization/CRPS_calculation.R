
shuffle_matrix_efficiently <- function(X) {
  n <- nrow(X)
  for (i in 1:(n-1)) {
    
    choices <- setdiff(1:n, i)
    j <- sample(choices, 1)
    
    X[c(i, j),] <- X[c(j, i),]
  }
  return(X)
}


reps = 50
n_Sex = 2
n_Age = 8
n_Time = 10
n_subsize = 100

# Case I (not depending in Y)
DEPEND_ON_Y = FALSE

# Case II (depending in Y)
DEPEND_ON_Y = TRUE

nCol = 2 * n_Age * n_Time
RW_CRPS = matrix(NA, nrow = reps, ncol = nCol)
Indep_CRPS = matrix(NA, nrow = reps, ncol = nCol)
Fixed_CRPS = matrix(NA, nrow = reps, ncol = nCol)
BL_CRPS = matrix(NA, nrow = reps, ncol = nCol)

for(i in 1:50){
  print(i)
  load(file = paste0("../simulation_fitted_model_itr_", i, "_age_",  n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))
  
  fitted_model_BL = models_fit$BL
  fitted_model_Fixed = models_fit$Fixed 
  fitted_model_Indep = models_fit$Indep
  fitted_model_RW = models_fit$RW
  sim_data = models_fit$sim_data 
  
  RW_sample_pis <- fitted_model_RW$poster_sample$sample_pis
  Indep_sample_pis <- fitted_model_Indep$poster_sample$sample_pis
  Fixed_sample_pis <- fitted_model_Fixed$poster_sample$sample_pis
  BL_sample_pis <- fitted_model_BL$poster_sample$sample_pis
  
  repeat {
    RW_sample_pis_prime <- shuffle_matrix_efficiently(RW_sample_pis)
    RW_is_not_same = all(sapply(1:nrow(RW_sample_pis), function(i) !identical(RW_sample_pis[i, ], RW_sample_pis_prime[i, ])))
    
    if (RW_is_not_same) {
      RW_CRPS[i, ] = - 0.5*colMeans(abs(RW_sample_pis_prime - RW_sample_pis)) + colMeans(abs(RW_sample_pis - sim_data$param$pi))
      break
    }
  }
  
  repeat {
    Indep_sample_pis_prime <- shuffle_matrix_efficiently(Indep_sample_pis)
    Indep_is_not_same = all(sapply(1:nrow(Indep_sample_pis), function(i) !identical(Indep_sample_pis[i, ], Indep_sample_pis_prime[i, ])))
    
    if (Indep_is_not_same) {
      Indep_CRPS[i, ] = - 0.5*colMeans(abs(Indep_sample_pis_prime - Indep_sample_pis)) + colMeans(abs(Indep_sample_pis - sim_data$param$pi))
      break
    }
  }
  
  repeat {
    Fixed_sample_pis_prime <- shuffle_matrix_efficiently(Fixed_sample_pis)
    Fixed_is_not_same = all(sapply(1:nrow(Fixed_sample_pis), function(i) !identical(Fixed_sample_pis[i, ], Fixed_sample_pis_prime[i, ])))
    
    if (Fixed_is_not_same) {
      Fixed_CRPS[i, ] = - 0.5*colMeans(abs(Fixed_sample_pis_prime - Fixed_sample_pis)) + colMeans(abs(Fixed_sample_pis - sim_data$param$pi))
      break
    }
  }
  
  repeat {
    BL_sample_pis_prime <- shuffle_matrix_efficiently(BL_sample_pis)
    BL_is_not_same = all(sapply(1:nrow(BL_sample_pis), function(i) !identical(BL_sample_pis[i, ], BL_sample_pis_prime[i, ])))
    
    if (BL_is_not_same) {
      BL_CRPS[i, ] = - 0.5*colMeans(abs(BL_sample_pis_prime - BL_sample_pis)) + colMeans(abs(BL_sample_pis - sim_data$param$pi))
      break
    }
  }
}
