#################################################################################
#####    Fit Model with Structured Prior with Aggregation to Time Level   ######
#################################################################################
library(BayesLogit)
library(gtools)
library(caret)

# bias 
#' @param estimate estimated data
#' @param true true data
get_bias = function(estimate, truth) {
  colMeans(estimate) - truth
}


#' @param sim_data simulated data set
#' @param K number of latent classes
#' @param E number of iterations
#' @param BURN_IN number of burn in
#' @param STRUCTURED_PRIOR Fixed, Indep, RW or BL as non-structured prior
fit_model_aggToTime <- function(sim_data, K = K, E = E, BURN_IN = BURN_IN, STRUCTURED_PRIOR = "Indep"){
  ##  ---------------  ##
  ##    observations
  ##  ---------------  ##
  # Ti is the time index (rename the variable as D)
  D = sim_data$data$Ti
  # Y.t is the indicator of cause-of-death, {0, 1, NA}
  Y.t = sim_data$data$Y.t
  # Y.t_true is the true cause-of-death
  Y.t_true = sim_data$data.truth$Y.t
  # X is the input data, n by q, {0, 1, NA}
  X = sim_data$data$X
  # Number of sex sub-population
  n_Sex = sim_data$n_Sex
  # Number of time sub-population
  n_Data = sim_data$n_Time
  # Number of age sub-population
  n_Age = sim_data$n_Age
  # Number of observation for each time
  n = sim_data$nt

  
  
  # Create dummy variables
   dummies <- dummyVars("~ .", data = data.frame(A =factor(sim_data$data$A)))
  # Transform data to get dummy variables
   encoded_Age <- predict(dummies, newdata = data.frame(A = factor(sim_data$data$A)))
   encoded_Age <- encoded_Age[,-1]
  
  X = cbind(X, sim_data$data$S, encoded_Age)
  
  # initialize H
  # H is the categorical latent class
  H = sample(seq(K), length(D), replace=TRUE, prob=rdirichlet(1, rep(1, K)))
  
  
  q = dim(X)[2]
  
  # parameters
  a_phi = 1; b_phi = 1
  a_omega = 8; b_omega = 2 
  
  
  if (STRUCTURED_PRIOR == "RW"){
    a_sigma = 0.5; b_sigma = 0.0009 # prior in Polson et al 2013
    sigma2_pi_init = 100
    sigma2_pi = 0.04
    mu_pi = 0 # non-informative
    ms = rep(NA, n_Data)
    ms[1] = rnorm(1, mu_pi, sd = sqrt(sigma2_pi_init))
    for (nData in 2:n_Data) {
      ms[nData] = rnorm(1, ms[nData-1], sd = sqrt(sigma2_pi))
    }
  }else if(STRUCTURED_PRIOR == "Indep"){
    a_sigma = 0.5; b_sigma = 0.0015
    sigma2_pi = 0.04
    mu_pi = 0
    ms = rnorm(n_Data, mu_pi, sd = sqrt(sigma2_pi))
  }else if(STRUCTURED_PRIOR == "Fixed"){
    sigma2_pi = 100
    mu_pi = 0 # non-informative
    ms = rnorm(n_Data, mu_pi, sd = sqrt(sigma2_pi))
  }else if(STRUCTURED_PRIOR == "BL"){
    a_pi = 1
    b_pi = 1
    pis = rbeta(n_Data, 1, 1) 
  }

  
  # initialize phi
  phi_0 = array(rbeta(n_Data*q*K, a_phi, b_phi), dim = c(n_Data, q, K))
  phi_1 = array(rbeta(n_Data*q*K, a_phi, b_phi), dim = c(n_Data, q, K))
  
  # initialize omega
  omega_0 = rgamma(n_Data, a_omega, b_omega)
  omega_1 = rgamma(n_Data, a_omega, b_omega)
  
  V_0 = matrix(rbeta(n_Data*K, 1, omega_0), nrow = n_Data)
  V_1 = matrix(rbeta(n_Data*K, 1, omega_1), nrow = n_Data)
  V_0[, K] = V_1[, K] = 1
  
  # initialize lambda
  lambda_0 = matrix(rep(NA, n_Data*K), nrow = n_Data)
  lambda_1 = matrix(rep(NA, n_Data*K), nrow = n_Data)
  
  for(nData in 1:n_Data){
    lambda_0[nData, 1] = V_0[nData, 1]
    lambda_1[nData, 1] = V_1[nData, 1]
    for (k in 2:(K-1)) {
      V0_prod = V1_prod = 1
      for (j in 1:(k-1)) {
        V0_prod = V0_prod*(1 - V_0[nData, j])
        V1_prod = V1_prod*(1 - V_1[nData, j])
      }
      lambda_0[nData, k] = V_0[nData, k]*V0_prod
      lambda_1[nData, k] = V_1[nData, k]*V1_prod
    }
    lambda_0[nData, K] = 1 - sum(lambda_0[nData, 1:(K-1)])
    lambda_1[nData, K] = 1 - sum(lambda_1[nData, 1:(K-1)])
  }
  lambda_0_prior = lambda_0
  lambda_1_prior = lambda_1
  

  omega_samp = rep(NA, n_Data)
  for (nData in 1:n_Data) {
    omega_samp[nData] = rpg(1, h = n[nData], z = 0)
  }
  
  
  # store the samples
  sample_pis = matrix(0, nrow = E-BURN_IN, ncol = n_Data)
  sample_omega = matrix(0, nrow = E-BURN_IN, ncol = n_Data)
  sample_lambda_0 = sample_lambda_1 = array(NA, dim = c(E-BURN_IN, n_Data, K))
  sample_omega_0 = sample_omega_1 = array(NA, dim = c(E-BURN_IN, n_Data))
  sample_phi_0 = sample_phi_1 = array(NA, dim = c(E-BURN_IN, n_Data, q, K))
  sample_sigma2_pi = sample_sigma2_pi_init = rep(0, E-BURN_IN)


  require("Rcpp")
  require("RcppArmadillo")
  
  if (STRUCTURED_PRIOR == "RW"){
    out <- fit_model_RW_aggToTime_internal(E = E, BURN_IN = BURN_IN, 
                                                    n_Data = n_Data,  a_sigma= a_sigma, b_sigma = b_sigma, 
                                                    a_phi = a_phi, b_phi = b_phi, a_omega = a_omega, b_omega = b_omega,
                                                    ms_r = ms, 
                                                    phi_0_r = phi_0, phi_1_r = phi_1, lambda_0_r = lambda_0, lambda_1_r = lambda_1,
                                                    V_0_r = V_0, V_1_r = V_1, omega_0_r = omega_0, omega_1_r = omega_1,
                                                    X_r = as.matrix(X), D_r = D - 1, Yt_r = Y.t, H_r = H - 1, K_r = K, mu_pi_r = mu_pi, sigma2_pi_r = sigma2_pi, sigma2_pi_init_r = sigma2_pi_init, sub_size_n_r = n)
  }else if (STRUCTURED_PRIOR == "Indep"){
    out <- fit_model_Indep_aggToTime_internal(E = E, BURN_IN = BURN_IN, 
                                                       n_Data = n_Data,  a_sigma= a_sigma, b_sigma = b_sigma, 
                                                       a_phi = a_phi, b_phi = b_phi, a_omega = a_omega, b_omega = b_omega,
                                                       ms_r = ms, 
                                                       phi_0_r = phi_0, phi_1_r = phi_1, lambda_0_r = lambda_0, lambda_1_r = lambda_1,
                                                       V_0_r = V_0, V_1_r = V_1, omega_0_r = omega_0, omega_1_r = omega_1,
                                                       X_r = as.matrix(X), D_r = D - 1, Yt_r = Y.t, H_r = H - 1, K_r = K, sigma2_pi_r = sigma2_pi, sub_size_n_r = n)
  }else if (STRUCTURED_PRIOR == "Fixed"){
    out <- fit_model_Fixed_aggToTime_internal(E = E, BURN_IN = BURN_IN, 
                                                       n_Data = n_Data, 
                                                       a_phi = a_phi, b_phi = b_phi, a_omega = a_omega, b_omega = b_omega,
                                                       ms_r = ms, 
                                                       phi_0_r = phi_0, phi_1_r = phi_1, lambda_0_r = lambda_0, lambda_1_r = lambda_1,
                                                       V_0_r = V_0, V_1_r = V_1, omega_0_r = omega_0, omega_1_r = omega_1,
                                                       X_r = as.matrix(X), D_r = D - 1, Yt_r = Y.t, H_r = H - 1, K_r = K, sigma2_pi_r = sigma2_pi, sub_size_n_r = n)
  }else if (STRUCTURED_PRIOR == "BL"){
    out <- fit_model_BL_aggToTime_internal(E = E, BURN_IN = BURN_IN, 
                                                    n_Data = n_Data, 
                                                    a_phi = a_phi, b_phi = b_phi, a_omega = a_omega, b_omega = b_omega,
                                                    a_pi = a_pi, b_pi = b_pi,
                                                    pis_r = pis, 
                                                    phi_0_r = phi_0, phi_1_r = phi_1, lambda_0_r = lambda_0, lambda_1_r = lambda_1,
                                                    V_0_r = V_0, V_1_r = V_1, omega_0_r = omega_0, omega_1_r = omega_1,
                                                    X_r = as.matrix(X), D_r = D - 1, Yt_r = Y.t, H_r = H - 1, K_r = K, sub_size_n_r = n)
  }

  sample_pis = out$pis
  sample_lambda_1 = out$lambda_1
  sample_lambda_0 = out$lambda_0

  
  bias = get_bias(sample_pis, sim_data$param$pi)

  
  poster_sample = list(sample_pis = sample_pis, sample_lambda_0 = sample_lambda_0, sample_lambda_1 = sample_lambda_1)
  return(list(bias = bias, fit = out, poster_sample = poster_sample))
}



