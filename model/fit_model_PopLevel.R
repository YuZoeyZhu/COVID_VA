################################################
#####    Fit Model Population Level   ######
################################################
expit <- function(x){exp(x)/(1 + exp(x))}
logit <- function(x){log(x/(1-x))}
library(gtools)
# bias 
#' @param estimate estimated data
#' @param true true data
get_bias = function(estimate, truth) {
  colMeans(estimate) - truth
}


# accuracy of individual-level classification
#' @param unobs_true true unobserved data
#' @param unobs_pred estimated unobserved data
get_acc <- function(unobs_true, unobs_pred) {
  n_unobs = length(unobs_true)
  n_pred_correct = n_unobs - sum(abs(unobs_pred  - unobs_true))
  acc = n_pred_correct / n_unobs
  return(acc)
}



#' @param sim_data simulated data set
#' @param E number of iterations
#' @param BURN_IN number of burn in
fit_model_PopLevel <- function(sim_data, K = K, E = E, BURN_IN = BURN_IN){
  require(BayesLogit)
  library(gtools)
  library(caret)
  ##  ---------------  ##
  ##    observations
  ##  ---------------  ##
  # Y.t is the indicator of cause-of-death, {0, 1, NA}
  Y.t = sim_data$data$Y.t
  # Y.t_true is the true cause-of-death
  Y.t_true = sim_data$data.truth$Y.t
  # X is the input data, n by q, {0, 1, NA}
  X = sim_data$data$X
  
  # Create dummy variables
  dummies <- dummyVars("~ .", data = data.frame(A =factor(sim_data$data$A)))
  # Transform data to get dummy variables
  encoded_Age <- predict(dummies, newdata = data.frame(A = factor(sim_data$data$A)))
  encoded_Age <- encoded_Age[,-1]
  
  X = cbind(X, sim_data$data$S, encoded_Age)

  
  # initialize H
  # H is the categorical latent class
  H = sample(seq(K), length(Y.t), replace=TRUE, prob=rdirichlet(1, rep(1, K)))
  
  Y.t_na_index = which(is.na(Y.t))
  

  q = dim(X)[2]
  
  # parameters set up
  a_phi = 1; b_phi = 1
  a_omega = 8; b_omega = 2 
  a_pi = 1; b_pi = 1
  
  
  phi_0 = array(rbeta(q*K, a_phi, b_phi), dim = c(q, K))
  phi_1 = array(rbeta(q*K, a_phi, b_phi), dim = c(q, K))
  
  
  omega_0 = rgamma(1, a_omega, b_omega)
  omega_1 = rgamma(1, a_omega, b_omega)
  
  V_0 = rbeta(K, 1, omega_0)
  V_1 = rbeta(K, 1, omega_1)
  V_0[K] = V_1[K] = 1
  
  lambda_0 = rep(NA, K)
  lambda_1 = rep(NA, K)
  
  # initialize lambda
    lambda_0[1] = V_0[1]
    lambda_1[1] = V_1[1]
    for (k in 2:(K-1)) {
      V0_prod = V1_prod = 1
      for (j in 1:(k-1)) {
        V0_prod = V0_prod*(1 - V_0[j])
        V1_prod = V1_prod*(1 - V_1[j])
      }
      lambda_0[ k] = V_0[k]*V0_prod
      lambda_1[ k] = V_1[ k]*V1_prod
    }
    lambda_0[ K] = 1 - sum(lambda_0[ 1:(K-1)])
    lambda_1[ K] = 1 - sum(lambda_1[ 1:(K-1)])

  lambda_0_prior = lambda_0
  lambda_1_prior = lambda_1
  
  
  pi = rbeta(1, 1, 1) 
  
  
  # Store the samples
  sample_pi = rep(NA, E-BURN_IN)
  sample_Y.t = matrix(0, nrow = E-BURN_IN, ncol = length(Y.t))

  

  require("Rcpp")
  require("RcppArmadillo")
  out <- fit_model_PopLevel_internal(E = E, BURN_IN = BURN_IN, 
                                                    a_phi = a_phi, b_phi = b_phi, a_omega = a_omega, b_omega = b_omega,
                                                    a_pi = a_pi, b_pi = b_pi,
                                                    pi = pi, 
                                                    phi_0_r = phi_0, phi_1_r = phi_1, lambda_0_r = lambda_0, lambda_1_r = lambda_1,
                                                    V_0_r = V_0, V_1_r = V_1, omega_0 = omega_0, omega_1 = omega_1,
                                                    X_r = as.matrix(X), Yt_r = Y.t, H_r = H - 1, K_r = K)
  sample_pi = out$pi
  sample_Y.t= out$Yt
  
  
  # Bias of estimated lambda
  bias = get_bias(sample_pi, sim_data$param$pi)
  
  poster_sample = list(sample_pi = sample_pi)
  
  return(list(bias = bias, poster_sample = poster_sample))
}


