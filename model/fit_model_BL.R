#################################################
#####   Fit Latent Class Model Baseline ######
#################################################
library(MASS)
library(mnormt)
library(BayesLogit)
library(gsignal)
library(gtools)
library(tidyr)


#' @param sim_data simulated data set
#' @param K number of latent classes
#' @param E number of iterations
#' @param BURN_IN number of burn in
fit_model_BL <- function(sim_data, K = 10, E = 10000, BURN_IN = 5000){
  ##  ---------------  ##
  ##    observations
  ##  ---------------  ##
  # S is the sex index
  S = sim_data$data$S
  # Ti is the time index
  Ti = sim_data$data$Ti
  # A is the age index
  A = sim_data$data$A
  # Y.t is the indicator of cause-of-death, {0, 1, NA}
  Y.t = sim_data$data$Y.t
  # Y.t_true is the true cause-of-death
  Y.t_true = sim_data$data.truth$Y.t
  # X is the input data, n by q, {0, 1, NA}
  X = sim_data$data$X
  # Number of sex sub-population
  n_Sex = sim_data$n_Sex
  # Number of time sub-population
  n_Time = sim_data$n_Time
  # Number of age sub-population
  n_Age = sim_data$n_Age
  # Number of observation for each sex, time and age
  nsta = sim_data$nsta

  
  # initialize H
  # H is the categorical latent class
  H = sample(seq(K), length(Ti), replace=TRUE, prob=rdirichlet(1, rep(1, K)))
  
  q = dim(X)[2]
  
  # parameters
  a_phi = 1; b_phi = 1
  a_omega = 8; b_omega = 2 
  a_pi = 1; b_pi = 1
  
  
  # initialize phi
  phi_0 = array(rbeta(n_Time*n_Age*q*K, a_phi, b_phi), dim = c(n_Time, n_Age, q, K))
  phi_1 = array(rbeta(n_Time*n_Age*q*K, a_phi, b_phi), dim = c(n_Time, n_Age, q, K))
  
  # initialize omega
  omega_0 = array(rgamma(n_Sex*n_Time*n_Age, a_omega, b_omega), dim = c(n_Sex, n_Time, n_Age))
  omega_1 = array(rgamma(n_Sex*n_Time*n_Age, a_omega, b_omega), dim = c(n_Sex, n_Time, n_Age))
  
  # initialize V
  V_0 = array(rbeta(n_Sex*n_Time*n_Age*K, 1, omega_0), dim = c(n_Sex, n_Time, n_Age, K))
  V_1 = array(rbeta(n_Sex*n_Time*n_Age*K, 1, omega_1), dim = c(n_Sex, n_Time, n_Age, K))
  V_0[, , , K] = V_1[, , , K] = 1
  
  # initialize lambda
  lambda_0 = array(rep(NA, n_Sex*n_Time*n_Age*K), dim = c(n_Sex, n_Time, n_Age, K))
  lambda_1 = array(rep(NA, n_Sex*n_Time*n_Age*K), dim = c(n_Sex, n_Time, n_Age, K))
  
  for(nSex in 1:n_Sex){
    for(nTime in 1:n_Time){
      for(nAge in 1:n_Age){
        lambda_0[nSex, nTime, nAge, 1] = V_0[nSex, nTime, nAge, 1]
        lambda_1[nSex, nTime, nAge, 1] = V_1[nSex, nTime, nAge, 1]
        for (k in 2:(K-1)) {
          V0_prod = V1_prod = 1
          for (j in 1:(k-1)) {
            V0_prod = V0_prod*(1 - V_0[nSex, nTime, nAge, j])
            V1_prod = V1_prod*(1 - V_1[nSex, nTime, nAge, j])
          }
          lambda_0[nSex, nTime, nAge, k] = V_0[nSex, nTime, nAge, k]*V0_prod
          lambda_1[nSex, nTime, nAge, k] = V_1[nSex, nTime, nAge, k]*V1_prod
        }
        lambda_0[nSex, nTime, nAge, K] = 1 - sum(lambda_0[nSex, nTime, nAge, 1:(K-1)])
        lambda_1[nSex, nTime, nAge, K] = 1 - sum(lambda_1[nSex, nTime, nAge, 1:(K-1)])
      }
    }
  }
  lambda_0_prior = lambda_0
  lambda_1_prior = lambda_1
  
  
  # initialize pi
  pis = rbeta(n_Sex*n_Time*n_Age, 1, 1)
  
  
  # store the samples
  sample_pis = matrix(0, nrow = E-BURN_IN, ncol = n_Sex*n_Time*n_Age)
  
  
  sample_lambda_0 = sample_lambda_1 = array(NA, dim = c(E-BURN_IN, n_Sex, n_Time, n_Age, K))
  sample_omega_0 = sample_omega_1 = array(NA, dim = c(E-BURN_IN, n_Sex, n_Time, n_Age))
  sample_phi_0 = sample_phi_1 = array(NA, dim = c(E-BURN_IN, q, K))  
  
  
  require("Rcpp")
  require("RcppArmadillo")
  
  out <- fit_model_BL_internal(E = E, BURN_IN = BURN_IN, n_Sex = n_Sex, n_Time = n_Time, n_Age = n_Age,
                                                            a_phi = a_phi, b_phi = b_phi, a_omega = a_omega, b_omega = b_omega,
                                                            a_pi = a_pi, b_pi = b_pi,
                                                            pis_r = pis, 
                                                            phi_0_r = phi_0, phi_1_r = phi_1, lambda_0_r = lambda_0, lambda_1_r = lambda_1,
                                                            V_0_r = V_0, V_1_r = V_1, omega_0_r = omega_0, omega_1_r = omega_1,
                                                            X_r = as.matrix(X), S_r = S, Ti_r = Ti - 1, A_r = A - 1, Yt_r = Y.t, H_r = H - 1, K_r = K,
                                                            sub_size_n_r = nsta)
  
  sample_pis = out$pis
  sample_Yt = out$Yt
  sample_H = out$H
  
  for(i in 1:(E - BURN_IN)){
    for(k in 1:K){
      for(nTime in 1:n_Time){
        for(nAge in 1:n_Age){
          for(nSex in 1:n_Sex){
            sample_lambda_1[i, nSex, nTime, nAge, k] = get_lambda(E - BURN_IN, n_Sex, n_Time, n_Age, K, out$lambda_1, i, nSex, nTime, nAge, k)
            sample_lambda_0[i, nSex, nTime, nAge, k] = get_lambda(E - BURN_IN, n_Sex, n_Time, n_Age, K, out$lambda_0, i, nSex, nTime, nAge, k)
          }
        }
      }
      for(j in 1:q){
        sample_phi_1[i, j, k] = get_phi(E - BURN_IN, n_Time, n_Age, q, K, out$phi_1, i, 1, 1, j, k)
        sample_phi_0[i, j, k] = get_phi(E - BURN_IN, n_Time, n_Age, q, K, out$phi_0, i, 1, 1, j, k)
      }
    }
  }
  
  bias = get_bias(sample_pis, sim_data$param$pi)
  
  poster_sample = list(sample_pis = sample_pis, sample_phi_1 = sample_phi_1, sample_phi_0 = sample_phi_0, sample_lambda_0 = sample_lambda_0, sample_lambda_1 = sample_lambda_1, 
                       sample_Yt = sample_Yt, sample_H = sample_H)
  
  return(list(bias = bias, poster_sample = poster_sample))
  
}



################# dependency functions #############################
expit <- function(x){exp(x)/(1 + exp(x))}
logit <- function(x){log(x/(1-x))}

get_bias = function(estimate, truth) {
  colMeans(estimate) - truth
}


get_phi <- function(dim1, dim2, dim3, dim4, dim5, array_input, index1, index2, index3, index4, index5){
  index <- index1 + (index2 - 1) * dim1 + (index3 - 1) * dim1 * dim2 + (index4 - 1) * dim1 * dim2 * dim3 + (index5 - 1) * dim1 * dim2 * dim3 * dim4
  
  return(array_input[index])
}

get_lambda <- function(dim1, dim2, dim3, dim4, dim5, array_input, index1, index2, index3, index4, index5){
  # Calculate the index in R-style indexing
  index <- index1 + (index2 - 1) * dim1 + (index3 - 1) * dim1 * dim2 + (index4 - 1) * dim1 * dim2 * dim3 + (index5 - 1) * dim1 * dim2 * dim3 * dim4
  return(array_input[index])
}


