library(gtools)
expit <- function(x){exp(x)/(1 + exp(x))}
logit <- function(x){log(x/(1-x))}


#' @param pi_true simulated prevalence
#' @param n_subsize sample size of each sub-population
#' @param q dimension of signs/symptoms X
#' @param n_Sex number of sex group
#' @param n_Time number of time group
#' @param n_Age number of age group
#' @param K number of latent classes
#' @param DEPEND_ON_Y bool, FALSE or TRUE, corresponding to case i and ii of whether sampling with dependency on Y
sim_with_sex_time_age <- function(pi_true, n_subsize, q = 10, n_Sex = 2, n_Time = 10, n_Age = 8, K = 10, DEPEND_ON_Y = FALSE){
  n = n_subsize * n_Sex * n_Time * n_Age
  S = c(rep(0, n_subsize * n_Time * n_Age), rep(1, n_subsize * n_Time * n_Age))
  Ti = rep(rep(1:n_Time, each = n_subsize*n_Age), n_Sex)
  A =  rep(rep(1:n_Age, each = n_subsize, n_Time), n_Sex)
  STA = rep(1:(n_Sex * n_Time * n_Age), each = n_subsize)
  
  
  a_phi = 1; b_phi = 1
  a_omega = 8; b_omega = 2 
  
  # generate phi
  phi_0 = array(rbeta(q*K, a_phi, b_phi), dim = c(q, K))
  phi_1 = array(rbeta(q*K, a_phi, b_phi), dim = c(q, K))
  
  # generate omega
  omega_0 = array(rgamma(n_Sex*n_Time*n_Age, a_omega, b_omega), dim = c(n_Sex, n_Time, n_Age))
  omega_1 = array(rgamma(n_Sex*n_Time*n_Age, a_omega, b_omega), dim = c(n_Sex, n_Time, n_Age))
  
  # generate V
  V_0 = array(rbeta(n_Sex*n_Time*n_Age*K, 1, omega_0), dim = c(n_Sex, n_Time, n_Age, K))
  V_1 = array(rbeta(n_Sex*n_Time*n_Age*K, 1, omega_1), dim = c(n_Sex, n_Time, n_Age, K))
  V_0[, , , K] = V_1[, , , K] = 1
  
  # generate lambda
  lambda_0 = array(rep(NA, n_Sex*n_Time*n_Age*K), dim = c(n_Sex, n_Time, n_Age, K))
  lambda_1 = array(rep(NA, n_Sex*n_Time*n_Age*K), dim = c(n_Sex, n_Time, n_Age, K))
  
  for(nSex in 1:n_Sex){
    for(nTime in 1:n_Time){
      for(nAge in 1:n_Age){
        lambda_0[nSex, nTime, nAge, ] = rdirichlet(1, rep(0.1, K))
        lambda_1[nSex, nTime, nAge, ] = rdirichlet(1, rep(0.1, K))
      }
    }
  }
  
  
  Y = rep(NA, n)
  Z = rep(NA, n)
  X = matrix(NA, nrow = n, ncol = q)
  
  
  for (nSex in c(0, 1)){
    for(nTime in 0:(n_Time-1)){
      for(nAge in 1:n_Age){
        i = nSex * n_Time * n_Age + nTime * n_Age + nAge
        indx = (i-1)*n_subsize+1
        sub_sample = rep(0, n_subsize)
        replace_1_indx = sample(1:n_subsize, pi_true[i]*n_subsize)
        sub_sample[replace_1_indx] = 1
        Y[indx:(indx+n_subsize-1)] = sub_sample
      }
    }
  }
  
  
  
  for(i in 1:n){
    if(Y[i] == 0){
      lambda = lambda_0[S[i]+1, Ti[i], A[i], ]
    }else{
      lambda = lambda_1[S[i]+1, Ti[i], A[i], ]
    }
    Z[i] = sample(1:K, 1, prob = lambda)
    for(j in 1:q){
      X[i, j] = rbinom(1, 1,  ifelse(Y[i] == 0, phi_0[Z[i], j], phi_1[Z[i], j]))
    }
  }

  
  sim_data_prep = list()
  sim_data_prep$n_Time = n_Time
  sim_data_prep$n_Age = n_Age
  sim_data_prep$n_Sex = n_Sex
  
  sim_data.truth = list()
  sim_data.truth$Y.t = Y
  sim_data.data = list()
  sim_data.data$Y.t = Y
  sim_data.data$X = X
  sim_data.data$S = S
  sim_data.data$Ti = Ti
  sim_data.data$A = A
  sim_data.data$STA = STA
  sim_data.data$K = K

  
  missing_props = array(NA, dim = c(n_Sex, n_Time, n_Age))
  alpha_T = c(1.2, rep(0.1, n_Time - 2), 1.2)
  alpha_A = c(0.4, 0.4, rep(-1.6, n_Age - 4), 0.4, 0.4)
  
  if(DEPEND_ON_Y == FALSE){
    for(i in 1:n_Time){
      beta = rep(0, dim(X)[2])
      indexes <- sample(1:dim(X)[2], 3, replace = FALSE)
      beta[indexes] <- rep(0.1, 3)
      for (s in c(0, 1)){
        for(j in 1:n_Age){
          subset_sample = which(S == s & Ti == i & A == j)
          X_modified <- X[subset_sample, ]
          X_modified[is.na(X_modified)] <- 0
          subset_sample_miss = expit(alpha_T[i] + alpha_A[j] + beta %*% t(X_modified))
          Z <- rbinom(subset_sample, 1, subset_sample_miss)
          sim_data.data$Y.t[subset_sample[which(Z == 0)]] <- NA
          missing_props[s+1, i, j] = mean(1-Z)
        }
      }
    }
  }else{
    gamma1 <- runif(1, min = -0.4, max = 0)
    gamma2 <- -gamma1
    
    for(i in 1:n_Time){
      beta = rep(0, dim(X)[2])
      indexes <- sample(1:dim(X)[2], 3, replace = FALSE)
      beta[indexes] <- rep(0.1, 3)
      for (s in c(0, 1)){
        for(j in 1:n_Age){
          subset_sample = which(S == s & Ti == i & A == j)
          X_modified <- X[subset_sample, ]
          X_modified[is.na(X_modified)] <- 0
          subset_sample_miss = expit(alpha_T[i] + alpha_A[j] + beta %*% t(X_modified) + gamma1 %*% Y[subset_sample] + gamma2 %*% (1-Y[subset_sample]))
          Z <- rbinom(subset_sample, 1, subset_sample_miss)
          sim_data.data$Y.t[subset_sample[which(Z == 0)]] <- NA
          missing_props[s+1, i, j] = mean(1-Z)
        }
      }
    }
  }
  
  
  sim_data.param = list()
  sim_data.param$pi = pi_true
  
  sim_data_prep$ns = unname(table(sim_data.data$S))
  sim_data_prep$nt = unname(table(sim_data.data$Ti))
  sim_data_prep$na = unname(table(sim_data.data$A))
  sim_data_prep$nsta =  unname(table(sim_data.data$STA))
  sim_data_prep$data = sim_data.data
  sim_data_prep$data.truth = sim_data.truth
  sim_data_prep$param = sim_data.param
  sim_data_prep$missing_props = missing_props
  
  return(sim_data_prep)
  
}
