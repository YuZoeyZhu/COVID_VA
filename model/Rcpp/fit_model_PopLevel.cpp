#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <stdio.h>
#include <math.h> 

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace RcppArmadillo;
using namespace arma;




// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

// FCN prototypes
double samplepg(double);
double exprnd(double);
double tinvgauss(double, double);
double truncgamma();
double randinvg(double);
double aterm(int, double, double);

// [[Rcpp::export]]
NumericVector rcpp_pgdraw(NumericVector b, NumericVector c)
{
  int m = b.size();
  int n = c.size();
  NumericVector y(n);
  // Setup
  int i, j, bi = 1;
  if (m == 1)
  {
    bi = b[0];
  }
  
  // Sample
  for (i = 0; i < n; i++)
  {
    if (m > 1)
    {
      bi = b[i];
    }
    
    // Sample
    y[i] = 0;
    for (j = 0; j < (int)bi; j++)
    {
      y[i] += samplepg(c[i]);
    }
  }
  return y;
}


// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg(double z)
{
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;
  
  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);
  
  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q); 
  
  double u, X;
  
  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1) 
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = tinvgauss(z, t);
    }
    
    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;
    
    while(1) 
    {
      Sn = Sn + asgn * aterm(i, X, t);
      
      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }
      
      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }
      
      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

// Generate exponential distribution random variates
double exprnd(double mu)
{
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
// Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm(int n, double x, double t)
{
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }    
  return (double)exp(f);
}

// Generate inverse gaussian random variates
double randinvg(double mu)
{
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {    
    out = mu*mu / out; 
  }    
  return out;
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables 
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma()
{
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;  
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss(double z, double t)
{
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma 
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }  
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }    
  return X;
}



// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// p is probability
// return integer from 0 to P-1
// [[Rcpp::export]]
int sample_prob(vec &x){
  int len = x.size();
  vec p = cumsum(x);
  
  double u = runif(1, 0.0, 1.0)(0);
  int i;
  for(i = 0; i < len; i++){
    if(u <= p(i)){
      return(i);
    }   
  }
  return(-1);
}

// [[Rcpp::export]]
double Polya_Gamma_Aug_sampling(double pi_c, double h_x, int sum_yt, int n, double Sigma){
  NumericVector n_input(1);
  n_input(0) = n;
  NumericVector pi_c_input(1);
  pi_c_input(0) = pi_c;
  NumericVector omega_cur = rcpp_pgdraw(n_input, pi_c_input);

  double var_pi = 1/(1/Sigma + omega_cur(0));

  double mean_pi = var_pi*(h_x/Sigma + sum_yt-n/2);
  
  double m_cur = Rcpp::rnorm(1, mean_pi, var_pi)(0);

  return(m_cur);
}

mat Polya_Gamma_Aug_joint_sampling(vec pi_c, colvec mu, vec sum_yt, vec n, mat inv_Sigma, int n_Data){
  int i;
  vec omega_cur(n_Data);
  NumericVector n_input(1);
  NumericVector pi_c_input(1);
  for(i=0; i<n_Data; i++){
    n_input(0) = n(i);
    pi_c_input(0) = pi_c(i);
    omega_cur(i) = rcpp_pgdraw(n_input, pi_c_input)(0);
  }
  mat diag_omega_cur = diagmat(omega_cur);
  mat inv_var_pi = inv_Sigma + diag_omega_cur;
  mat var_pi = inv(inv_var_pi);
  colvec sum_yt_minus_n(n_Data);
  for(i=0; i<n_Data; i++){
    sum_yt_minus_n(i) = sum_yt(i) - n(i)/2;
  }
  vec mean_pi = var_pi * (inv_Sigma * mu + sum_yt_minus_n);
  mat m_cur = mvnrnd(mean_pi, var_pi, 1);
  return(m_cur);
}

mat get_inv_Sigma(double sigma2_pi_init, double sigma2_pi, int n_Data){
  int i,j;
  mat A(n_Data-2, n_Data-2, fill::eye);
  mat B = {{-1, 2, -1}};
  mat C_1 = conv2(A, B);
  // C = conv2(a = diag(1, n_Data-2), b = t(c(-1, 2, -1)))
  vec D(n_Data);
  D(0) = 1;
  D(1) = -1;
  for(i = 2; i < n_Data; i++ ){
    D(i) = 0;
  }
  mat C_2(n_Data-1, n_Data);
  for(j = 0; j < n_Data; j++){
    C_2(0, j) = D(j);
  }
  for(i = 1; i < n_Data-1; i++){
    for(j = 0; j < n_Data; j++){
      C_2(i, j) = C_1(i-1, j);
    }
  }
  
  vec E(n_Data);
  for(i = 0; i < n_Data-2; i++ ){
    E(i) = 0;
  }
  E(n_Data-2) = -1;
  E(n_Data-1) = 1;
  
  mat C_3(n_Data, n_Data);
  for(i = 0; i < n_Data-1; i++){
    for(j = 0; j < n_Data; j++){
      C_3(i, j) = C_2(i, j);
    }
  }
  for(j = 0; j < n_Data; j++){
    C_3(n_Data-1, j) = E(j);
  }
  
  mat F(n_Data, n_Data, fill::zeros);
  F(0, 0) = 1;
  
  mat inv_Sigma_matrix = 1/sigma2_pi_init * F + 1/sigma2_pi * C_3;
  return(inv_Sigma_matrix);
}

// [[Rcpp::export]]
mat create_inv_Sigma_pi(double sigma2_pi, int n_Data) {
  mat inv_Sigma_pi(n_Data, n_Data);
  
  inv_Sigma_pi.fill(0.0);
  
  for (int i = 0; i < n_Data; i++) {
    inv_Sigma_pi(i, i ) = 1.0 / sigma2_pi;
  }
  
  return inv_Sigma_pi;
}

// [[Rcpp::export]]
double expit(double x){
  double y = exp(x)/(1 + exp(x));
  return(y);
}


// [[Rcpp::export]]
SEXP fit_model_PopLevel_internal(int E, int BURN_IN, double a_phi,double b_phi, double a_omega, double b_omega, double a_pi, double b_pi,
                                              double pi, SEXP phi_0_r, SEXP phi_1_r, SEXP lambda_0_r, SEXP lambda_1_r, SEXP V_0_r, SEXP V_1_r, double omega_0, double omega_1,
                                              SEXP X_r, SEXP Yt_r, SEXP H_r, int K_r) {
  mat phi_0 = as<mat>(phi_0_r);
  mat phi_1 = as<mat>(phi_1_r);
  vec lambda_0 = as<vec>(lambda_0_r);
  vec lambda_1 = as<vec>(lambda_1_r);
  vec V_0 = as<vec>(V_0_r);
  vec V_1 = as<vec>(V_1_r);
  mat X = as<mat>(X_r);
  vec Yt = as<vec>(Yt_r);
  vec H = as<vec>(H_r);
  int K = K_r;
   
  int n = X.n_rows;
  int q = X.n_cols;
  vec Yt_pred(n);
  mat Yt_out(E - BURN_IN, n);
  vec pi_out(E - BURN_IN);
  int count_Yt_0;
  int count_Yt_1;
  vec count_H_k(K);
  
  
  int i, j, k, r, itr;
  double first_term, second_term, first_term_H, tmp0, tmp1;
  double y_1_h_k_x_1_sum, y_1_h_k_x_0_sum, y_0_h_k_x_0_sum, y_0_h_k_x_1_sum;
  
  vec numerator_H(K);
  double denominator_H;
  vec sample_mult(4);
  vec sample_mult_H(K);
  IntegerVector seq_K(K);
  double H_sum;
  double h_k;
  double V0_prod;
  double V1_prod;
  double prod_V0K;
  double prod_V1K;
  double sum_lambda_0;
  double sum_lambda_1;
  
  vec lambda(K);
  mat phi(q, K);
  
  for(itr = 0; itr < E; itr++){
    
    if(itr % 100 == 0) Rcout << "."; 
    count_Yt_1 = 0;
    count_Yt_0 = 0;
    
    count_H_k.zeros();
    
    
    for(i = 0; i < n; i++){
      // Sample Yt_i under Z == 1, S == 0, notice the index is respect to the full dataset
      if(ISNA(Yt(i))){
        k = H(i);
        first_term = pi*lambda_1(H(i));
        second_term = (1 - pi)*lambda_0(H(i));
        for(j = 0; j < q; j++){
          if(!ISNA(X(i, j))){
            first_term *= pow(phi_1(j, H(i)), X(i, j));
            first_term *= pow(1 - phi_1(j, H(i)), 1 - X(i, j));
            second_term *= pow(phi_0(j, H(i)), X(i, j));
            second_term *= pow(1 - phi_0(j, H(i)), 1 - X(i, j));
          }
        }
        Yt_pred(i) = Rcpp::rbinom(1, 1, first_term / (first_term + second_term))(0);
      }
      
      // add the known obs of Yt 
      if(!ISNA(Yt(i))){
        Yt_pred(i) = Yt(i);
      }
      
      
      // count Yt == 0/1
      count_Yt_1 += Yt_pred(i);
      count_Yt_0 += 1 - Yt_pred(i);
      
      
      // sample H_i
      numerator_H.zeros();
      denominator_H = 0;
      if (Yt_pred(i) == 0){
        lambda = lambda_0;
        phi = phi_0;
      }else{
        lambda = lambda_1;
        phi = phi_1;
      }
      for(k = 0; k < K; k++){
        first_term_H = lambda(k);
        for(j = 0; j < q; j++){ 
          if(!ISNA(X(i, j))){
            first_term_H *= pow(phi(j, k), X(i, j));
            first_term_H *= pow(1 - phi(j, k), 1 - X(i, j));
          }
        }
        numerator_H(k) = first_term_H;
        denominator_H += first_term_H;
      }
      
      for (k = 0; k < K; k++){
        sample_mult_H(k) = numerator_H(k) / denominator_H;
      }
      
      seq_K = seq_len(K);
      
      H(i) =  Rcpp::RcppArmadillo::sample(seq_K, 1, false, sample_mult_H)(0) - 1;
      
      for(k = 0; k < K; k++) {
        if(H(i) == k){
          count_H_k(k) += 1;
        }
      }
    }
    
    
    
    // sample V
      for(k = 0; k < (K-1); k++) {
        H_sum = 0;
        for(r = k+1; r < K; r++){
          H_sum += count_H_k(r);
        }
        V_0(k) = Rcpp::rbeta(1, 1 + count_H_k(k), omega_0 + H_sum)(0);
        V_1(k) = Rcpp::rbeta(1, 1 + count_H_k(k), omega_1 + H_sum)(0);
      }
    
    
    // get lambda
      lambda_0(0) = V_0(0);
      lambda_1(0) = V_1(0);
      sum_lambda_0 = lambda_0(0);
      sum_lambda_1 = lambda_1(0);
      for (k = 1; k < (K-1); k++) {
        V0_prod = 1;
        V1_prod = 1;
        for (j = 0; j < k; j++) {
          V0_prod = V0_prod*(1 - V_0(j));
          V1_prod = V1_prod*(1 - V_1(j));
        }
        lambda_0(k) = V_0(k)*V0_prod;
        lambda_1(k) = V_1(k)*V1_prod;
        
        sum_lambda_0 += lambda_0(k);
        sum_lambda_1 += lambda_1(k);
      }
      lambda_0(K-1) = 1 - sum_lambda_0;
      lambda_1(K-1) = 1 - sum_lambda_1;
  
  
    
    // sample omega
      prod_V0K = 1;
      prod_V1K = 1;
      for (k = 0; k < (K-1); k++) {
        prod_V0K = prod_V0K*(1 - V_0(k));
        prod_V1K = prod_V1K*(1 - V_1(k));
      }
      omega_0 = Rcpp::rgamma(1, a_omega + K - 1, 1/(b_omega - log(prod_V0K)))(0);
      omega_1 = Rcpp::rgamma(1, a_omega + K - 1, 1/(b_omega - log(prod_V1K)))(0);
    
    
    // sample pi
      pi = Rcpp::rbeta(1, count_Yt_1 + a_pi, n - count_Yt_1 + b_pi)(0);
  
    // Sample phi
    for(j = 0; j < q; j++){
      for(k = 0; k < K; k++){
        y_1_h_k_x_1_sum = 0;
        y_1_h_k_x_0_sum = 0;
        y_0_h_k_x_1_sum = 0;
        y_0_h_k_x_0_sum = 0;
        
        for(i = 0; i < n; i++){
          if(H(i) == k){
            h_k = 1;
          }else{
            h_k = 0;
          }
          if(!ISNA(X(i, j))){
            y_1_h_k_x_1_sum += Yt_pred(i) * h_k * X(i, j);
            y_1_h_k_x_0_sum += Yt_pred(i) * h_k * (1 - X(i, j));
            y_0_h_k_x_1_sum += (1 - Yt_pred(i)) * h_k *X(i, j);
            y_0_h_k_x_0_sum += (1 - Yt_pred(i)) * h_k *(1 - X(i, j));
          }
        }
   
        tmp0 = Rcpp::rbeta(1, a_phi + y_0_h_k_x_1_sum, b_phi + y_0_h_k_x_0_sum)(0);
        tmp1 = Rcpp::rbeta(1, a_phi + y_1_h_k_x_1_sum, b_phi + y_1_h_k_x_0_sum)(0);
        phi_0(j, k) = tmp0;
        phi_1(j, k) = tmp1;
      }
    }
    
    
    // store Yt_pred and pi
    if(itr >= BURN_IN){
      for(i = 0; i < n; i++){
        Yt_out(itr - BURN_IN, i) = Yt_pred(i);
      }
      pi_out(itr - BURN_IN) = pi;
    }
  }
  
  List out;
  out("pi") = pi_out;
  out("Yt") = Yt_out;
  return out;
}





