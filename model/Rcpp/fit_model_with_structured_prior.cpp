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
  // Rcout << "omega_cur " << omega_cur << "\n";
  double var_pi = 1/(1/Sigma + omega_cur(0));
  // Rcout << "var_pi " << var_pi << "\n";
  double mean_pi = var_pi*(h_x/Sigma + sum_yt-n/2);
  // Rcout << "mean_pi " << mean_pi << "\n";
  double m_cur = Rcpp::rnorm(1, mean_pi, var_pi)(0);
  // Rcout << "m_cur " << m_cur << "\n";
  return(m_cur);
}

vec Polya_Gamma_Aug_joint_sampling(vec etas, colvec mu, vec sum_yt, vec n, mat inv_Sigma, int nSTA, mat map_STA){
  vec omega_cur(nSTA);
  
  for (int x = 0; x < nSTA; x++) {
    rowvec row = map_STA.row(x);
    double result = dot(row, etas);
    NumericVector n_input(1);
    n_input(0) = n(x);
    NumericVector pi_c_input(1);
    pi_c_input(0) = result;
    omega_cur(x) = rcpp_pgdraw(n_input, pi_c_input)(0);
  }
  mat var_pi = inv(inv_Sigma +  map_STA.t() * diagmat(omega_cur) * map_STA);
  colvec sum_yt_minus_n(nSTA);
  for(int i=0; i<nSTA; i++){
    sum_yt_minus_n(i) = sum_yt(i) - n(i)/2;
  }
  vec mean_pi = var_pi * (inv_Sigma * mu + map_STA.t() * sum_yt_minus_n);
  // Rcout << "omega_cur " << omega_cur << "\n";
  // Rcout << "var_pi " << var_pi << "\n";
  // Rcout << "mean_pi " << mean_pi << "\n";
  vec etas_cur = mvnrnd(mean_pi, var_pi, 1);
  // Rcout << "omega_cur " << omega_cur << "\n";
  // Rcout << "diag_omega_cur " << diag_omega_cur << "\n";
  // Rcout << "var_pi " << var_pi << "\n";
  // Rcout << "mean_pi " << mean_pi << "\n";
  return(etas_cur);
}

mat get_inv_Sigma(double sigma2_pi_init, double sigma2_pi, int n_sub){
  int i,j;
  mat A(n_sub-2, n_sub-2, fill::eye);
  mat B = {{-1, 2, -1}};
  mat C_1 = conv2(A, B);
  // C = conv2(a = diag(1, n_sub-2), b = t(c(-1, 2, -1)))
  vec D(n_sub);
  D(0) = 1;
  D(1) = -1;
  for(i = 2; i < n_sub; i++ ){
    D(i) = 0;
  }
  mat C_2(n_sub-1, n_sub);
  for(j = 0; j < n_sub; j++){
    C_2(0, j) = D(j);
  }
  for(i = 1; i < n_sub-1; i++){
    for(j = 0; j < n_sub; j++){
      C_2(i, j) = C_1(i-1, j);
    }
  }
  
  vec E(n_sub);
  for(i = 0; i < n_sub-2; i++ ){
    E(i) = 0;
  }
  E(n_sub-2) = -1;
  E(n_sub-1) = 1;
  
  mat C_3(n_sub, n_sub);
  for(i = 0; i < n_sub-1; i++){
    for(j = 0; j < n_sub; j++){
      C_3(i, j) = C_2(i, j);
    }
  }
  for(j = 0; j < n_sub; j++){
    C_3(n_sub-1, j) = E(j);
  }
  
  mat F(n_sub, n_sub, fill::zeros);
  F(0, 0) = 1;
  
  mat inv_Sigma_matrix = 1/sigma2_pi_init * F + 1/sigma2_pi * C_3;
  return(inv_Sigma_matrix);
}

// [[Rcpp::export]]
mat invSigmaPi(mat inv_Sigma_pi_T, mat inv_Sigma_pi_A, double sigma2_pi_0, double sigma2_pi_S, int n_Time, int n_Age) {
  mat inv_Sigma_pi(n_Time + n_Age + 1 + 1, n_Time + n_Age + 1 + 1);
  inv_Sigma_pi.zeros();
  inv_Sigma_pi.submat(1+1, 1+1, n_Time+1, n_Time+1) = inv_Sigma_pi_T;
  inv_Sigma_pi.submat(n_Time + 1+1, n_Time + 1+1, n_Time + n_Age+1, n_Time + n_Age+1) = inv_Sigma_pi_A;
  // inv_Sigma_pi.submat(1+1, n_Time + 1+1, n_Time+1, n_Time + n_Age+1) = zeros(n_Time, n_Age);
  // inv_Sigma_pi.submat(n_Time + 1+1, 1+1, n_Time + n_Age+1, n_Time+1) = zeros(n_Age, n_Time);
  inv_Sigma_pi(0, 0) = 1.0 / sigma2_pi_0;
  inv_Sigma_pi(1, 1) = 1.0 / sigma2_pi_S;
  // inv_Sigma_pi.submat(0, 1, 0, n_Time + n_Age) = zeros(1, n_Time + n_Age);
  // inv_Sigma_pi.submat(1, 0, n_Time + n_Age, 0) = zeros(n_Time + n_Age, 1);
  
  return inv_Sigma_pi;
}


// [[Rcpp::export]]
mat invSigmaPiInteraction(mat inv_Sigma_pi_T, mat inv_Sigma_pi_A, double sigma2_pi_0, double sigma2_pi_S, double sigma2_pi_STA, int n_Time, int n_Age, int n_Sex) {
  int n = n_Time + n_Age + 1 + 1 + n_Sex*n_Time*n_Age;
  mat inv_Sigma_pi(n, n, fill::zeros);
  
  // set P[1,1] to 1/100
  inv_Sigma_pi(0, 0) = 1.0 / sigma2_pi_0;
  
  inv_Sigma_pi(1, 1) = 1.0 / sigma2_pi_S;
  
  // set diagonal elements of T*T matrix to 1
  inv_Sigma_pi.submat(1+1, 1+1, n_Time+1, n_Time+1) = inv_Sigma_pi_T;
  
  // set diagonal elements of A*A matrix to 2
  inv_Sigma_pi.submat(n_Time + 1+1, n_Time + 1+1, n_Time + n_Age+1, n_Time + n_Age+1) = inv_Sigma_pi_A;
  
  // set diagonal elements of (T*A)*(T*A) matrix to 4
  inv_Sigma_pi.submat(n_Time + n_Age + 1+1, n_Time + n_Age + 1+1, n - 1, n - 1) = (1.0 / sigma2_pi_STA)* mat(n_Sex*n_Time*n_Age, n_Sex*n_Time*n_Age, fill::eye);
  
  
  return inv_Sigma_pi;
}

// [[Rcpp::export]]
mat create_inv_Sigma_pi(double sigma2_pi_0, double sigma2_pi_T, double sigma2_pi_A, int n_Time, int n_Age) {
  int n = 1 + n_Time + n_Age;
  mat inv_Sigma_pi(n, n);
  
  inv_Sigma_pi.fill(0.0);
  inv_Sigma_pi(0, 0) = 1.0 / sigma2_pi_0;
  
  for (int i = 0; i < n_Time; i++) {
    inv_Sigma_pi(i + 1, i + 1) = 1.0 / sigma2_pi_T;
  }
  
  for (int i = 0; i < n_Age; i++) {
    inv_Sigma_pi(n_Time + i + 1, n_Time + i + 1) = 1.0 / sigma2_pi_A;
  }
  
  return inv_Sigma_pi;
}

// [[Rcpp::export]]
mat create_inv_Sigma_pi_interaction(double sigma2_pi_0, double sigma2_pi_S, double sigma2_pi_T, double sigma2_pi_A, double sigma2_pi_STA, int n_Time, int n_Age, int n_Sex) {
  int n = n_Time + n_Age + 1 + 1 + n_Sex*n_Time*n_Age;
  mat inv_Sigma_pi(n, n);
  
  inv_Sigma_pi.fill(0.0);
  inv_Sigma_pi(0, 0) = 1.0 / sigma2_pi_0;
  inv_Sigma_pi(1, 1) = 1.0 / sigma2_pi_S;
  
  for (int i = 0; i < n_Time; i++) {
    inv_Sigma_pi(i + 1 + 1, i + 1 + 1) = 1.0 / sigma2_pi_T;
  }
  
  for (int i = 0; i < n_Age; i++) {
    inv_Sigma_pi(n_Time + i + 1 + 1, n_Time + i + 1 + 1) = 1.0 / sigma2_pi_A;
  }
  
  // set diagonal elements of (T*A)*(T*A) matrix to 4
  inv_Sigma_pi.submat(n_Time + n_Age + 1 + 1, n_Time + n_Age + 1 + 1, n - 1, n - 1) = (1.0 / sigma2_pi_STA)* mat(n_Sex*n_Time*n_Age, n_Sex*n_Time*n_Age, fill::eye);
  
  
  return inv_Sigma_pi;
}



// [[Rcpp::export]]
vec zerosRow(int n_Time, int n_Age) {
  vec zeros_row(n_Time + n_Age + 1);
  zeros_row.zeros();
  
  return zeros_row;
}

// [[Rcpp::export]]
vec zerosRowInteraction(int n_Time, int n_Age, int n_Sex) {
  vec zeros_row(n_Sex*n_Time * n_Age + n_Time + n_Age + 1 + 1);
  zeros_row.zeros();
//zeros_row.fill(1.5);
  
  return zeros_row;
}

// [[Rcpp::export]]
vec vectorize_mat_row_by_row(mat matrx) {
  int n = matrx.n_rows;
  int p = matrx.n_cols;
  vec vect(n * p);
  
  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      vect(i * p + j) = matrx(i, j);
    }
  }
  
  return vect;
}

// [[Rcpp::export]]
vec vectorize_cube_row_by_row(cube matrx) {
  int n = matrx.n_rows; // sex
  int p = matrx.n_cols; // time
  int q = matrx.n_slices; // age
  vec vect(n * p * q);
  
  for (int s = 0; s < n; s++) {
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < q; j++){
        vect(s * p * q  + i * q + j) = matrx(s, i, j);
      }
    }
  }

  
  return vect;
}

// [[Rcpp::export]]
double expit(double x){
  double y = exp(x)/(1 + exp(x));
  return(y);
}

// [[Rcpp::export]]
field<cube> read_4d_array(SEXP arr, int n_Time, int n_Age, int qq, int kk) {
  Rcpp::NumericVector arr_vec(arr);
  field<cube> arr_field(n_Time);
  
  for (int nTime = 0; nTime < n_Time; nTime++) {
    arr_field(nTime) = arma::cube(n_Age, qq, kk, arma::fill::zeros);
    for (int nAge = 0; nAge < n_Age; nAge++) {
      for (int q = 0; q < qq; q++) {
        for (int k = 0; k < kk; k++) {
          arr_field(nTime)(nAge,q,k) = arr_vec(nTime * n_Age * qq * kk + nAge * qq * kk + q * kk + k);
        }
      }
    }
  }
  return arr_field;
}

// Function 2: Return the element at a specific index in the 6-dimensional array
// Input: NumericVector representing the 6-dimensional array, index1, index2, index3, index4, index5, index6
// Output: Element at the specified index
// [[Rcpp::export]]
double get5DElementAtIndex(NumericVector array_input, int dim1, int dim2, int dim3, int dim4, int dim5, int index1, int index2, int index3, int index4, int index5) {
  
  // int index = index5 + dim5 * (index4 + dim4 * (index3 + dim3 * (index2 + dim2 * index1)));
  int index = 0;
  index = index1 + index2 * dim1 + index3 * dim1 * dim2 + index4 * dim1 * dim2 * dim3 + index5 * dim1 * dim2 * dim3 * dim4;
  
  return array_input(index);
}

// Function 3: Update the element at a specific index in the 6-dimensional array
// Input: NumericVector representing the 6-dimensional array, new_value, index1, index2, index3, index4, index5, index6
// Output: None (array_input is modified in place)
// [[Rcpp::export]]
void update5DElementAtIndex(NumericVector array_input, int dim1, int dim2, int dim3, int dim4, int dim5, double new_value, int index1, int index2, int index3, int index4, int index5) {
  
  // int index = index5 + dim5 * (index4 + dim4 * (index3 + dim3 * (index2 + dim2 * index1)));
  int index = 0;
  index = index1 + index2 * dim1 + index3 * dim1 * dim2 + index4 * dim1 * dim2 * dim3 + index5 * dim1 * dim2 * dim3 * dim4;
  array_input(index) = new_value;
}

// Function 3: Update the element at a specific index in the 6-dimensional array
// Input: NumericVector representing the 6-dimensional array, new_value, index1, index2, index3, index4, index5, index6
// Output: None (array_input is modified in place)
// [[Rcpp::export]]
void update4DElementAtIndex(NumericVector array_input, int dim1, int dim2, int dim3, int dim4, double new_value, int index1, int index2, int index3, int index4) {
  
  // int index = index5 + dim5 * (index4 + dim4 * (index3 + dim3 * (index2 + dim2 * index1)));
  int index = 0;
  index = index1 + index2 * dim1 + index3 * dim1 * dim2 + index4 * dim1 * dim2 * dim3;
  array_input(index) = new_value;
}

// [[Rcpp::export]]
SEXP fit_model_with_structured_prior_internal(int E, int BURN_IN, int n_Sex, int n_Time, int n_Age, 
                                                                          double a_phi, double b_phi, double a_omega, double b_omega, 
                                                                          double a_sigma_T, double b_sigma_T, double a_sigma_A, double b_sigma_A, double a_sigma_STA, double b_sigma_STA,
                                                                          SEXP ms_r, SEXP phi_0_r, SEXP phi_1_r, SEXP lambda_0_r, SEXP lambda_1_r, SEXP V_0_r, SEXP V_1_r, SEXP omega_0_r, SEXP omega_1_r,
                                                                          SEXP X_r, SEXP Ti_r, SEXP A_r, SEXP S_r, SEXP Yt_r, SEXP H_r, int K_r, SEXP map_STA_r,
                                                                          SEXP etas_r,  double sigma2_pi_0, double sigma2_pi_S, double sigma2_pi_T, double sigma2_pi_A, double sigma2_pi_init_T, double sigma2_pi_init_A, double sigma2_pi_STA, SEXP sub_size_n_r, 
                                                                          std::string STRUCTURED_PRIOR) {
  
  
  int K = K_r;
  field<cube> lambda_0(n_Sex);
  field<cube> lambda_1(n_Sex);
  lambda_0 = read_4d_array(lambda_0_r, n_Sex, n_Time, n_Age, K);
  lambda_1 = read_4d_array(lambda_1_r, n_Sex, n_Time, n_Age, K);
  
  field<cube> V_0(n_Sex);
  field<cube> V_1(n_Sex);
  V_0 = read_4d_array(V_0_r, n_Sex, n_Time, n_Age, K);
  V_1 = read_4d_array(V_1_r, n_Sex, n_Time, n_Age, K);

  cube omega_0 = as<cube>(omega_0_r);
  cube omega_1 = as<cube>(omega_1_r);
  
  vec ms = as<vec>(ms_r);
  mat X = as<mat>(X_r);
  vec S = as<vec>(S_r);
  vec Ti = as<vec>(Ti_r);
  vec A = as<vec>(A_r);
  vec Yt = as<vec>(Yt_r);
  vec H = as<vec>(H_r);

  
  vec etas = as<vec>(etas_r);
  mat map_STA = as<mat>(map_STA_r);
  vec sub_size_n = as<vec>(sub_size_n_r);
  int n = X.n_rows;
  int q = X.n_cols;
  vec Yt_pred(n);
  vec pis(n_Sex*n_Time*n_Age);
  cube pis_matrix(n_Sex, n_Time, n_Age);
  mat Yt_out(E - BURN_IN, n);
  mat H_out(E - BURN_IN, n);
  mat pis_out(E - BURN_IN, n_Sex*n_Time*n_Age);
  mat etas_interaction_out(E - BURN_IN,  n_Sex*n_Time*n_Age + n_Time + n_Age + 1 + 1);
  
  vec sigma2_pi_T_out(E - BURN_IN);
  vec sigma2_pi_A_out(E - BURN_IN);
  vec sigma2_pi_STA_out(E - BURN_IN);
  cube count_Yt_0_mat(n_Sex, n_Time, n_Age);
  cube count_Yt_1_mat(n_Sex, n_Time, n_Age);
  vec count_Yt_0(n_Sex * n_Time * n_Age);
  vec count_Yt_1(n_Sex * n_Time * n_Age);
  

  field<cube> count_H_k(n_Sex);
  field<cube> count_H_k_Y_1(n_Sex);
  field<cube> count_H_k_Y_0(n_Sex);
  cube count_init(n_Time, n_Age, K);
  for(int nSex = 0; nSex<n_Sex; nSex++){
    count_H_k(nSex) = count_init;
    count_H_k_Y_1(nSex) = count_init;
    count_H_k_Y_0(nSex) = count_init;
  }
  double eta_0;
  double eta_S;
  vec eta_T(n_Time);
  vec eta_A(n_Age);
  vec eta_STA(n_Sex*n_Time*n_Age);
  
  
  int nSex, nTime, nAge, s, i, j, k, r, itr;
  double first_term, second_term, first_term_k, second_term_k, first_term_H, tmp0, tmp1;
  cube y_1_h_k_x_1(n_Sex, n_Time, n_Age);
  cube y_1_h_k_x_0(n_Sex, n_Time, n_Age);
  cube y_0_h_k_x_0(n_Sex, n_Time, n_Age);
  cube y_0_h_k_x_1(n_Sex, n_Time, n_Age);
  vec y_1_c_1(n_Age); 
  vec y_1_c_0(n_Age);
  vec y_0_c_0(n_Age);
  vec y_0_c_1(n_Age);
  double y_1_h_k_x_1_sum, y_1_h_k_x_0_sum, y_0_h_k_x_0_sum, y_0_h_k_x_1_sum;
  
  mat inv_Sigma_pi(n_Time + n_Age + 1 + 1, n_Time + n_Age + 1 + 1);
  colvec mu_m_interaction(n_Sex*n_Time * n_Age + n_Time + n_Age + 1 + 1);
  mat inv_Sigma_pi_interaction(n_Time + n_Age + 1 + 1 + n_Sex*n_Time*n_Age, n_Time + n_Age + 1 + 1 + n_Sex*n_Time*n_Age);
  
  
  vec numerator_H(K);
  double denominator_H;
  vec sample_mult(4);
  vec sample_mult_H(K);
  IntegerVector seq_K(K);
  double H_Y1_sum, H_Y0_sum;
  double h_k;
  double V0_prod;
  double V1_prod;
  double prod_V0K;
  double prod_V1K;
  double sum_lambda_0;
  double sum_lambda_1;
  double tol = 0.00001;
  
  int nSTA = n_Sex*n_Time*n_Age;
  field<cube> lambda(n_Sex);
  cube lambda_init(n_Time, n_Age, K);
  for(nSex = 0; nSex<n_Sex; nSex++){
    lambda(nSex) = lambda_init;
  }
  
  field<cube> phi_0(n_Time);
  field<cube> phi_1(n_Time);
  phi_0 = read_4d_array(phi_0_r, n_Time, n_Age, q, K);
  phi_1 = read_4d_array(phi_1_r, n_Time, n_Age, q, K);
  field<cube> phi(n_Time);
  cube phi_init(n_Age, q, K);
  for(nTime = 0; nTime<n_Time; nTime++){
    phi(nTime) = phi_init;
  }
  
  NumericVector lambda_1_out((E - BURN_IN)*n_Sex*n_Time*n_Age*K);
  NumericVector lambda_0_out((E - BURN_IN)*n_Sex*n_Time*n_Age*K);
  NumericVector phi_1_out((E - BURN_IN)*n_Time*n_Age*q*K);
  NumericVector phi_0_out((E - BURN_IN)*n_Time*n_Age*q*K);
  
  
  // start posterior sampling
  for(itr = 0; itr < E; itr++){
    
    if(itr % 100 == 0) Rcout << "."; 
    // initialization counts
    count_Yt_1.zeros();
    count_Yt_0.zeros();
    count_Yt_1_mat.zeros();
    count_Yt_0_mat.zeros();
    count_init.zeros();
    for(nSex = 0; nSex<n_Sex; nSex++){
      count_H_k(nSex) = count_init;
      count_H_k_Y_1(nSex) = count_init;
      count_H_k_Y_0(nSex) = count_init;
    }
    
    for(i = 0; i < nSTA; i++){
      pis(i) = expit(ms(i));
    }
    
    for (s = 0; s < n_Sex; s++) {
      for (int i = 0; i < n_Time; i++) {
        for (int j = 0; j < n_Age; j++) {
          pis_matrix(s, i, j) = pis(s * n_Time * n_Age + i * n_Age + j);
        }
      }
    }
    
    for(i = 0; i < n; i++){
      // Sample Yt_i under Z == 1, S == 0, notice the index is respect to the full dataset
      first_term = 0;
      second_term = 0;
      if(ISNA(Yt(i))){
        for(k = 0; k < K; k++){
          first_term_k = lambda_1(S(i))(Ti(i), A(i), H(i));
          second_term_k = lambda_0(S(i))(Ti(i), A(i), H(i));
          for(j = 0; j < q; j++){
            if(!ISNA(X(i, j))){
              first_term_k *= pow(phi_1(Ti(i))(A(i), j, k), X(i, j));
              first_term_k *= pow(1 - phi_1(Ti(i))(A(i), j, k), 1 - X(i, j));
              second_term_k *= pow(phi_0(Ti(i))(A(i), j, k), X(i, j));
              second_term_k *= pow(1 - phi_0(Ti(i))(A(i), j, k), 1 - X(i, j));
            }
          }
          first_term += first_term_k;
          second_term += second_term_k;
        }
        first_term *= pis_matrix(S(i), Ti(i), A(i));
        second_term *= 1 - pis_matrix(S(i), Ti(i), A(i));
        Yt_pred(i) = Rcpp::rbinom(1, 1, first_term / (first_term + second_term))(0);
      }
   
      // add the known obs of Yt 
      if(!ISNA(Yt(i))){
        Yt_pred(i) = Yt(i);
      }

      // count Yt == 0/1
      count_Yt_1_mat(S(i), Ti(i), A(i)) += Yt_pred(i);
      count_Yt_0_mat(S(i), Ti(i), A(i)) += 1 - Yt_pred(i);
      
      
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
        first_term_H = lambda(S(i))(Ti(i), A(i), k);
        for(j = 0; j < q; j++){ 
          if(!ISNA(X(i, j))){
            first_term_H *= pow(phi(Ti(i))(A(i), j, k), X(i, j));
            first_term_H *= pow(1 - phi(Ti(i))(A(i), j, k), 1 - X(i, j));
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
          count_H_k(S(i))(Ti(i), A(i), k) += 1;
          if(Yt_pred(i) == 1){
            count_H_k_Y_1(S(i))(Ti(i), A(i), k) += 1;
          }else{
            count_H_k_Y_0(S(i))(Ti(i), A(i), k) += 1;
          }
        }
      }
    }
    // transfer count_Yt_1_mat to vector count_Yt_1 after loop completed
    count_Yt_1 = vectorize_cube_row_by_row(count_Yt_1_mat); // require: complete the function
    count_Yt_0 = vectorize_cube_row_by_row(count_Yt_0_mat);

    
    for(nSex = 0; nSex < n_Sex; nSex++){
      for(nTime = 0; nTime < n_Time; nTime++){
        for(nAge = 0; nAge < n_Age; nAge++){
          
          // sample V
          for(k = 0; k < (K-1); k++) {
            H_Y1_sum = 0;
            H_Y0_sum = 0;
            for(r = k+1; r < K; r++){
              H_Y1_sum += count_H_k_Y_1(nSex)(nTime, nAge, r);
              H_Y0_sum += count_H_k_Y_0(nSex)(nTime, nAge, r);
            }
            V_0(nSex)(nTime, nAge, k) = Rcpp::rbeta(1, 1 + count_H_k_Y_0(nSex)(nTime, nAge, k), omega_0(nSex, nTime, nAge) + H_Y0_sum)(0);
            V_1(nSex)(nTime, nAge, k) = Rcpp::rbeta(1, 1 + count_H_k_Y_1(nSex)(nTime, nAge, k), omega_1(nSex, nTime, nAge) + H_Y1_sum)(0);
          }
          
          // get lambda
          lambda_0(nSex)(nTime, nAge, 0) = V_0(nSex)(nTime, nAge, 0);
          lambda_1(nSex)(nTime, nAge, 0) = V_1(nSex)(nTime, nAge, 0);
          sum_lambda_0 = lambda_0(nSex)(nTime, nAge, 0);
          sum_lambda_1 = lambda_1(nSex)(nTime, nAge, 0);
          for (k = 1; k < (K-1); k++) {
            V0_prod = 1;
            V1_prod = 1;
            for (j = 0; j < k; j++) {
              V0_prod = V0_prod*(1 - V_0(nSex)(nTime, nAge, j));
              V1_prod = V1_prod*(1 - V_1(nSex)(nTime, nAge, j));
            }
            lambda_0(nSex)(nTime, nAge, k) = V_0(nSex)(nTime, nAge, k)*V0_prod;
            lambda_1(nSex)(nTime, nAge, k) = V_1(nSex)(nTime, nAge, k)*V1_prod;
            
            sum_lambda_0 += lambda_0(nSex)(nTime, nAge, k);
            sum_lambda_1 += lambda_1(nSex)(nTime, nAge, k);
          }
          if(sum_lambda_0 > 1 - tol){
            lambda_0(nSex)(nTime, nAge, K-1) = tol;
          }else{
            lambda_0(nSex)(nTime, nAge, K-1) = 1 - sum_lambda_0;
          }
          if(sum_lambda_1 > 1 - tol){
            lambda_1(nSex)(nTime, nAge, K-1) = tol;
          }else{
            lambda_1(nSex)(nTime, nAge, K-1) = 1 - sum_lambda_1;
          }
          
          // sample omega
          prod_V0K = 1;
          prod_V1K = 1;
          for (k = 0; k < (K-1); k++) { //
            prod_V0K = prod_V0K*(1 - V_0(nSex)(nTime, nAge, k));
            prod_V1K = prod_V1K*(1 - V_1(nSex)(nTime, nAge, k));
          }
          omega_0(nSex, nTime, nAge) = Rcpp::rgamma(1, a_omega + K - 1, 1/(b_omega - log(prod_V0K)))(0); 
          omega_1(nSex, nTime, nAge) = Rcpp::rgamma(1, a_omega + K - 1, 1/(b_omega - log(prod_V1K)))(0); 
        }
      }
    }
    

    
    
    
    
    // sample m
    // get inv_Sigma_pi
    if(STRUCTURED_PRIOR == "RW"){
      mat inv_Sigma_pi_T = get_inv_Sigma(sigma2_pi_init_T, sigma2_pi_T, n_Time);
      mat inv_Sigma_pi_A = get_inv_Sigma(sigma2_pi_init_A, sigma2_pi_A, n_Age);
      inv_Sigma_pi_interaction =  invSigmaPiInteraction(inv_Sigma_pi_T, inv_Sigma_pi_A, sigma2_pi_0, sigma2_pi_S, sigma2_pi_STA, n_Time, n_Age, n_Sex);
      // get mu_m
      mu_m_interaction = zerosRowInteraction(n_Time, n_Age, n_Sex);
      // get etas and m
      etas = Polya_Gamma_Aug_joint_sampling(etas, mu_m_interaction, count_Yt_1, sub_size_n, inv_Sigma_pi_interaction, nSTA, map_STA);
    }else if((STRUCTURED_PRIOR == "Indep") || (STRUCTURED_PRIOR == "Fixed")){
      inv_Sigma_pi_interaction =  create_inv_Sigma_pi_interaction(sigma2_pi_0, sigma2_pi_S, sigma2_pi_T, sigma2_pi_A, sigma2_pi_STA, n_Time, n_Age, n_Sex);
      // get mu_m
      mu_m_interaction = zerosRowInteraction(n_Time, n_Age, n_Sex);
      // get etas and m
      etas = Polya_Gamma_Aug_joint_sampling(etas, mu_m_interaction, count_Yt_1, sub_size_n, inv_Sigma_pi_interaction, nSTA, map_STA);
    }
    

    
    // transfer centered etas to m
    ms = map_STA * etas;
    
    eta_0 = etas(0);
    eta_S = etas(1);
    eta_T = etas.subvec(2, n_Time+1);
    eta_A = etas.subvec(n_Time+1+1, n_Time+n_Age+1);
    
    eta_STA = etas.subvec(n_Time+n_Age+1+1, n_Time+n_Age + n_Sex*n_Time*n_Age+1);

    
    // Sample phi
    for(j = 0; j < q; j++){
      for(k = 0; k < K; k++){
        y_1_h_k_x_1.zeros();
        y_1_h_k_x_0.zeros();
        y_0_h_k_x_1.zeros();
        y_0_h_k_x_0.zeros();
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
            y_1_h_k_x_1(S(i), Ti(i), A(i)) += Yt_pred(i) * h_k * X(i, j);
            y_1_h_k_x_0(S(i), Ti(i), A(i)) += Yt_pred(i) * h_k * (1 - X(i, j));
            y_0_h_k_x_1(S(i), Ti(i), A(i)) += (1 - Yt_pred(i)) * h_k *X(i, j);
            y_0_h_k_x_0(S(i), Ti(i), A(i)) += (1 - Yt_pred(i)) * h_k *(1 - X(i, j));
          }
        }
        for(nSex = 0; nSex < n_Sex; nSex++){
          for(nTime = 0; nTime < n_Time; nTime++){
            for(nAge = 0; nAge < n_Age; nAge++){
              y_1_h_k_x_1_sum += y_1_h_k_x_1(nSex, nTime, nAge);
              y_1_h_k_x_0_sum += y_1_h_k_x_0(nSex, nTime, nAge);
              y_0_h_k_x_1_sum += y_0_h_k_x_1(nSex, nTime, nAge);
              y_0_h_k_x_0_sum += y_0_h_k_x_0(nSex, nTime, nAge);
            }
          }
        }
        tmp0 = Rcpp::rbeta(1, a_phi + y_0_h_k_x_1_sum, b_phi + y_0_h_k_x_0_sum)(0);
        tmp1 = Rcpp::rbeta(1, a_phi + y_1_h_k_x_1_sum, b_phi + y_1_h_k_x_0_sum)(0);
        for(nTime = 0; nTime < n_Time; nTime++){
          for(nAge = 0; nAge < n_Age; nAge++){
            phi_0(nTime)(nAge, j, k) = tmp0;
            phi_1(nTime)(nAge, j, k) = tmp1;
          }
        }
      }
    }
    
    
    // sample sigma2_pi_0, sigma2_pi_T and sum_pis_sq_A
    if(STRUCTURED_PRIOR == "RW"){
      double sum_pis_sq_T = 0.0;
      for (int i = 1; i < n_Time; i++) {
        sum_pis_sq_T += pow(eta_T(i) - eta_T(i - 1), 2);
      }
      sigma2_pi_T =  1 / Rcpp::rgamma(1, a_sigma_T + (n_Time-1)/2, 1/(b_sigma_T + sum_pis_sq_T/2))(0);
      
      
      double sum_pis_sq_A = 0.0;
      for (int i = 1; i < n_Age; i++) {
        sum_pis_sq_A += pow(eta_A(i) - eta_A(i - 1), 2);
      }
      sigma2_pi_A =  1 / Rcpp::rgamma(1, a_sigma_A + (n_Age-1)/2, 1/(b_sigma_A + sum_pis_sq_A/2))(0);
      
      double sum_pis_sq_STA = 0.0;
      for (int i = 0; i < n_Sex*n_Time*n_Age; i++) {
        sum_pis_sq_STA += pow(eta_STA(i), 2);
      }
      sigma2_pi_STA =  1 / Rcpp::rgamma(1, a_sigma_STA + (n_Sex*n_Time*n_Age)/2, 1/(b_sigma_STA + sum_pis_sq_STA/2))(0);
      
    }else if(STRUCTURED_PRIOR == "Indep"){
      double sum_pis_sq_T = 0.0;
      for (int i = 0; i < n_Time; i++) { //
        sum_pis_sq_T += pow(eta_T(i), 2);
      }
      sigma2_pi_T = 1 / Rcpp::rgamma(1, a_sigma_T + n_Time/2,  1 / (b_sigma_T + sum_pis_sq_T/2))(0);
      
      double sum_pis_sq_A = 0.0;
      for (int i = 0; i < n_Age; i++) { //
        sum_pis_sq_A += pow(eta_A(i), 2);
      }
      sigma2_pi_A = 1 / Rcpp::rgamma(1, a_sigma_A + n_Age/2,  1 / (b_sigma_A + sum_pis_sq_A/2))(0);
      
      double sum_pis_sq_STA = 0.0;
      for (int i = 0; i < n_Sex*n_Time*n_Age; i++) {
        sum_pis_sq_STA += pow(eta_STA(i), 2);
      }
      sigma2_pi_STA =  1 / Rcpp::rgamma(1, a_sigma_STA + (n_Sex*n_Time*n_Age)/2, 1 / (b_sigma_STA + sum_pis_sq_STA/2))(0);
      
    }else if(STRUCTURED_PRIOR == "Fixed"){
      double sum_pis_sq_STA = 0.0;
      for (int i = 0; i < n_Sex*n_Time*n_Age; i++) {
        sum_pis_sq_STA += pow(eta_STA(i), 2);
      }
      sigma2_pi_STA =  1 / Rcpp::rgamma(1, a_sigma_STA + (n_Sex*n_Time*n_Age)/2, 1/(b_sigma_STA + sum_pis_sq_STA/2))(0);
    }
    
    
    
    // store Yt_pred and pi
    if(itr >= BURN_IN){
      for(i = 0; i < n; i++){
        Yt_out(itr - BURN_IN, i) = Yt_pred(i);
        H_out(itr - BURN_IN, i) = H(i) + 1;
      }
      for(i = 0; i < nSTA; i++){
        pis_out(itr - BURN_IN, i) = pis(i);
      }
      for(nTime = 0; nTime < n_Time; nTime++){
        for(nAge = 0; nAge < n_Age; nAge++){
          for(k = 0; k < K; k++){
            for(nSex = 0; nSex < n_Sex; nSex++){
              update5DElementAtIndex(lambda_1_out, E - BURN_IN, n_Sex, n_Time, n_Age, K, lambda_1(nSex)(nTime, nAge, k), itr - BURN_IN, nSex, nTime, nAge, k);
              update5DElementAtIndex(lambda_0_out, E - BURN_IN, n_Sex, n_Time, n_Age, K, lambda_0(nSex)(nTime, nAge, k), itr - BURN_IN, nSex, nTime, nAge, k);
            } 
            for(j = 0; j < q; j++){
              update5DElementAtIndex(phi_1_out, E - BURN_IN, n_Time, n_Age, q, K, phi_1(nTime)(nAge, j, k), itr - BURN_IN, nTime, nAge, j, k);
              update5DElementAtIndex(phi_0_out, E - BURN_IN, n_Time, n_Age, q, K, phi_0(nTime)(nAge, j, k), itr - BURN_IN, nTime, nAge, j, k);
            }
          }
        }
      }
      
      for(i = 0; i < n_Sex*n_Time*n_Age + n_Time + n_Age + 1 + 1; i++){
        etas_interaction_out(itr - BURN_IN, i) = etas(i);
      }
      sigma2_pi_T_out(itr - BURN_IN) = sigma2_pi_T;
      sigma2_pi_A_out(itr - BURN_IN) = sigma2_pi_A;
      sigma2_pi_STA_out(itr - BURN_IN) = sigma2_pi_STA;
    }
  }
  List out;
  out("pis") = pis_out;
  out("Yt") = Yt_out;
  out("H") = H_out;
  out("etas") = etas_interaction_out;
  out("sigma2_pi_T") = sigma2_pi_T_out;
  out("sigma2_pi_A") = sigma2_pi_A_out;
  out("sigma2_pi_STA") = sigma2_pi_STA_out;
  out("phi_1") = phi_1_out;
  out("phi_0") = phi_0_out;
  out("lambda_1") = lambda_1_out;
  out("lambda_0") = lambda_0_out;
  return out;
}





