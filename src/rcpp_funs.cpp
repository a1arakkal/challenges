
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat chol_inv(arma::mat& x){
  int n = x.n_rows;
  int p = x.n_cols;
  arma::mat L(n, p);
  
  for(int j = 0; j < p; j++){
    double sum = 0;
    for(int k = 0; k < j; k++){
      sum += L(j,k) * L(j,k);
    }
    L(j,j) = sqrt(x(j,j) - sum);
    
    for(int i = j + 1; i < p; i++){
      sum = 0;
      for (int k = 0; k < j; k++){
        sum += L(i,k) * L (j,k);
      }
      L(i,j) = (1.0 / L(j,j) * (x(i,j) - sum));
    }
  }
  arma::mat L_inv = inv(L);
  return L_inv.t()*L_inv;
}


// [[Rcpp::export]]
double log_like_cpp(arma::colvec& params, arma::mat& X, arma::colvec& y, 
                    arma::colvec& tpts){
  int p = X.n_cols;
  arma::colvec beta = params(arma::span(0, p-1));
  arma::colvec gamma = params(arma::span(p, p+1));
  double b = std::exp(params(p+2));
  double sigma2 = std::exp(params(p+3));
  arma::colvec X_beta = X*beta;
  arma::colvec X_gamma = X*gamma;
  arma::colvec a = X_beta.transform( [](double val) { return std::exp(val); } );
  arma::colvec c = X_gamma.transform( [](double val) { return std::exp(val); } );

  
  arma::colvec C_T = -c % tpts;
  arma::colvec mu_1 = -b * C_T.transform( [](double val) { return std::exp(val); } );
  arma::colvec mu_2 = mu_1.transform( [](double val) { return std::exp(val); } );
  arma::colvec mu = a % mu_2;
  
  arma::colvec y_mu = y-mu;
  arma::colvec y_mu_2 = y_mu.transform( [](double val) { return pow(val, 2); } );
  arma::colvec out_vec = (-1/ double(2))*std::log(sigma2) - (1/double((2*sigma2)))*y_mu_2;
  double out = accu(out_vec);
  return out;
}

// [[Rcpp::export]]
double log_post_cpp(arma::colvec& params, arma::mat& X, arma::colvec& y, 
                    arma::colvec& tpts){

  int p = X.n_cols; 
  arma::colvec beta_gamma = params(arma::span(0, p+1));
  double b = std::exp(params(p+2));
  double sigma2 = std::exp(params(p+3));

  arma::mat sigma_inv = arma::diagmat(arma::ones(p*2));
  arma::colvec sigma_inv_mean = sigma_inv*beta_gamma;
  double p2_prior = -0.5*dot(beta_gamma, sigma_inv_mean);
  double out = log_like_cpp(params, X, y, tpts) + p2_prior + (R::dgamma(b, 1, 1, 1)) + (R::dgamma(sigma2, 1, 1, 1)) + params(p+2) + params(p+3);

  return out;
}


// [[Rcpp::export]]

int MH(arma::mat& MH_draws, arma::mat& proposal_cov,
       arma::mat& X, arma::colvec& y,arma::colvec& tpts){

  int n = MH_draws.n_rows;
  arma::colvec beta_proposal = arma::zeros<arma::colvec>(MH_draws.n_cols);
  int acc_count = 0;
  
  // Obtain environment containing function
  Rcpp::Environment package_env("package:MASS"); 

  // Make function callable from C++
  Rcpp::Function mvtnomrfun = package_env["mvrnorm"]; 
  
  for(int j = 1; j < n; j++){
    MH_draws.row(j) = MH_draws.row(j-1);
    arma::colvec temp = MH_draws.row(j).t();
    beta_proposal = Rcpp::as<arma::colvec>(mvtnomrfun(1, temp, proposal_cov));
    double unif = R::runif(0,1);
    double acc_prob = 0.0;
    int any_neg = 0;

    for (int k = 0; k == 3; k++){
      if(beta_proposal(k) < 0){
        any_neg = any_neg + 1;
      }
    }

    if(any_neg > 0){
      acc_prob = R_NegInf;
    } else {
      acc_prob = std::exp(log_post_cpp(beta_proposal, X, y, tpts) - 
                          log_post_cpp(temp, X, y, tpts));
    }
    
    if(unif < acc_prob){
      MH_draws.row(j) = beta_proposal.t();
      acc_count =  acc_count + 1;
    }
    
  }
    
  return acc_count;
}  