
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

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