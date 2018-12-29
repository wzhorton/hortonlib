#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat covprec, bool is_prec){
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  if (is_prec) {
    return sum((x_cen * covprec) % x_cen, 1);
  } else {
    return sum((x_cen * covprec.i()) % x_cen, 1);
  }
}

// [[Rcpp::export]]
arma::vec dmvnorm(arma::mat x,  arma::rowvec mu,  arma::mat covprec, bool log = false, bool is_prec = false) {
  arma::vec distval = Mahalanobis(x,  mu, covprec, is_prec);
  double logdet = sum(arma::log(arma::eig_sym(covprec)));
  if(is_prec) {
    logdet = -1*logdet;
  }
  double log2pi = std::log(2.0 * M_PI);
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2 );

  if (log) {
    return(logretval);
  } else {
    return(exp(logretval));
  }
}
