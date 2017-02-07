// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double gaussCondLogLik(const arma::mat & Y, const arma::mat & X, 
    const arma::mat & ThetaYY, const arma::mat & ThetaXY) {
  
  double n = (double)Y.n_rows;
  double logDetThetaYY;
  double sign;
  arma::log_det(logDetThetaYY, sign, ThetaYY);
  
  //ignoring constant of log likelihood
  double term1 = -0.5 * arma::trace( Y * ThetaYY * Y.t() );
  double term2 = -1*arma::trace( X * ThetaXY * Y.t() );
  double term3 = 0.5*n*logDetThetaYY;
  arma::mat InvThetaYY = arma::zeros<arma::mat>(Y.n_cols, Y.n_cols);
  bool invertible = arma::inv(InvThetaYY, ThetaYY);
  double result;
  if (invertible) {
    double term4 = -0.5*arma::trace( X * ThetaXY * InvThetaYY * ThetaXY.t() * X.t() );
    result = term1 + term2 + term3 + term4;
  } else {
    //Make the model unfavorable to select.
    result = arma::datum::inf;
  }
  return result; 
}
