// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include<math.h>
using namespace Rcpp;

//DESCRIPTION: 
//Function(1): Calculate the log likelihood of the multivariate gaussian 

// [[Rcpp::export]]
double gaussLogLik(const arma::mat & Theta, const arma::mat & Y) {
  //Convert to concentration matrix
  //arma::mat Theta = parcor2conc(R,d);
  //sample size
  //double n = (double)Y.n_rows;
  
  //log determinant of concentration matrix
  double logDetTheta;
  double sign;
  arma::log_det(logDetTheta, sign, Theta);
  //log likelihood (up to a constant)
  double logLik =   0.5*logDetTheta - 0.5*arma::trace ( Theta * cov(Y) );   
  return logLik;
}
 
