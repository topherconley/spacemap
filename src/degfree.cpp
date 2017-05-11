#include <Rcpp.h>
#include <cmath> 
using namespace Rcpp;

//' @title Count the number of x-y edges
//' 
//' @description Count the number of x-y edges after fitting spaceMap model. 
//' 
//' @param Gamma The estimated regression coefficient matrix, \eqn{\hat\Gamma}, 
//' output from \code{\link{spacemap}}.
//' @param aszero Positive numeric value indicating at what point to consider extremely 
//' small parameter estimates of \eqn{\Gamma} as zero. 
//' @return Number of x-y edges.
//' @export
// [[Rcpp::export]]
unsigned int nonZeroWhole(const NumericMatrix & Gamma, const double aszero) {

  const unsigned int P = Gamma.nrow();
  const unsigned int Q = Gamma.ncol();
  //non-zero count
  unsigned int nzCount = 0; 
  for(unsigned int i = 0; i < P; i++) {
    for(unsigned int j = 0; j < Q; j++) {
      if (std::abs(Gamma(i,j)) > aszero)  
        nzCount = nzCount + 1;
    }
  }
  return nzCount;
}

//' @title Count the number of y-y edges
//' 
//' @description Count the number of y-y edges after fitting spaceMap model. 
//' 
//' @param ParCor The estimated partial correlation matrix, \eqn{\hat\rho}, 
//' output from \code{\link{spacemap}}. Assumed to be symmetric.
//' @param aszero Positive numeric value indicating at what point to consider extremely 
//' small parameter estimates of \eqn{\rho} as zero. 
//' @return Number of y-y edges.
//' @export 
// [[Rcpp::export]]
unsigned int nonZeroUpper(const NumericMatrix & ParCor, const double aszero) {

  const unsigned int Q = ParCor.ncol();
  //non-zero count
  unsigned int nzCount = 0; 
  //Symmetric matrix because ParCor(j,k) = Parcor(k,j). Iter over upper triangle.
  for(unsigned int j = 0; j < Q; j++) {
    for(unsigned int k = j + 1; k < Q; k++) {
      if (std::abs(ParCor(j,k)) > aszero)  
        nzCount = nzCount + 1;
    }
  }
  return nzCount;
} 
