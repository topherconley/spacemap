#include <Rcpp.h>
#include <cmath> 
using namespace Rcpp;

//DESCRIPTION: 

 
// [[Rcpp::export]]
unsigned int nonZeroWhole(const NumericMatrix & Gamma, const double tol) {

  const unsigned int P = Gamma.nrow();
  const unsigned int Q = Gamma.ncol();
  //non-zero count
  unsigned int nzCount = 0; 
  for(unsigned int i = 0; i < P; i++) {
    for(unsigned int j = 0; j < Q; j++) {
      if (std::abs(Gamma(i,j)) > tol)  
        nzCount = nzCount + 1;
    }
  }
  return nzCount;
}

// [[Rcpp::export]]
unsigned int nonZeroUpper(const NumericMatrix & ParCor, const double tol) {

  const unsigned int Q = ParCor.ncol();
  //non-zero count
  unsigned int nzCount = 0; 
  //Symmetric matrix because ParCor(j,k) = Parcor(k,j). Iter over upper triangle.
  for(unsigned int j = 0; j < Q; j++) {
    for(unsigned int k = j + 1; k < Q; k++) {
      if (std::abs(ParCor(j,k)) > tol)  
        nzCount = nzCount + 1;
    }
  }
  return nzCount;
} 
