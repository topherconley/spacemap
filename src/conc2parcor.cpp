#include <Rcpp.h>
using namespace Rcpp;

// Convert a concentration matrix to a partial correlation matrix. 

// [[Rcpp::export]]
NumericMatrix conc2parcor(const NumericMatrix & D) {
  int M = D.nrow();
  NumericMatrix Rho(M, M); 
  for(int i = 0; i < M; i++) {
    for(int j = 0; j < M; j++) {
      if (D(i,i) == 0 || D(j,j) == 0) {
        Rho(i,j) = NA_REAL;
      } else {
        Rho(i,j) = -1 * D(i,j) / sqrt(D(i,i) * D(j,j));
      }
    }
  }
  return Rho;
}
