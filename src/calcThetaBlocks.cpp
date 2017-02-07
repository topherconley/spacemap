#include <Rcpp.h>
using namespace Rcpp;

//DESCRIPTION: Estimate the blocks XY, YY of the concentration matrix of Y|X ~ Normal

/*
   //calcThetaXY: 
   Following similar index notation of RemMap

   Predictor variables, X, indices: i = 1,...,P
   Response variable, Y, indices: j,k = 1,...,Q 

   Gamma is a P X Q coefficient matrix from remMap(...) estimation.
   
   Using Lemma 1 from space paper, see section 2.1 first and second paragraphs.
   Considering for a moment that X's are a response like the Y's in the joint (X,Y) 
   partial correlation estimation of space. 

   Gamma(i,j) = - ThetaXY(i,j) / ThetaYY(j,j);

   Therefore: 

   ThetaXY(i,j) = - Gamma(i,j) * ThetaYY(j,j);

 */

// [[Rcpp::export]]
NumericMatrix calcThetaXY(const NumericMatrix & Gamma, const NumericVector & diagThetaYY) {

  const int P = Gamma.nrow();
  const int Q = Gamma.ncol();

  NumericMatrix ThetaXY(P,Q);
  for(int i = 0; i < P; i++) {
    for(int j = 0; j < Q; j++) {
      ThetaXY(i,j) = -Gamma(i,j) * diagThetaYY(j);
    }
  }
  return ThetaXY;
}


/*
   //calcThetaYY: 

   1. ThetaYY is the Q X Q YY block of the concentration matrix Y|X. 
   2. rho is the partial correlation between Y_j and Y_k give the rest of 
   the Y's and X's.  
   
   Response variable, Y, indices: j,k = 1,...,Q 

   By lemma 1 from space paper: 

   Beta(j,k) = - ThetaYY(j,k) / ThetaYY(j,j) = rho(j,k) * sqrt ( ThetaYY(k,k) / ThetaYY(j,j) );
   
   Therefore: 
   
   ThetaYY(j,k) = - Beta(j,k) * ThetaYY(j,j); 
   
   Substituting for Beta ...
   ThetaYY(j,k) = - rho(j,k) * sqrt ( ThetaYY(k,k) / ThetaYY(j,j) ) * ThetaYY(j,j); 
   ThetaYY(j,k) = - rho(j,k) * sqrt( ThetaYY(k,k) * ThetaYY(j,j) );

*/

// [[Rcpp::export]]
NumericMatrix calcThetaYY(const NumericMatrix & ParCor, const NumericVector & diagThetaYY) {

  const int Q = ParCor.ncol();
  NumericMatrix ThetaYY(Q,Q);

  //Symmetric matrix because ParCor(j,k) = Parcor(k,j). Iter over diagonal and upper triangle.
  for(int j = 0; j < Q; j++) {
    for(int k = j; k < Q; k++) {
      if(j == k) 
        ThetaYY(j,k) = diagThetaYY(j);
      else {
        ThetaYY(j,k) = -ParCor(j,k) * sqrt( diagThetaYY(k) * diagThetaYY(j) );
        //symmetry
        ThetaYY(k,j) = ThetaYY(j,k);
      }
    }
  }
  return ThetaYY;
}
