// [[Rcpp::depends(RcppArmadillo)]]
#include <algorithm>
#include <vector>
#include <RcppArmadillo.h>
using namespace Rcpp;

struct IncGenerator {
    int current_;
    IncGenerator (int start) : current_(start) {}
    int operator() () { return current_++; }
};

//Lederer et al. (2013) "Trust but verify...." methodology. 
bool cLSLasso(const arma::mat & Xnz, const arma::mat & coef, double numCoef, double cthresh) { 
  double Fs = arma::accu ( arma::sign(coef) == arma::sign( arma::inv( Xnz.t()*Xnz ) * coef ) ) / numCoef;
  return Fs < cthresh; 
}

// [[Rcpp::export]]
bool olsRefitSpacemap(const arma::mat & Y, const arma::mat & X,
    arma::mat & ParCor, arma::mat & Gamma, arma::colvec & sigma, arma::colvec & RSS, 
    double tol, double ridge = 0.1) {
  
  //reset RSS to zero 
  RSS.zeros(); 
  const unsigned int Q = ParCor.n_cols;
  const unsigned int N = Y.n_rows;
  arma::mat Beta = arma::zeros<arma::mat>(Q,Q);
  arma::umat olsId = arma::zeros<arma::umat>(Q,Q);
  bool sparse = true; 

  for(unsigned int j = 0; j < Q; j++) {
    //do not regress jth variable on jth variable. 
    ParCor(j,j) = 0;
    //Use only lasso-selected covariates (i.e. those not shrunck to zero).
    arma::uvec nzIdY = find(arma::abs(ParCor.col(j)) > tol);
    unsigned int numNzY = nzIdY.n_elem;
    arma::uvec nzIdX = find(arma::abs(Gamma.col(j)) > tol);
    unsigned int numNzX = nzIdX.n_elem;

    arma::uvec jth(1);  jth.fill(j);
    arma::colvec resid = arma::zeros<arma::colvec>(N);
    arma::colvec coef = arma::zeros<arma::colvec>(N);

    //no predictors selected 
    if(numNzY == 0 && numNzX == 0)  { 
      resid = Y.col(j);
      sigma(j) = 1 / arma::as_scalar( arma::cov( resid ) ) ;
    } //only Y predictors selected  
    else if (numNzY > 0 && numNzX == 0) {
      if (numNzY > N) {
        sparse = false;
        break;
      } 
      //Solving equation: Y(,-j) \beta = Y(,j) 
      arma::mat Ynz = Y.cols(nzIdY);
      arma::mat Ynzt = Ynz.t();
      sparse = arma::solve(coef, Ynzt*Ynz + ridge*arma::eye<arma::mat>(numNzY,numNzY),  Ynzt * Y.col(j) ); 
      if (!sparse)  { 
        break;
      }

/*
      if(!cLSLasso(Ynz, coef, (double) numNzY, cthresh))
        continue;
*/

      resid = Y.col(j) - Ynz*coef;
      sigma(j) = 1 / arma::as_scalar(arma::cov(resid));
     //store OLS-based Beta estimates;
      Beta.submat(nzIdY,jth) = coef;
      //keep track of cases where beta is updated
      olsId.submat(nzIdY,jth) = arma::ones<arma::umat>(numNzY, 1); 
    }
    // only X predictors selected 
    else if (numNzX > 0 && numNzY == 0) {
      if (numNzX > N) {
        sparse = false;
        break;
      } 
      //Solving equation: X(nonzeros,j) \beta = Y(,j) 
      arma::mat Xnz = X.cols(nzIdX);
      arma::mat Xnzt = Xnz.t();
      sparse = arma::solve(coef, Xnzt*Xnz + ridge*arma::eye<arma::mat>(numNzX,numNzX), Xnzt * Y.col(j) ); 
      if (!sparse)  { 
        break;
      }

/*
      if(!cLSLasso(Xnz, coef, (double) numNzX, cthresh))
        continue;
*/
      resid = Y.col(j) - Xnz*coef;
      sigma(j) = 1 / arma::as_scalar(arma::cov(resid));
      //store OLS-based Gamma estimates;
      Gamma.submat(nzIdX,jth) = coef;
    } // Both Y and X selected  
    else if(numNzY > 0 && numNzX > 0 ) {
      if (numNzY + numNzX > N) {
        sparse = false;
        break;
      }
      //combine (Y(-j),X) lasso-selected predictors
      arma::mat nzJointYX = join_horiz( Y.cols(nzIdY), X.cols(nzIdX) );
      arma::mat nzJointYXt = nzJointYX.t();
      //Solving equation: X \beta = Y 
      sparse = arma::solve(coef, 
          nzJointYXt*nzJointYX + ridge*arma::eye<arma::mat>(nzJointYX.n_cols, nzJointYX.n_cols),
          nzJointYXt * Y.col(j) ); 
      if (!sparse)  { 
        break;
      }
/*
      if(!cLSLasso(nzJointYX, coef, (double) numNzY + numNzX, cthresh))
        continue;
*/
      resid = Y.col(j) - nzJointYX*coef;
      // \sigma^ii = 1 / var(e_{i}
      sigma(j) = 1 / arma::as_scalar(arma::cov(resid));

      //index y coefficients with (0,...,numNzY - 1)
      IncGenerator gY (0);
      arma::uvec yCoefId(numNzY);
      std::generate(yCoefId.begin(), yCoefId.end(), gY); // Fill with the result of calling g() repeatedly.

      //index x coefficients with  (numNzY, ..., numNzY + numNzX - 1)
      IncGenerator gX (numNzY);
      arma::uvec xCoefId(numNzX);
      std::generate(xCoefId.begin(), xCoefId.end(), gX); // Fill with the result of calling g() repeatedly.

      //store OLS-based regression coefficient estimates;
      Beta.submat(nzIdY,jth) = coef(yCoefId);
      //keep track of cases where beta is updated
      olsId.submat(nzIdY,jth) = arma::ones<arma::umat>(numNzY,1); 
      Gamma.submat(nzIdX,jth) = coef(xCoefId);
    } 
    RSS(0) += arma::sum(arma::square(resid));
  }
  
  //Now convert Beta-OLS updated to Partial Correlation:
  //Temporarily store the ols-based ParCor in Beta so we can subset from a matrix.
  Beta = sign(Beta) % sqrt( abs ( Beta % Beta.t()  ) );
  arma::uvec olsIdUvec = arma::find(olsId);
  ParCor.elem(olsIdUvec) = Beta.elem(olsIdUvec);
  //Partial correlation with itself is always one 1.
  ParCor.diag() = arma::ones<arma::vec>(Q);
  return sparse;
}


// [[Rcpp::export]]
bool olsRefitSpace(const arma::mat & Y, arma::mat & ParCor, arma::colvec & sigma, arma::colvec & RSS, 
    double tol, double ridge = 0.1) {
  
  //reset RSS to zero 
  RSS.zeros(); 
  //Rcpp::Rcout << "Starting Refit: " << std::endl;
  const unsigned int Q = ParCor.n_cols;
  const unsigned int N = Y.n_rows;
  arma::mat Beta = arma::zeros<arma::mat>(Q,Q);
  //Rcpp::Rcout << "Setting up olsId: " << std::endl;
  arma::umat olsId = arma::zeros<arma::umat>(Q,Q);
  bool sparse = true; 

  for(unsigned int j = 0; j < Q; j++) {
    //Rcpp::Rcout << "Jth: " << j << std::endl;
    //do not regress jth variable on jth variable. 
    ParCor(j,j) = 0;
    //Use only lasso-selected covariates (i.e. those not shrunck to zero).
    arma::uvec nzIdY = find(arma::abs(ParCor.col(j)) > tol);
    unsigned int numNzY = nzIdY.n_elem;
    
    //define resid here to update RSS
    arma::colvec resid = arma::zeros<arma::colvec>(N);
    //define coef here for error handling of no solutio of solve(...)
    arma::colvec coef = arma::zeros<arma::colvec>(N);
    //Rcpp::Rcout << "numNzY: " << numNzY << std::endl;
    //nzIdY.print("nzIdY:");
    //no predictors selected 
    if(numNzY == 0)  { 
     //Rcpp::Rcout << "Zeros: " << std::endl;
      resid = Y.col(j);
      sigma(j) = 1 / arma::as_scalar( arma::cov( resid ) ) ;
      Beta(j,j) = arma::as_scalar(sigma(j));
    } //Y predictors selected  
    else if (numNzY > 0) {
      if (numNzY > N) {
        sparse = false;
        break;
      } 
      //Solving equation: Y(,-j) \beta = Y(,j) 

      //Rcpp::Rcout << "Solving: " << j << std::endl;
      arma::mat Ynz = Y.cols(nzIdY);
      arma::mat Ynzt = Ynz.t();
      sparse = arma::solve(coef, Ynzt*Ynz + ridge*arma::eye<arma::mat>(numNzY,numNzY),  Ynzt * Y.col(j) ); 
      if (!sparse)  { 
        break;
      }
/*
      if(!cLSLasso(Ynz, coef, (double) numNzY, cthresh))
        continue;
*/
      //Rcpp::Rcout << "Residual: " << j << std::endl;
      resid = Y.col(j) - Ynz*coef;
      sigma(j) = 1 / arma::as_scalar(arma::cov(resid));
      //store OLS-based Beta estimates;
      arma::uvec jth(1);  jth.fill(j);
      //Rcpp::Rcout << "Beta.submat call: " << j << std::endl;
      Beta.submat(nzIdY,jth) = coef;
      Beta(j,j) = arma::as_scalar(sigma(j));
      //Beta.submat(nzIdY,jth).print("Submatrix contents.");
      //Rcpp::Rcout << "numNzY: " << numNzY << std::endl;
      //coef.print("Coef:");
      //keep track of cases where beta is updated
      olsId.submat(nzIdY,jth) = arma::ones<arma::umat>(numNzY, 1); 
    }
    RSS(0) += arma::sum(arma::square(resid));
  }
   

  //Rcpp::Rcout << "Exited loop: " << std::endl;
  //Now convert Beta-OLS updated to Partial Correlation:
  //Temporarily store the ols-based ParCor in Beta so we can subset from a matrix.
  Beta = sign(Beta) % sqrt( abs ( Beta % Beta.t()  ) );
  arma::uvec olsIdUvec = arma::find(olsId);
  ParCor.elem(olsIdUvec) = Beta.elem(olsIdUvec);
  //Partial correlation with itself is always one 1.
  ParCor.diag() = arma::ones<arma::vec>(Q);
  
  return sparse;
}

/*
// [[Rcpp::export]]
List fastLm(const arma::vec & y, const arma::mat & X) {

    int n = X.n_rows, k = X.n_cols;
   
    arma::colvec coef = arma::solve(X, y); 
    arma::colvec resid = y - X*coef; 
   
    double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
    arma::colvec stderrest = 
        arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
   
   return List::create(Named("coefficients") = coef,
                       Named("stderr")       = stderrest);
}
*/


