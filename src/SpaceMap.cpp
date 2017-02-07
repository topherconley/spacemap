#include <Rcpp.h>
#include "RemMap.h"
#include "Space.h"
#include <iostream>
#include <math.h>

using namespace Rcpp;

double residualSS(const NumericMatrix & Em) {

  int N = Em.nrow();
  int Q = Em.ncol();
  double rss = 0;
  for (int n = 0; n < N; n++) {
    for (int q = 0; q < Q; q++) {
      rss = rss + Em(n,q)*Em(n,q);
    }
  }
  return rss;
}
 
// [[Rcpp::export]]
List doRemMap(const NumericMatrix & Ym, const NumericMatrix & Xm,
    const NumericMatrix & Wm, const double rlasso, const double rgroup,
    const double tol, const int maxIter) {

  RemMap rmap(Xm.nrow(), Xm.ncol(), Ym.ncol(), rlasso, rgroup, tol);
  rmap.estimateInitial(Xm, Ym, Wm);
  //residual 
  NumericMatrix Em(Ym.nrow(), Ym.ncol()); 
  //predicted value (N X Q) 
  NumericMatrix Yhat(Ym.nrow(), Ym.ncol()); 
  //Efficiently predict respective parameters and update common residual
  rmap.predict(Yhat, Em, Xm, Ym);

  bool active = true;
  int nIter = 0;
  double maxDiff ,aMaxDiff;
  maxDiff = 100; 
  while (maxDiff > tol && nIter < maxIter) {
    //get active set 
    rmap.updateActiveList();
    aMaxDiff = 100;
    while (aMaxDiff > tol && nIter < maxIter) {
      aMaxDiff = rmap.updateCoef(active, Em, Xm, Wm, nIter); 
      //Rcout << "aMaxDiff: " << aMaxDiff << std::endl;
      //Rcout << "nIter: " << nIter << std::endl;
    }
    maxDiff = rmap.updateCoef(!active, Em, Xm, Wm, nIter); 
  }

  List ret;
  //Issue warning for lack of convergence. 
  if (maxDiff > tol) {
    ret["convergence"] = wrap(false);
    Rcout << "------------WARNING--------------" << std::endl;
    Rcout << "No RemMap convergence." << std::endl; 
    Rcout << "rlasso: " << rlasso << "rgroup: " << rgroup << std::endl;
    Rcout << "cd_iter: " << maxIter << std:: endl; 
    Rcout << "Consider raising penalty, convergence tolerance, or maximum iterations." << std::endl;
  } else {
    ret["convergence"] = wrap(true);
  }

  ret["Gamma"] = rmap.reportCoef(); 
  ret["rss"] = residualSS(Em);
  ret["deltaMax"] = maxDiff;
  return ret;
 
}

// [[Rcpp::export]]
List doSpace(const NumericMatrix & Ym, NumericVector & sigma_sr, const double lasso, 
    const double ridge, const double tol, const int maxIter, 
    const NumericVector & beta_init, const bool init) {

  Space space(Ym.nrow(), Ym.ncol(), lasso, ridge, tol);
  if (init) { 
    std::vector<double> beta_init2 = as< std::vector<double> >(beta_init);
    space.initBeta(Ym, sigma_sr, beta_init2);
  } else { 
    space.estimateInitial(Ym, sigma_sr);
  }
  
  //residual 
  NumericMatrix Em(Ym.nrow(), Ym.ncol()); 
  //predicted value (N X Q) 
  NumericMatrix Yhat(Ym.nrow(), Ym.ncol()); 
  //Efficiently predict respective parameters and update common residual
  space.predict(Yhat, Em, Ym);

  //Expect to converge if penalty is so high that no edges detected. 
  bool converged = space.updateOneBeta(Em,Ym);
  double maxDiff = tol  + 100;
  if (!converged)  { 
    bool active = true;
    int nIter = 0;
    while (maxDiff > tol && nIter < maxIter) {
      //get active set 
      space.updateActiveList();
      double aMaxDiff = 100;
      while (aMaxDiff > tol && nIter < maxIter) {
        //Rcout << "Active List Size: " << space.getActiveSize() << std::endl;  
        aMaxDiff = space.updateCoef(active, Em, Ym, nIter); 
      }
      maxDiff = space.updateCoef(!active, Em, Ym, nIter); 
    }
  } else {
    Rcout << "Warning: Converged at initial step. Likely too high of a penalty parameter." << std::endl;  
    maxDiff = tol - 0.01;
  }

  List ret;
  //Issue warning for lack of convergence. 
  if (maxDiff > tol) {
    ret["convergence"] = wrap(false);
    Rcout << "------------WARNING--------------" << std::endl;
    Rcout << "No Space convergence." << std::endl; 
    Rcout << "slasso: " << lasso << std::endl;
    Rcout << "cd_iter: " << maxIter << std:: endl; 
    Rcout << "Consider raising penalty, convergence tolerance, or maximum iterations." << std::endl;
  } else {
    ret["convergence"] = wrap(true);
  }

  ret["rho"] = wrap(space.reportCoef()); 
  ret["rss"] = residualSS(Em);
  ret["deltaMax"] = maxDiff;
  return ret;
}

// [[Rcpp::export]]
List doSpaceMap(const NumericMatrix & Ym, const NumericMatrix & Xm,
    const NumericMatrix & Wm, NumericVector & sigma_sr, const double slasso, const double sridge, 
    const double rlasso, const double rgroup, const double tol, const int maxIter) {

  RemMap rmap(Xm.nrow(), Xm.ncol(), Ym.ncol(), rlasso, rgroup, tol);
  Space space(Ym.nrow(), Ym.ncol(), slasso, sridge, tol);

  rmap.estimateInitial(Xm, Ym, Wm);
  space.estimateInitial(Ym, sigma_sr);

  //residual 
  NumericMatrix Em(Ym.nrow(), Ym.ncol()); 
  //predicted value (N X Q) 
  NumericMatrix Yhat(Ym.nrow(), Ym.ncol()); 
  //Efficiently predict respective parameters and update common residual
  rmap.predict(Yhat, Em, Xm, Ym);
  space.predict(Yhat, Em, Ym);
  
  /*Identify initial active predictors. 
      o  If active predictors exist
            initialize coordinate descent. 
         else 
            Report as converged. 
  */
  bool converged = space.updateOneBeta(Em,Ym);

  //initialize to be larger than convergence tolerance. 
  double maxDiff = tol + 100;

  if (!converged) {
    bool active = true;
    int nIter = 0;
    double aMaxDiff, rMaxDiff, sMaxDiff;
    while (maxDiff > tol && nIter < maxIter) {
      //get active set 
      rmap.updateActiveList();
      space.updateActiveList();
      aMaxDiff = 100;
      while (aMaxDiff > tol && nIter < maxIter) {
        double arMaxDiff = rmap.updateCoef(active, Em, Xm, Wm, nIter); 
        double asMaxDiff = space.updateCoef(active, Em, Ym, nIter); 
        aMaxDiff = std::max(arMaxDiff, asMaxDiff);

      }
      rMaxDiff = rmap.updateCoef(!active, Em, Xm, Wm, nIter); 
      sMaxDiff = space.updateCoef(!active, Em, Ym, nIter); 
      maxDiff = std::max(rMaxDiff, sMaxDiff);
    }
  } else {
    Rcout << "Warning: Converged at initial step. Likely too high of a penalty parameter." << std::endl;  
    maxDiff = tol - 0.01;
  }

  List ret;
  //Issue warning for lack of convergence. 
  if (maxDiff > tol) {
    ret["convergence"] = wrap(false);
    ret["convergence"] = wrap(false);
    Rcout << "------------WARNING--------------" << std::endl;
    Rcout << "No Spacemap convergence." << std::endl; 
    Rcout << "slasso: " << slasso << "rlasso: " << rlasso << "rgroup: " << rgroup << std::endl;
    Rcout << "cd_iter: " << maxIter << std:: endl; 
    Rcout << "Consider raising penalty, convergence tolerance, or maximum iterations." << std::endl;
  } else {
    ret["convergence"] = wrap(true);
  }

  ret["Gamma"] = rmap.reportCoef(); 
  ret["rho"] =  space.reportCoef();
  ret["rss"] = residualSS(Em);
  ret["deltaMax"] = maxDiff;
  return ret; 
}
