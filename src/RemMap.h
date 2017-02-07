//
//  RemMap.h
//  
//
//  Created by Christopher Conley on 11/10/2014.
//
//


#ifndef ____RemMap__
#define ____RemMap__

#include <iostream>
#include <vector>
#include <list> 
#include <Rcpp.h>
#include <math.h> 

class RemMap {
  private:
    //INSTANCE VARIABLES 
    //Sample size
    const int N;
    //Number of dimensions in predictor space (X)
    const int P;
    //Number of dimensions in response space (Y)
    const int Q; 
    //Lasso Penalty
    const double lasso;
    //Group Penalty
    const double group;
    //Convergence Tolerance
    const double tol;
    //scaled relative convergence criteria
    double scaledConv;
    //normed rows of the design matrix (P x 1 vector) 
    std::vector<double> Xnorm2;
    //The lasso solution for Beta (previously Beta, remMap v0.1). 
    Rcpp::NumericMatrix BetaLasso;
    //The norm of Beta(intermediate lasso, and final) 
    std::vector<double> Bnorm;
    //The current estimate of the target variable
    Rcpp::NumericMatrix phi;
    //The active set of variables
    std::list<int> activeList;
    //The non-active set of variables (needed to update non-penalized)
    std::list<int> nonactiveList;
    //All the predictors in a list format
    std::list<int> allPredictorList;
    
    //FUNCTIONS
    double softShrink(double y, double lam);
    void ell2Norm(const int p, const Rcpp::NumericMatrix & Wm, const bool lassoCase);
    double maxDiffCoef(const Rcpp::NumericMatrix & phi_last);
    void updateGeneral(Rcpp::NumericMatrix & Em, const Rcpp::NumericMatrix & Xm, 
        const Rcpp::NumericMatrix & Wm, Rcpp::NumericMatrix & phi_old, int & nIter, 
        const std::list<int> & predictorList);
    void updateNonpenalized(Rcpp::NumericMatrix & Em, const Rcpp::NumericMatrix & Xm, 
        const Rcpp::NumericMatrix & Wm, Rcpp::NumericMatrix & phi_old, int & nIter, 
        const std::list<int> & predictorList);


  public:
    //constructor
    RemMap(const int _N, const int _P, const int _Q, const double _lasso, const double _group,
        const double _tol);
    //destructor
    ~RemMap();

    void estimateInitial(const Rcpp::NumericMatrix & Xm, 
        const Rcpp::NumericMatrix & Ym, const Rcpp::NumericMatrix & Wm);

    void predict(Rcpp::NumericMatrix & Yhat, Rcpp::NumericMatrix & Em,
        const Rcpp::NumericMatrix & Xm, const Rcpp::NumericMatrix & Ym);

    void updateActiveList();

    double updateCoef(bool active, Rcpp::NumericMatrix & Em, 
        const Rcpp::NumericMatrix & Xm,
        const Rcpp::NumericMatrix & Wm, int & nIter);

    Rcpp::NumericMatrix reportCoef();

    double reportScaledConv();
};


#endif /* defined(____RemMap__) */
