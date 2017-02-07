//
//  Space.h
//  
//
//  Created by Christopher Conley on 11/10/2014.
//
//


#ifndef ____Space__
#define ____Space__

#include <iostream>
#include <vector>  //std::vector
#include <utility> //std::pair
#include <list>    //std::list
#include <Rcpp.h>
#include <math.h> 

class Space {
  private:
    //INSTANCE VARIABLES 
    //Sample size
    const int N;
    //dimension in response space (Y)
    const int Q; 
    //Lasso Penalty
    const double lasso;
    //Ridge Penalty
    const double ridge;
    //Convergence Tolerance
    const double tol;
    //Scaled convergence
    double scaledConv;

    Rcpp::NumericMatrix Bs;
    Rcpp::NumericMatrix B;
    std::vector<double> beta_new; 
    std::vector<double> beta_old;
    std::vector<double> beta_last;

    int cur_i,cur_j;
    int change_i, change_j;
    double beta_change;
    

    //The active set of variables
    std::list<std::pair <int,int> > activeList;
    //All the predictors in a list format
    std::list<std::pair <int,int> > allPredictorList;

    //FUNCTIONS
    void updateCurrentBeta(Rcpp::NumericMatrix & Em,
        const Rcpp::NumericMatrix & Ym); 

    double maxDiffCoef();

  public:
    //constructor
    Space(const int _N, const int _Q, const double _lasso, const double _ridge,
        const double _tol);
    //destructor
    ~Space();

    void estimateInitial(const Rcpp::NumericMatrix & Ym, 
        const Rcpp::NumericVector & sigma_sr);

    void initBeta(const Rcpp::NumericMatrix & Ym, 
        const Rcpp::NumericVector & sigma_sr, const std::vector<double> & beta_init);


    void predict(Rcpp::NumericMatrix & Yhat, Rcpp::NumericMatrix & Em,
        const Rcpp::NumericMatrix & Ym);

    bool updateOneBeta(Rcpp::NumericMatrix & Em,
        const Rcpp::NumericMatrix & Ym);

    void updateActiveList();

    int getActiveSize();
    
    double updateCoef(bool active, Rcpp::NumericMatrix & Em, 
        const Rcpp::NumericMatrix & Ym, int & nIter);

    Rcpp::NumericVector reportCoef();
    double reportScaledConv();

};


#endif /* defined(____Space__) */
