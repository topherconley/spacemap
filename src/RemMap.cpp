#include "RemMap.h"

using namespace Rcpp; //need this for sugar functionality

//constructor
RemMap::RemMap(const int _N, const int _P, const int _Q, 
    const double _lasso, const double _group, const double _tol) :
  N(_N), P(_P), Q(_Q), lasso(_lasso), group(_group), tol(_tol), 
  scaledConv(100), Xnorm2(_P), BetaLasso(_P, _Q), 
  Bnorm(_P, 0), phi(_P, _Q) {
    
    for(int p = 0; p < P; p++) 
      allPredictorList.push_back(p);
  }
//destructor
RemMap::~RemMap() {}

//std::vector<double> norm2(const NumericMatrix X) {

//}

double RemMap::softShrink(double y, double lam) {
  double temp;
  double result;
  if(y>0)
    temp=y-lam;
  else
    temp=-y-lam;
  if(temp<=0)
    result=0;
  else
  {
    result=temp;
    if(y<0)
      result=temp*(-1);
  }
  return result;
}

void RemMap::ell2Norm(const int p, const NumericMatrix & Wm, const bool lassoCase) {
  int q;
  Bnorm.at(p)=0;
  if(lassoCase) {
    for(q=0;q<Q;q++) {
      if(Wm(p,q) > 0)
        Bnorm.at(p) = Bnorm.at(p) + BetaLasso(p,q)*BetaLasso(p,q);
    }
  } else {
    for(q=0;q<Q;q++) {
      if(Wm(p,q) > 0)
        Bnorm.at(p) = Bnorm.at(p) + phi(p,q)*phi(p,q);
    }
  }
  Bnorm.at(p) = std::sqrt(Bnorm.at(p));
  return;
}

double RemMap::maxDiffCoef(const NumericMatrix & phi_last) {
  
  //scaled convergence variables
  double oldNorm = 0;
  double deltaNorm = 0;
  double result = 0;

  for(int p=0;p<P;p++) {
    for(int q=0;q<Q;q++) {
      double delta = phi_last(p,q) - phi(p,q);
      deltaNorm = deltaNorm + delta*delta;
      oldNorm = oldNorm + phi_last(p,q)*phi_last(p,q);
      double temp = std::fabs(delta);
      if(temp > result)
        result = temp;
    }
  }
  //class variable 
  if (oldNorm > 0) 
    scaledConv = std::sqrt(deltaNorm / oldNorm);
  else 
    scaledConv = 1e10;

  return(result);
}

void RemMap::estimateInitial(const NumericMatrix & Xm, const NumericMatrix & Ym,
    const NumericMatrix & Wm) {

  //row-wise norm of Xm required for lasso solution of beta (compute only once).
  int p,q,n;
  double temp, temp1;

  for(p=0;p<P;p++) {
    for(n=0;n<N;n++) {
      Xnorm2.at(p) = Xnorm2.at(p) + Xm(n,p)*Xm(n,p);
    }
  }

  //calculate lasso solution of Beta and its norm.
  for(p=0;p<P;p++) {
    for(q=0;q<Q;q++) {
      ///infinite weight -> not considered in the model
      if(std::isinf(Wm(p,q)))         
        BetaLasso(p,q) = 0;
      else {
        temp= 0;
        for(n=0;n<N;n++)
          temp = temp + Xm(n,p)*Ym(n,q);
        //penalized
        if(Wm(p,q) > 0) {
          BetaLasso(p,q) = softShrink(temp, lasso*Wm(p,q))  /  Xnorm2.at(p);
          //Active set variable to only penalized parameters.
          //Bnorm.at(p) = Bnorm.at(p) + BetaLasso(p,q)*BetaLasso(p,q);
        }
        else 
          BetaLasso(p,q) = temp/Xnorm2.at(p);
      }
    }
    //Bnorm.at(p) = std::sqrt(Bnorm.at(p));
    ell2Norm(p, Wm, true);
  }

  //calculate group lasso solution of BetaLasso.
  for(p=0;p<P;p++) {
    for(q=0;q<Q;q++) {
      // not considered in the model
      if(std::isinf(Wm(p,q))) 
        phi(p,q) = 0;
      /// penalized
      else if(Wm(p,q) > 0 && Bnorm.at(p)> tol) { 	
        temp = Xnorm2.at(p)*Bnorm.at(p);
        temp1 = softShrink(temp, group);
        phi(p,q) = temp1*BetaLasso(p,q)/temp;
      }
      // not penalized
      else
        phi(p,q) = BetaLasso(p,q); 
    }
    //reset Bnorm based on phi
    ell2Norm(p, Wm, false);
  }

  return;
}

void RemMap::predict(Rcpp::NumericMatrix & Yhat, Rcpp::NumericMatrix & Em,
    const Rcpp::NumericMatrix & Xm, const Rcpp::NumericMatrix & Ym){ 
  int n,q,p;
  double temp;
  for(n=0;n<N;n++) {
    for(q=0;q<Q;q++) {
      temp = 0;
      for(p=0;p<P;p++)
        temp = temp + phi(p,q)*Xm(n,p);
      //add whatever was passed in (zero in this case) 
      Yhat(n,q) = Yhat(n,q) + temp;
      Em(n,q) = Ym(n,q) - Yhat(n,q);
    }
  }
  return;
}

void RemMap::updateActiveList() {
  activeList.clear();
  nonactiveList.clear();
  for (int p = 0; p < P; p++) {
    if (Bnorm.at(p) > tol)
      activeList.push_back(p);
    else 
      nonactiveList.push_back(p);
  }
  return;
}


void RemMap::updateGeneral(NumericMatrix & Em, const NumericMatrix & Xm, 
    const NumericMatrix & Wm, NumericMatrix & phi_old, int & nIter, 
    const std::list<int> & predictorList) {

  std::list<int>::const_iterator it;
  std::list<int>::const_iterator lastIndex = predictorList.end();
  int n,p,q;
  double temp, temp1;

  //update penalized coefficients 
  for (it = predictorList.begin(); it != lastIndex; it++) {
    nIter += 1;
    //dereference only once, instead of many times
    p = *it;

    // calculate lasso solution Beta
    for(q=0;q<Q;q++) {
      ///infinite weight -> not considered in the model
      if(std::isinf(Wm(p,q))) {         
        BetaLasso(p,q) = 0;
      } else {
        temp=0;
        for(n=0;n<N;n++)
          temp = temp + Em(n,q)*Xm(n,p);
        temp1 = temp + phi_old(p,q)*Xnorm2.at(p);
        //positive weight -> penalized
        if(Wm(p,q) > 0)
          BetaLasso(p,q) = softShrink(temp1, lasso*Wm(p,q)) / Xnorm2.at(p);
        else
          BetaLasso(p,q) = temp1 / Xnorm2.at(p);
      }
    }

    //reset norm based on lasso solution Beta
    ell2Norm(p, Wm, true);   

    //update phi where it depends on group penalization (BetaNormLasso)
    for(q=0;q<Q;q++) {
      if(std::isinf(Wm(p,q))) {
        phi(p,q) = 0;
      }
      else if(Wm(p,q) > 0 && Bnorm.at(p) > tol) {
        temp = Xnorm2.at(p)*Bnorm.at(p);
        phi(p,q) = BetaLasso(p,q)*softShrink(temp, group)/temp;
      } else 
        phi(p,q) = BetaLasso(p,q);
    }    

    //update residue
    for(q=0;q<Q;q++) {
      for(n=0;n<N;n++) {
        Em(n,q) = Em(n,q) + (phi_old(p,q) - phi(p,q))*Xm(n,p);
      }
    }

    //update phi_old
    ///i am not sure why????????????
    for(q=0;q<Q;q++)
      phi_old(p,q) = phi(p,q);

    //reset norm based on phi (for active set) 
    ell2Norm(p, Wm, false);   
  }
  return;
}

void RemMap::updateNonpenalized(NumericMatrix & Em, const NumericMatrix & Xm, 
    const NumericMatrix & Wm, NumericMatrix & phi_old, int & nIter, 
    const std::list<int> & predictorList) {

  int n,p,q;
  double temp;
 
  std::list<int>::const_iterator it;
  std::list<int>::const_iterator lastIndex = predictorList.end();

  //update coefficients not penalized within non-active list
  for (it = predictorList.begin(); it != lastIndex; it++) {
    //dereference only once, instead of many times
    p=*it;
    for(q=0; q<Q; q++) {
      //non-penalized coefficients
      if(Wm(p,q) == 0) {
        temp=0;
        for(n=0;n<N;n++)
          temp = temp + Em(n,q)*Xm(n,p);
        phi(p,q) = temp/Xnorm2.at(p) + phi_old(p,q);	
        //update residue
        for(n=0; n<N;n++)
          Em(n,q) = Em(n,q) + (phi_old(p,q) - phi(p,q))*Xm(n,p);
        //update phi_old
        phi_old(p,q) = phi(p,q);                 
        nIter += 1;
      } /// end if
    }
  }
  return;
}

double RemMap::updateCoef(bool active, NumericMatrix & Em, const NumericMatrix & Xm,
    const NumericMatrix & Wm, int & nIter) {

  NumericMatrix phi_last(clone(phi));
  NumericMatrix phi_old(clone(phi));

  if(active) {
    updateGeneral(Em, Xm, Wm, phi_old, nIter, activeList);
    updateNonpenalized(Em, Xm, Wm, phi_old, nIter, nonactiveList);
  } else {
    updateGeneral(Em, Xm, Wm, phi_old, nIter, allPredictorList);
  }
  return maxDiffCoef(phi_last);
}

NumericMatrix RemMap::reportCoef() {
  return phi;
}

double RemMap::reportScaledConv() {
  return scaledConv;
}
