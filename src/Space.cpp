#include "Space.h"

using namespace Rcpp; //need this for sugar functionality

//constructor
Space::Space(const int _N, const int _Q, 
    const double _lasso, const double _ridge, const double _tol) :
  N(_N), Q(_Q), lasso(_lasso), ridge(_ridge), tol(_tol), scaledConv(100), 
  Bs(_Q, _Q), B(_Q, _Q), beta_new(_Q*_Q, 0), 
  beta_old(_Q*_Q, 0), beta_last(_Q*_Q, 0) {
     
    //initialize all \rho indices (upper triangle)   
    for(int i = 0; i < Q-1 ; i++)
      for (int j = i + 1; j < Q ; j++)
        allPredictorList.push_back(std::make_pair(i,j));
}
//destructor
Space::~Space() {}

void Space::estimateInitial(const NumericMatrix & Ym, const NumericVector & sigma_sr) { 
  int i,j,k;
  double temp, temp1;
  //Only needed for initialization
  NumericMatrix A(Q,Q);
  //normed rows of the response matrix (Q x 1 vector) 
  std::vector<double> normx(Q,0);

  //column wise 
  for(j=0;j<Q;j++)  
    for(i=0;i<N;i++) 
      normx.at(j) = normx.at(j)+ Ym(i,j)*Ym(i,j);

  for(i=0;i<Q;i++) 
    for(j=0;j<Q;j++)
      B(i,j) = (sigma_sr[i]/sigma_sr[j]);

  for(i = 0;i<Q;i++) {
    for(j=0;j<(i+1);j++) {
      Bs(i,j) = B(i,j)*B(i,j)*normx.at(i) +  B(j,i)*B(j,i)*normx.at(j);
      Bs(j,i) = Bs(i,j);
    }
  }

  for(i = 0;i<Q;i++) {
    for(j=0;j<(i+1);j++) {
      temp=0;
      for(k = 0; k < N ; k++ ) 
        temp = temp + Ym(k,i)*Ym(k,j);
      A(i,j)=temp*B(j,i);
      A(j,i)=temp*B(i,j);
    }
  }

  for(i = 0;i<Q;i++) {
    for(j=0;j<(i+1);j++) {
      temp1 =   A(i,j) + A(j,i);
      if ( temp1 > 0 )
        temp = temp1 - lasso;
      else
        temp = - temp1 - lasso;
      if(temp < 0 )
        beta_new.at(i*Q+j) = 0;
      else {
        beta_new.at(i*Q+j) =  temp /(Bs(i,j)*(1+ridge));
        if ( temp1 < 0)
          beta_new.at(i*Q+j) = (-1) * beta_new.at(i*Q+j);
      }
      beta_new.at(j*Q+i)=beta_new.at(i*Q+j);
    }
  }

  //not modeling self-regulation
  for(i = 0;i<Q;i++) beta_new.at(i*Q+i) = 0 ; 

  return;
}

void Space::initBeta(const NumericMatrix & Ym, const NumericVector & sigma_sr, 
    const std::vector<double> & beta_init) {

  int i,j;
  //normed rows of the response matrix (Q x 1 vector) 
  std::vector<double> normx(Q,0);

  //column wise 
  for(j=0;j<Q;j++)  
    for(i=0;i<N;i++) 
      normx.at(j) = normx.at(j)+ Ym(i,j)*Ym(i,j);

  for(i=0;i<Q;i++) 
    for(j=0;j<Q;j++)
      B(i,j) = (sigma_sr[i]/sigma_sr[j]);

  for(i = 0;i<Q;i++) {
    for(j=0;j<(i+1);j++) {
      Bs(i,j) = B(i,j)*B(i,j)*normx.at(i) +  B(j,i)*B(j,i)*normx.at(j);
      Bs(j,i) = Bs(i,j);
    }
  }

  beta_new = beta_init;
  //not modeling self-regulation
  for(i = 0;i<Q;i++) beta_new.at(i*Q+i) = 0 ; 

  return;
}


void Space::predict(Rcpp::NumericMatrix & Yhat, Rcpp::NumericMatrix & Em,
    const Rcpp::NumericMatrix & Ym) {

  int k,j,i;

  for(k=0; k<N; k++) {
    for(j=0; j<Q; j++) {
      for(i=0; i<Q; i++)
        Yhat(k,j) = Yhat(k,j) + Ym(k,i)*beta_new.at(i*Q+j)*B(i,j);
      Em(k,j) = Ym(k,j) - Yhat(k,j);
    }
  }
  return;
}

bool Space::updateOneBeta(Rcpp::NumericMatrix & Em,
        const Rcpp::NumericMatrix & Ym) {

  int j,i; 

  //Indicate if converged at the initial step
  bool converged = false;

  activeList.clear();
  //find first coefficient in active set.
  for(j = Q-1; j>=1; j--) {
    for(i = j-1; i>=0; i--) {
      if( beta_new.at(i*Q+j) > tol || beta_new.at(i*Q+j) < -tol ) {
        activeList.push_back(std::make_pair(i,j));
        break;
      }
    }
    if (activeList.size() > 0)
      break;
  }

  //converge at initial step?
  if (activeList.size() == 0) {
    converged = true;
    return converged;
  }

  //clone beta
  beta_old = beta_new;
  
  cur_i = activeList.front().first;
  cur_j = activeList.front().second;
  updateCurrentBeta(Em,Ym);
  
  return converged;
}

void Space::updateCurrentBeta(Rcpp::NumericMatrix & Em,
        const Rcpp::NumericMatrix & Ym) {
  int k;
  double temp, temp1, beta_next;
  double Aij = 0;
  double Aji = 0;

  for(k=0; k<N; k++) {
    Aij = Aij + Em(k,cur_j)*Ym(k,cur_i);
    Aji = Aji + Em(k,cur_i)*Ym(k,cur_j);
  }
  Aij = Aij*B(cur_i,cur_j);
  Aji = Aji*B(cur_j,cur_i);

  beta_next = (Aij + Aji)/Bs(cur_i,cur_j) + beta_old.at(cur_i*Q+cur_j);

  ///shrink beta_next
  temp1 = beta_next;
  if ( beta_next > 0 )
    temp = beta_next - lasso/Bs(cur_i,cur_j);
  else
    temp = - beta_next - lasso/Bs(cur_i,cur_j);
  if(temp < 0)
    temp = 0;
  else
  {
    temp =  temp /(1+ridge);
    if(temp1 < 0 )
      temp = (-1) * temp;
  }

  beta_change = beta_old.at(cur_i*Q+cur_j) - temp;

  beta_new.at(cur_i*Q+cur_j) = temp;
  beta_new.at(cur_j*Q+cur_i) = temp;

  change_i=cur_i;
  change_j=cur_j;

  return;
}

void Space::updateActiveList() {
 
  //clear of previous active list 
  activeList.clear();

  int j,i;
  for(j = Q-1; j>=1; j--) {
    for(i = j-1; i>=0; i--) {
      if( beta_new.at(i*Q+j) > tol || beta_new.at(i*Q+j) < -tol ) {
        activeList.push_back(std::make_pair(i,j));
      }
    }
  }
  return;
}

int Space::getActiveSize() {
  return activeList.size();
}

double Space::maxDiffCoef() {

 //scaled convergence variables
  double oldNorm = 0;
  double deltaNorm = 0;
  double maxdif=-100;

  for(int i = 0;i<Q;i++) {
    for(int j=0;j<Q;j++) {
      double delta = beta_last.at(i*Q+j) - beta_new.at(i*Q+j);
      deltaNorm = deltaNorm + delta*delta;
      oldNorm = oldNorm + beta_last.at(i*Q+j)*beta_last.at(i*Q+j);
      double temp = std::fabs(delta);
      if( temp > maxdif )
        maxdif = temp;
    }
  }
 //class variable 
  if (oldNorm > 0) 
    scaledConv = std::sqrt(deltaNorm / oldNorm);
  else 
    scaledConv = 1e10;


  return maxdif;
}

double Space::updateCoef(bool active, NumericMatrix & Em, 
    const NumericMatrix & Ym, int & nIter) {

    //required for finding max difference in update
  beta_last = beta_new;

  std::list<std::pair <int,int> >::iterator it;
  std::list<std::pair <int,int> >::iterator firstIndex;
  std::list<std::pair <int,int> >::iterator lastIndex;

//a reference to either full parameters or active set
  if (active) {
    firstIndex = activeList.begin();
    lastIndex = activeList.end();
  } else {
    firstIndex = allPredictorList.begin();
    lastIndex = allPredictorList.end();
  }

  int k;
  for (it = firstIndex; it != lastIndex; it++) {

    cur_i = it->first;
    cur_j = it->second;

    beta_old.at(change_i*Q+change_j) = beta_new.at(change_i*Q+change_j);
    beta_old.at(change_j*Q+change_i) = beta_new.at(change_j*Q+change_i);

    // Update Residue   
    for(k=0; k<N; k++) {
      Em(k,change_i) = Em(k,change_i) + Ym(k,change_j)*beta_change*B(change_j,change_i);
      Em(k,change_j) = Em(k,change_j) + Ym(k,change_i)*beta_change*B(change_i,change_j);
    }

    updateCurrentBeta(Em,Ym);
    nIter += 1;
  }

  return maxDiffCoef();
}

NumericVector Space::reportCoef() {
  return wrap(beta_new);
}

double Space::reportScaledConv() {
  return scaledConv;
}
