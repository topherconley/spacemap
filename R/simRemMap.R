#'@title Simulate remMap scenario.
#'@description Generate data suitable for remMap. 
#'@details Wrapper around example from remMap documentation.
simRemMap <- function(n = 100, p = 300, q = 300, seed = 1) {

  set.seed(seed)
  ### generate X matrix
  rho=0.5; Sig<-matrix(0,p,p)
  for(i in 2:p){ for(j in 1: (i-1)){
    Sig[i,j]<-rho^(abs(i-j))
    Sig[j,i]<-Sig[i,j]
  }}
  diag(Sig)<-1
  R<-chol(Sig)
  X.m<-matrix(rnorm(n*p),n,p)
  X.m<-X.m%*%R
  
  ### generate coefficient
  coef.m<-matrix(0,p,q)
  hub.n=20
  hub.index=sample(1:p, hub.n)
  for(i in 1:q){
    cur=sample(1:3,1)
    temp=sample(hub.index, cur)
    coef.m[temp,i]<-runif(length(temp), min=2, max=3)
  }
  
  ### generate responses
  E.m<-matrix(rnorm(n*q),n,q)
  Y.m<-X.m%*%coef.m+E.m
  
  list(X.m = X.m, Y.m = Y.m, E.m = E.m, 
       coef.m = coef.m)
}


