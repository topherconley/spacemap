# 
# ####BIC
#   temp.size<-sum(abs(coef)>1e-6)
# Sigma.fit.inv<-(-ParCor.fit)*matrix(sqrt(sig.fit),p,p,byrow=FALSE)*matrix(sqrt(sig.fit),p,p,byrow=TRUE) 
#   diag(Sigma.fit.inv)=sig.fit
# bic<-n*(-log(det(Sigma.fit.inv))+sum(diag(Sigma.fit.inv%*%s)))+temp.size*log(n)
# 
# 
# ####BIC with ols refit
# 
#   Sigma.fit.inv<-OLS.SigInv(Y,try.adj)
# bic.ols<-n*(-log(det(Sigma.fit.inv))+sum(diag(Sigma.fit.inv%*%s)))+temp.size*log(n)
# 
#   cur.list<-list("parcor"=temp.cor,"sig"=sig.fit, "summary"=summary.JSRM.fit,"time"=end-start,"bic"=bic,"bic.ols"=bic.ols) # Pei 
#   lasso.joint.ol.deweight.c[[l]]<-cur.list    
# 
# ###


##
OLS.SigInv<-function(Y,try.adj){
  p<-ncol(Y)
  beta<-apply(matrix(1:p),1, OLS.Sol,Y=Y,try.adj=try.adj)
  beta.a<-abs(beta)
  ParCor.fit<-sign(beta)*sqrt(beta.a*t(beta.a))
  diag(ParCor.fit)<-1
  sig.fit<-diag(beta) 
  Sigma.fit.inv<-(-ParCor.fit)*matrix(sqrt(sig.fit),p,p,byrow=FALSE)*matrix(sqrt(sig.fit),p,p,byrow=TRUE) 
  diag(Sigma.fit.inv)=sig.fit
  return(Sigma.fit.inv) 
  
}


####OLS solution after model selection
OLS.Sol<-function(Y,try.adj,i){
  ##print(i)
  p<-ncol(Y)
  n<-nrow(Y)
  
  diag(try.adj)<-0
  Y.c<-matrix(Y[,i],ncol=1)
  index.c<-try.adj[i,]
  
  result<-numeric(p)
  
  if(sum(index.c)==0){
    
    e.cur<-Y.c
    sig.ii<-var(e.cur)
    result[i]<-1/sig.ii 
    return(result)
  }
  
  if(sum(index.c)>n){
    temp<-numeric(p)
    temp[(index.c==1)][1:n]<-1
    index.c<-temp
    
  }
  
  X.c<-matrix(Y[,index.c==1],ncol=sum(index.c))
  beta.c<-solve(t(X.c)%*%X.c)%*%t(X.c)%*%Y.c
  
  result[index.c==1]<-beta.c
  
  e.cur<-Y.c-X.c%*%beta.c
  sig.ii<-var(e.cur)
  result[i]<-1/sig.ii
  
  return(result)
}

rolsRefit <- function(spaceObj, XY, tol) {
  beta <- apply(matrix(1:nrow(spaceObj$ParCor)),1, OLS.Sol,Y=XY,try.adj=abs(spaceObj$ParCor) > tol)
  beta.a<-abs(beta)
  ParCor.fit<-sign(beta)*sqrt(beta.a*t(beta.a))
  diag(ParCor.fit) <- 1
  #calculate rss
  diag(beta) <- 0
  esti.XY=XY%*%beta
  residue=XY-esti.XY
  rss <- sum(residue*residue)
  
  list(ParCor = ParCor.fit, sig.fit = diag(beta), rss = rss)
}

####OLS spacemap solution after model selection
OLS.Sol2<-function(Y, X, ParCor, Gamma, sigma.fit, tol){
  
  Q <- ncol(ParCor)
  N <- nrow(Y)
  P <- nrow(Gamma)
  Beta <- matrix(0.0, Q,Q)
  sparse <- TRUE
  
  for (j in seq_len(Q)) {
    Y.j <- matrix(Y[,j],ncol=1)
    ParCor[j,j] <- 0
    nzIdY <- which(abs(ParCor[,j]) > tol)
    nzIdX <- which(abs(Gamma[,j]) > tol)
    numNzY <- length(nzIdY)
    numNzX <- length(nzIdX)
    #no predictors selected
    if(numNzY == 0 & numNzX == 0) {
      sigma.fit[j] <- 1/var(Y.j)
    } else if (numNzY > 0 & numNzX == 0) {
      if(numNzY > N) {
        sparse <- FALSE
        break;
      }
      Ynz <- Y[,nzIdY]
      coef <- fastLmPure(as.matrix(Ynz),Y.j)$coefficients
      resid <- Y.j - Ynz%*%coef
      sigma.fit[j] <- 1/var(resid)
      Beta[nzIdY,j] <- coef
    } else if (numNzY == 0 & numNzX >  0) { 
      if(numNzX > N) {
        sparse <- FALSE
        break;
      }
      Xnz <- X[,nzIdX]
      coef <- fastLmPure(as.matrix(Xnz),Y.j)$coefficients
      resid <- Y.j - Xnz%*%coef
      sigma.fit[j] <- 1/var(resid)
      Gamma[nzIdX,j] <- coef
    } else if (numNzY > 0 & numNzX >  0) {
      if(numNzX + numNzY > N) {
        sparse <- FALSE
        break;
      }
      nzJointYX <- cbind(Y[,nzIdY], X[,nzIdX])
      coef <- fastLmPure(as.matrix(nzJointYX),Y.j)$coefficients
      resid <- Y.j - nzJointYX%*%coef
      sigma.fit[j] <- 1/var(resid)
      Beta[nzIdY,j] <- coef[1:numNzY]
      Gamma[nzIdX,j] <- coef[(numNzY + 1):length(coef)]
    }
  }
  
  if(!sparse) {
    ParCor <- matrix(NA, Q, Q)
    Gamma <- matrix(NA, P, Q)
    sigma.fit <- rep(NA, times = Q)
  } else {
    ParCor <- sign(Beta) * sqrt( abs( Beta * t(Beta) )  )
    diag(ParCor) <- 1
  }
  return(list(ParCor=ParCor, Gamma=Gamma, sig.fit=sigma.fit))
}


###
RSS.OLS<-function(Y,try.adj){
  p<-ncol(Y)
  beta<-apply(matrix(1:p),1, OLS.Sol,Y=Y,try.adj=try.adj)
  rss<-1/diag(beta)*(n-1)
  return(rss)
}

