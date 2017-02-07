# #' Obtain Robust Networks Through Spacemap
# #'
# #' We highlight key functions for learning a robust network from \pkg{spacemap}.
# #'
# #' There are several steps for learning robust genetic regulatory networks
# #' from \pkg{spacemap}. First, tune the model with the CV.Vote procedure found in
# #' \code{\link{crossValidation}} function. Second, consider further tuning the
# #' spacemap model under CV-selected tuning parameters with Boot.Vote procedures
# #' found in \code{\link{reprodEdge}}. Third, derive biological insight from the
# #' inferred networks from steps 1,2 with the network analysis toolkit.
# #' @rdname spacemap-package
# #' @references
# #' Christopher J. Conley, Pei Wang,  and Jie Peng. Characterizing major functional consequences of DNA copy number alterations
# #' in breast tumors through a new conditional graphical model. Submitted to ISMB-ECCB 2017.
# "_PACKAGE"
# #> [1] "_PACKAGE"


#'@title Fit spacemap model
#'@description Estimate partial correlations between response variables (Y)
#' and perform model selection of predictors (X) based on the spacemap model 
#' fitted to multivariate data (Y,X).
#'@param Y.m A numeric matrix (N x Q) representing the response of Q variables on 
#'N independent observations.
#'@param X.m A numeric matrix (N X P) representing P predictor variables on N 
#'independent observations.  
#'@param slasso A numeric value, the Space-lasso penalty parameter (l_1 norm). 
#'@param sridge A numeric value (optional). If not specified, lasso regression 
#'is used in the Joint Sparse Regression Model (JSRM). Otherwise, elastic net 
#'regression is used in JSRM and sridge serves as the squared l_2 norm penalty parameter.
#'@param rlasso The MAP-lasso penalty. 
#'@param rgroup The MAP-group penalty.
#'@param sig A numeric vector with length being the same as the number of 
#'columns of Y.m. It is the vector of \eqn{\sigma^{ii}} (the diagonal of 
#'the inverse covariance matrix). If not specified, σ^{ii} will be 
#'estimated during the model fitting with initial values \code{rep(1,p)} and  
#'the iterations of the model fitting (\code{iter}) \eqn{\ge 2}. 
#'Note, the scale of \code{sig} does not matter.
#'@param weight A numeric value or vector (Q x 1) specifying the weights or the type of
#' weights used for each regression in JSRM. The default value is NULL, 
#' which means all regressions will be weighted equally in the joint model. 
#' If \code{weight=1}, residue variances will be used for weights. 
#' If \code{weight=2}, the estimated degree of each variable will be used for weights. 
#' Otherwise, it should be a positive numeric vector, whose length is equal to 
#' the number of columns of Y.m.
#' @param remWeight optional weighting scheme on [fill in details].
#' @param iter An integer specifying the iterations in JSRM for estimating 
#' \eqn{\sigma^{ii}} and partial correlations. When \code{sig=NULL} and/or 
#' (\code{weight=NULL} or \code{weight=2}), then \code{iter} \eqn{\ge 2}.
#' 
#' @return A list containing 
#'\enumerate{
#'   \item \code{ParCor} The estimated partial correlation matrix. 
#'   \item \code{sig.fit} The estimated diagonal \eqn{\sigma^{ii}}.
#'   \item \code{Gamma} The estimated coefficients of Y on X.
#'   \item \code{rss} The residual sums of squares from the model fit. 
#'   \item \code{deltaMax} The maximum change in parameter values between the penultimate and ultimate iteration  
#' }
#' 
#' @examples
#' data(sim1)
#' net <- spacemap(Y.m = sim1$Y, X.m = sim1$X, slasso = 70, rlasso = 28.8, rgroup = 12.38)
spacemap <-function(Y.m, X.m, slasso, sridge=0, rlasso, rgroup, sig=NULL, 
                    weight=NULL, remWeight = NULL, iter=3, tol = 1e-6, cd_iter = 1e7L, 
                    verbose = FALSE) {
  
  ####################### return value
  ## A list: the estimated \{\rho^{ij}\}: p by p matrix, and $\{\sigma^{ii}\}$: p by 1 vector
  #######################
  
  n=nrow(Y.m)
  p=ncol(Y.m)
  ITER=iter
  
  #check for missing data
  if (any(is.na(Y.m))) { 
    stop("Missing values found in input matrix Y.m; imputation is required prior to Spacemap fitting.")
  } 
  if (any(is.na(X.m))) { 
    stop("Missing values found in input matrix X.m; imputation is required prior to Spacemap fitting.")
  }
  
  #Standardize the response vector.
  Y.s <- scale(Y.m)
  X.s <- scale(X.m)
  #Keep the same scale as user input:
  # W = diag(apply(Y, 2, sd))
  # R = inv(W) \Sigma inv(W)
  #...Therefore
  # \Sigma = inv(W) inv(R) inv(W)
  Y.stddev <- attr(Y.s, "scaled:scale")
  
  #RemMap Weight 
  if(is.null(remWeight)) {
    remWeight = matrix(1, nrow = ncol(X.s), ncol = ncol(Y.s))
  } else {
    stopifnot(is.matrix(remWeight))
  }
  
  ################### preparation
  if(!is.null(sig))
  { #### if the user specify "sig", do not need to update sig.
    SIG.update=F
    SIG=sig
  } else { #### otherwise, need to update sig in each iteration
    SIG.update=T
    SIG=rep(1, p) 
  } 
  
  if(length(weight)==0 | (length(weight)>1 & length(weight)<p)) ### weight==NULL 
  {
    WEIGHT.tag=0 ### equal weight 
    WEIGHT=rep(1, p)
    WEIGHT.update=F
  } 
  if(length(weight)==1)
  {
    if(weight==1)
    {
      WEIGHT.tag=1 ### sig based weight
      WEIGHT=SIG
      WEIGHT.update=T
      ITER=max(2, iter)
    } else {
      WEIGHT.tag=2 ### degree based weight
      WEIGHT=rep(1,p)
      WEIGHT.update=T
      ITER=max(2, iter)
    } 
  }
  if(length(weight)==p)  
  {
    WEIGHT.tag=3 ## prespecified weight
    WEIGHT=weight
    WEIGHT.update=F
  }
  ################## begin to iterate
  if (verbose)
    cat("iter: ")
  for(i in 1:iter) {
    if (verbose)
      cat(i, ", "); flush.console();
    Y.u <- Y.s*matrix(sqrt(WEIGHT),n,p,byrow=TRUE)
    X.u <- X.s  #*matrix(sqrt(WEIGHT),n,ncol(X.s),byrow=TRUE)
    sig.u<-SIG/WEIGHT
    
    mod.fit <- doSpaceMap(Y.u, X.u, remWeight, sig.u^0.5, 
		    slasso, sridge, rlasso, rgroup, tol, cd_iter)
    ParCor.fit <- matrix(mod.fit$rho, p,p, byrow=T)
    
    diag(ParCor.fit)<-1
    
    coef<-ParCor.fit[upper.tri(ParCor.fit)]
    beta.cur<-Beta.coef(coef,SIG) 
    
    
    if(!WEIGHT.update & !SIG.update)
    {
      break
    } else { ## update sig or weight
      if(SIG.update)
      {
        SIG<-gInvSig.diag.new(Y.s, X.s, beta.cur, mod.fit$Gamma)   
      }
      if(WEIGHT.update)
      {
        if(WEIGHT.tag==1)
        {        #### sig based
          WEIGHT=SIG
        } else { #### degree based
          temp.w<-apply(abs(ParCor.fit) > tol,1,sum)
          temp.w<-temp.w+max(temp.w)     
          WEIGHT<-temp.w/sum(temp.w)*p    
        }
      } 
    } ### end updating WEIGHT and SIG
    #Do not continue estimating if it didn't converge in the first iteration. 
    if (!mod.fit$convergence) break;
  } ### end iteration
  
  #Keep the same scale as user input:
  SIG <- SIG / (Y.stddev * Y.stddev) 
  
  result<-list(ParCor=ParCor.fit,sig.fit=SIG, Gamma = mod.fit$Gamma, 
               rss = mod.fit$rss, convergence = mod.fit$convergence,
	       deltaMax = mod.fit$deltaMax)
  return(result)  
}


########################################
##### Estimate diagonal sigma
##use the fact that 1/sig^{ii} is the residual variance in one versus all others setting


gInvSig.diag.new<-function(Y, X, Beta, Gamma){
  ################### parameters
  ### Y:    n by p data matrix; should be standardized to sd 1;
  ### Beta: beta matrix for the regression model
  p=ncol(Y)
  Beta.temp=Beta
  diag(Beta.temp)=0
  esti.Y=Y%*%Beta.temp + X %*% Gamma
  residue=Y-esti.Y
  result=apply(residue^2, 2, mean)
  return(1/result)
}


