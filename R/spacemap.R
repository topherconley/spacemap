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


#'@title Fit the spacemap model.
#'@description The conditional graphical model called spacemap learns a
#'sparse network encoding the interactions between two high dimensional data inputs denoted as response and predictor variables. 
#'@details The two data types input are response random vector 
#'\eqn{\textbf{y} = (y_1, \dots, y_Q)} and predictor random vector \eqn{\textbf{x} = (x_1, \dots, x_P)}. 
#'Provided \eqn{N} iid samples from \eqn{(\textbf{x},\textbf{y})}, spacemap learns a sparse network, where the degree
#' of sparsity is determined by tuning penalties \eqn{\lambda_1, \lambda_2, \lambda_3}.
#'When we fit the spaceMap model, we seek to learn the network comprised of node set \eqn{\textbf{V}= (\textbf{x},\textbf{y})}
#'and edge set \eqn{\textbf{E}= \{ x_p - y_q \} \cup \{ y_q - y_l : q \neq l \}}. 
#'@param Y.m Numeric matrix \eqn{(N \times Q)} containing N iid samples of the response vector \eqn{\textbf{y}}. 
#'@param X.m Numeric matrix \eqn{(N \times P)} containing N iid samples of the predictor vector \eqn{\textbf{x}}. 
#'@param slasso Non-negative numeric value, the space-lasso penalty corresponding to \eqn{\lambda_1}, 
#'which subjects the partial correlation vector, \eqn{\bf \rho_{yy}}, to  the \eqn{l_1} norm. 
#'It induces overall sparsity of \eqn{\{ y_q - y_l : q \neq l \}} edges. 
#'@param sridge Non-negative numeric value defaults to 0, which is the recommended value. Otherwise, 
#' a positive value applies an elastic net penalty to \eqn{\bf \rho_{yy}}, 
#'where \code{sridge} serves as the squared \eqn{l_2} norm penalty parameter.
#'@param rlasso Non-negative numeric value, the MAP-lasso penalty corresponding to \eqn{\lambda_2}, 
#'which subjects the regression coefficients in \eqn{ \bf \Gamma}, to the  \eqn{l_1} norm.  
#'It induces overall sparsity of \eqn{\{ x_p - y_q \}} edges. 
#'@param rgroup Non-negative numeric value, the MAP-group penalty corresponding to \eqn{\lambda_3},
#'which subjects the regression coefficients in the \eqn{p}th row of \eqn{ \bf \Gamma_{P \times Q}}, to the \eqn{l_2} norm. 
#'It encourages selection of predictor variables with many edges to response variables, in other words, hub nodes. 
#'@param sig Positive numeric vector (\eqn{p \times 1}) representing the estimate of \eqn{\sigma^{ii}}, the diagonal of 
#'the inverse covariance matrix. It defaults to NULL and
#'and will be estimated \code{iter} times during the model fitting with initial 
#'values set to the ones vector.    
#'@param weight Numeric value or vector (Q x 1) specifying the weights or the type of
#' weights used for each regression in JSRM. The default value is NULL, 
#' which means all regressions will be weighted equally in the joint model. 
#' If \code{weight=1}, residue variances will be used for weights. 
#' If \code{weight=2}, the estimated degree of each variable will be used for weights. 
#' Otherwise, it should be a positive numeric vector, whose length is equal to 
#' the number of columns of Y.m.
#' @param remWeight Optional weighting scheme [fill in the details].
#' @param iter Positive integer specifying the number of iterations for estimating 
#' \eqn{\sigma^{ii}} and partial correlations. Defaults to 3. When \code{sig=NULL} and/or 
#' (\code{weight=NULL} or \code{weight=2}), then \code{iter} \eqn{\ge 2}.
#' @param tol Positive numeric value specifying the convergence tolerance of the coordinate descent algorithm; 
#' in other words, it is criterion that stops parameter estimation when no parameter changes value exceeding \code{tol} 
#' between iterations. \code{tol} defaults to 1e-6, but may be lowered (e.g. 1e-4) to speed up network learning. 
#' @param cd_iter Positive integer specifiying the maximum number of parameter updates allowed before reporting the algorithm
#' as having failed to converge. Default may need to be increased for inferring very large-scale networks (i.e. \eqn{p,q > 1000}).
#' @param verbose Logical indicating to print out the current status of \code{iter}. Default is FALSE. 
#' @param iscale Logical indicating to standardize the input data internally. Defaults to TRUE if just running \code{spacemap}. If running 
#' \code{\link{crossValidation}} or \code{\link{bootEnsemble}}, user should set to FALSE.  
#' @return list containing 
#'\enumerate{
#'   \item \code{ParCor} The estimated partial correlation matrix (\eqn{P \times P}), 
#'   where off-diagonals  \eqn{ |\hat \rho^{p,q}_{yy}| > } \code{tol} encode the edges \eqn{\{ y_q - y_l : q \neq l \} }
#'   and the diagonals are 1's. 
#'   \item \code{sig.fit} The estimated diagonal \eqn{\hat \sigma^{ii}}.
#'   \item \code{Gamma} The estimated regression coefficients matrix  (\eqn{P \times Q}), 
#'   where elements \eqn{ |\hat \gamma_{pq}| > } \code{tol} encode the edges \eqn{\{ x_p - y_q \}}. 
#'   \item \code{rss} The residual sums of squares from the model fit. 
#'   \item convergence logical: true for successful convergence, otherwise failed to converge. Failure can be 
#'   mitigated by increasing \code{tol} and/or \code{cd_iter}. 
#'   \item \code{deltaMax} The maximum change in parameter values between the penultimate and ultimate iteration. 
#' }
#' @seealso  \code{\link{crossValidation}}, \code{\link{bootVote}}
#' @examples
#' data(sim1)
#' net <- spacemap(Y.m = sim1$Y, X.m = sim1$X, slasso = 70, rlasso = 28.8, rgroup = 12.38)
#' @export
spacemap <-function(Y.m, X.m, slasso, sridge=0, rlasso, rgroup, sig=NULL, rho = NULL,
                    weight=NULL, remWeight = NULL, iter=3, tol = 1e-6, cd_iter = 1e7L, 
                    verbose = FALSE, iscale = TRUE) {
  
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
  
  if(iscale) {
    #Standardize the response vector.
    Y.s <- scale(Y.m)
    X.s <- scale(X.m)
    #Keep the same scale as user input:
    # W = diag(apply(Y, 2, sd))
    # R = inv(W) \Sigma inv(W)
    #...Therefore
    # \Sigma = inv(W) inv(R) inv(W)
    Y.stddev <- attr(Y.s, "scaled:scale")
  } else { 
    Y.s <- Y.m
    X.s <- X.m
    Y.stddev <- rep(1,p)
  }
  
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
    
    if (!is.null(rho) & is.matrix(rho)) {
      stopifnot(nrow(rho) == p, ncol(rho) == p)
      mod.fit <- doSpaceMap(Y.u, X.u, remWeight, sig.u^0.5, 
                            slasso, sridge, rlasso, rgroup, tol, cd_iter,
                            as.vector(rho), TRUE)
    } else {
      mod.fit <- doSpaceMap(Y.u, X.u, remWeight, sig.u^0.5, 
                            slasso, sridge, rlasso, rgroup, tol, cd_iter,
                            as.vector(matrix(0, p, p)), FALSE)
    }
    
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


