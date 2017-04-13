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


#pass defaults if not specified by user
varySpaceMapParam <- function(...) {
  opt <- list(...)
  
  #l_2 penalty paramter
  if (is.null(opt$sridge)) { opt$sridge <- 0.0 }
  
  # Numeric value or vector (Q x 1) specifying the weights or the type of
  #weights used in the loss function of spac. The default value is NULL, 
  #which means all regressions will be weighted equally in the space step of the model. 
  #If \code{weight=1}, residue variances will be used for weights. 
  #If \code{weight=2}, the estimated degree of each variable will be used for weights. 
  #Otherwise, it should be a positive numeric vector, whose length is equal to 
  #the number of columns of Y.
  if (is.null(opt$weight)) { opt$weight <- NULL }
  
  #Needs further exploration 
  if (is.null(opt$remWeight)) { opt$remWeight <- NULL }
  opt
}


#' @title Fit the spacemap model.
#' @description The conditional graphical model called spacemap learns a
#' sparse network encoding the interactions between two high dimensional data inputs denoted as response and predictor variables. 
#' @details The two data types input are response random vector 
#' \eqn{\textbf{y} = (y_1, \dots, y_Q)} and predictor random vector \eqn{\textbf{x} = (x_1, \dots, x_P)}. 
#' Provided \eqn{N} iid samples from \eqn{(\textbf{x},\textbf{y})}, spacemap learns a sparse network, where the degree
#' of sparsity is determined by tuning penalties \eqn{\lambda_1, \lambda_2, \lambda_3}.
#' When we fit the spaceMap model, we seek to learn the network comprised of node set \eqn{\textbf{V}= (\textbf{x},\textbf{y})}
#' and edge set \eqn{\textbf{E}= \{ x_p - y_q \} \cup \{ y_q - y_l : q \neq l \}}. 
#' @param Y Numeric matrix \eqn{(N \times Q)} containing N iid samples of the response vector \eqn{\textbf{y}}. 
#' @param X Numeric matrix \eqn{(N \times P)} containing N iid samples of the predictor vector \eqn{\textbf{x}}. 
#' @param lam1 Non-negative numeric value, the space-lasso penalty corresponding to \eqn{\lambda_1}, 
#' which subjects the partial correlation vector, \eqn{\bf \rho_{yy}}, to  the \eqn{l_1} norm. 
#' It induces overall sparsity of \eqn{\{ y_q - y_l : q \neq l \}} edges. .
#' @param lam2 Non-negative numeric value, the MAP-lasso penalty corresponding to \eqn{\lambda_2}, 
#' which subjects the regression coefficients in \eqn{ \bf \Gamma}, to the  \eqn{l_1} norm.  
#' It induces overall sparsity of \eqn{\{ x_p - y_q \}} edges. 
#' @param lam3 Non-negative numeric value, the MAP-group penalty corresponding to \eqn{\lambda_3},
#' which subjects the regression coefficients in the \eqn{p}th row of \eqn{ \bf \Gamma_{P \times Q}}, to the \eqn{l_2} norm. 
#' It encourages selection of predictor variables with many edges to response variables, in other words, hub nodes. 
#' @param sig Positive numeric vector (\eqn{p \times 1}) representing the estimate of \eqn{\sigma^{ii}}, the diagonal of 
#' the inverse covariance matrix. It defaults to NULL and
#' and will be estimated \code{iter} times during the model fitting with initial 
#' values set to the ones vector.    
#' @param rho Numeric matrix (\eqn{p \times p}) representing the estimate of \eqn{\hat\rho}, the partial correlation matrix. 
#' It defaults to NULL and and will be estimated \code{iter} times during the model fitting. 
#' @param iter Positive integer specifying the number of iterations for estimating 
#' \eqn{\sigma^{ii}} and \eqn{\rho}. Defaults to 3. 
#' @param tol Positive numeric value specifying the convergence tolerance of the coordinate descent algorithm; 
#' in other words, it is criterion that stops parameter estimation when no parameter changes value exceeding \code{tol} 
#' between iterations. \code{tol} defaults to 1e-6, but may be lowered (e.g. 1e-4) to speed
#' up network learning. 
#' @param iscale Logical indicating to standardize the whole input data. Defaults to TRUE. 
#' See \code{\link{base::scale(x, center = TRUE, scale = TRUE)}} for details of standardization. 
#' @param cdmax Positive integer specifiying the maximum number of parameter updates allowed before reporting the algorithm
#' as having failed to converge. Default may need to be increased for inferring very large-scale networks (i.e. \eqn{p,q > 1000}).
#' @return A list containing 
#'\enumerate{
#'   \item \code{ParCor} The estimated partial correlation matrix (\eqn{P \times P}), 
#'   where off-diagonals  \eqn{ |\hat \rho^{p,q}_{yy}| > 1e-6}  encode the edges \eqn{\{ y_q - y_l : q \neq l \} }
#'   and the diagonals are 1's. 
#'   \item \code{sig.fit} The estimated diagonal \eqn{\hat \sigma^{ii}}.
#'   \item \code{Gamma} The estimated regression coefficients matrix  (\eqn{P \times Q}), 
#'   where elements \eqn{ |\hat \gamma_{pq}| > 1e-6} encode the edges \eqn{\{ x_p - y_q \}}. 
#'   \item \code{rss} The residual sums of squares from the model fit. 
#'   \item convergence logical: true for successful convergence, otherwise failed to converge. Failure can be 
#'   mitigated by increasing \code{tol} and/or \code{cdmax}. 
#'   \item \code{deltaMax} The maximum change in parameter values between the penultimate and ultimate iteration. 
#'   If \code{spacemap} does not converge, \code{deltaMax} provides some measure of how far away it was from converging
#'   when compared to \code{tol}. 
#' }
#' @seealso  \code{\link{cvVote}}, \code{\link{bootEnsemble}}, \code{\link{bootVote}}
#' @examples
#' data(sim1)
#' net <- spacemap(Y = sim1$Y, X = sim1$X, lam1 = 70, lam2 = 28.8, lam3 = 12.38)
#' #adjacency matrix of y-y and x-y edges. 
#' adjnet <- adjacency(net)
#' @export
spacemap <- function(Y, X, lam1, lam2, lam3, sig=NULL, rho = NULL,
                     iter=3, tol = 1e-6, iscale = TRUE, cdmax = 1e7L, ...) {
  
  #experimental/optional parameters
  opt <- varySpaceMapParam(...)
  
  n=nrow(Y)
  p=ncol(Y)
  ITER=iter
  
  #check for missing data
  if (any(is.na(Y))) { 
    stop("Missing values found in input matrix Y; imputation is required prior to Spacemap fitting.")
  } 
  if (any(is.na(X))) { 
    stop("Missing values found in input matrix X; imputation is required prior to Spacemap fitting.")
  }
  
  if(iscale) {
    #Standardize the response vector.
    Ys <- scale(Y)
    Xs <- scale(X)
    #Keep the same scale as user input:
    # W = diag(apply(Y, 2, sd))
    # R = inv(W) \Sigma inv(W)
    #...Therefore
    # \Sigma = inv(W) inv(R) inv(W)
    Y.stddev <- attr(Ys, "scaled:scale")
  } else { 
    Ys <- Y
    Xs <- X
    Y.stddev <- rep(1,p)
  }
  
  #RemMap Weight 
  if(is.null(opt$remWeight)) {
    opt$remWeight = matrix(1, nrow = ncol(Xs), ncol = ncol(Ys))
  } else {
    stopifnot(is.matrix(opt$remWeight))
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
  
  if(length(opt$weight)==0 | (length(opt$weight)>1 & length(opt$weight)<p)) ### opt$weight==NULL 
  {
    WEIGHT.tag=0 ### equal weight 
    WEIGHT=rep(1, p)
    WEIGHT.update=F
  } 
  if(length(opt$weight)==1)
  {
    if(opt$weight==1)
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
  if(length(opt$weight)==p)  
  {
    WEIGHT.tag=3 ## prespecified weight
    WEIGHT=opt$weight
    WEIGHT.update=F
  }
  ################## begin to iterate
  for(i in 1:iter) {
    Yu <- Ys*matrix(sqrt(WEIGHT),n,p,byrow=TRUE)
    Xu <- Xs  #*matrix(sqrt(WEIGHT),n,ncol(Xs),byrow=TRUE)
    sig.u<-SIG/WEIGHT
    
    if (!is.null(rho) & is.matrix(rho)) {
      stopifnot(nrow(rho) == p, ncol(rho) == p)
      mod.fit <- doSpaceMap(Yu, Xu, opt$remWeight, sig.u^0.5, 
                            lam1, opt$sridge, lam2, lam3, tol, cdmax,
                            as.vector(rho), TRUE)
    } else {
      mod.fit <- doSpaceMap(Yu, Xu, opt$remWeight, sig.u^0.5, 
                            lam1, opt$sridge, lam2, lam3, tol, cdmax,
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
        SIG<-gInvSig.diag.new(Ys, Xs, beta.cur, mod.fit$Gamma)   
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

#' Adjacency matrix from spacemap and spacemap models
#' 
#' Output from \code{\link{spacemap}} and \code{\link{space.joint}} is 
#' transformed into an adjacency matrix encoding a network.
#' 
#' @param net Object output from \code{\link{spacemap}} and \code{\link{space.joint}}
#' @param aszero Positive numeric value (defaults to 1e-6) indicating at what point to consider extremely 
#' small parameter estimates of \eqn{\Gamma} and \eqn{\rho} as zero. 
#' @return A list containing 
#' \itemize{ 
#'  \item  \code{yy} A \eqn{p \times p} adjancency matrix 
#'  indicates an edge\eqn{y_q - y_l} when \code{yy[q,l] == 1} and 0 otherwise. 
#'  \item  \code{xy} A \eqn{p \times q} adjancency matrix 
#'  indicates an edge\eqn{x_p - y_q} when \code{xy[p,q] == 1} and 0 otherwise. 
#' }
#' @export
adjacency <- function(net,  aszero  = 1e-6) { 
  
  #set the diagonal to zero to avoid self-loops
  diag(net$ParCor) <- 0
  if (!is.null(net$Gamma)) { 
    out <- list(yy = (abs(net$ParCor) > aszero) + 0L, 
                xy = (abs(net$Gamma) > aszero) + 0L)
  } else  { 
    out <- list(yy = (abs(net$ParCor) > aszero) + 0L)
  }
  out
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


