#' Ensemble networks through boostrapping. 
#' 
#' Fits spacemap model to B bootstrapped data replicates to create an ensemble network. 
#' 
#' @inheritParams spacemap
#' @param tune List with names being the tuning penalties
#' \code{lam1, lam2, lam3} and each element being a single numeric value that is input to \code{spacemap}. One might  
#' consider the output from \code{\link{cvVote}} (i.e. list \code{minTune}) as input for this parameter. 
#' @param method Character value specifying to use \code{\link{spacemap}} or 
#' \code{\link{space}} in fitting ensemble network. 
#' @param boot Logical. Default is \code{boot=TRUE} and implies sampling with replacement (i.e. boostrap). 
#' \code{boot=FALSE} will resample without replacement. 
#' @param B Positive integer denoting the number of model fits making up the ensemble. Default is 1000, but can 
#' be lowered to save time. 
#' @param aszero Positive numeric value (defaults to 1e-6) indicating at what point to consider extremely 
#' small parameter estimates of \eqn{\Gamma} and \eqn{\rho} as zero. 
#' @param seed Positive integer allowing user to set the random seed for reproducing the network ensemble. 
#' @param p0 Positive numeric not exceeding the default of 1, which represents the proportion of the original 
#' samples that will be down-sampled for each model fitting iteration of the ensemble. 
#' @param ... Additional arguments for \code{\link{spacemap}} or  \code{\link{space}}.
#' @return A list up to length \code{B} of convergent and sparse model fits from 
#' either \code{\link{spacemap}} or \code{\link{space}}. The list object should typically
#' not be modified by the user, but passed on to the  \code{\link{bootVote}} function.
#' Model fits that are not convergent or produce networks that are not sparse are not reported. 
#' @importFrom rngtools RNGseq setRNG
#' @importFrom Matrix Matrix
#' @import RcppArmadillo
#' @export 
#' @examples 
#' 
#' #Load simulation
#' library(spacemap)
#' data(sim1)
#' 
#' #Boostrap Ensemble (B = 10) for illustration only
#' tune <- data.frame(lam1 = 70, lam2 = 28, lam3 = 17.5)
#' ens <- bootEnsemble(Y = sim1$Y, X = sim1$X, tune = tune, 
#'                     method = "spacemap", B = 10)
#' @seealso \code{\link{bootVote}}, \code{\link{cvVote}}
#' 
#' 
bootEnsemble <- function(Y, X = NULL, tune, method = c("spacemap", "space"), 
                         boot = TRUE, B =  1000L, iscale = TRUE, aszero = 1e-6,  
                         seed = sample.int(n = 1e6, size = 1), 
                         p0 = 1.0, ...) {
  
  method <- match.arg(method)
  requireNamespace("rngtools")
  requireNamespace("spacemap")
  requireNamespace("foreach")
  requireNamespace("Matrix")
  rng <- RNGseq(B, seed)
  N <- nrow(Y)
  P <- ncol(X)
  Q <- ncol(Y)
  rN <- floor(p0*N)
  
  #variable number of parameters
  opt <- variTrainParam(...)
  
  if(!is.matrix(Y)) stop("Y is not a matrix.")
  givenX <- !is.null(X) 
  Xindex <- NULL 
  Yindex <- NULL
  if(method == "space" & givenX)  {
    if(!is.matrix(X)) stop("X is not a matrix.")
    if(iscale) { 
      X <- scale(X)
      Y <- scale(Y)
    }
    XY = cbind(X, Y)
    Xindex <- seq_len(ncol(X))
    Yindex <- (ncol(X) + 1):(ncol(X) + ncol(Y))
  } else if(method == "space" & !givenX){ 
    if(iscale) { 
      Y <- scale(Y)
    }
    XY = Y
  } else if (method == "spacemap") { 
    if(!is.matrix(X)) stop("X is not a matrix.")
    if(iscale) { 
      X <- scale(X)
      Y <- scale(Y)
    }
  }
  
  if (method == "spacemap") { 
    ens <- foreach(j = seq_len(B),  r = rng[seq_len(B)], .packages = c("Matrix", "rngtools")) %dopar% {
      rngtools::setRNG(r)
      sidx <- sample.int(n = N, size = rN, replace = boot)
      fit <- spacemap::spacemap(Y = Y[sidx,], X = X[sidx,], 
                                lam1 = tune$lam1[1], sridge = opt$sridge, 
                                lam2 = tune$lam2[1], lam3 = tune$lam3[1],
                                sig = opt$sig, rho = opt$rho, 
                                weight = opt$weight, remWeight = opt$remWeight, 
                                iter = opt$iter, tol = opt$tol, cdmax = opt$cdmax,
                                iscale = FALSE)
      #zero out those below tolerance
      fit$ParCor[abs(fit$ParCor) <= aszero] <- 0.0
      fit$Gamma[abs(fit$Gamma) <= aszero] <- 0.0
      if (!fit$convergence) { 
        fit$sparse <- FALSE
        return(fit[c("convergence", "sparse")])
      }
      fit$sparse <- olsRefitSpacemap(Y = Y[sidx,], X = X[sidx,], 
                                     ParCor = fit$ParCor, Gamma = fit$Gamma, 
                                     sigma = fit$sig.fit, RSS = fit$rss, tol = opt$tol, 
                                     ridge = opt$refitRidge)
      if (!fit$sparse) { 
        return(fit[c("convergence", "sparse")])
      } 
      fit$xy <- sparseAdjMat(fit$Gamma)
      fit$yy <- sparseAdjMat(fit$ParCor)
      fit[c("xy", "yy", "convergence", "sparse")]
    }
  } else if (method == "space") { 
    ens <- foreach(j = seq_len(B),  r = rng[seq_len(B)], .packages = c("Matrix", "rngtools")) %dopar% {
      rngtools::setRNG(r)
      sidx <- sample.int(n = N, size = rN, replace = boot)
      fit <- spacemap::space(Y = XY[sidx,], lam1 = tune$lam1[1], 
                                   sridge = opt$sridge, iter = opt$iter, 
                                   cdmax = opt$cdmax, tol = opt$tol, 
                                   iscale = FALSE)
      #zero out those below tolerance. 
      fit$ParCor[abs(fit$ParCor) <= aszero] <- 0.0
      if (!fit$convergence) { 
        fit$sparse <- FALSE
        return(fit[c("convergence", "sparse")])
      }
      
      fit$sparse <- olsRefitSpace(Y = XY[sidx,],
                                  ParCor = fit$ParCor, 
                                  sigma = fit$sig.fit, RSS = fit$rss, tol = opt$tol,
                                  ridge = opt$refitRidge)
      if (!fit$sparse) { 
        return(fit[c("convergence", "sparse")])
      } 
      
      if (givenX) { 
        fit$xy <- sparseAdjMat(fit$ParCor[Xindex,Yindex])
        fit$yy <- sparseAdjMat(fit$ParCor[Yindex,Yindex])
        return(fit[c("xy", "yy", "convergence", "sparse")])
      } else {
        fit$yy <- sparseAdjMat(fit$ParCor)
        return(fit[c("yy", "convergence", "sparse")])
      }
    }
  }
  
  goodIndex <- which(sapply(ens, function(x) x$convergence & x$sparse))
  ngi <- length(goodIndex)
  if (ngi < B) { 
    message(paste("Out of ", B, "replicate fits, ", ngi, "converged and were sufficiently sparse."))
  }
  good <- ens[goodIndex]
  attr(good, which = "method") <- method
  good
}  

#Make a sparse adjacency matrix
sparseAdjMat <- function(x) { 
  Matrix((abs(x) > 0) + 0)
}

#' Model aggregation through Boot.Vote
#' 
#' Aggregate boostrap replicates of spacemap into a final Boot.Vote model. 
#' @param bfits List of fitted spacemap models returned from \code{\link{bootEnsemble}}.
#' @param thresh Positive numeric threshold for the minimum proportion of bootstrap 
#' replicate model fits with a particular edge such that the edge is included in the Boot.Vote model. 
#' @param givenX Logical. Defaults to FALSE. Should be set to TRUE when 
#' \code{attr(bfits, "method) == "space} and \code{space} was used
#' to infer (x--x, x--y, y--y) edges  but only reported (x--y, y--y) edges. 
#' @return Returns a list of lists. 
#' First list is \code{bv}, which encodes the edges in two logical adjacency matrices.
#' 
#' \enumerate{
#'   \item \code{yy} Adjacency matrix where 1 for the (q,l) off-diagonals element indicate an edge 
#'   between the qth and lth response variables, and 0 otherwise. 
#'   \item \code{xy}  Adjacency matrix where 1 for the (p,q) element indicate an edge 
#'   between the pth predictor and qth response variable, and 0 otherwise. 
#' }
#' 
#' Second list is \code{bdeg}, which contains the degree distribution for each bootstrap replicate fit. 
#'
#' \enumerate{
#'   \item \code{yy} Integer matrix (\eqn{B \times Q} where the (q,b) off-diagonals element indicates
#'    the out-degree of the qth response variable for the bth converged model based on the bth bootstrap replicate. 
#'   \item \code{xy} Integer matrix (\eqn{B \times P} where the (p,b) element indicates
#'    the out-degree of the pth predictor variable for the bth converged model based on the bth bootstrap replicate. 
#' }
#' 
#' 
#'  Third list is \code{bc}, which stores several additional statistics on the ensemble network fits. 
#' \enumerate{ 
#'    \item \code{yy} Integer matrix containing the y--y edge selection frequency out of B replicates.
#'    \item \code{xy}  Integer matrix containing the x--y edge selection frequency out of B replicates.
#'    \item \code{dfyy} Integer vector containing the total number of y--y edges for each fit.
#'    \item \code{dfxy} Integer vector containing the total number of x--y edges for each fit.
#' }
#' 
#' Note: If \code{method == "space" & givenX == FALSE}, 
#' then no \code{xy, dfxy} elements will be reported in the above lists. 
#' 
#' @seealso  \code{\link{bootEnsemble}}
#' 
#' @examples 
#' 
#' 
#' #Load simulation
#' library(spacemap)
#' data(sim1)
#' 
#' #Boostrap Ensemble (B = 10) for illustration only
#' tune <- data.frame(lam1 = 70, lam2 = 28, lam3 = 17.5)
#' 
#' #suppress warnings because parallel backend not set up. 
#' ens <- suppressWarnings(bootEnsemble(Y = sim1$Y, X = sim1$X, tune = tune, 
#'                     method = "spacemap", B = 10))
#'                     
#' bv <- suppressWarnings(bootVote(ens))         
#' @importFrom foreach %dopar%
#' @export
bootVote <- function(bfits, thresh = 0.5, givenX = FALSE) {
  
  if(!is.list(bfits)) { 
    stop("bfits is not a list")
  }

  method <- attr(bfits, "method")
  
  #requireNamespace("spacemap")
  requireNamespace("foreach")
  
  if (method == "spacemap" | (method == "space" & givenX)) { 
    
    addmods <- function(m1, m2) { 
      list(xy = m1$xy + m2$xy, yy = m1$yy + m2$yy, 
           dfxy = c(m1$dfxy, m2$dfxy), dfyy = c(m1$dfyy, m2$dfyy),
           ngood = c(m1$ngood, m2$ngood))
    }
    
    bout <- foreach(bfit = bfits, .combine = addmods, .packages = c("spacemap")) %dopar% { 
      bxy <- as.matrix(bfit$xy)
      byy <- as.matrix(bfit$yy)
      list(xy = bxy, yy = byy, 
           dfxy = nonZeroWhole(bxy, 0), 
           dfyy = nonZeroUpper(byy, 0), ngood = bfit$convergence)
    }
    diag(bout$yy) <- 0
    #adjust for the number that actually converged
    effB <- sum(bout$ngood)
    #vote based on effective bootstraps replicates. 
    bv <- list(yy = (bout$yy > thresh*effB) + 0,
               xy = (bout$xy > thresh*effB) + 0)
    
    #degree distributions of ensembles
    degxy <- foreach(bfit = bfits, .combine = 'rbind') %dopar% { 
      rowSums(as.matrix(bfit$xy))
    }
    degyy <- foreach(bfit = bfits, .combine = 'rbind') %dopar% { 
      rowSums(as.matrix(bfit$yy)) + colSums(as.matrix(bfit$xy))
    }
    bdeg <- list(yy = degyy, xy = degxy)
    return(list(bv = bv,
                bdeg = bdeg,
                bc = bout[!(names(bout) %in% "ngood")]))
  } else if (method == "space" & !givenX) { 
    
    addmodsyy <- function(m1, m2) { 
      list(yy = m1$yy + m2$yy, 
           dfyy = c(m1$dfyy, m2$dfyy),
           ngood = c(m1$ngood, m2$ngood))
    }
    
    bout <- foreach(bfit = bfits, .combine = addmodsyy, .packages = c("spacemap")) %dopar% { 
      byy <- as.matrix(bfit$yy)
      list(yy = byy, 
           dfyy = nonZeroUpper(byy, 0), ngood = bfit$convergence)
    }
    diag(bout$yy) <- 0
    #adjust for the number that actually converged
    effB <- sum(bout$ngood)
    
    #vote based on effective bootstraps replicates. 
    bv <- (bout$yy > thresh*effB) + 0
    
    #degree distributions of ensembles
    bdeg <- foreach(bfit = bfits, .combine = 'rbind') %dopar% { 
      rowSums(as.matrix(bfit$yy))
    }
    
    return(list(bv = list(yy = bv),
                bdeg = bdeg,
                bc = bout[!(names(bout) %in% "ngood")]))
  } else { 
    stop("Wrong method specified.")
  }
}

######################################################
#OLD code
######################################################

# sampIndPairs <- function(reps = 10) { 
#   pool <- combn(x = reps, 2)
#   begin <- match(1:(reps - 1),pool[1,])
#   last <- c(begin[2:length(begin)] - 1,tail(begin,1))
#   blocks <- lapply(seq_along(begin), function(i) begin[i]:last[i])
#   blocks <- blocks[sample.int(n = length(blocks))]
#   indpairs <- matrix(NA, nrow = 2, ncol = reps/2)
#   pairmembers <- c()
#   paircounter <- 1
#   for(bl in blocks) { 
#     if(any(pool[1,bl[1]] %in% pairmembers)) next;
#     candidx <- sample(x = bl, size = 1)
#     candpair <- pool[,candidx]
#     go_next <- FALSE
#     while (any(candpair %in% pairmembers)) {
#       bl <- setdiff(bl, candidx)
#       if(length(bl) == 0) { 
#         go_next <- TRUE
#         break;
#       } 
#       candidx <- sample(x = bl, size = 1)
#       candpair <- pool[,candidx]
#     }
#     if(go_next) next; 
#     pairmembers <- c(pairmembers, candpair)
#     indpairs[,paircounter] <- candpair
#     paircounter <- paircounter + 1
#   }
#   indpairs
# }
# 
# commonPairwiseBoot <- function(boots, P, tol, mat, pwi = NULL) { 
#   
#   if(is.null(pwi)) {
#     pwi <- combn(seq_along(boots), 2)
#   }
#   bpw <- vector(mode = "numeric", length = P)
#   prev1 <- 1
#   bspmap1 <- boots[[prev1]]
#   b1 <- bspmap1[[mat]] > tol + 0
#   check1 <- all(bspmap1 == "Not Sparse") | all(bspmap1 == "No Convergence")
#   cur1 <- 1
#   prev2 <- 2
#   bspmap2 <- boots[[prev2]]
#   b2 <- bspmap2[[mat]] > tol + 0
#   check2 <- all(bspmap2 == "Not Sparse") | all(bspmap2 == "No Convergence")
#   cur2 <- 2
#   nna <- 0
#   
#   pb <- txtProgressBar(min = 1, max = ncol(pwi), style = 3)
#   bpw <- matrix(nrow = ncol(pwi), ncol = P)
#   for(j in 1:ncol(pwi)) {
#     
#     cur1 <- pwi[1,j]
#     if(cur1 != prev1) { 
#       prev1 <- cur1
#       bspmap1 <- boots[[cur1]]
#       check1 <- !(all(bspmap1 == "Not Sparse") | all(bspmap1 == "No Convergence"))
#       if(check1) { 
#         b1 <- bspmap1[[mat]] > tol + 0
#       } else { 
#         nna <- nna + 1
#         next;
#       }
#       
#     }
#     
#     cur2 <- pwi[2,j]
#     if(cur2 != prev2) { 
#       prev2 <- cur2
#       bspmap2 <- boots[[cur2]]
#       check2 <- !(all(bspmap2 == "Not Sparse") | all(bspmap2 == "No Convergence"))
#       if(check2) { 
#         b2 <- bspmap2[[mat]] > tol + 0
#       } else { 
#         nna <- nna + 1
#         next;
#       }
#     } 
#     bpw[j,] <- rowSums(as.matrix(b1*b2))
#     setTxtProgressBar(pb, j)
#   }  
#   close(pb)
#   list(bootpw = bpw, nna = nna)
# }