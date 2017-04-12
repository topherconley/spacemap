#' Ensemble networks through boostrapping. 
#' 
#' Fits spacemap model to B bootstrapped data replicates to create an ensemble network. 
#' 
#' @inheritParams spacemap
#' @param tune List with names being the tuning penalties
#' \code{lam1, lam2, lam3} and each element being a single numeric value that is input to \code{spacemap}. One might  
#' consider the output from \code{\link{cvVote}} (i.e. list \code{minTune}) as input for this parameter. 
#' @param boot Logical. Default is \code{boot=TRUE} and implies sampling with replacement (i.e. boostrap). 
#' \code{boot=FALSE} will resample without replacement. 
#' @param B Positive integer denoting the number of model fits making up the ensemble. Default is 1000, but can 
#' be lowered to save time. 
#' @param seed Positive integer allowing user to set the random seed for reproducing the network ensemble. 
#' @param p0 Positive numeric not exceeding the default of 1, which represents the proportion of the original 
#' samples that will be down-sampled for each model fitting iteration of the ensemble. 
#' @return A list up to length \code{B} of convergent and sparse model fits from 
#' either \code{\link{spacemap}} or \code{\link{space.joint}}. The list object should typically
#' not be modified by the user, but passed on to the  \code{\link{bootVote}} function.
#' Model fits that are not convergent or produce networks that are not sparse are not reported. 
#' 
#' @examples 
#' #Load simulation
#' data(sim1)
#' #Boostrap Ensemble (B = 10) for illustration only
#' tune <- data.frame(lam1 = 70, lam2 = 28.8, lam3 = 14)
#' out <- bootEnsemble(Y = sim1$Y, X = sim1$X, tune = tune, B = 10, iter = 3,
#'                     cdmax = 1e7, tol = 1e-6)
#' @export
#' @seealso bootVote
#' 
#' 
bootEnsemble <- function(Y, X = NULL, tune, method = c("spacemap", "space"), 
                         boot = TRUE, B =  1000L, iscale = TRUE,  
                         seed = sample.int(n = 1, size = 1e6), 
                         p0 = 1.0, ...) {
  
  method <- match.arg(method)
  library(doRNG)
  library(Matrix)
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
    ens <- foreach(j = seq_len(B),  r = rng[seq_len(B)]) %dopar% {
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
      fit$ParCor[abs(fit$ParCor) <= tol] <- 0.0
      fit$Gamma[abs(fit$Gamma) <= tol] <- 0.0
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
    ens <- foreach(j = seq_len(B),  r = rng[seq_len(B)]) %dopar% {
      rngtools::setRNG(r)
      sidx <- sample.int(n = N, size = rN, replace = boot)
      fit <- spacemap::space.joint(Y.m = XY[sidx,], lam1 = tune$lam1[1], 
                                   lam2 = opt$sridge, iter = opt$iter, 
                                   cdmax = opt$cdmax, tol = opt$tol, 
                                   iscale = FALSE)
      #zero out those below tolerance. 
      fit$ParCor[abs(fit$ParCor) <= tol] <- 0.0
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
#' @param givenX Logical. Defaults to FALSE. Should be set to TRUE when 
#' \code{attr(bfits, "method) == "space} and \code{space.joint} was used
#' to infer (x--x, x--y, y--y) edges  but only reported (x--y, y--y) edges. 
#' @param thresh Positive numeric threshold for the minimum proportion of bootstrap 
#' replicate model fits with a particular edge such that the edge is included in the Boot.Vote model. 
#' @return List containing two lists. 
#' First list is \code{bv}, which encodes the edges in two logical adjacency matrices.
#' 
#' \enumerate{
#'   \item \code{yy} Logical matrix where TRUE for the (q,l) off-diagonals element indicate an edge 
#'   between the qth and lth response variables. 
#'   \item \code{xy}  Logical matrix where TRUE for the (p,q) element indicate an edge 
#'   between the pth predictor and qth response variable. 
#' }
#' 

#'  Second list is \code{bc}, which stores several additional statistics on the ensemble network fits. 
#' \enumerate{ 
#'    \item \code{yy} Integer matrix containing the y--y edge selection frequency out of B replicates.
#'    \item \code{xy}  Integer matrix containing the x--y edge selection frequency out of B replicates.
#'    \item \code{dfyy} Integer vector containing the total number of y--y edges for each fit.
#'    \item \code{dfxy} Integer vector containing the total number of x--y edges for each fit.
#' }
#' 
#' @examples 
#' #Load simulation
#' data(sim1)
#' #Boostrap Ensemble (B = 10) for illustration only
#' tune <- data.frame(lam1 = 70, lam2 = 28.8, lam3 = 14)
#' out <- bootEnsemble(Y = sim1$Y, X = sim1$X, tune = tune, B = 10, iter = 3,
#'                     cdmax = 1e7, tol = 1e-6, refitRidge = 0)
#' bv <- bootVote(fits = out[[1]], P = ncol(sim1$X), Q = ncol(sim1$Y))
#' @export
bootVote <- function(bfits, givenX = FALSE, thresh = 0.5) {
  
  if(!is.list(bfits)) { 
    stop("bfits is not a list")
  }

  method <- attr(bfits, "method")
  
  library(spacemap)
  library(foreach)
  
  if (method == "spacemap" | (method == "space" & givenX)) { 
    bout <- foreach(bfit = bfits, .combine = spacemap::addmods) %dopar% { 
      list(xy = as.matrix(bfit$xy), yy = as.matrix(bfit$yy), 
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
      rowSums(as.matrix(bfit$yy))
    }
    bdeg <- list(xy = degxy, yy = degyy)
    return(list(bv = bv,
                bdeg = bdeg,
                bc = bout))
  } else if (method == "space" & !givenX) { 
    bout <- foreach(bfit = bfits, .combine = spacemap::addmodsyy) %dopar% { 
      list(yy = as.matrix(bfit$yy), 
           dfyy = nonZeroUpper(byy, 0), ngood = bfit$convergence)
    }
    diag(bout$yy) <- 0
    #adjust for the number that actually converged
    effB <- sum(bout$ngood)
    #should be equal
    stopifnot(effB == length(bout))
    
    #vote based on effective bootstraps replicates. 
    bv <- (bout$yy > thresh*effB) + 0

    #degree distributions of ensembles
    bdeg <- foreach(bfit = bfits, .combine = 'rbind') %dopar% { 
      rowSums(as.matrix(bfit$yy))
    }
    
    return(list(bv = bv,
                bdeg = bdeg,
                bc = bout))
  } else { 
    stop("Wrong method specified.")
  }
}

#combining reps
addmods <- function(m1, m2) { 
  list(xy = m1$xy + m2$xy, yy = m1$yy + m2$yy, 
       dfxy = c(m1$dfxy, m2$dfxy), dfyy = c(m1$dfyy, m2$dfyy),
       ngood = c(m1$ngood, m2$ngood))
}

addmodsyy <- function(m1, m2) { 
  list(yy = m1$yy + m2$yy, 
       dfyy = c(m1$dfyy, m2$dfyy),
       ngood = c(m1$ngood, m2$ngood))
}

######################################################
#OLD code
######################################################

sampIndPairs <- function(reps = 10) { 
  pool <- combn(x = reps, 2)
  begin <- match(1:(reps - 1),pool[1,])
  last <- c(begin[2:length(begin)] - 1,tail(begin,1))
  blocks <- lapply(seq_along(begin), function(i) begin[i]:last[i])
  blocks <- blocks[sample.int(n = length(blocks))]
  indpairs <- matrix(NA, nrow = 2, ncol = reps/2)
  pairmembers <- c()
  paircounter <- 1
  for(bl in blocks) { 
    if(any(pool[1,bl[1]] %in% pairmembers)) next;
    candidx <- sample(x = bl, size = 1)
    candpair <- pool[,candidx]
    go_next <- FALSE
    while (any(candpair %in% pairmembers)) {
      bl <- setdiff(bl, candidx)
      if(length(bl) == 0) { 
        go_next <- TRUE
        break;
      } 
      candidx <- sample(x = bl, size = 1)
      candpair <- pool[,candidx]
    }
    if(go_next) next; 
    pairmembers <- c(pairmembers, candpair)
    indpairs[,paircounter] <- candpair
    paircounter <- paircounter + 1
  }
  indpairs
}

commonPairwiseBoot <- function(boots, P, tol, mat, pwi = NULL) { 
  
  if(is.null(pwi)) {
    pwi <- combn(seq_along(boots), 2)
  }
  bpw <- vector(mode = "numeric", length = P)
  prev1 <- 1
  bspmap1 <- boots[[prev1]]
  b1 <- bspmap1[[mat]] > tol + 0
  check1 <- all(bspmap1 == "Not Sparse") | all(bspmap1 == "No Convergence")
  cur1 <- 1
  prev2 <- 2
  bspmap2 <- boots[[prev2]]
  b2 <- bspmap2[[mat]] > tol + 0
  check2 <- all(bspmap2 == "Not Sparse") | all(bspmap2 == "No Convergence")
  cur2 <- 2
  nna <- 0
  
  pb <- txtProgressBar(min = 1, max = ncol(pwi), style = 3)
  bpw <- matrix(nrow = ncol(pwi), ncol = P)
  for(j in 1:ncol(pwi)) {
    
    cur1 <- pwi[1,j]
    if(cur1 != prev1) { 
      prev1 <- cur1
      bspmap1 <- boots[[cur1]]
      check1 <- !(all(bspmap1 == "Not Sparse") | all(bspmap1 == "No Convergence"))
      if(check1) { 
        b1 <- bspmap1[[mat]] > tol + 0
      } else { 
        nna <- nna + 1
        next;
      }
      
    }
    
    cur2 <- pwi[2,j]
    if(cur2 != prev2) { 
      prev2 <- cur2
      bspmap2 <- boots[[cur2]]
      check2 <- !(all(bspmap2 == "Not Sparse") | all(bspmap2 == "No Convergence"))
      if(check2) { 
        b2 <- bspmap2[[mat]] > tol + 0
      } else { 
        nna <- nna + 1
        next;
      }
    } 
    bpw[j,] <- rowSums(as.matrix(b1*b2))
    setTxtProgressBar(pb, j)
  }  
  close(pb)
  list(bootpw = bpw, nna = nna)
}