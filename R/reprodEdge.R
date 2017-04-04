#' Ensemble networks through boostrapping. 
#' 
#' Fits spacemap model to B bootstrapped data replicates to create an ensemble network. 
#' 
#' @inheritParams spacemap
#' @param tune Data.frame with three columns requiring names to be the tuning penalties
#' \code{slasso, rlasso, rgroup} that are input to \code{spacemap}. Each row is a tuning penalty set 
#' that indexes an ensemble network of \code{spacemap} fits to \code{B} boostrap replicates of \code{Y.m, X.m}. 
#' @param boot Logical. Default is \code{boot=TRUE} and implies sampling with replacement (i.e. boostrap). 
#' \code{boot=FALSE} will resample without replacement. 
#' @param B Positive integer denoting the number of model fits making up the ensemble. Default is 1000, but can 
#' be lowered to save time. 
#' @param seed Positive integer allowing user to set the random seed for reproducing the network ensemble. 
#' @param p0 Positive numeric not exceeding the default of 1, which represents the proportion of the original 
#' samples that will be down-sampled for each model fitting iteration of the ensemble. 
#' @param refitRidge Positive numeric value specifying a ridge penalty, which is applied when 
#' refitting the model to the data after having selected the network with \code{spacemap}. Defaults to 0.1, 
#' but can be set to zero for ordinary least squares refit. 
#' @return a list of length B, where the bth element contains the refit of \code{spacemap} model. The model output is identical to 
#' \code{\link{spacemap}}, but the partial correlations matrix \code{ParCor} and the regression coefficients matrix \code{Gamma}
#' are stored in a sparse matrix format, (dgCMatrix format), from the \code{\link{Matrix}} package. If the bth model is not 
#' sufficiently sparse, a character vector will be returned as "Not Sparse". If the bth model did not converge, a 
#' character vector will be returned as "No Convergence". 
#' 
#' @examples 
#' #Load simulation
#' data(sim1)
#' #Boostrap Ensemble (B = 10) for illustration only
#' tune <- data.frame(slasso = 70, rlasso = 28.8, rgroup = 14)
#' out <- bootEnsemble(Y.m = sim1$Y, X.m = sim1$X, tune = tune, B = 10, iter = 3,
#'                     cd_iter = 1e7, tol = 1e-6, refitRidge = 0)
#' bv <- bootVote(fits = out[[1]], P = ncol(sim1$X), Q = ncol(sim1$Y))
#' @export
bootEnsemble <- function(Y.m, X.m, tune, boot = TRUE, B =  1000L, seed = 55139, p0 = 1.0, iter = 3, 
                       tol = 1e-4, cd_iter = 100e7, refitRidge = 0.1, iscale = FALSE) {
  
  nt <- nrow(tune)
  library(doRNG)
  library(Matrix)
  rng <- RNGseq(B, seed)
  N <- nrow(Y.m)
  P <- ncol(X.m)
  Q <- ncol(Y.m)
  rN <- floor(p0*N)
  
  foreach(i = seq_len(nt)) %:% foreach(j = seq_len(B),  r = rng[seq_len(B)]) %dopar% {
    # set RNG seed
    rngtools::setRNG(r)
    sidx <- sample.int(n = N, size = rN, replace = boot)
    nuniq <- length(unique(sort(sidx)))
    fit <- spacemap::spacemap(Y.m = Y.m[sidx,], X.m = X.m[sidx,], 
                              slasso = tune$slasso[i], rlasso = tune$rlasso[i], rgroup = tune$rgroup[i], iter = iter, 
                              tol = tol, cd_iter = cd_iter, iscale = iscale)
    #zero out those below tolerance
    fit$ParCor[abs(fit$ParCor) <= tol] <- 0.0
    fit$Gamma[abs(fit$Gamma) <= tol] <- 0.0
    if (!fit$convergence) return("No Convergence")
    sparseFit <- spacemap::olsRefitSpacemap(Y = Y.m[sidx,], X = X.m[sidx,], 
                                            ParCor = fit$ParCor, Gamma = fit$Gamma, 
                                            sigma = fit$sig.fit, RSS = fit$rss, tol = tol, 
                                            ridge = refitRidge)
    if (!sparseFit) return("Not Sparse")
    
    #Sparsify output
    fit$Gamma <- Matrix::Matrix(fit$Gamma, sparse = TRUE)
    fit$ParCor <- Matrix::Matrix(fit$ParCor, sparse = TRUE)
    fit
  }
}  


reprodEdgeSpace <- function(Y.m, tune, boot = TRUE, B =  1000L, seed = 55139, p0 = 1.0, iter = 3, 
                       tol = 1e-4, cd_iter = 100e7, refitRidge = 0.1, sridge = 0) {
  
  nt <- nrow(tune)
  library(doRNG)
  library(Matrix)
  rng <- RNGseq(B, seed)
  N <- nrow(Y.m)
  Q <- ncol(Y.m)
  rN <- floor(p0*N)
  
  foreach(i = seq_len(nt)) %:% foreach(j = seq_len(B),  r = rng[seq_len(B)]) %dopar% {
    # set RNG seed
    rngtools::setRNG(r)
    sidx <- sample.int(n = N, size = rN, replace = boot)
    nuniq <- length(unique(sort(sidx)))
    
    fit <- spacemap::space.joint(Y.m = Y.m[sidx,], lam1 = tune$slasso[i], lam2 = sridge, iter = iter, 
                                 cd_iter = cd_iter, tol = tol)
    #zero out those below tolerance. 
    fit$ParCor[abs(fit$ParCor) <= tol] <- 0.0
    if (!fit$convergence) return("No Convergence")
    
    sparseFit <- spacemap::olsRefitSpace(Y = Y.m[sidx,],
                                         ParCor = fit$ParCor, 
                                         sigma = fit$sig.fit, RSS = fit$rss, tol = tol,
                                         ridge = refitRidge)
    if (!sparseFit) return("Not Sparse")

    #Sparsify output
    fit$ParCor <- Matrix::Matrix(fit$ParCor, sparse = TRUE)
    fit
  }
}  

############################################
#Evaluate Statistics


#combining reps
addmods <- function(m1, m2) { 
  list(Gamma = m1$Gamma + m2$Gamma, ParCor = m1$ParCor + m2$ParCor, 
       dfGamma = c(m1$dfGamma, m2$dfGamma), dfParCor = c(m1$dfParCor, m2$dfParCor),
       nuniq = c(m1$nuniq, m2$nuniq))
}


evalReprodEdge <- function(Y.m, X.m, tune, boot = TRUE, B =  1000L, seed = 55139, p0 = 1.0, iter = 3, 
                       tol = 1e-4, cd_iter = 100e7, refitRidge = 0.1) {
  
  nt <- nrow(tune)
  library(doRNG)
  rng <- RNGseq(B, seed)
  N <- nrow(Y.m)
  P <- ncol(X.m)
  Q <- ncol(Y.m)
  rN <- floor(p0*N)
  
  foreach(i = seq_len(nt)) %:% foreach(j = seq_len(B),  r = rng[seq_len(B)], .combine = addmods) %dopar% {
    # set RNG seed
    rngtools::setRNG(r)
    sidx <- sample.int(n = N, size = rN, replace = boot)
    nuniq <- length(unique(sort(sidx)))
    fit <- spacemap::spacemap(Y.m = Y.m[sidx,], X.m = X.m[sidx,], 
                              slasso = tune$slasso[i], rlasso = tune$rlasso[i], rgroup = tune$rgroup[i], iter = iter, 
                              tol = tol, cd_iter = cd_iter)
    #zero out those below tolerance
    fit$ParCor[abs(fit$ParCor) <= tol] <- 0.0
    fit$Gamma[abs(fit$Gamma) <= tol] <- 0.0
    if (!fit$convergence) return(list("Gamma" = matrix(0.0,P,Q), "ParCor" = matrix(0.0,Q,Q),  
                                      "dfGamma" = "No Convergence", "dfParCor" = "No Convergence", "nuniq" = nuniq))
    sparseFit <- spacemap::olsRefitSpacemap(Y = Y.m[sidx,], X = X.m[sidx,], 
                                            ParCor = fit$ParCor, Gamma = fit$Gamma, 
                                            sigma = fit$sig.fit, RSS = fit$rss, tol = tol, 
                                            ridge = refitRidge)
    if (!sparseFit) return(list("Gamma" = matrix(0.0,P,Q), "ParCor" = matrix(0.0,Q,Q),  
                                "dfGamma" = "Not Sparse", "dfParCor" = "Not Sparse", "nuniq" = nuniq))
    list("Gamma" = abs(fit$Gamma) > 0.0, "ParCor" = abs(fit$ParCor) > 0.0, 
         "dfGamma" = spacemap::nonZeroWhole(fit$Gamma, 0.0), 
         "dfParCor" = spacemap::nonZeroUpper(fit$ParCor, 0.0), "nuniq" = nuniq)
  }
}

bootDegree <- function(boots, P, Q, type = c("ParCor", "Gamma"), tol) { 
  library(foreach)
  
  D <- NULL
  if (type == "ParCor") { 
    D <- Q
  } else if (type == "Gamma") { 
    D <- P
  } else { 
    return(message("Did not specify correct type."))
  }
  
  bdeg <- foreach(bfit = boots, .combine = 'rbind') %dopar% { 
    library(Matrix)
    if (all(bfit == "Not Sparse")) { 
      rep(NA, times = D)  
    } else if (bfit$convergence) { 
      xy <- as.matrix(bfit[[type]])
      if (type == "ParCor") { 
        diag(xy) <- 0
      } 
      rowSums(abs(xy) > tol)
    } else { 
      rep(NA, times = D)  
    }
  }
  bdeg
}


#' Model aggregation through Boot.Vote
#' 
#' Aggregate boostrap replicates of spacemap into a final Boot.Vote model. 
#' 
#' @param fits List of fitted spacemap models returned from \code{\link{bootEnsemble}}.
#' @param P Integer for the number of predictor variables
#' @param Q Integer for the number of response variables
#' @param vote_thresh  Positive numeric threshold for the minimum proportion of bootstrap 
#' replicate model fits with a particular edge such that the edge is included in the Boot.Vote model. 
#' @param btol The numeric lower bound that indicates parameter estimates equal to or less than \code{btol} ought 
#' to not be considered as an edge in a model fit. Do not change from zero if used in tandem 
#' with \code{\link{bootEnsemble}}.
#' @return List containing two lists. 
#' First list is \code{bv}, which encodes the edges in two logical adjacency matrices.
#' 
#' \enumerate{
#'   \item \code{ParCor} Logical matrix where TRUE for the (q,l) off-diagonals element indicate an edge 
#'   between the qth and lth response variables. 
#'   \item \code{Gamma}  Logical matrix where TRUE for the (p,q) element indicate an edge 
#'   between the pth predictor and qth response variable. 
#' }
#' 

#'  Second list is \code{bc}, which stores several additional statistics on the ensemble network fits. 
#' \enumerate{ 
#'    \item \code{ParCor} Integer matrix containing the y--y edge selection frequency out of B replicates.
#'    \item \code{Gamma}  Integer matrix containing the x--y edge selection frequency out of B replicates.
#'    \item \code{dfParCor} Integer vector containing the total number of y--y edges for each fit.
#'    \item \code{dfGamma} Integer vector containing the total number of x--y edges for each fit.
#'    \item \code{nuniq} Logical vector where TRUE indicates the model fit converged and was sparse and FALSE otherwise. 
#' }
#' 
#' @examples 
#' #Load simulation
#' data(sim1)
#' #Boostrap Ensemble (B = 10) for illustration only
#' tune <- data.frame(slasso = 70, rlasso = 28.8, rgroup = 14)
#' out <- bootEnsemble(Y.m = sim1$Y, X.m = sim1$X, tune = tune, B = 10, iter = 3,
#'                     cd_iter = 1e7, tol = 1e-6, refitRidge = 0)
#' bv <- bootVote(fits = out[[1]], P = ncol(sim1$X), Q = ncol(sim1$Y))
#' @export
bootVote <- function(fits, P, Q, vote_thresh = 0.5, btol = 0.0) {
  
  library(spacemap)
  library(foreach)
  boot_spmap <- foreach(bspmap = fits, .combine = spacemap::addmods) %do% { 
    if (all(bspmap == "Not Sparse") | all(bspmap == "No Convergence")) { 
      list("Gamma" = matrix(0.0,P,Q), "ParCor" = matrix(0.0,Q,Q),  
           "dfGamma" = "No Convergence", "dfParCor" = "No Convergence", "nuniq" = FALSE)
    } else if (bspmap$convergence) { 
      bxy <- abs(as.matrix(bspmap$Gamma)) > btol
      byy <- abs(as.matrix(bspmap$ParCor)) > btol
      list(Gamma = bxy, ParCor = byy, 
           dfGamma = nonZeroWhole(bxy, btol), 
           dfParCor <- nonZeroUpper(byy, btol), nuniq = bspmap$convergence)
    } else { 
      list("Gamma" = matrix(0.0,P,Q), "ParCor" = matrix(0.0,Q,Q),  
           "dfGamma" = "No Convergence", "dfParCor" = "No Convergence", "nuniq" = FALSE)
    }
  }
  diag(boot_spmap$ParCor) <- 0
  
  #adjust for the number that actually converged
  effB <- sum(boot_spmap$nuniq)
  
  #vote based on effective bootstraps replicates. 
  bv <- list(ParCor = boot_spmap$ParCor > vote_thresh*effB,
             Gamma = boot_spmap$Gamma > vote_thresh*effB)
  
  list(bv = bv,
       bc = boot_spmap)
}



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