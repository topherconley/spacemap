#' Identify scale of tuning penalties
#' 
#' When initially choosing tuning penalties, it can be 
#' challenging to find the appropriate scale. This function fits either 
#' \code{\link{spacemap}} or \code{\link{space}} once for a specified 
#' tuning grid on the whole data. It reports the corresponding number of
#' \eqn{y-y} edges and \eqn{x-y} edges for each tuning penalty set. 
#' Having a prior understanding of how sparse the network ought to 
#' be can help narrow the scale of the tuning grid 
#' based on the output of this function.   
#' 
#' 
#' @inheritParams cvVote
#' @return If \code{method=="spacemap"} or (\code{method=="space"} and \code{X!=NULL}),
#'  return a  data.frame where the first column \code{nyy} reports the 
#'  number of \eqn{y-y} edges  and the second column \code{nxy} reports 
#'  the number of \eqn{x-y} edges. Rows of the data.frame correspond to 
#'  the input parameter \code{tuneGrid}. 
#'  
#'  If \code{method=="space"} and \code{X==NULL}, return a vector 
#'  of the number of \eqn{y-y} edges. 
#' @seealso \code{\link{cvVote}}
initFit <- function(Y, X = NULL, tuneGrid, method = c("spacemap", "space"), iscale = TRUE, aszero  = 1e-6, ...) {
  #additional arguments
  library(foreach)
  
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

    nedges <- foreach(l = seq_len(nrow(tuneGrid)), .combine = 'rbind') %dopar% {
      fit <- spacemap(Y = Y, X = X, 
                      lam1 = tuneGrid$lam1[l], 
                      lam2 = tuneGrid$lam2[l], 
                      lam3 = tuneGrid$lam3[l], 
                      sig=opt$sig, rho = opt$rho, 
                      iter= opt$iter,
                      tol = opt$tol,
                      iscale = FALSE,
                      cdmax = opt$cdmax, ...)
      if (fit$convergence) {
        net <- adjacency(net = fit, aszero = aszero)
        matrix(c(nyy = nonZeroUpper(net$yy,0.0), 
                 nxy = nonZeroWhole(net$xy,0.0)),
               nrow = 1, ncol = 2)
      } else {
        matrix(NA, nrow = 1, ncol = 2)
      } # if else convergence
      
    } # end foreach over tuneGrid
    nedges <- as.data.frame(nedges)
    names(nedges) <- c("nyy", "nxy")
  } else if (method == "space") {
    
    combf <- ifelse(givenX, 'rbind', 'c')
    nedges <- foreach(l = seq_len(nrow(tuneGrid)), .combine = combf) %dopar% {
      
      fit <- spacemap::space(Y = XY, lam1 = tuneGrid$lam1[l], sridge = opt$sridge, 
                                   sig = opt$sig, rho = opt$rho, iter = opt$iter, 
                                   tol = opt$tol, cd_iter = opt$cd_iter, iscale = FALSE)
      if (fit$convergence) {
        net <- adjacency(net = fit, aszero = aszero)
        if(givenX) { 
          matrix(c(nyy = nonZeroUpper(net$yy[Yindex, Yindex],0.0), 
                   nxy = nonZeroWhole(net$yy[Xindex, Yindex],0.0)),
                 nrow = 1, ncol = 2)
        } else { 
          nonZeroUpper(net$yy,0.0)  
        }
      } else {
        if(givenX) { 
          matrix(NA, nrow = 1, ncol = 2)
        } else { 
          NA
        }
      }  # if else convergence
    }  # end foreach over tuneGrid
    if (givenX) { 
      nedges <- as.data.frame(nedges)
      names(nedges) <- c("nyy", "nxy")
    }
  } else {
    stop("Incorrect method arguments.")
  }
  nedges
}

#
