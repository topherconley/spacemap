cvVoteRDS <- function(mod_lambda, tol, major_thresh = 0.5, 
                   method = c("spacemap", "space", "scggm")) { 
  fold_fits <- lapply(mod_lambda, readRDS)
  npotential_fits <- length(fold_fits)
  nfits <- 0
  p <- nrow(as.matrix(fold_fits[[1]]$ThetaXY))
  q <- ncol(as.matrix(fold_fits[[1]]$ThetaXY))
  vote_xy <- matrix(0, nrow = p, ncol = q)
  vote_yy <- matrix(0, nrow = q, ncol = q)
  deg_xy <- matrix(0, nrow = npotential_fits, ncol = p)
  deg_yy <- matrix(0, nrow = npotential_fits, ncol = q)
  for (i in seq_along(fold_fits)) {
    if (method == "scggm") { 
      valid <- TRUE
    } else { 
      valid <- fold_fits[[i]]$convergence & fold_fits[[i]]$sparseFit  
    }
    
    if (!valid) { 
      next;
    }
    xy <- abs(as.matrix(fold_fits[[i]]$ThetaXY)) > tol
    yy <- abs(as.matrix(fold_fits[[i]]$ThetaYY)) > tol
    deg_xy[i,] <- rowSums(xy)
    deg_yy[i,] <- rowSums(yy)
    nfits <- nfits + 1
    vote_xy <- vote_xy + xy
    vote_yy <- vote_yy + yy
  }
  
  if (nfits < ceiling(major_thresh*length(fold_fits))) { 
    return(NA)
  } else { 
    cutoff <- floor(major_thresh*nfits)
    vote_xy <- matrix(as.numeric(vote_xy > cutoff), nrow = p, ncol = q)
    vote_yy <- matrix(as.numeric(vote_yy > cutoff), nrow = q, ncol = q)
  }
  list(xy = vote_xy, yy = vote_yy, nfits = nfits, deg_xy = deg_xy, deg_yy = deg_yy)
}

