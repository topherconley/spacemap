cvVoteRDS <- function(files, tol, thresh = 0.5, 
                   method = c("spacemap", "space", "scggm"), 
                   givenX = TRUE) { 
  method <- match.arg(method)
  fold_fits <- lapply(files, readRDS)
  npotential_fits <- length(fold_fits)
  nfits <- 0
  
  #vanilla space (no conditioning on x)
  part_xy <- !(method == "space" & !givenX)

  q <- ncol(as.matrix(fold_fits[[1]]$ThetaYY))
  vote_yy <- matrix(0, nrow = q, ncol = q)
  deg_yy <- matrix(0, nrow = npotential_fits, ncol = q)
  if(part_xy) { 
    p <- nrow(as.matrix(fold_fits[[1]]$ThetaXY))
    vote_xy <- matrix(0, nrow = p, ncol = q)
    deg_xy <- matrix(0, nrow = npotential_fits, ncol = p)
  }
  
  for (i in seq_along(fold_fits)) {
    if (method == "scggm") { 
      valid <- TRUE
    } else { 
      valid <- fold_fits[[i]]$convergence & fold_fits[[i]]$sparseFit  
    }
    
    if (!valid) { 
      next;
    }
    
    yy <- abs(as.matrix(fold_fits[[i]]$ThetaYY)) > tol
    deg_yy[i,] <- rowSums(yy)
    nfits <- nfits + 1
    vote_yy <- vote_yy + yy
    
    if (part_xy) { 
      xy <- abs(as.matrix(fold_fits[[i]]$ThetaXY)) > tol
      deg_xy[i,] <- rowSums(xy)
      vote_xy <- vote_xy + xy
    }
  }
  
  if (nfits < ceiling(thresh*length(fold_fits))) { 
    return(NA)
  } else { 
    cutoff <- floor(thresh*nfits)
    vote_yy <- (vote_yy > cutoff) + 0
    if (part_xy) { 
      vote_xy <- (vote_xy > cutoff) + 0
    }
  }
  if (part_xy) { 
    ret <- list(xy = vote_xy, yy = vote_yy, nfits = nfits, deg_xy = deg_xy, deg_yy = deg_yy)
  } else { 
    ret <- list(yy = vote_yy, nfits = nfits, deg_yy = deg_yy)
  }
  ret
}

