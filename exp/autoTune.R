

autoTuneSpace <- function(maxHours = 0.5, maxTuneIter = 2, ngrid = 25) { 
  
  #total hours computed
  th <- 0.0
  #seconds in an hour
  sih <- 60^2
  #current tune iteration
  cti <- 1
  
  
  ltunes <- list()
  
  #initial space tuning parameter
  lam1start <- function(n, q, alpha) { 
    sqrt(n) * qnorm(1 - (alpha/ (2*q^2)))
  }
  lam0 <- lam1start(n = floor(N - N*.10), q = Q, alpha = 0.01)
  tsp <- expand.grid(slasso = seq(lam0*(1 - eps.10), 2.5*lam0, length = ngrid))
  
  while(TRUE) {
    
    
    it <- system.time({cvSpaceYY <- spacemap::crossValidation(dat = datyy, fold = K, 
                                                              trainIds = trainSetIds, testIds = testSetIds, 
                                                              method = "space", tuneGrid = tsp,
                                                              sridge = sridge, tol = tol, iter = iter, cd_iter = cd_iter, 
                                                              residConditional = FALSE, 
                                                              refitRidge = refitRidge)})[3]
    th <- th +  it
    ltunes[[cti]]$cvobj <- cvSpaceYY
    ltunes[[cti]]$iter_time <- it/sih
    ltunes[[cti]]$tg <- tsp
    #is the tuning on the boundary? 
    stp <- cvSpaceYY$tuneGrid
    max_slasso <- max(tsp$slasso)
    min_slasso <- min(tsp$slasso)
    #define tuning grid based on past tuning iteration
    
    #upper boundary
    if(stp == max_slasso) { 
      tsp <- expand.grid(slasso = seq(max_slasso + 0.5, 1.5*max_slasso, length = ngrid))
    } else if (stp == min_slasso){ 
      tsp <- expand.grid(slasso = seq(0.75*min_slasso, min_slasso - 0.5, length = ngrid))
    } else { 
      tsp <- expand.grid(slasso = seq(stp*0, stp, length = ngrid))
    }
    
   
    cti <- cti + 1
    if ((th/sih) > maxHours | cti >  maxTuneIter) break;
  }
  
}